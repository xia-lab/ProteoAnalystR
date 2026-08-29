## =========================================================================
## analysis_summary.txt  --  consolidated, human-readable analysis report
## -------------------------------------------------------------------------
## Writes a single plain-text summary into the session working directory. The
## Downloads page lists every *.txt in that directory (DownloadBean.OnlyExt),
## so this file is offered for download with no Java/JSF change beyond the call
## that triggers it.
##
## The file opens with the dataset "Text Summary" block (the same fields shown
## on the data-input page) and then carries one section per analysis that has
## actually been run. Every section is reconstructed from persisted state
## (paramSet / msgSet / analSet / dataSet) on each call, so the file always
## reflects whatever the user has done so far and is regenerated whenever the
## Downloads page is opened.
##
## Reporting fields requested in the reviewer revision are stated explicitly:
## number of features tested, raw and adjusted p-value columns, FDR threshold,
## fold-change threshold, enrichment multiple-testing correction, and the
## background/universe feature set used for enrichment.
##
## Every accessor is defensive: a missing field yields "NA" or an omitted
## section and never aborts the calling analysis.
## =========================================================================

# Safe field getter for a list/environment-like object.
.pa_g <- function(obj, field, default = "NA") {
  tryCatch({
    if (is.null(obj)) return(default)
    v <- obj[[field]]
    if (is.null(v) || length(v) == 0) return(default)
    if (is.character(v) && length(v) == 1 && !nzchar(v)) return(default)
    v
  }, error = function(e) default)
}

# "Label:                         value"  (label padded to a fixed width).
.pa_line <- function(label, value) {
  value <- tryCatch(paste(as.character(value), collapse = "; "),
                    error = function(e) "NA")
  if (is.null(value) || length(value) == 0 || !nzchar(value)) value <- "NA"
  sprintf("%-34s%s", paste0(label, ":"), value)
}

# nrow() that tolerates matrices, data.frames and NULL.
.pa_nrow <- function(x) {
  tryCatch({
    if (is.null(x)) return(NA_integer_)
    d <- dim(x)
    if (!is.null(d)) return(d[1])
    length(x)
  }, error = function(e) NA_integer_)
}

# Readable label for an ID-type code (the code is what R stores; the friendly
# label otherwise lives only in the Java UI). Falls back to the raw code.
.pa_idtype_label <- function(code) {
  code <- tolower(as.character(code))
  if (length(code) == 0 || !nzchar(code) || identical(code, "na")) return("NA")
  map <- c(uniprot = "UniProt Accession", entrez = "Entrez Gene ID",
           symbol = "Official Gene Symbol", genesymbol = "Official Gene Symbol",
           refseq = "RefSeq", embl = "EMBL", genbank = "GenBank",
           ensembl = "Ensembl Gene", ensp = "Ensembl Protein",
           ensg = "Ensembl Gene", enst = "Ensembl Transcript",
           kegg = "KEGG", none = "NA")
  if (!is.na(map[code])) unname(map[code]) else code
}

# Readable label for a PPI database code (NetController.getPpiOpts()).
.pa_ppidb_label <- function(code) {
  code <- tolower(as.character(code))
  if (length(code) == 0 || !nzchar(code) || identical(code, "na")) return("NA")
  map <- c(string = "STRING", innate = "InnateDB", innatedb = "InnateDB",
           rolland = "Rolland", huri = "HuRI", intact = "IntAct",
           irefinx = "iRefIndex", biogrid = "BioGRID",
           tissue = "Tissue-specific", `cross-species` = "Host-parasite/microbiome")
  if (!is.na(map[code])) unname(map[code]) else code
}

## ---- Dataset block ------------------------------------------------------
.pa_dataset_section <- function(paramSet, msgSet, dataSet) {
  sv <- .pa_g(msgSet, "summaryVec", NULL)
  # summaryVec (1-based): 1 matched, 2 percent, 3 total, 4 unmatched,
  # 5 samples, 6 factors, 7 totalCount, 8 avg, 9 min, 10 max, 11 groups,
  # 12 missing-count.  (See data_idanot_utils.R.)
  gsv <- function(i) if (!is.null(sv) && length(sv) >= i) sv[[i]] else "NA"

  # dataSet$type is the reliable per-dataset type ("prot"/"count"/"array");
  # paramSet$data.type is only set for some flows (e.g. "phospho").
  data.type <- as.character(.pa_g(dataSet, "type", .pa_g(paramSet, "data.type")))
  data.type.readable <- switch(data.type,
                               prot    = "Proteomics data",
                               count   = "RNA-seq count data",
                               array   = "Microarray data",
                               phospho = "Phosphoproteomics data",
                               data.type)
  data.format <- tolower(paste(as.character(.pa_g(paramSet, "data.format", "")), collapse = ""))
  quant <- if (grepl("phospho", data.format)) "Phosphosite-level" else
           if (grepl("peptide", data.format)) "Peptide-level" else "Protein-level"

  miss.pct <- "NA"
  tryCatch({
    mc <- suppressWarnings(as.numeric(gsv(12)))
    fc <- suppressWarnings(as.numeric(gsv(3)))
    sc <- suppressWarnings(as.numeric(gsv(5)))
    if (!is.na(mc) && !is.na(fc) && !is.na(sc) && (fc * sc) > 0)
      miss.pct <- sprintf("%.2f%%", 100 * mc / (fc * sc))
  }, error = function(e) NULL)

  c("== Dataset ==",
    .pa_line("Data type", data.type.readable),
    .pa_line("Filename", .pa_g(dataSet, "name")),
    .pa_line("Quantification level", quant),
    .pa_line("Organism", .pa_g(paramSet, "data.org", .pa_g(paramSet, "org"))),
    .pa_line("ID type", .pa_idtype_label(.pa_g(paramSet, "data.idType",
                                               .pa_g(dataSet, "id.current")))),
    .pa_line("Total feature number", gsv(3)),
    .pa_line("Matched feature number", gsv(1)),
    .pa_line("Unmatched feature number", gsv(4)),
    .pa_line("Percent matched", gsv(2)),
    .pa_line("Sample number", gsv(5)),
    .pa_line("Number of experimental factors", gsv(6)),
    .pa_line("Group names", gsv(11)),
    .pa_line("Total missing percentage", miss.pct))
}

## ---- Data filtering + normalization ------------------------------------
## Parameters persisted by PerformNormalization() (norm_utils.R): the two
## percentile filter knobs (var.perc / abun.perc), the valid-value filter
## (remove.missing + missing.percent), the matrix normalization (norm.opt),
## and any sample-wise normalization (sample.norm.opt). Reported so the summary
## states exactly how the matrix was filtered and normalized before analysis.
.pa_filternorm_section <- function(paramSet, msgSet) {
  norm.opt  <- .pa_g(paramSet, "norm.opt", NULL)
  var.perc  <- .pa_g(paramSet, "var.perc", NULL)
  abun.perc <- .pa_g(paramSet, "abun.perc", NULL)
  rm.miss   <- .pa_g(paramSet, "remove.missing", NULL)
  # Emit nothing until normalization/filtering has actually been run.
  if (is.null(norm.opt) && is.null(var.perc) && is.null(abun.perc) && is.null(rm.miss))
    return(NULL)

  pct.filter <- function(v, by) {
    v <- suppressWarnings(as.numeric(v))
    if (length(v) == 0 || is.na(v)) return("NA")
    if (v <= 0) "not applied" else sprintf("bottom %g%% by %s removed", v, by)
  }
  rm.on   <- tryCatch(isTRUE(as.logical(rm.miss)), error = function(e) FALSE)
  miss.pc <- suppressWarnings(as.numeric(.pa_g(paramSet, "missing.percent", NA)))
  valid.line <- if (rm.on && is.finite(miss.pc))
                  sprintf("features with > %g%% missing values removed", miss.pc)
                else "not applied"

  filt <- c("== Data Filtering ==",
    .pa_line("Low-abundance filter", pct.filter(abun.perc, "mean abundance")),
    .pa_line("Low-variance filter",  pct.filter(var.perc, "variance")),
    .pa_line("Missing-value filter", valid.line))
  pf <- .pa_g(msgSet, "prefilter.msg", NULL)
  if (!is.null(pf) && !identical(as.character(pf)[1], "NA"))
    filt <- c(filt, .pa_line("Additional filtering", paste(pf, collapse = " ")))

  norm.readable <- switch(as.character(norm.opt),
    none      = "none",
    log       = "Log2 transformation",
    logMedian = "Log2 transformation + sample-median centering",
    vsn       = "Variance-stabilizing normalization (VSN)",
    quantile  = "Quantile normalization",
    MORlog    = "DESeq2 median-of-ratios size factors + log2(x+1)",
    msstats   = "MSstats equalizeMedians + Tukey median-polish",
    Rlr       = "Robust linear regression (log2)",
    Loess     = "Cyclic LOESS (log2)",
    as.character(norm.opt))
  sn <- as.character(.pa_g(paramSet, "sample.norm.opt", "none"))
  sn.readable <- if (!nzchar(sn) || sn %in% c("none", "NA")) "none" else switch(sn,
    SumNorm    = "Normalization to constant sum",
    MedianNorm = "Normalization to sample median",
    ProbNorm   = "Probabilistic Quotient Normalization (PQN)",
    SpecNorm   = "Reference-feature (specific) normalization",
    sn)
  imp <- as.character(.pa_g(paramSet, "impute.opt", "not run"))
  imp.readable <- switch(imp,
    min       = "MNAR limit-of-detection replacement (1/5 feature minimum)",
    mindet    = "MinDet (imputeLCMD)",
    minprob   = "MinProb (imputeLCMD)",
    qrilc     = "QRILC left-censored (imputeLCMD)",
    mean      = "Feature mean replacement",
    median    = "Feature median replacement",
    colmin    = "Feature minimum replacement",
    knn_var   = "KNN by feature (impute::impute.knn)",
    knn_smp   = "KNN by sample (impute::impute.knn)",
    seqknn    = "Sequential KNN",
    bpca      = "Bayesian PCA",
    ppca      = "Probabilistic PCA",
    svdImpute = "SVD imputation",
    impseq    = "Sequential covariance imputation",
    exclude   = "Exclude features containing missing values",
    imp)

  norm <- c("== Normalization ==",
    .pa_line("Normalization / transform", norm.readable),
    .pa_line("Sample-wise normalization", sn.readable),
    .pa_line("Protein/feature-level missing-value imputation", imp.readable))

  c(filt, "", norm)
}

## ---- Peptide-to-protein summarization + imputation ---------------------
.pa_summarization_section <- function(paramSet) {
  m <- .pa_g(paramSet, "summ.method", NULL)
  if (is.null(m)) return(NULL)  # not a peptide-level dataset / summarization not run
  method.readable <- switch(as.character(m),
    tukey         = "Tukey median polish (stats::medpolish per protein)",
    median_polish = "Median sweep (fast Tukey approximation)",
    top_n         = paste0("Top-N (Hi3: mean of the ",
                           .pa_g(paramSet, "summ.topn", 3),
                           " most-intense peptides; Silva et al. 2006)"),
    mean          = "Mean of log2 peptide intensities",
    median        = "Median of log2 peptide intensities",
    sum           = "Summed linear intensity (iBAQ-style)",
    as.character(m))
  imp <- as.character(.pa_g(paramSet, "summ.impute", "none"))
  imp.readable <- switch(imp,
    none    = "none (missingness handled by the differential model)",
    qrilc   = "QRILC left-censored (imputeLCMD)",
    minprob = "MinProb (imputeLCMD)",
    mindet  = "MinDet (imputeLCMD)",
    mnar    = "MNAR-adaptive (imputeLCMD impute.MAR.MNAR; MAR=KNN, MNAR=QRILC)",
    imp)
  c("== Peptide-to-Protein Summarization ==",
    .pa_line("Summarization method", method.readable),
    .pa_line("Minimum peptides per protein", .pa_g(paramSet, "summ.minpep")),
    .pa_line("Peptide-level imputation", imp.readable))
}

## ---- Differential expression (proteomics + phospho) --------------------
.pa_de_section <- function(paramSet) {
  dataName <- as.character(.pa_g(paramSet, "dataName", ""))
  if (!nzchar(dataName) || identical(dataName, "NA")) return(NULL)
  dataSet <- tryCatch(readDataset(dataName), error = function(e) NULL)
  if (is.null(dataSet)) return(NULL)

  comp <- .pa_g(dataSet, "comp.res", NULL)
  sig  <- .pa_g(dataSet, "sig.mat", NULL)
  if (is.null(comp) && is.null(sig)) return(NULL)  # DE has not been run

  n.tested <- .pa_nrow(comp)
  n.sig    <- .pa_nrow(sig)

  cn <- tryCatch(colnames(as.data.frame(comp)), error = function(e) character())
  pcol  <- intersect(c("P.Value", "PValue", "pvalue", "p.value"), cn)
  apcol <- intersect(c("adj.P.Val", "padj", "FDR", "adj.pval"), cn)
  pcol  <- if (length(pcol))  pcol[1]  else "P.Value"
  apcol <- if (length(apcol)) apcol[1] else "adj.P.Val"

  fdr    <- .pa_g(dataSet, "pval",   .pa_g(paramSet, "pvalu"))
  fc     <- .pa_g(dataSet, "fc.lvl", .pa_g(paramSet, "fc.thresh"))
  method.requested <- as.character(.pa_g(
    dataSet, "de.method.requested", .pa_g(dataSet, "de.method")
  ))
  method.effective <- as.character(.pa_g(
    dataSet, "de.method.effective", method.requested
  ))
  method.display <- if (!identical(method.requested, method.effective)) {
    paste0(method.requested, " requested; ", method.effective, " used")
  } else {
    method.effective
  }

  # Comparison / contrast being tested. Multi-group ("default") designs are an
  # overall F-test across all groups rather than a single pairwise contrast.
  comp.type <- as.character(.pa_g(dataSet, "comp.type", ""))
  comparison <- if (identical(comp.type, "default"))
                  "overall F-test across groups"
                else .pa_g(dataSet, "active.comp.nm")

  # Whether significance is called on the BH-adjusted p or the raw p is the
  # user's FDR toggle, persisted by GetSigfeatures() as paramSet$use.fdr.
  use.fdr <- tryCatch(isTRUE(as.logical(.pa_g(paramSet, "use.fdr", TRUE))),
                      error = function(e) TRUE)
  sig.line <- if (use.fdr) paste0(apcol, " (BH-FDR) < ", fdr)
              else            paste0(pcol, " (raw) < ", fdr)

  # up/down split for a single-contrast logFC (omitted for multi-group omnibus).
  sig.count <- n.sig
  updown <- tryCatch({
    sdf <- as.data.frame(sig)
    lf  <- intersect(c("logFC", "log2FC"), colnames(sdf))
    if (length(lf) == 1 && !identical(comp.type, "default")) {
      v <- suppressWarnings(as.numeric(sdf[[lf[1]]]))
      paste0(sig.count, "  (up ", sum(v > 0, na.rm = TRUE),
             " / down ", sum(v < 0, na.rm = TRUE), ")")
    } else as.character(sig.count)
  }, error = function(e) as.character(sig.count))

  c("== Differential Expression ==",
    .pa_line("Method", method.display),
    if (!is.null(.pa_g(dataSet, "deqms.count.source", NULL)))
      .pa_line("DEqMS count source", .pa_g(dataSet, "deqms.count.source")),
    if (!is.null(.pa_g(dataSet, "de.fallback.reason", NULL)) &&
        !is.na(.pa_g(dataSet, "de.fallback.reason", NA_character_)))
      .pa_line("Fallback reason", .pa_g(dataSet, "de.fallback.reason")),
    .pa_line("Comparison", comparison),
    .pa_line("Number of features tested", n.tested),
    .pa_line("Significance", sig.line),
    .pa_line("Fold-change threshold (log2)", fc),
    .pa_line("p-value columns in results", paste0(pcol, " (raw), ", apcol, " (BH-FDR)")),
    .pa_line("Significant features", updown))
}

## ---- Functional enrichment (ORA) ---------------------------------------
.pa_enrich_section <- function(paramSet, analSet) {
  enr <- .pa_g(analSet, "enr.mat", NULL)
  if (is.null(enr)) return(NULL)
  cn <- tryCatch(colnames(as.data.frame(enr)), error = function(e) character())
  pcol  <- intersect(c("Pval", "P.Value", "pvalue"), cn)
  apcol <- intersect(c("FDR", "adj.P.Val", "padj"), cn)

  bg <- .pa_g(paramSet, "universe.opt.readable", .pa_g(paramSet, "universe.opt"))
  bg.size <- tryCatch({
    u <- .pa_g(paramSet, "backgroundUniverse", NULL)
    if (!is.null(u) && length(u) > 1) length(u) else "NA"
  }, error = function(e) "NA")

  c("== Functional Enrichment (ORA) ==",
    .pa_line("Gene-set library", .pa_g(paramSet, "init.lib")),
    .pa_line("Background / universe", bg),
    .pa_line("Background feature count", bg.size),
    .pa_line("Number of terms tested", .pa_nrow(enr)),
    .pa_line("Raw p-value column", if (length(pcol)) pcol[1] else "Pval"),
    .pa_line("Adjusted p-value column", if (length(apcol)) apcol[1] else "FDR"),
    .pa_line("Multiple-testing correction", "Benjamini-Hochberg (FDR)"))
}

## ---- Gene-set enrichment (GSEA / fgsea) --------------------------------
.pa_gsea_section <- function(paramSet) {
  lib <- .pa_g(paramSet, "gsea.lib", NULL)
  if (is.null(lib)) return(NULL)
  c("== Gene-Set Enrichment (GSEA) ==",
    .pa_line("Gene-set library", lib),
    .pa_line("Ranking statistic", .pa_g(paramSet, "gsea.rank.opt",
                                        .pa_g(paramSet, "gseaRankOpt"))),
    .pa_line("Number of gene sets tested", .pa_g(paramSet, "gsea.n.tested")),
    .pa_line("Enrichment engine", "fgsea (permutation-based)"),
    .pa_line("Statistics reported", "NES, pval (raw), padj (BH-FDR)"),
    .pa_line("Multiple-testing correction", "Benjamini-Hochberg (FDR)"))
}

## ---- PTM / phosphosite occupancy ---------------------------------------
.pa_ptm_section <- function(paramSet) {
  n <- .pa_g(paramSet, "ptm.occ.n.sites", NULL)
  if (is.null(n)) return(NULL)
  c("== PTM / Protein-adjusted phosphosite change (relative; not stoichiometry) ==",
    .pa_line("Comparison", .pa_g(paramSet, "ptm.occ.contrast")),
    .pa_line("Adjustment model", "site log2FC - parent-protein log2FC (delta method)"),
    .pa_line("Degrees of freedom", "Satterthwaite approximation"),
    .pa_line("Sites analyzed", n),
    .pa_line("Columns in results", "Delta.Occupancy = protein-adjusted site log2FC (relative), Occ.Pvalue (raw), Occ.FDR (BH-FDR)"),
    .pa_line("Multiple-testing correction", "Benjamini-Hochberg (FDR)"))
}

## ---- Kinase enrichment (phospho) ---------------------------------------
.pa_kinase_section <- function(analSet) {
  ke <- .pa_g(analSet, "kinase.enrich", NULL)
  if (is.null(ke)) return(NULL)
  cn <- tryCatch(colnames(as.data.frame(ke)), error = function(e) character())
  pcol  <- intersect(c("P_value", "Pval", "P.Value", "pvalue"), cn)
  apcol <- intersect(c("FDR", "adj.P.Val", "padj"), cn)

  ked <- tryCatch(as.data.frame(ke), error = function(e) NULL)
  ksea.lines <- NULL
  if (!is.null(ked) && "KSEA_Z" %in% cn) {
    okz <- is.finite(suppressWarnings(as.numeric(ked$KSEA_Z)))
    if (any(okz)) {
      src <- if ("KSEA_Input" %in% cn) as.character(ked$KSEA_Input[which(okz)[1]]) else "site-level"
      ksea.lines <- c(
        .pa_line("Signed KSEA z-score", "Casado et al. 2013 formula; two-sided normal p, BH-adjusted"),
        .pa_line("KSEA input fold changes", paste0(src, " log2FC")),
        .pa_line("Kinases with KSEA z-score", sum(okz)))
    }
  }

  c("== Kinase Enrichment ==",
    .pa_line("Number of kinases tested", .pa_nrow(ke)),
    .pa_line("Raw p-value column", if (length(pcol)) pcol[1] else "P_value"),
    .pa_line("Adjusted p-value column", if (length(apcol)) apcol[1] else "FDR"),
    .pa_line("Multiple-testing correction", "Benjamini-Hochberg (FDR)"),
    .pa_line("Background set", "all quantified phosphosites"),
    ksea.lines)
}

## ---- Biomarker / ROC ----------------------------------------------------
.pa_biomarker_section <- function(analSet) {
  roc   <- .pa_g(analSet, "ROCtest", NULL)
  multi <- .pa_g(analSet, "multiROC", NULL)
  frank <- .pa_g(analSet, "feat.rank.mat", NULL)
  if (is.null(roc) && is.null(multi) && is.null(frank)) return(NULL)

  c("== Biomarker / ROC Analysis ==",
    .pa_line("Model built", if (!is.null(roc) || !is.null(multi)) "yes" else "no"),
    .pa_line("Ranked candidate features", .pa_nrow(frank)),
    .pa_line("Performance", "cross-validated AUC (see roc/biomarker CSVs)"),
    .pa_line("Covariate/batch adjustment", "re-fit within each CV fold (leakage-free)"))
}

## ---- Network ------------------------------------------------------------
.pa_network_section <- function(paramSet, analSet) {
  comps <- .pa_g(analSet, "ppi.comps", NULL)
  stats <- .pa_g(analSet, "net.stats", NULL)
  db    <- .pa_g(paramSet, "ppi.db.name", NULL)
  if (is.null(comps) && is.null(stats) && is.null(db)) return(NULL)

  n.sub <- tryCatch(if (is.list(comps)) length(comps) else .pa_nrow(comps),
                    error = function(e) "NA")

  lines <- "== Network Analysis =="
  if (!is.null(db)) {
    req.exp <- tryCatch(isTRUE(as.logical(.pa_g(paramSet, "ppi.require.exp", FALSE))),
                        error = function(e) FALSE)
    lines <- c(lines,
      .pa_line("PPI database", .pa_ppidb_label(db)),
      .pa_line("Interaction confidence cutoff", .pa_g(paramSet, "ppi.conf")),
      .pa_line("Experimental evidence required", if (req.exp) "yes" else "no"))
  }
  c(lines, .pa_line("Number of subnetworks", n.sub))
}

## ---- Co-expression network (WGCNA / CEMiTool) --------------------------
.pa_coexp_section <- function(paramSet) {
  cm <- .pa_g(paramSet, "coexp.cor.method", NULL)
  nm <- .pa_g(paramSet, "coexp.n.modules", NULL)
  if (is.null(cm) && is.null(nm)) return(NULL)

  power <- .pa_g(paramSet, "coexp.soft.power")
  c("== Co-expression Network (CEMiTool / WGCNA) ==",
    .pa_line("Correlation method", cm),
    .pa_line("Soft-thresholding power (beta)", power),
    .pa_line("Minimum module size (genes)", .pa_g(paramSet, "coexp.min.ngen")),
    .pa_line("Number of modules", nm))
}

## ---- Meta-analysis ------------------------------------------------------
.pa_meta_section <- function(paramSet, analSet) {
  mm <- .pa_g(analSet, "meta.mat", .pa_g(analSet, "meta.mat.all", NULL))
  if (is.null(mm)) return(NULL)
  c("== Meta-analysis ==",
    .pa_line("Method", .pa_g(paramSet, "meta.method", .pa_g(paramSet, "selectedMetaMethod"))),
    .pa_line("Number of features in meta result", .pa_nrow(mm)),
    .pa_line("Multiple-testing correction", "Benjamini-Hochberg (FDR)"))
}

## ---- Top-level entry point ---------------------------------------------
#' Write the consolidated analysis_summary.txt for the Downloads page.
#'
#' Reconstructs the dataset block and one section per analysis that has been
#' run from persisted state, and writes them to `out.file` in the working
#' (session) directory. Best-effort: on any error it still writes a short
#' placeholder rather than propagating the failure to the caller.
#' @param out.file output file name (relative to the session working dir).
#' @return the output file name, invisibly.
WriteAnalysisSummary <- function(out.file = "analysis_summary.txt") {
  res <- tryCatch({
    paramSet <- tryCatch(readSet(NULL, "paramSet"), error = function(e) NULL)
    msgSet   <- tryCatch(readSet(NULL, "msgSet"),   error = function(e) NULL)
    analSet  <- tryCatch(readSet(NULL, "analSet"),  error = function(e) NULL)

    # Nothing loaded yet: don't leave an all-NA file behind.
    if (is.null(paramSet) && is.null(msgSet)) return(invisible(out.file))

    dataName <- as.character(.pa_g(paramSet, "dataName", ""))
    dataSet  <- if (nzchar(dataName) && !identical(dataName, "NA"))
                  tryCatch(readDataset(dataName), error = function(e) NULL) else NULL

    add <- function(acc, sec) if (is.null(sec)) acc else c(acc, sec, "")
    body <- NULL
    body <- add(body, .pa_dataset_section(paramSet, msgSet, dataSet))
    body <- add(body, .pa_filternorm_section(paramSet, msgSet))
    body <- add(body, .pa_summarization_section(paramSet))
    body <- add(body, .pa_de_section(paramSet))
    body <- add(body, .pa_ptm_section(paramSet))
    body <- add(body, .pa_enrich_section(paramSet, analSet))
    body <- add(body, .pa_gsea_section(paramSet))
    body <- add(body, .pa_kinase_section(analSet))
    body <- add(body, .pa_biomarker_section(analSet))
    body <- add(body, .pa_network_section(paramSet, analSet))
    body <- add(body, .pa_coexp_section(paramSet))
    body <- add(body, .pa_meta_section(paramSet, analSet))

    header <- c("ProteoAnalyst - Analysis Summary",
                paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                strrep("=", 60), "")
    footer <- c(strrep("-", 60),
                "Notes:",
                " - Adjusted p-values use the Benjamini-Hochberg (FDR) procedure.",
                " - 'Background / universe' is the feature set against which enrichment",
                "   is tested; the default is the measured proteome/phosphoproteome.",
                " - Thresholds shown are the values applied in the current session.")

    writeLines(c(header, body, footer), con = out.file, useBytes = TRUE)
    out.file
  }, error = function(e) {
    tryCatch(writeLines(paste("Analysis summary unavailable:", conditionMessage(e)),
                        out.file), error = function(e2) NULL)
    out.file
  })
  invisible(res)
}
