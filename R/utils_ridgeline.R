##################################################
## R script for ProteoAnalyst
## Description: Compute Ridgeline plot
## Authors: 
## Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
## Jessica Ewald, jessica.ewald@mail.mcgill.ca
##################################################

compute.ridgeline <- function(dataSet, imgNm = "abc", dpi=96, format="png", fun.type = "kegg", ridgeType = "ora", ridgeColor = "teal",rankOpt="fc", sigLevel = 0.05, pwNum=20, inx = 1){
  
  #save.image("ridge.RData");
  paramSet <- readSet(paramSet, "paramSet");
  msgSet <- readSet(msgSet, "msgSet");
  analSet <- readSet(analSet, "analSet");
  imageName <- paste0(imgNm, "dpi" , dpi, ".", format);
  jsonNm <- paste0(imgNm, ".json");
  anal.type <- paramSet$anal.type;
  require("dplyr");
  require("fgsea");
  require("reshape");
  require("ggplot2");
  require("ggridges");
  safe_qread <- function(path) {
    if (!file.exists(path)) {
      AddErrMsg(paste0("Required file missing: ", path)); return(0)
    }
    tryCatch(
      ov_qs_read(path),
      error = function(e) {
        AddErrMsg(paste0("Failed to read ", path)); return(0)
      }
    )
  }
  # process colors
  if(ridgeColor == "teal"){
    high.col = "#3C7C60";
    low.col = "#C9E1D6";
  } else if (ridgeColor == "orange"){
    high.col = "#ffa34d";
    low.col = "#994a00";    
  } else {
    high.col = "#00da00";
    low.col = "#005000";   
  }
  
  #get pw library
  setres <- .loadEnrichLib(fun.type, paramSet);
  current.featureset <- setres$current.featureset;
  # get DEGs
  if(anal.type == "proteinlist"){
    if(paramSet$numOfLists > 1){
      dataSet <- readDataset(paramSet$selDataNm);
    }
    sigmat <- as.data.frame(dataSet$prot.mat)

    sigmat$entrez <- rownames(sigmat);
    universe <- unique(unlist(current.featureset));
    expr.vec <- sigmat[,1];

    if(sum(expr.vec) == 0){
      msgSet$current.msg <- "Uploaded gene list needs to contain fold-change value to perform Ridgeline analysis!";
      saveSet(msgSet, "msgSet");
      return(-2);
    }
  }else if(anal.type == "onedata"){
    if(ridgeType == "ora"){
      sigmat <- dataSet$sig.mat;
      # Check if there are any significant features
      if(is.null(sigmat) || nrow(sigmat) == 0){
        msgSet <- readSet(msgSet, "msgSet");
        msgSet$current.msg <- "No significant features found for enrichment analysis. Please relax the differential expression threshold (e.g., increase p-value cutoff or decrease fold-change threshold) to obtain more features for pathway enrichment.";
        saveSet(msgSet, "msgSet");
        return(0);
      }
    }else{
      sigmat <- dataSet$comp.res;
      allmat <- dataSet$comp.res;
    }
    sigmat$entrez <- rownames(sigmat);
    universe <- rownames(dataSet$data.norm);
  }else{
    meta.avgFC <- analSet$meta.avgFC;
    inx <- 1;
    if(paramSet$selDataNm == "meta_default"){
        if(ridgeType == "ora"){
          sigmat <- analSet$meta.mat;
        # Check if there are any significant features
        if(is.null(sigmat) || nrow(sigmat) == 0){
          msgSet <- readSet(msgSet, "msgSet");
          msgSet$current.msg <- "No significant features found for enrichment analysis. Please relax the differential expression threshold (e.g., increase p-value cutoff or decrease fold-change threshold) to obtain more features for pathway enrichment.";
          saveSet(msgSet, "msgSet");
          return(0);
        }
        allmat <- safe_qread("meta.resTable.qs");
        sigmat <- cbind(unname(meta.avgFC[rownames(sigmat)]), sigmat);

      }else{
        allmat <- safe_qread("meta.resTable.qs");
        sigmat <- allmat;
        sigmat <- cbind(unname(meta.avgFC[rownames(sigmat)]), sigmat);

      }
      allmat$logFC <- unname(meta.avgFC[rownames(allmat)]);
      universe <- rownames(allmat);
    }else{
      dataSet <- readDataset(paramSet$selDataNm);
      if(ridgeType == "ora"){
        sigmat <- dataSet$sig.mat;
        # Check if there are any significant features
        if(is.null(sigmat) || nrow(sigmat) == 0){
          msgSet <- readSet(msgSet, "msgSet");
          msgSet$current.msg <- "No significant features found for enrichment analysis. Please relax the differential expression threshold (e.g., increase p-value cutoff or decrease fold-change threshold) to obtain more features for pathway enrichment.";
          saveSet(msgSet, "msgSet");
          return(0);
        }
      }else{
        sigmat <- dataSet$comp.res;
        allmat <- dataSet$comp.res;
      }
      sigmat <- as.data.frame(sigmat);
      sigmat$entrez <- rownames(sigmat);
      universe <- rownames(dataSet$data.norm);
    }
  }
  if(ridgeType == "ora"){
    gene.vec <- rownames(sigmat);

    # Check if phospho data and use appropriate functions
    is_phospho <- (!is.null(paramSet$data.type) && paramSet$data.type == "phospho")


    # Check if there are any significant genes
    if(length(gene.vec) == 0){
      msgSet <- readSet(msgSet, "msgSet");
      msgSet$current.msg <- "No significant features found for enrichment analysis. Please relax the differential expression threshold (e.g., increase p-value cutoff or decrease fold-change threshold) to obtain more features for pathway enrichment.";
      saveSet(msgSet, "msgSet");
      return(0);
    }

    # When the caller (the multi-library enrichment driver) already computed the
    # ORA for this fun.type over the SAME query, reuse the enr.mat.qs /
    # hits_query.qs it wrote instead of recomputing — the recompute here was
    # ~half the multi-library step's runtime. Gated on the ov.enrich.reuse option
    # (set only by .ai_run_enrich_libs); default-off keeps manual / matrix /
    # refine callers recomputing, where the ridgeline query (sig.mat) differs.
    .reuse.enr <- isTRUE(getOption("ov.enrich.reuse", FALSE)) &&
                  file.exists("enr.mat.qs") && file.exists("hits_query.qs");
    if (.reuse.enr) {
      message("[compute.ridgeline] reusing precomputed enr.mat.qs for ", fun.type,
              " (skipped ORA recompute)");
    } else if (is_phospho) {
      # For phospho data, gene.vec contains phosphosite IDs (uniprot+site)
      .performEnrichAnalysisPhospho(dataSet, imgNm, fun.type, gene.vec, "ridgeline")
    } else {
      # Regular data - use standard enrichment
      sym.vec <- doEntrez2SymbolMapping(gene.vec, paramSet$data.org, paramSet$data.idType);
      names(gene.vec) <- sym.vec;
      .performEnrichAnalysis(dataSet, imgNm, fun.type, gene.vec, "ridgeline")
    }


    # Check if file exists before reading
    if (!file.exists("enr.mat.qs")) {
      AddErrMsg("Enrichment analysis did not produce results. Please try different parameters."); return(0)
    }

    res <- safe_qread("enr.mat.qs");
    colnames(res) <- c("size", "expected", "overlap", "pval", "padj");

    res <- res[,c(4,5,3,1,2)]
    res <- as.data.frame(cbind(pathway=rownames(res), res));
    res$padj <- as.numeric(res$padj)
    res$pval <- as.numeric(res$pval)

    # Load hits.query saved by enrichment analysis (contains UniProt IDs after conversion)
    if (!file.exists("hits_query.qs")) {
      AddErrMsg("Enrichment results file not found. Please re-run enrichment analysis."); return(0)
    }
    saved.hits.query <- safe_qread("hits_query.qs");
  } else {
    # GSEA branch

    rankedVec<- ComputeRankedVec(dataSet, rankOpt, paramSet$selectedFactorInx);

    # For proteomics data, names(rankedVec) are UniProt IDs, need to convert to Entrez
    # Universe for proteinlist is Entrez IDs from library, for onedata is UniProt IDs from data
    if(anal.type == "proteinlist"){
      # Universe already contains Entrez IDs from current.featureset
      gene.vec <- universe;
    } else {
      # Universe contains UniProt IDs, need to convert to Entrez

      # Normalize UniProt IDs (remove phosphosite annotations, isoforms)
      normalized.universe <- sub("_[A-Z]_\\d+$", "", universe)  # Remove phosphosites
      normalized.universe <- sub("-\\d+$", "", normalized.universe)  # Remove isoforms

      uniprot.map <- queryGeneDB("entrez_uniprot", paramSet$data.org)
      hit.inx <- match(normalized.universe, uniprot.map[, "accession"])
      entrez.vec <- uniprot.map[hit.inx, "gene_id"]
      gene.vec <- entrez.vec[!is.na(entrez.vec)]
    }

    # Convert Entrez IDs to symbols for fgsea matching
    sym.vec <- doEntrez2SymbolMapping(gene.vec, paramSet$data.org, "entrez");
    gene.nms <- sym.vec;

    # Create symbol-based gene sets from Entrez-based current.featureset
    current.featureset.symb <- lapply(current.featureset,
                         function(x) {
                           gene.nms[gene.vec%in%unlist(x)];
    }
    );

    # Convert rankedVec names from UniProt to symbols (via Entrez)

    # Normalize UniProt IDs (remove phosphosite annotations, isoforms)
    normalized.names <- sub("_[A-Z]_\\d+$", "", names(rankedVec))  # Remove phosphosites
    normalized.names <- sub("-\\d+$", "", normalized.names)        # Remove isoforms

    uniprot.map <- queryGeneDB("entrez_uniprot", paramSet$data.org)
    hit.inx <- match(normalized.names, uniprot.map[, "accession"])
    entrez.vec <- uniprot.map[hit.inx, "gene_id"]
    sym.vec <- doEntrez2SymbolMapping(entrez.vec, paramSet$data.org, "entrez")
    names(rankedVec) <- sym.vec
    # Remove NAs
    rankedVec <- rankedVec[!is.na(names(rankedVec)) & nzchar(names(rankedVec))]
    # Deduplicate: multiple isoforms/UniProt IDs can map to the same symbol;
    # keep the entry with the largest absolute value (strongest signal)
    if (any(duplicated(names(rankedVec)))) {
      rankedVec <- rankedVec[order(-abs(rankedVec))]
      rankedVec <- rankedVec[!duplicated(names(rankedVec))]
    }


    if(fun.type %in% c("go_bp", "go_mf", "go_cc")){
      res <- fgsea::fgsea(pathways = current.featureset.symb, 
                          stats    = rankedVec,
                          minSize  = 5,
                          maxSize = 500,
                          scoreType = "std",
                          nperm=10000)    
    }else{
      res <- fgsea::fgsea(pathways = current.featureset.symb, 
                          stats    = rankedVec,
                          minSize  = 5,
                          maxSize = 500,
                          scoreType = "std")   
      
    }


  }
  res <- .signif_df(res, 4);
  res <- res[order(res$pval),];
  resTable <- res;
  # process results;
  res <- res[,c(1,2,3)];
  colnames(res) <- c("name", "pval", "adj.pval");
  res.sig <- res;
  totalSigPws <- dim(res[res$pval < sigLevel, ])[1];
  
  if(pwNum != -1){
    if(dim(res.sig)[1] > pwNum){ # limit size if too many sig results
      res.sig <- res.sig[1:pwNum, ]
    }
  }

  ridge.pathway.compartments <- .getRidgePathwayCompartments(
    pathways = as.character(resTable$pathway),
    fun.type = fun.type,
    current.featureset = current.featureset,
    current.setids = setres$current.setids,
    paramSet = paramSet
  )
  # Get hits per pathway early (needed for gs.plot creation)
  # For ORA, use the saved hits.query from enrichment analysis (contains UniProt IDs)
  # For GSEA, compute from current.featureset and convert to UniProt
  if(ridgeType == "ora"){
    # Use saved hits.query (already has UniProt IDs after conversion)
    hits.query <- saved.hits.query[resTable$pathway];
  } else {
    # GSEA: compute hits.query from current.featureset (Entrez IDs)
    # Use gene.vec (Entrez) as the universe for matching
    hits.query <- lapply(current.featureset,
                         function(x) {
                           x[x %in% gene.vec];
                         }
    );
    hits.query <- hits.query[resTable$pathway];

    # Convert Entrez IDs back to UniProt IDs for consistency
    uniprot.map <- queryGeneDB("entrez_uniprot", paramSet$data.org)

    hits.query <- lapply(hits.query, function(entrez.list) {
      # Find UniProt IDs for these Entrez IDs
      hit.inx <- match(entrez.list, uniprot.map[, "gene_id"])
      uniprot.list <- uniprot.map[hit.inx, "accession"]
      uniprot.list[!is.na(uniprot.list)]
    })

  }

  # ── Full-distribution ORA ridgeline ──────────────────────────────────────
  # Originally the ORA ridge plotted only the SIGNIFICANT proteins overlapping
  # each pathway (dataSet$sig.mat). When few proteins pass the DE cut that
  # leaves 0-1 points per pathway and the ridges render empty. Mirror the GSEA
  # ridge instead: plot every pathway member that has a fold change (all
  # measured proteins), and flag the significant ones so the renderer can
  # highlight them while the rest are drawn grey. Pathway selection and the
  # p-value fill stay ORA-derived. Scoped to onedata (sig.mat / comp.res share
  # the data's UniProt id space); other anal.types keep the prior behaviour.
  ridge.sig.ids <- NULL
  if (ridgeType == "ora" && anal.type == "onedata" && !is.null(dataSet$comp.res)) {
    ridge.sig.ids <- rownames(sigmat)                 # significant proteins (sig.mat)
    # Fold changes come from the primary comparison of interest (a single-contrast
    # table); for multi-group designs comp.res is the overall F-test with no single
    # logFC, so fall back through the resolver. col 1 = logFC.
    primary.fc <- if (exists(".ov_primary_comp")) .ov_primary_comp(dataSet) else dataSet$comp.res
    full.fc  <- as.data.frame(primary.fc)
    univ.ids <- rownames(full.fc)
    umap <- tryCatch(queryGeneDB("entrez_uniprot", paramSet$data.org),
                     error = function(e) NULL)
    if (!is.null(umap)) {
      norm.ids <- sub("-\\d+$", "", sub("_[A-Z]_\\d+$", "", univ.ids))
      univ.ent <- umap[match(norm.ids, umap[, "accession"]), "gene_id"]
      full.hits <- lapply(current.featureset,
                          function(es) univ.ids[as.character(univ.ent) %in% as.character(es)])
      full.hits <- full.hits[resTable$pathway]
      if (any(vapply(full.hits, length, integer(1)) > 0)) {
        hits.query <- full.hits                       # all pathway members in data
        sigmat <- full.fc                             # fold-change source = all proteins
      }
    }
  }

  # prepare data for plotting
  # For phospho data, rownames(sigmat) are phosphosite IDs that need to be converted to entrez
  is_phospho <- (!is.null(paramSet$data.type) && paramSet$data.type == "phospho")

  if (is_phospho && !is.null(paramSet$phospho.mapping)) {
    # Convert phosphosite IDs to entrez IDs using the mapping
    phospho.to.entrez <- paramSet$phospho.mapping$phospho.to.entrez
    phosphosite.ids <- rownames(sigmat)

    # Map phosphosites to entrez (keep only mapped ones)
    entrez.ids <- phospho.to.entrez[phosphosite.ids]
    entrez.ids <- entrez.ids[!is.na(entrez.ids)]

    # Get fold changes for mapped phosphosites
    phosphosite.with.entrez <- names(entrez.ids)
    fc.values <- sigmat[phosphosite.with.entrez, inx]

    # Create data frame with entrez IDs and fold changes
    degs.plot <- data.frame(entrez = as.character(entrez.ids),
                           log2FC = fc.values,
                           stringsAsFactors = FALSE)
    degs.plot <- reshape::melt(degs.plot);
    colnames(degs.plot)[1] <- "entrez";
  } else {
    # Regular data - use rownames as-is
    degs.plot <- data.frame(entrez = rownames(sigmat), log2FC = sigmat[,inx]);
    degs.plot <- reshape::melt(degs.plot);
    colnames(degs.plot)[1] <- "entrez";
  }


  # Build a separate Entrez-based merge key. The saved hits_query.qs is a
  # display list of gene symbols, so it must not be used as the merge key.
  degs.plot$merge_id <- as.character(degs.plot$entrez);
  sample_ids <- degs.plot$merge_id[!is.na(degs.plot$merge_id) & nzchar(degs.plot$merge_id)];
  sample_ids <- head(sample_ids, 20);
  degs_are_entrez <- length(sample_ids) > 0 && mean(grepl("^[0-9]+$", sample_ids)) >= 0.8;

  if (!is_phospho && !degs_are_entrez) {
    normalized.ids <- sub("_[A-Z]_\\d+$", "", degs.plot$merge_id);
    normalized.ids <- sub("-\\d+$", "", normalized.ids);
    uniprot.map <- queryGeneDB("entrez_uniprot", paramSet$data.org);
    hit.inx <- match(normalized.ids, uniprot.map[, "accession"]);
    degs.plot$merge_id <- as.character(uniprot.map[hit.inx, "gene_id"]);
    degs.plot <- degs.plot[!is.na(degs.plot$merge_id), ];
    msg("[Ridgeline] Converted degs.plot IDs to Entrez merge IDs")
  } else if (!is_phospho) {
    msg("[Ridgeline] degs.plot IDs look like Entrez; using them directly for merge")
  }

  if (ridgeType == "ora" && !is_phospho) {
    hits.merge <- lapply(current.featureset, function(x) {
      ids <- as.character(x);
      ids[ids %in% degs.plot$merge_id];
    });
    hits.merge <- hits.merge[resTable$pathway];
    gs.plot <- reshape::melt(hits.merge);
    colnames(gs.plot) <- c("merge_id", "name");
  } else {
    gs.plot <- reshape::melt(hits.query);
    colnames(gs.plot) <- c("merge_id", "name");
  }

  if (!is_phospho) {
    gs.sample.ids <- as.character(gs.plot$merge_id);
    gs.sample.ids <- gs.sample.ids[!is.na(gs.sample.ids) & nzchar(gs.sample.ids)];
    gs.sample.ids <- head(gs.sample.ids, 20);
    gs.are.entrez <- length(gs.sample.ids) > 0 && mean(grepl("^[0-9]+$", gs.sample.ids)) >= 0.8;
    if (!gs.are.entrez) {
      normalized.gs.ids <- sub("_[A-Z]_\\d+$", "", as.character(gs.plot$merge_id));
      normalized.gs.ids <- sub("-\\d+$", "", normalized.gs.ids);
      uniprot.map <- queryGeneDB("entrez_uniprot", paramSet$data.org);
      hit.inx <- match(normalized.gs.ids, uniprot.map[, "accession"]);
      gs.plot$merge_id <- as.character(uniprot.map[hit.inx, "gene_id"]);
      gs.plot <- gs.plot[!is.na(gs.plot$merge_id), ];
    }
  }

  # For phospho data in ORA mode, gs.plot$merge_id contains phosphosite IDs
  # Need to convert to Entrez IDs to match degs.plot
  if (is_phospho && ridgeType == "ora") {

    # Strip phosphosite IDs to UniProt IDs and remove isoforms
    uniprot_ids <- sapply(strsplit(as.character(gs.plot$merge_id), "_"), function(x) x[1])
    uniprot_ids <- sub("-\\d+$", "", uniprot_ids)

    # Convert UniProt to Entrez
    uniprot.map <- queryGeneDB("entrez_uniprot", paramSet$data.org)
    hit.inx <- match(uniprot_ids, uniprot.map[, "accession"])
    entrez.ids <- uniprot.map[hit.inx, "gene_id"]

    # Update gs.plot with Entrez IDs
    gs.plot$merge_id <- as.character(entrez.ids)

    # Remove rows with NA Entrez IDs
    gs.plot <- gs.plot[!is.na(gs.plot$merge_id), ]

  }

  # For phospho data in GSEA mode, gs.plot$merge_id has UniProt IDs (converted from Entrez
  # at lines 283-288), but degs.plot$merge_id has Entrez IDs (from phospho.mapping) or
  # phosphosite IDs (else branch). Align them.
  if (is_phospho && ridgeType == "gsea") {
    if (!is.null(paramSet$phospho.mapping)) {
      # degs.plot$merge_id = Entrez IDs; convert gs.plot UniProt → Entrez to match
      uniprot.map <- queryGeneDB("entrez_uniprot", paramSet$data.org)
      gs.uniprot <- sub("-\\d+$", "", as.character(gs.plot$merge_id))
      hit.inx <- match(gs.uniprot, uniprot.map[, "accession"])
      gs.plot$merge_id <- as.character(uniprot.map[hit.inx, "gene_id"])
      gs.plot <- gs.plot[!is.na(gs.plot$merge_id), ]
    } else {
      # degs.plot$merge_id = phosphosite IDs; normalize to plain UniProt IDs to match gs.plot
      degs.plot$merge_id <- sub("_[A-Z]_\\d+$", "", degs.plot$merge_id)
      degs.plot$merge_id <- sub("-\\d+$", "", degs.plot$merge_id)
      gs.plot$merge_id <- sub("-\\d+$", "", as.character(gs.plot$merge_id))
    }
  }

  df <- merge(res.sig, gs.plot, by = "name", all.x = TRUE, all.y = FALSE);
  df <- merge(df, degs.plot, by = "merge_id", all.x = TRUE, all.y = FALSE);
  df <- na.omit(df)

  # Flag the significant proteins among the plotted pathway members so the
  # renderer highlights them in red (others grey). For GSEA (and any path that
  # didn't set ridge.sig.ids) fall back to the DE-significant set so both ORA and
  # GSEA mark the same proteins.
  if (is.null(ridge.sig.ids) && !is.null(dataSet$sig.mat) && nrow(dataSet$sig.mat) > 0)
    ridge.sig.ids <- rownames(dataSet$sig.mat)
  sig.flag <- if (!is.null(ridge.sig.ids)) df$entrez %in% ridge.sig.ids else rep(FALSE, nrow(df))
  df$sig <- factor(ifelse(sig.flag, "significant", "other"),
                   levels = c("significant", "other"))

  # Check if we have any data to plot
  if (nrow(df) == 0) {
    msgSet <- readSet(msgSet, "msgSet");
    msgSet$current.msg <- "No data available for ridgeline plot. Enrichment results may not have overlapping genes with expression data.";
    saveSet(msgSet, "msgSet");
    if(file.exists(jsonNm)) {
      unlink(jsonNm);
    }
    return(0);
  }

  # calculate the mean fold change to order the pathways in the plot
  means <- aggregate(df$value, by = list(df$name), mean);
  means <- means[order(means$x, decreasing = FALSE), ];
  df$name <- factor(df$name, levels = means$Group.1);
  
  # The density ridge is built from ALL pathway members, drawn LOW and semi-
  # transparent so it sits behind the member fold-changes (the "hits") instead of
  # covering them — significant members in red, the rest grey. A low scale also
  # stops a dense pathway's ridge from towering over and hiding the sparse rows
  # next to it. fill = adj. pval is each method's own p (ORA vs GSEA stay distinct).
  # The hits are a SEPARATE point layer (always rendered) over a low, semi-
  # transparent density. geom_density_ridges drops a group entirely — density AND
  # its jittered points — when it can't estimate a density (pathways with only a
  # couple of members), which is what left sparse rows blank; a standalone
  # geom_point keeps every member's fold-change visible regardless. Significant
  # members in red, the rest grey; fill = each method's own adj. pval.
  rp <- ggplot(df, aes(x = value, y = name)) +
    geom_density_ridges(
      aes(fill = adj.pval),
      alpha = 0.45, color = "white", scale = 0.85, rel_min_height = 0.005, size = 0.25) +
    geom_point(aes(colour = sig), shape = 124, size = 3.8, alpha = 1,
               position = position_jitter(height = 0.06, width = 0)) +
    geom_vline(xintercept = 0, color = "red") +
    scale_y_discrete(expand = expansion(add = c(0.3, 0.9)), name = "Gene Set",
                     labels = function(x) ifelse(nchar(x) > 45L,
                                                 paste0(substr(x, 1L, 42L), "..."), x)) +
    scale_x_continuous(expand = c(0, 0), name = "log2FC") +
    scale_fill_gradient("adj. pval", low = high.col, high = low.col) +
    ggplot2::scale_colour_manual("Protein",
      values = c(significant = "#D7263D", other = "#9aa0a6"), drop = FALSE) +
    coord_cartesian(clip = "off") +
    theme_ridges(center = TRUE) +
    theme(legend.position = "right",
          text = element_text(size=12, color = "black"),
          axis.title = element_text(size=12, face = "bold"),
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(size=12,color = "black"))
  
  msg("[Ridge] printing plot")
  Cairo::Cairo(file=imageName, width=10, height=8, type=format, bg="white", dpi=dpi, unit="in");
  print(rp);
  dev.off();
  msg("[Ridge] plot done")

  ##interative ridge json data
  ridge_bw <- rp$layers[[1]]$computed_stat_params$bandwidth;
  msg("[Ridge] ridge_bw=", if(is.null(ridge_bw)) "NULL" else as.character(ridge_bw))

  # Map IDs to symbols. Prefer Entrez merge IDs where available, while keeping
  # the original UniProt/phosphosite fallback for legacy rows.
  msg("[Ridge] df cols=", paste(colnames(df), collapse=","), " nrow=", nrow(df))
  msg("[Ridge] df$merge_id sample=", paste(head(as.character(df$merge_id), 5), collapse=","))
  if(all(grepl("^[0-9]+$", as.character(df$merge_id)))) {
    msg("[Ridge] symbol mapping via Entrez")
    symb <- doEntrez2SymbolMapping(as.character(df$merge_id), paramSet$data.org, "entrez");
  } else {
    msg("[Ridge] symbol mapping via UniProt")
    symb <- convert.uniprot.to.symbols(df$entrez, paramSet$data.org);
  }
  symb[is.na(symb)] <- df$entrez[is.na(symb)];
  df$symbol <- symb;
  msg("[Ridge] symbols done; NA count=", sum(is.na(symb)))

  # Add compartment annotation per gene-pathway row
  msg("[Ridge] computing compartments")
  df$compartment <- .getRidgeCompartments(df$merge_id, paramSet)
  df$pathway_compartment <- unname(ridge.pathway.compartments$primary[as.character(df$name)])
  df$pathway_compartment[is.na(df$pathway_compartment) | !nzchar(df$pathway_compartment)] <- "Unknown"
  msg("[Ridge] compartments done")

  data.list <- list();
  gene.list <- list();
  pval.list <- list();
  col.list <- list();
  pathway.compartment.list <- list();
  pathway.compartment.tooltip.list <- list();

  msg("[Ridge] building per-pathway lists; n levels=", length(levels(df$name)))
  for(i in 1:length(levels(df$name))){
    nm <- as.character(levels(df$name)[i]);
    data.list[[ nm ]] <- as.vector(unlist(df[which(df$name == nm), "value"]));
    gene.list[[ nm ]] <- as.vector(unlist(df[which(df$name == nm), "symbol"]));
    pval.list[[ nm ]] <- unname(unlist(res[which(res$name == nm), "pval"]));
    pathway.compartment.list[[ nm ]] <- unname(ridge.pathway.compartments$primary[nm]);
    pathway.compartment.tooltip.list[[ nm ]] <- unname(ridge.pathway.compartments$tooltip[nm]);
  }
  msg("[Ridge] per-pathway lists done")

  minFc <- min(df$value);
  maxFc <- max(df$value);
  minPval <- min(df$pval);
  maxPval <- max(df$pval);
  msg("[Ridge] minFc=", minFc, " maxFc=", maxFc, " minPval=", minPval, " maxPval=", maxPval)

  # enr result for table display (hits.query already assigned earlier)
  msg("[Ridge] building enr.res; ridgeType=", ridgeType)
  msg("[Ridge] resTable class=", class(resTable)[1], " ncol=", ncol(resTable), " nrow=", nrow(resTable))
  msg("[Ridge] resTable colnames=", paste(colnames(resTable), collapse=","))
  fun.anot <- if(ridgeType == "ora" && exists("saved.hits.query")) saved.hits.query[resTable$pathway] else hits.query
  msg("[Ridge] fun.anot length=", length(fun.anot))
  if(ridgeType == "ora"){
    total <- resTable[,5]; if(length(total) ==1) { total <- matrix(total) };
    fun.pval <- resTable[,"pval"]; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval) };
    fun.padj <- resTable[,"padj"]; if(length(fun.padj) ==1) { fun.padj <- matrix(fun.padj) };
    if(ridgeType == "ora"){
      hit.num <- resTable[,4]; if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
    }else{
      hit.num <- resTable[,"size"]; if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
    }
  }else{
    total <- as.list(resTable[,5])[[1]]; if(length(total) ==1) { total <- matrix(total) };
    fun.pval <- as.list(resTable[,"pval"])[[1]]; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval) };
    fun.padj <- as.list(resTable[,"padj"])[[1]]; if(length(fun.padj) ==1) { fun.padj <- matrix(fun.padj) };
    if(ridgeType == "ora"){
      hit.num <- as.list(resTable[,4])[[1]]; if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
    }else{
      hit.num <- as.list(resTable[,"size"])[[1]]; if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
    }

  }
  msg("[Ridge] fun.pval class=", class(fun.pval), " length=", length(fun.pval))
  msg("[Ridge] fun.padj class=", class(fun.padj), " length=", length(fun.padj))
  msg("[Ridge] hit.num class=", class(hit.num), " length=", length(hit.num))

  # Get pathway IDs - handle NULL cases
  msg("[Ridge] building fun.ids")
  if(!is.null(setres$current.setids) && length(names(fun.anot)) > 0) {
    fun.ids <- as.vector(setres$current.setids[names(fun.anot)]);
    # Replace NAs with empty strings
    fun.ids[is.na(fun.ids)] <- ""
  } else {
    # Fallback: use pathway names as IDs
    fun.ids <- names(fun.anot)
  }
  if(length(fun.ids) ==1) { fun.ids <- matrix(fun.ids) };
  msg("[Ridge] fun.ids length=", length(fun.ids))

  msg("[Ridge] building enr.res list")
  enr.res <- list(
    fun.anot = fun.anot,
    fun.ids = fun.ids,
    fun.pval = fun.pval,
    fun.padj = fun.padj,
    hit.num = hit.num,
    total= total,
    pathway.compartments = unname(ridge.pathway.compartments$primary[as.character(resTable$pathway)]),
    pathway.compartment.method = unname(ridge.pathway.compartments$method[as.character(resTable$pathway)]),
    pathway.compartment.score = unname(ridge.pathway.compartments$score[as.character(resTable$pathway)])
  );
  msg("[Ridge] enr.res built")

  if(ridgeType == "gsea"){
    msg("[Ridge] adding ES to enr.res")
    enr.res[["ES"]] <- unname(unlist(resTable[,"ES"]));
    msg("[Ridge] ES added; length=", length(enr.res[["ES"]]))
  }

  msg("[Ridge] building res.list")
  msg("[Ridge] df subset cols check: name=", "name" %in% colnames(df),
      " entrez=", "entrez" %in% colnames(df),
      " merge_id=", "merge_id" %in% colnames(df),
      " value=", "value" %in% colnames(df),
      " pval=", "pval" %in% colnames(df),
      " symbol=", "symbol" %in% colnames(df),
      " compartment=", "compartment" %in% colnames(df),
      " pathway_compartment=", "pathway_compartment" %in% colnames(df))
  res.list <- list(data=data.list,
                   proteinlist=gene.list,
                   df=df[, c("name","entrez","merge_id","value","pval","symbol","compartment","pathway_compartment"), drop=FALSE],
                   minPval=minPval,
                   maxPval=maxPval,
                   min=minFc,
                   max=maxFc,
                   minPval = min(res.sig$pval),
                   maxPval = max(res.sig$pval),
                   bandwidth=ridge_bw,
                   pathwayPvals = pval.list,
                   pathwayCols = col.list,
                   pathwayCompartments = pathway.compartment.list,
                   pathwayCompartmentTooltips = pathway.compartment.tooltip.list,
                   enrRes = enr.res,
                   dat.opt = paramSet$selDataNm,
                   naviString="ridge");
  msg("[Ridge] res.list built")

  if(ridgeType == "gsea"){
    msg("[Ridge] writing GSEA CSV; resTable leadingEdge class=", class(resTable$leadingEdge))
    resTable_csv <- resTable[, setdiff(colnames(resTable), "leadingEdge"), drop=FALSE]
    csv.nm <- paste0(imgNm, ".csv");
    fast.write(resTable_csv, file=csv.nm);
    msg("[Ridge] CSV written")
  }

  msg("[Ridge] saving analSet")
  analSet$ridgeline <- res.list;
  saveSet(analSet, "analSet");
  msg("[Ridge] analSet saved")

  msg("[Ridge] converting to JSON")
  json.obj <- rjson::toJSON(res.list);
  msg("[Ridge] JSON done; writing to", jsonNm)
  sink(jsonNm);
  cat(json.obj);
  sink();
  msg("[Ridge] JSON written")

  
  #for link sharing
  paramSet$jsonNms$ridge <- jsonNm
  paramSet$partialToBeSaved <- c( paramSet$partialToBeSaved, c(jsonNm));
  saveSet(paramSet, "paramSet");
  
  imgSet <- readSet(imgSet, "imgSet");
  rownames(resTable) <- NULL;
  
  imgSet$compute.ridgeline <- imageName;
    saveSet(imgSet, "imgSet");
  
  # Return the number of displayed pathways so Java gets a positive value on success
  # even when no pathways reach the significance threshold (totalSigPws == 0).
  return(max(totalSigPws, nrow(res.sig)))
}

# Look up compartment for each gene in a ridgeline df.
# merge_id is either all-numeric Entrez IDs or UniProt/phosphosite IDs.
.getRidgeCompartments <- function(merge.ids, paramSet) {
  n <- length(merge.ids)
  comp.vec <- rep("Unknown", n)
  org <- paramSet$data.org
  if (is.null(org) || org == "omk" || org == "NA") return(comp.vec)

  lib.path <- if (exists("api.lib.path")) api.lib.path else paramSet$lib.path
  loc.path <- paste0(lib.path, org, "/", org, "_localization.qs")
  if (!file.exists(loc.path)) loc.path <- paste0(lib.path, org, "/localization.qs")
  if (!file.exists(loc.path)) return(comp.vec)

  loc.data <- try(ov_qs_read(loc.path), silent = TRUE)
  if (inherits(loc.data, "try-error") || is.null(loc.data)) return(comp.vec)
  if (!all(c("EntrezID", "Broad.category") %in% colnames(loc.data))) return(comp.vec)

  ids <- as.character(merge.ids)
  if (all(grepl("^[0-9]+$", ids))) {
    # Already Entrez IDs
    entrez.ids <- ids
  } else {
    # UniProt or phosphosite IDs — bridge via entrez_uniprot table
    uniprot.db <- try(queryGeneDB("entrez_uniprot", org), silent = TRUE)
    if (!inherits(uniprot.db, "try-error") && !is.null(uniprot.db) &&
        "accession" %in% colnames(uniprot.db) && "gene_id" %in% colnames(uniprot.db)) {
      up2entrez <- setNames(as.character(uniprot.db$gene_id), as.character(uniprot.db$accession))
      clean.ids <- sub("_.*$", "", ids)
      entrez.ids <- up2entrez[clean.ids]
    } else {
      return(comp.vec)
    }
  }

  loc.entrez <- as.character(loc.data$EntrezID)
  hits <- which(!is.na(entrez.ids) & entrez.ids %in% loc.entrez)
  for (i in hits) {
    loc.rows <- loc.data[loc.entrez == as.character(entrez.ids[i]), , drop = FALSE]
    broad <- paste(unique(as.character(loc.rows$Broad.category)), collapse = "; ")
    main <- if ("Main.location" %in% colnames(loc.rows)) paste(unique(as.character(loc.rows$Main.location)), collapse = "; ") else NA_character_
    comp.vec[i] <- .paPrimaryCompartment(broad, main)$primary
  }
  comp.vec
}

.getRidgePathwayCompartments <- function(pathways, fun.type, current.featureset, current.setids, paramSet) {
  pathways <- as.character(pathways)
  empty <- list(
    primary = setNames(rep("Unknown", length(pathways)), pathways),
    tooltip = setNames(rep("Compartment: Unknown", length(pathways)), pathways),
    method = setNames(rep("inferred", length(pathways)), pathways),
    score = setNames(rep(0, length(pathways)), pathways)
  )
  if (length(pathways) == 0 || is.null(current.featureset) || length(current.featureset) == 0) return(empty)

  lib.path <- if (exists("api.lib.path")) api.lib.path else paramSet$lib.path
  org <- paramSet$data.org
  kegg.ids <- rep("", length(pathways))
  if (!is.null(current.setids) && length(current.setids) > 0) {
    kegg.ids <- as.vector(current.setids[pathways])
    kegg.ids[is.na(kegg.ids)] <- ""
  }

  if (identical(fun.type, "kegg") && exists(".annotateKeggPathwayCompartments", mode = "function")) {
    ann <- try(.annotateKeggPathwayCompartments(pathways, kegg.ids, current.featureset, org, lib.path), silent = TRUE)
    if (!inherits(ann, "try-error") && !is.null(ann) && nrow(ann) == length(pathways)) {
      primary <- setNames(as.character(ann$Primary.Compartment), pathways)
      tooltip <- setNames(paste0(
        "Compartment: ", ann$Primary.Compartment,
        "; method: ", ann$Compartment.Method,
        "; score: ", ann$Compartment.Score,
        "; distribution: ", ann$Compartment.Distribution
      ), pathways)
      return(list(
        primary = primary,
        tooltip = tooltip,
        method = setNames(as.character(ann$Compartment.Method), pathways),
        score = setNames(as.numeric(ann$Compartment.Score), pathways)
      ))
    }
  }

  loc.data <- NULL
  if (exists(".loadPathwayLocalization", mode = "function")) {
    loc.data <- .loadPathwayLocalization(org, lib.path)
  }
  if (is.null(loc.data)) return(empty)

  rows <- lapply(pathways, function(pw) {
    genes <- if (pw %in% names(current.featureset)) current.featureset[[pw]] else character(0)
    ann <- if (exists(".inferKeggPathwayCompartment", mode = "function")) {
      .inferKeggPathwayCompartment(genes, loc.data)
    } else {
      list(primary = "Unknown", all = "Unknown", method = "inferred", score = 0, distribution = "Unknown")
    }
    ann
  })
  primary <- setNames(vapply(rows, `[[`, character(1), "primary"), pathways)
  method <- setNames(vapply(rows, `[[`, character(1), "method"), pathways)
  score <- setNames(vapply(rows, function(x) as.numeric(x$score), numeric(1)), pathways)
  tooltip <- setNames(vapply(rows, function(x) {
    paste0("Compartment: ", x$primary, "; method: ", x$method,
           "; score: ", round(as.numeric(x$score), 3),
           "; distribution: ", x$distribution)
  }, character(1)), pathways)
  list(primary = primary, tooltip = tooltip, method = method, score = score)
}
