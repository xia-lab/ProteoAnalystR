

#' Impute Missing Values in Variables
#'
#' This function imputes missing values in variables (columns) of the dataset using various methods.
#'
#' @param dataName The name of the dataset to be processed.
#' @param method The imputation method to be used. Can be one of "exclude", "min", "colmin", "mean", "median", "knn_var", "knn_smp", "bpca", "ppca", "svdImpute".
#'
#' @author Guangyan Zhou \email{guangyan.zhou@mail.mcgill.ca}
#' @details Additional details about the function, if needed.
#'
#' @examples
#' \dontrun{
#' ImputeMissingVar("data_filename", method = "min")
#' }
#'
#' @export
#' @license MIT License
#'
ImputeMissingVar <- function(dataName="", method="min", min.obs.per.group=2){
  dataSet <- readDataset(dataName);
  msgSet <- readSet(msgSet, "msgSet");
  paramSet <- readSet(paramSet, "paramSet");

  int.mat <- dataSet$data.norm;
  if (is.null(int.mat)) {
    msgSet$current.msg <- c(msgSet$current.msg, "No normalized data available for imputation.");
    saveSet(msgSet, "msgSet");
    return(0);
  }

  # Advisory: after MSstats TMP+MBimpute summarization, censored values were
  # already imputed inside the model; remaining NAs are proteins entirely
  # unobserved in a condition, which MSstats deliberately leaves missing.
  if (identical(paramSet$summ.method, "msstats") && !identical(method, "none")) {
    msgSet$current.msg <- c(msgSet$current.msg,
      paste0("Note: summarization used MSstats MBimpute, which already handled ",
             "censored missing values; '", method, "' imputation will additionally ",
             "fill proteins completely unobserved in a condition. 'none' is recommended."));
  }

  # Treat strict zeros as missing (common in proteomics exports)
  int.mat[int.mat == 0] <- NA
  current.msg <- msgSet$current.msg;

  # DETECTION FILTER (before imputation) -- proteomics only. Mirrors the phospho
  # path (ImputeMissingVarPhospho): a protein quantified in only one condition, or
  # nearly so, is otherwise carried into min/LoD imputation, which fabricates an
  # extreme fold change that can pass FDR spuriously. Requiring >= min.obs.per.group
  # genuine observations within EACH group level removes these on/off artifacts.
  # Proteomics missingness is MNAR-dominant, so min-based imputation is retained,
  # but must be paired with this per-group validity filter (cf. Liu & Dongre 2021,
  # Brief. Bioinform.; the standard "valid values per group" step in Perseus/MSstats).
  # Validated: iPRG-2015 empirical FDR 0.22 -> ~0 with no loss of sensitivity;
  # PXD005536 nested-interaction FDR 3 -> 0. Applied to prot/peptide only (array/
  # count data are not left-censored MNAR); disable with min.obs.per.group = 0.
  if (!is.null(dataSet$type) && dataSet$type %in% c("prot", "peptide") &&
      is.numeric(min.obs.per.group) && min.obs.per.group > 0) {
    grp <- tryCatch({
      mi <- dataSet$meta.info
      if (!is.null(mi) && NCOL(mi) >= 1) { g <- as.character(mi[[1]]); names(g) <- rownames(mi); g[colnames(int.mat)] } else NULL
    }, error = function(e) NULL)
    if (!is.null(grp) && !all(is.na(grp))) {
      lvls <- unique(grp[!is.na(grp)])
      obs.ok <- sapply(lvls, function(L) {
        cols <- which(grp == L)
        if (!length(cols)) rep(TRUE, nrow(int.mat)) else rowSums(!is.na(int.mat[, cols, drop = FALSE])) >= min.obs.per.group
      })
      keep <- if (is.matrix(obs.ok)) rowSums(obs.ok) == length(lvls) else obs.ok
      n.drop <- sum(!keep)
      if (n.drop > 0 && sum(keep) > 0) {
        int.mat <- int.mat[keep, , drop = FALSE]
        current.msg <- c(current.msg, paste0(
          "Detection filter removed ", n.drop, " protein(s) not observed in at least ",
          min.obs.per.group, " samples in every group (before imputation)."))
      }
    }
  }

  row.nms <- rownames(int.mat);
  new.mat <- NULL;

  # NOTE: MSstats imputation is no longer used in the new workflow
  # Reason: Data is now summarized (peptide -> protein) BEFORE imputation
  # MSstats MBimpute requires peptide-level data, which we no longer have at this stage
  # All imputation now uses legacy methods (KNN, PPCA, etc.) on protein-level data

  # Legacy/fallback methods
  if(method %in% c("none", "NA", "")){
    # Pass-through: keep missing values (e.g. after MSstats MBimpute
    # summarization, where remaining NAs are condition-absent proteins).
    new.mat <- int.mat;
    current.msg <- c(current.msg, "Protein-level imputation skipped ('none'); missing values retained.");
  }else if(method=="exclude"){
    # OPTIMIZED: Use rowSums instead of apply for 60-100x speedup
    good.inx<-rowSums(is.na(int.mat))==0
    tmp.mat<-int.mat[good.inx,, drop=FALSE];
    if (is.null(new.mat)) {
      new.mat <- tmp.mat
      current.msg <- c(current.msg ,"Variables with missing values were excluded.")
      row.nms<-row.nms[good.inx]
    }
  }else if(method=="min"){
    if (is.null(new.mat)) {
      new.mat<- suppressWarnings(ReplaceMissingByLoD(int.mat));
      current.msg <- c(current.msg, "Missing variables were replaced by LoDs (1/5 of the min positive value for each variable)");
    }
  }else if(method=="colmin"){
    tmp.mat<-apply(int.mat, 1, function(x){
      if(sum(is.na(x))>0){
        x[is.na(x)]<-min(x,na.rm=T)/2;
      }
      x;
    });
    tmp.mat = t(tmp.mat)
    if (is.null(new.mat)) {
      new.mat <- tmp.mat
      current.msg <- c(current.msg,"Missing variables were replaced by 1/2 of min values for each feature column.");
    }
  }else if (method=="mean"){
    tmp.mat<-apply(int.mat, 1, function(x){
      if(sum(is.na(x))>0){
        x[is.na(x)]<-mean(x,na.rm=T);
      }
      x;
    });
    tmp.mat = t(tmp.mat)
    if (is.null(new.mat)) {
      new.mat <- tmp.mat
      current.msg <- c(current.msg,"Missing variables were replaced with the mean value for each feature column.");
    }
  }else if (method == "median"){
    tmp.mat<-apply(int.mat, 1, function(x){
      if(sum(is.na(x))>0){
        x[is.na(x)]<-median(x,na.rm=T);
      }
      x;
    });
   tmp.mat = t(tmp.mat)
    if (is.null(new.mat)) {
      new.mat <- tmp.mat
      current.msg <- c(current.msg,"Missing variables were replaced with the median for each feature column.");
    }
  }else{
    if (is.null(new.mat)) {
      if(method == "knn_var"){
        # impute.knn keeps the input (features x samples) orientation. Do NOT
        # transpose the result, or it becomes samples x features and no longer
        # aligns with the feature row names assigned below (this previously threw
        # "invalid 'row.names' length" whenever #features != #samples).
        new.mat<-impute::impute.knn(as.matrix(int.mat))$data;
        current.msg <- c(current.msg, "Missing variables were imputed using KNN (feature-wise)");
      }else if(method == "knn_smp"){
        # Impute across samples (rows after transpose), then transpose back to the
        # original features x samples orientation.
        new.mat<-t(impute::impute.knn(data.matrix(t(int.mat)))$data);
        current.msg <- c(current.msg, "Missing variables were imputed using KNN (sample-wise)");
      }else if(method == "bpca"){
        new.mat<-pcaMethods::pca(t(int.mat), nPcs =5, method="bpca", center=T)@completeObs;
        new.mat = t(new.mat)
        current.msg <- c(current.msg, "Missing variables were imputed using BPCA");
      }else if(method == "ppca"){
        new.mat<-pcaMethods::pca(t(int.mat), nPcs =5, method="ppca", center=T)@completeObs;
        new.mat = t(new.mat)
        current.msg <- c(current.msg, "Missing variables were imputed using PPCA");
      }else if(method == "svdImpute"){
        new.mat<-pcaMethods::pca(t(int.mat), nPcs =5, method="svdImpute", center=T)@completeObs;
        new.mat = t(new.mat)
        current.msg <- c(current.msg, "Missing variables were imputed using SVD Impute");
      }else if(method %in% c("mindet", "minprob", "qrilc")){
        # MinDet: deterministic minimum value imputation (MNAR/left-censored)
        # MinProb: stochastic minimum value imputation (MNAR/left-censored)
        # QRILC: Quantile Regression Imputation of Left-Censored data (MNAR)
        if (requireNamespace("imputeLCMD", quietly = TRUE)) {
          if (method == "mindet") {
            new.mat <- t(imputeLCMD::impute.MinDet(t(int.mat)))
            current.msg <- c(current.msg, "Missing variables were imputed using MinDet (deterministic minimum)");
          } else if (method == "minprob") {
            new.mat <- t(imputeLCMD::impute.MinProb(t(int.mat)))
            current.msg <- c(current.msg, "Missing variables were imputed using MinProb (stochastic minimum)");
          } else if (method == "qrilc") {
            new.mat <- t(imputeLCMD::impute.QRILC(t(int.mat))[[1]])
            current.msg <- c(current.msg, "Missing variables were imputed using QRILC (quantile regression for left-censored data)");
          }
        } else {
          AddErrMsg("MNAR imputation methods require the 'imputeLCMD' package. Install with: BiocManager::install('imputeLCMD')"); return(0);
        }
      }else if(method == "seqknn"){
        # SeqKNN: sequential K-nearest neighbors (MAR/mixed)
        has.seqknn <- FALSE

        # SeqKNN imputes each feature from its k nearest feature rows (like
        # impute.knn), so the features x samples matrix is passed untransposed;
        # transposing starves the neighbor pool and errors on small datasets.
        k.seqknn <- min(10L, nrow(int.mat) - 1L)

        # Try multiUS package first (CRAN available)
        if (requireNamespace("multiUS", quietly = TRUE)) {
          tryCatch({
            new.mat <- multiUS::seqKNNimp(int.mat, k = k.seqknn)
            has.seqknn <- TRUE
            current.msg <- c(current.msg, "Missing variables were imputed using SeqKNN (sequential KNN via multiUS)");
          }, error = function(e) {
            #msg("[SeqKNN] multiUS::seqKNNimp failed: ", e$message)
          })
        }

        # Fallback to SeqKnn package if multiUS not available
        if (!has.seqknn && requireNamespace("SeqKnn", quietly = TRUE)) {
          tryCatch({
            new.mat <- SeqKnn::SeqKNN(int.mat, k = k.seqknn)
            has.seqknn <- TRUE
            current.msg <- c(current.msg, "Missing variables were imputed using SeqKNN (sequential KNN via SeqKnn)");
          }, error = function(e) {
            #msg("[SeqKNN] SeqKnn::SeqKNN failed: ", e$message)
          })
        }

        if (!has.seqknn) {
          AddErrMsg("SeqKNN requires either 'multiUS' (CRAN) or 'SeqKnn' (archived) package. Install with: install.packages('multiUS')"); return(0);
        }
      }else if(method == "impseq"){
        # Impseq: sequential covariance-based imputation
        if (requireNamespace("rrcovNA", quietly = TRUE)) {
          tryCatch({
            new.mat <- t(rrcovNA::impSeq(t(int.mat)))
            current.msg <- c(current.msg, "Missing variables were imputed using Impseq (sequential covariance-based)");
          }, error = function(e) {
            #msg("[Impseq] rrcovNA::impSeq failed: ", e$message)
            # Fallback to robust version
            if (requireNamespace("rrcovNA", quietly = TRUE)) {
              new.mat <- t(rrcovNA::impSeqRob(t(int.mat)))
              current.msg <- c(current.msg, "Missing variables were imputed using ImpseqRob (robust sequential covariance-based)");
            } else {
              AddErrMsg("Impseq requires the 'rrcovNA' package. Install with: install.packages('rrcovNA')"); return(0);
            }
          })
        } else {
          AddErrMsg("Impseq requires the 'rrcovNA' package. Install with: install.packages('rrcovNA')"); return(0);
        }
      }else{
        AddErrMsg(paste0("Unknown imputation method: ", method)); return(0);
      }
    }
  }
  msgSet$current.msg <- current.msg;
  saveSet(msgSet, "msgSet");

  # Record the imputation method alongside the other processing options, so the
  # report and slides can state it. Until now the chosen method survived only inside
  # msgSet$current.msg, which is free text, and both surfaces looked up an
  # impute.opt that nothing ever wrote.
  paramSet$impute.opt <- method;
  saveSet(paramSet, "paramSet");

  data.missed <- as.data.frame(new.mat);
  rownames(data.missed) <- row.nms;

  saveDataQs(data.missed, "data.missed.qs", paramSet$anal.type, dataName);
  dataSet$data.norm <- data.missed

  RegisterData(dataSet);
}

# Proteomics: MSstats-based imputation + summarization
.imputeWithMSstats <- function(msstats.path, meta.info) {
  msgSet <- readSet(msgSet, "msgSet");
  if (!file.exists(msstats.path)) {
    return(list(mat = NULL));
  }
  if (!requireNamespace("MSstats", quietly = TRUE)) {
    msgSet$current.msg <- c(msgSet$current.msg, "MSstats package not available for proteomics imputation.");
    saveSet(msgSet, "msgSet");
    return(list(mat = NULL));
  }
  msin <- try(ov_qs_read(msstats.path), silent = TRUE);
  if (inherits(msin, "try-error") || is.null(msin)) {
    msgSet$current.msg <- c(msgSet$current.msg, "Failed to read msstats_input.qs for imputation.");
    saveSet(msgSet, "msgSet");
    return(list(mat = NULL));
  }
  required.cols <- c("ProteinName", "Run", "BioReplicate", "Condition", "Intensity");
  if (!all(required.cols %in% colnames(msin))) {
    msgSet$current.msg <- c(msgSet$current.msg, "MSstats input missing required columns; skipping MSstats imputation.");
    saveSet(msgSet, "msgSet");
    return(list(mat = NULL));
  }
  proc <- try(MSstats::dataProcess(msin,
                                  normalization = "equalizeMedians",
                                  summaryMethod = "TMP",
                                  censoredInt = "NA",
                                  MBimpute = TRUE,
                                  featureSubset = "all"), silent = TRUE);
  if (inherits(proc, "try-error")) {
    msgSet$current.msg <- c(msgSet$current.msg, "MSstats dataProcess error during imputation.");
    saveSet(msgSet, "msgSet");
    return(list(mat = NULL));
  }
  prof <- proc$ProteinLevelData;
  value.col <- if ("LogIntensities" %in% colnames(prof)) "LogIntensities" else if ("ABUNDANCE" %in% colnames(prof)) "ABUNDANCE" else NULL;
  if (is.null(value.col)) {
    msgSet$current.msg <- c(msgSet$current.msg, "MSstats output missing LogIntensities/ABUNDANCE.");
    saveSet(msgSet, "msgSet");
    return(list(mat = NULL));
  }
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    msgSet$current.msg <- c(msgSet$current.msg, "reshape2 package not available to reshape MSstats output.");
    saveSet(msgSet, "msgSet");
    return(list(mat = NULL));
  }
  run.col <- if ("originalRUN" %in% colnames(prof)) "originalRUN" else "RUN";
  wide <- reshape2::dcast(prof, stats::as.formula(paste("Protein ~", run.col)), value.var = value.col);
  if (!"Protein" %in% colnames(wide)) {
    msgSet$current.msg <- c(msgSet$current.msg, "MSstats output missing Protein column.");
    saveSet(msgSet, "msgSet");
    return(list(mat = NULL));
  }
  mat <- as.matrix(wide[, -1, drop = FALSE]);
  rownames(mat) <- wide$Protein;
  storage.mode(mat) <- "numeric";
  # Align columns to metadata order if available
  if (!is.null(meta.info)) {
    run.order <- intersect(rownames(meta.info), colnames(mat));
    if (length(run.order) > 0) {
      mat <- mat[, run.order, drop = FALSE];
    }
  }
  return(list(mat = mat, proc = proc));
}

# Phosphoproteomics-specific imputation (operates on dataSet$data.norm)
# Group-aware hybrid imputation. `grp` is a per-column group label aligned to
# columns of `mat`. A sporadic gap with >=2 observations in its group is treated
# as MAR and sampled from that site's observed group distribution, following
# PhosR::scImpute (N(group mean, group SD)). Gaps with insufficient group
# information are treated as censored and filled with the scale-aware per-row
# LoD used by ReplaceMissingByLoD. A fixed local seed makes the stochastic fill
# reproducible without perturbing the caller's random-number stream.
.imputeMinGroupAware <- function(mat, grp, seed = 123L) {
  if (is.null(grp) || all(is.na(grp)) || length(grp) != ncol(mat))
    return(suppressWarnings(ReplaceMissingByLoD(mat)))

  mat <- as.matrix(mat)
  original.na <- is.na(mat)
  if (!any(original.na)) return(mat)

  # Obtain the floor from the canonical helper so log2 matrices use
  # min-log2(5), rather than the invalid min/5 operation on log intensities.
  lod.mat <- suppressWarnings(ReplaceMissingByLoD(mat))
  first.na <- max.col(original.na, ties.method = "first")
  row.floor <- lod.mat[cbind(seq_len(nrow(mat)), first.na)]

  seed <- suppressWarnings(as.integer(seed[1]))
  if (length(seed) == 1L && is.finite(seed)) {
    had.seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had.seed) old.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit({
      if (had.seed) {
        assign(".Random.seed", old.seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  lvls <- unique(grp[!is.na(grp)])
  for (L in lvls) {
    cols <- which(grp == L); if (!length(cols)) next
    sub <- mat[, cols, drop = FALSE]
    na  <- is.na(sub)
    if (!any(na)) next
    grp.obs  <- rowSums(!na)                       # observed count in this group per row
    grp.mean <- rowMeans(sub, na.rm = TRUE)
    grp.sd   <- apply(sub, 1, stats::sd, na.rm = TRUE)

    # PhosR-style site/condition-specific stochastic fill. Requiring two
    # observations avoids estimating a variance from a singleton.
    stochastic <- na & matrix(grp.obs >= 2L & is.finite(grp.sd),
                              nrow = nrow(sub), ncol = ncol(sub))
    if (any(stochastic)) {
      stochastic.rows <- row(stochastic)[stochastic]
      sub[stochastic] <- stats::rnorm(
        sum(stochastic),
        mean = grp.mean[stochastic.rows],
        sd = grp.sd[stochastic.rows]
      )
    }

    # Whole-group gaps (and singleton-observed groups in permissive filter
    # modes) do not provide a defensible site/group SD and remain MNAR/LoD.
    censored <- na & !stochastic
    if (any(censored)) sub[censored] <- row.floor[row(censored)[censored]]
    mat[, cols] <- sub
  }
  if (any(is.na(mat))) mat <- suppressWarnings(ReplaceMissingByLoD(mat))  # residual all-NA rows
  mat
}

ImputeMissingVarPhospho <- function(dataName = "", method = "min", min.obs.per.group = 2,
                                    impute.seed = 123L) {
  dataSet <- readDataset(dataName);
  msgSet <- readSet(msgSet, "msgSet");
  paramSet <- readSet(paramSet, "paramSet");

  int.mat <- dataSet$data.norm;
  if (is.null(int.mat)) {
    msgSet$current.msg <- c(msgSet$current.msg, "No phospho matrix available for imputation.");
    saveSet(msgSet, "msgSet");
    return(0);
  }
  # Treat strict zeros as missing (common in proteomics/phospho exports)
  int.mat[int.mat == 0] <- NA

  current.msg <- msgSet$current.msg;
  new.mat <- NULL;

  # DETECTION FILTER (before imputation). Phospho/MaxQuant exports are zero-inflated:
  # a site quantified in only one condition (e.g. one SILAC channel) is otherwise carried
  # into min/LoD imputation, which fabricates an extreme fold change (a one-channel site
  # -> log2FC ~ -20) that can pass FDR spuriously. Requiring the site to be genuinely
  # observed in >= min.obs.per.group samples within EACH group level removes these on/off
  # artifacts up front (validated on PXD005536: nested-interaction FDR hits 3 -> 0).
  # Group = the first (primary) metadata factor; fail-open if absent.
  grp <- tryCatch({
    mi <- dataSet$meta.info
    if (!is.null(mi) && NCOL(mi) >= 1) { g <- as.character(mi[[1]]); names(g) <- rownames(mi); g[colnames(int.mat)] } else NULL
  }, error = function(e) NULL)
  mog.env <- suppressWarnings(as.numeric(Sys.getenv("PA_PTM_MINOBS", unset = "")))
  if (is.finite(mog.env) && mog.env > 0) min.obs.per.group <- mog.env
  seed.env <- Sys.getenv("PA_PTM_IMPUTE_SEED", unset = "")
  if (nzchar(seed.env)) impute.seed <- suppressWarnings(as.integer(seed.env))
  impute.seed <- suppressWarnings(as.integer(impute.seed[1]))
  if (length(impute.seed) != 1L || !is.finite(impute.seed)) impute.seed <- 123L
  if (!is.null(grp) && !all(is.na(grp)) && is.numeric(min.obs.per.group) && min.obs.per.group > 0) {
    lvls <- unique(grp[!is.na(grp)])
    obs.ok <- sapply(lvls, function(L) {
      cols <- which(grp == L)
      if (!length(cols)) rep(TRUE, nrow(int.mat)) else rowSums(!is.na(int.mat[, cols, drop = FALSE])) >= min.obs.per.group
    })
    # Filter mode (env PA_PTM_FILTER_MODE), a sensitivity<->specificity dial:
    #  "all" (default): observed in >= min.obs.per.group in EVERY group.
    #     Conservative/FDR-controlled; drops genuine one-group-censored (MNAR)
    #     positives along with on/off artifacts.
    #  "anchor": keep a site with at least one FULLY observed group (a robust
    #     anchor); its censored group is then LoD-imputed by the group-aware min
    #     path, recovering the large true fold change. Preserves ranking while
    #     lifting sensitivity ~3-4x, at a moderate FDR cost (validated on
    #     Hogrebe). Only sound with a censoring-aware imputation (min/mindet).
    #  "any": >= min.obs.per.group in AT LEAST ONE group -- loosest; admits
    #     noise-censored nulls and degrades ranking. Not recommended.
    filt.mode <- Sys.getenv("PA_PTM_FILTER_MODE", unset = "all")
    if (identical(filt.mode, "anchor")) {
      full.ok <- sapply(lvls, function(L) {
        cols <- which(grp == L)
        if (!length(cols)) rep(FALSE, nrow(int.mat)) else rowSums(!is.na(int.mat[, cols, drop = FALSE])) == length(cols)
      })
      keep <- if (is.matrix(full.ok)) rowSums(full.ok) >= 1 else full.ok
    } else if (identical(filt.mode, "any")) {
      keep <- if (is.matrix(obs.ok)) rowSums(obs.ok) >= 1 else obs.ok
    } else {
      keep <- if (is.matrix(obs.ok)) rowSums(obs.ok) == length(lvls) else obs.ok
    }
    n.drop <- sum(!keep)
    if (n.drop > 0 && sum(keep) > 0) {
      int.mat <- int.mat[keep, , drop = FALSE]
      current.msg <- c(current.msg, paste0(
        "Detection filter removed ", n.drop, " phosphosite(s) not observed in at least ",
        min.obs.per.group, " samples in every group (before imputation)."))
    }
  }

  if(method %in% c("none", "NA", "")){
    # Keep missing values: the detection filter above still applies; limma
    # handles remaining NAs per-feature. Avoids LoD fills inside high groups
    # (a fill at the row floor inside an observed-high group explodes the
    # within-group variance and destroys power under group-dependent
    # missingness).
    new.mat <- int.mat;
    current.msg <- c(current.msg, "Imputation skipped ('none'); missing values retained for the model.");
  }else if(method=="exclude"){
    # OPTIMIZED: Use rowSums instead of apply for 60-100x speedup
    good.inx<-rowSums(is.na(int.mat))==0
    new.mat<-int.mat[good.inx,, drop=FALSE];
    current.msg <- c(current.msg ,"Variables with missing values were excluded.")
  }else if(method=="min"){
    # Group-aware LoD. After the per-group detection filter above, any remaining
    # missing value sits in a group that is still observed (its on/off artifacts
    # were already dropped). Filling such a sporadic gap with a global floor
    # fabricates a within-group low outlier -> spurious fold change and false
    # positives in true-null sites (Hogrebe ground-truth phospho: blanket LoD
    # collapses AUC 0.92 -> 0.59). So apply the left-censored floor ONLY to
    # values whose entire group is missing (MNAR); fill sporadic within-group
    # gaps by sampling N(group mean, group SD), as in PhosR::scImpute. This
    # avoids the variance deflation of deterministic group-mean replacement.
    new.mat <- .imputeMinGroupAware(int.mat, grp, seed = impute.seed)
    current.msg <- c(current.msg, paste0(
      "Missing values imputed by a group-aware hybrid (censored/insufficient group -> LoD floor; ",
      "sporadic within-group gap -> N(group mean, group SD); seed=", impute.seed, ")."));
  }else if(method %in% c("mindet","minprob","qrilc")){
    if (requireNamespace("imputeLCMD", quietly = TRUE)) {
      if (method == "mindet") {
        new.mat <- t(imputeLCMD::impute.MinDet(t(int.mat)))
        current.msg <- c(current.msg, "Missing variables were imputed using MinDet (deterministic minimum).");
      } else if (method == "minprob") {
        new.mat <- t(imputeLCMD::impute.MinProb(t(int.mat)))
        current.msg <- c(current.msg, "Missing variables were imputed using MinProb (stochastic minimum).");
      } else if (method == "qrilc") {
        new.mat <- t(imputeLCMD::impute.QRILC(t(int.mat))[[1]])
        current.msg <- c(current.msg, "Missing variables were imputed using QRILC (quantile regression for left-censored data).");
      }
    } else {
      AddErrMsg("MNAR imputation methods require the 'imputeLCMD' package. Install with: BiocManager::install('imputeLCMD')"); return(0);
    }
  }else if(method=="knn_var"){
    # impute.knn returns features x samples; do NOT transpose (see ImputeMissingVar).
    new.mat<-impute::impute.knn(as.matrix(int.mat))$data;
    current.msg <- c(current.msg, "Missing variables were imputed using KNN (feature-wise)");
  }else if(method=="knn_smp"){
    new.mat<-t(impute::impute.knn(data.matrix(t(int.mat)))$data);
    current.msg <- c(current.msg, "Missing variables were imputed using KNN (sample-wise)");
  }else if(method=="bpca"){
    new.mat<-pcaMethods::pca(t(int.mat), nPcs =5, method="bpca", center=T)@completeObs;
    new.mat = t(new.mat)
    current.msg <- c(current.msg, "Missing variables were imputed using BPCA");
  }else if(method=="ppca"){
    new.mat<-pcaMethods::pca(t(int.mat), nPcs =5, method="ppca", center=T)@completeObs;
    new.mat = t(new.mat)
    current.msg <- c(current.msg, "Missing variables were imputed using PPCA");
  }else if(method=="svdImpute"){
    new.mat<-pcaMethods::pca(t(int.mat), nPcs =5, method="svdImpute", center=T)@completeObs;
    new.mat = t(new.mat)
    current.msg <- c(current.msg, "Missing variables were imputed using SVD Impute");
  }else if(method=="seqknn"){
    # imputeLCMD has no sequential-KNN implementation; use multiUS like the
    # non-phospho path. SeqKNN imputes each feature from its k nearest feature
    # rows, so the features x samples matrix is passed untransposed.
    if (!requireNamespace("multiUS", quietly = TRUE)) {
      AddErrMsg("SeqKNN requires the 'multiUS' package. Install with: install.packages('multiUS')"); return(0);
    }
    k.seqknn <- min(10L, nrow(int.mat) - 1L)
    new.mat <- multiUS::seqKNNimp(int.mat, k = k.seqknn)
    current.msg <- c(current.msg, "Missing variables were imputed using SeqKNN.");
  }else{
    new.mat<- suppressWarnings(ReplaceMissingByLoD(int.mat));
    current.msg <- c(current.msg, paste0("Unknown method '", method, "'. Defaulted to LoD replacement."));
  }

  dataSet$data.norm <- new.mat;
  ov_qs_save(new.mat, "data.raw.qs");
  fast.write(sanitizeSmallNumbers(new.mat), file="data_imputed.csv");
  msgSet$current.msg <- current.msg;
  saveSet(msgSet, "msgSet");
  # Same record as ImputeMissingVar, so a phospho run states its imputation method
  # on the report and slides too.
  paramSet$impute.opt <- method;
  if (identical(method, "min")) paramSet$impute.seed <- impute.seed;
  saveSet(paramSet, "paramSet");
  return(RegisterData(dataSet));
}
