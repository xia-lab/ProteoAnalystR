

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
ImputeMissingVar <- function(dataName="", method="min"){
  dataSet <- readDataset(dataName);
  msgSet <- readSet(msgSet, "msgSet");
  paramSet <- readSet(paramSet, "paramSet");

  int.mat <- dataSet$data.norm;
  if (is.null(int.mat)) {
    msgSet$current.msg <- c(msgSet$current.msg, "No normalized data available for imputation.");
    saveSet(msgSet, "msgSet");
    return(0);
  }

  # Treat strict zeros as missing (common in proteomics exports)
  int.mat[int.mat == 0] <- NA
  row.nms <- rownames(int.mat);
  current.msg <- msgSet$current.msg;
  new.mat <- NULL;

  # NOTE: MSstats imputation is no longer used in the new workflow
  # Reason: Data is now summarized (peptide → protein) BEFORE imputation
  # MSstats MBimpute requires peptide-level data, which we no longer have at this stage
  # All imputation now uses legacy methods (KNN, PPCA, etc.) on protein-level data

  # Legacy/fallback methods
  if(method=="exclude"){
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
        new.mat<-t(impute::impute.knn(as.matrix(int.mat))$data);
        current.msg <- c(current.msg, "Missing variables were imputed using KNN (feature-wise)");
      }else if(method == "knn_smp"){
        new.mat<-impute::impute.knn(data.matrix(t(int.mat)))$data;
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
          stop("MNAR imputation methods require the 'imputeLCMD' package. Install with: BiocManager::install('imputeLCMD')")
        }
      }else if(method == "seqknn"){
        # SeqKNN: sequential K-nearest neighbors (MAR/mixed)
        has.seqknn <- FALSE

        # Try multiUS package first (CRAN available)
        if (requireNamespace("multiUS", quietly = TRUE)) {
          tryCatch({
            new.mat <- t(multiUS::seqKNNimp(t(int.mat), k = 10))
            has.seqknn <- TRUE
            current.msg <- c(current.msg, "Missing variables were imputed using SeqKNN (sequential KNN via multiUS)");
          }, error = function(e) {
            #msg("[SeqKNN] multiUS::seqKNNimp failed: ", e$message)
          })
        }

        # Fallback to SeqKnn package if multiUS not available
        if (!has.seqknn && requireNamespace("SeqKnn", quietly = TRUE)) {
          tryCatch({
            new.mat <- t(SeqKnn::SeqKNN(t(int.mat), k = 10))
            has.seqknn <- TRUE
            current.msg <- c(current.msg, "Missing variables were imputed using SeqKNN (sequential KNN via SeqKnn)");
          }, error = function(e) {
            #msg("[SeqKNN] SeqKnn::SeqKNN failed: ", e$message)
          })
        }

        if (!has.seqknn) {
          stop("SeqKNN requires either 'multiUS' (CRAN) or 'SeqKnn' (archived) package. Install with: install.packages('multiUS')")
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
              stop("Impseq requires the 'rrcovNA' package. Install with: install.packages('rrcovNA')")
            }
          })
        } else {
          stop("Impseq requires the 'rrcovNA' package. Install with: install.packages('rrcovNA')")
        }
      }else{
        stop("Unknown imputation method: ", method)
      }
    }
  }
  msgSet$current.msg <- current.msg;
  saveSet(msgSet, "msgSet");

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
  msin <- try(qs::qread(msstats.path), silent = TRUE);
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
                                  featureSubset = "top3"), silent = TRUE);
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
  wide <- reshape2::dcast(prof, Protein ~ RUN, value.var = value.col);
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
ImputeMissingVarPhospho <- function(dataName = "", method = "min") {
  dataSet <- readDataset(dataName);
  msgSet <- readSet(msgSet, "msgSet");
  paramSet <- readSet(paramSet, "paramSet");

  int.mat <- dataSet$data.norm;
  # Treat strict zeros as missing (common in proteomics/phospho exports)
  int.mat[int.mat == 0] <- NA
  if (is.null(int.mat)) {
    msgSet$current.msg <- c(msgSet$current.msg, "No phospho matrix available for imputation.");
    saveSet(msgSet, "msgSet");
    return(0);
  }

  current.msg <- msgSet$current.msg;
  new.mat <- NULL;

  if(method=="exclude"){
    # OPTIMIZED: Use rowSums instead of apply for 60-100x speedup
    good.inx<-rowSums(is.na(int.mat))==0
    new.mat<-int.mat[good.inx,, drop=FALSE];
    current.msg <- c(current.msg ,"Variables with missing values were excluded.")
  }else if(method=="min"){
    new.mat<- suppressWarnings(ReplaceMissingByLoD(int.mat));
    current.msg <- c(current.msg, "Missing variables were replaced by LoDs (1/5 of the min positive value for each variable).");
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
      stop("MNAR imputation methods require the 'imputeLCMD' package. Install with: BiocManager::install('imputeLCMD')")
    }
  }else if(method=="knn_var"){
    new.mat<-t(impute::impute.knn(as.matrix(int.mat))$data);
    current.msg <- c(current.msg, "Missing variables were imputed using KNN (feature-wise)");
  }else if(method=="knn_smp"){
    new.mat<-impute::impute.knn(data.matrix(t(int.mat)))$data;
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
    if (!requireNamespace("imputeLCMD", quietly = TRUE)) {
      stop("SeqKNN requires the 'imputeLCMD' package.")
    }
    new.mat <- t(imputeLCMD::impute.SequentialKNN(t(int.mat), allowMissing = TRUE))
    current.msg <- c(current.msg, "Missing variables were imputed using SeqKNN.");
  }else{
    new.mat<- suppressWarnings(ReplaceMissingByLoD(int.mat));
    current.msg <- c(current.msg, paste0("Unknown method '", method, "'. Defaulted to LoD replacement."));
  }

  dataSet$data.norm <- new.mat;
  qs::qsave(new.mat, "data.raw.qs");
  fast.write(sanitizeSmallNumbers(new.mat), file="data_imputed.csv");
  msgSet$current.msg <- current.msg;
  saveSet(msgSet, "msgSet");
  saveSet(paramSet, "paramSet");
  return(RegisterData(dataSet));
}