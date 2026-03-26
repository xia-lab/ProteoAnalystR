##################################################
## R script for ProteoAnalyst
## Description: functions for quality check boxplot
## Authors: 
## Jeff Xia, jeff.xia@mcgill.ca
## Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
###################################################


#'Perform data normalization
#'@description Filtering and Normalizing gene expression data
#'@param norm.opt Normalization method to be used
#'@param var.thresh Variance threshold
#'@param abundance Relative abundance threshold
#'@param count.thresh Count threshold for RNA-seq data and abundance threshold for microarray data
#'@param filterUnmapped, boolean, whether to filter unmapped features or not
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: MIT
#'@export
#'
PerformNormalization <- function(dataName, norm.opt, var.thresh, count.thresh, filterUnmapped,
                                 islog = "false", countOpt = "sum", removeMissing = "false", missingPercent = 50,
                                 sampleNormOpt = "none", sampleNormParam = "__manual__") {
  paramSet <- readSet(paramSet, "paramSet");
  msgSet   <- readSet(msgSet, "msgSet");
  dataSet  <- readDataset(dataName);
  msg <- ""
  msgSet$sample.norm.msg <- character(0)

  #msg("[Norm] PerformNormalization start data=", dataName, " type=", dataSet$type, " norm.opt=", norm.opt,
  #        " var.thresh=", var.thresh, " count.thresh=", count.thresh, " filterUnmapped=", filterUnmapped)

  # Read the current normalization input, which is rebuilt from the immutable
  # annotated baseline on every normalization submit.
  if(file.exists("norm.input.anot.qs")){
    ds <- qs::qread("norm.input.anot.qs");
    #msg("[Norm] Using norm.input.anot.qs (current normalization input)")
  } else if(file.exists("orig.data.anot.qs")){
    ds <- qs::qread("orig.data.anot.qs");
    #msg("[Norm] Using orig.data.anot.qs (annotated baseline)")
  } else {
    # Fallback for backward compatibility (older data that doesn't have orig.data.anot.qs)
    ds <- qs::qread("data.anot.qs");
    #msg("[Norm] Fallback: using data.anot.qs")
  }
  #msg("[Norm] Loaded data head rows=", paste(utils::head(rownames(ds), 5), collapse=", "),
  #        " cols=", paste(utils::head(colnames(ds), 5), collapse=", "))
  # Ensure the dataset carries the annotated matrix before filtering/normalization
  if (isTRUE(dataSet$annotated)) {
    dataSet$data.norm <- ds
    #msg("[Norm] data.norm refreshed from annotated matrix; head IDs=", paste(utils::head(rownames(dataSet$data.norm), 5), collapse=", "))
  }
  msgSet$current.msg <- c(msgSet$current.msg, paste0("Diagnostic before summarization: samples=", ncol(ds),
                                                    " features=", nrow(ds),
                                                    " zeros=", sum(ds==0),
                                                    " min=", min(ds, na.rm=T),
                                                    " max=", max(ds, na.rm=T)))
  saveSet(msgSet,"msgSet")
    

  data <- PerformFiltering(dataSet, var.thresh, count.thresh, filterUnmapped, countOpt,
                           removeMissing = removeMissing, missingPercent = missingPercent)

  # Check if we have data after filtering
  if (nrow(data) == 0 || ncol(data) == 0) {
    cause.msg <- if (ncol(data) == 0) {
      "No columns matched between data and metadata sample names."
    } else {
      "All features were removed by filtering thresholds (abundance/missing/variance)."
    }
    error.msg <- paste0("ERROR: No data remaining after filtering! ",
                       "rows=", nrow(data), " cols=", ncol(data), ". ",
                       cause.msg, " Please check sample names and/or relax filtering thresholds.")
    msgSet$current.msg <- c(msgSet$current.msg, error.msg)
    saveSet(msgSet, "msgSet")
    stop(error.msg)
  }

  .save.annotated.data(data)
  msg <- paste(filt.msg, msg)

  data <- sanitizeSmallNumbers(data)
  diag.filtered <- paste0("Diagnostic after filtering: rows=", nrow(data),
                          " cols=", ncol(data),
                          " zeros=", sum(data == 0, na.rm=T))
  msgSet$current.msg <- c(msgSet$current.msg, diag.filtered)
  #msg("[Norm] After filtering head IDs=", paste(utils::head(rownames(data), 5), collapse=", "),
  #        " head samples=", paste(utils::head(colnames(data), 5), collapse=", "))
  saveSet(msgSet, "msgSet")

  # Apply sample normalization (row-wise) if specified
  if (!is.null(sampleNormOpt) && sampleNormOpt != "none" && sampleNormOpt != "NA") {
    data <- ApplySampleNormalization(data, sampleNormOpt, sampleNormParam)
    paramSet$sample.norm.opt <- sampleNormOpt
    paramSet$sample.norm.param <- sampleNormParam
    saveSet(paramSet, "paramSet")
  }

  if (dataSet$type == "prot") {
    diag.prot <- paste0("Diagnostic before normalization (prot): norm.opt=", norm.opt,
                        " rows=", nrow(data), " cols=", ncol(data))
    msgSet$current.msg <- c(msgSet$current.msg, diag.prot)
    saveSet(msgSet, "msgSet")

    # MSstats-based normalization + summarization (uses saved msstats_input.qs from data loading)
    if (identical(norm.opt, "msstats")) {
      msstats.path <- "msstats_input.qs"
      if (file.exists(msstats.path)) {
        #msg("[Norm] Using MSstats normalization + summarization (equalizeMedians + TMP)...")
        ms.out <- .normalizeWithMSstats(msstats.path, dataSet$meta.info)
        if (!is.null(ms.out$mat)) {
          data <- ms.out$mat
          msg <- paste("[MSstats] equalizeMedians + TMP summarization applied.", msg)
          #msg("[Norm] MSstats complete: ", nrow(data), " proteins, ", ncol(data), " samples")
        } else {
          msgSet$current.msg <- c(msgSet$current.msg, "MSstats normalization failed; keeping filtered matrix.")
          saveSet(msgSet, "msgSet")
        }
      } else {
        msgSet$current.msg <- c(msgSet$current.msg, "MSstats normalization skipped: msstats_input.qs not found (data may already be protein-level).")
        saveSet(msgSet, "msgSet")
      }
    }

    if (islog == "true" || norm.opt == "Rlr" || norm.opt == "Loess") {
      data <- NormalizeData(data, "log", "NA", "NA")
      if (length(data) == 1 && identical(as.numeric(data), 0)) {
        return(0)
      }
      msg  <- paste(norm.msg, msg)
    }
  }

  paramSet$norm.opt   <- norm.opt
  paramSet$var.perc   <- var.thresh
  paramSet$abun.perc  <- count.thresh

  if (identical(norm.opt, "MORlog")) {

    if (!requireNamespace("DESeq2", quietly = TRUE)) {
      AddErrMsg("MORlog normalization requires the 'DESeq2' package. Please install it.")
      return(0)
    }

    m <- as.matrix(data)

    # basic checks for counts
    if (any(m < 0, na.rm = TRUE)) {
      AddErrMsg("MORlog expects non-negative count data.")
      return(0)
    }
    # ensure integer-like counts for DESeq2
    if (!is.integer(m)) {
      m <- round(m)
    }

    # minimal colData (no outcome needed for size factors)
    # rownames must match sample names (columns of count matrix)
    cd <- S4Vectors::DataFrame(row.names = colnames(m))

    dds <- DESeq2::DESeqDataSetFromMatrix(countData = m, colData = cd, design = ~ 1)
    dds <- DESeq2::estimateSizeFactors(dds)
    norm_counts <- DESeq2::counts(dds, normalized = TRUE)
    data <- log2(norm_counts + 1)

    msg <- paste("[MORlog] Applied DESeq2 median-of-ratios size-factor normalization and log2(x+1).", msg)
  } else if (identical(norm.opt, "msstats")) {
    # Already normalized via MSstats above (for proteomics only)
    msg <- msg
  } else {
    # ---- existing generic normalization ----
    data <- NormalizeData(data, norm.opt, "NA", "NA")
    if (length(data) == 1 && identical(as.numeric(data), 0)) {
      return(0)
    }
    msg  <- paste(norm.msg, msg)
  }


  if (paramSet$oneDataAnalType == "dose" && min(data) < 0) {
    add.val <- abs(min(data)) + 0.05 * abs(min(data))
    data <- data + add.val
  }

  # Final safeguard: remap rownames to annotated IDs if available
  if (file.exists("annotation.qs")) {
    anot.id <- qs::qread("annotation.qs")
    if (!is.null(anot.id)) {
      if (!is.null(names(anot.id))) {
        mapped <- anot.id[rownames(data)]
      } else {
        mapped <- anot.id
      }
      hit <- !is.na(mapped)
      if (any(hit)) {
        rownames(data)[hit] <- mapped[hit]
        #msg("[Norm] Rownames remapped to annotated IDs; head now: ",
        #        paste(utils::head(rownames(data), 5), collapse=", "))
        # If remapping introduced duplicates, aggregate them
        if (any(duplicated(rownames(data)))) {
          dup.n <- sum(duplicated(rownames(data)))
          #msg("[Norm] Found ", dup.n, " duplicate IDs after remap; aggregating by mean")

          # Safety check: only aggregate if we have columns
          if (ncol(data) > 0) {
            df <- as.data.frame(data)
            df$Protein <- rownames(data)
            df <- stats::aggregate(. ~ Protein, data = df, FUN = function(x) mean(x, na.rm = TRUE))
            rownames(df) <- df$Protein
            df$Protein <- NULL
            data <- as.matrix(df)
          } else {
            #msg("[Norm] WARNING: Cannot aggregate duplicates - data has 0 columns!")
          }
        }
      }
    }
  }

  dataSet$data.norm <- data
  fast.write(sanitizeSmallNumbers(data), file = "data_normalized.csv")
  qs::qsave(data, file = "data.stat.qs")

  msgSet <- readSet(msgSet, "msgSet")
  sample_norm_msg <- msgSet$sample.norm.msg
  if (is.null(sample_norm_msg)) {
    sample_norm_msg <- character(0)
  }
  msgSet$current.msg <- c(sample_norm_msg, msg)
  saveSet(msgSet,   "msgSet")
  saveSet(paramSet, "paramSet")
  return(RegisterData(dataSet))
}



#' Apply format-specific filtering based on data format
#' Called before general filtering to apply MaxQuant/DIA-NN/FragPipe/Spectronaut specific filters
#' @param dataName Dataset name
#' @param dataFormat Format type (maxquant, diann, fragpipe, spectronaut)
#' @param formatOpts List of format-specific options
ApplyFormatSpecificFiltering <- function(dataName, dataFormat, formatOpts = list()) {
  paramSet <- readSet(paramSet, "paramSet")
  msgSet <- readSet(msgSet, "msgSet")

  # Always rebuild from the immutable annotated baseline so repeated
  # normalization submits do not accumulate the previous filter state.
  if (file.exists("orig.data.anot.qs")) {
    data <- qs::qread("orig.data.anot.qs")
  } else {
    msgSet$current.msg <- c(msgSet$current.msg, "[Format Filter] No orig.data.anot.qs found, skipping format-specific filtering")
    saveSet(msgSet, "msgSet")
    return(1L)
  }

  original_count <- nrow(data)
  result <- NULL
  stats <- NULL

  # Apply format-specific filtering based on dataFormat
  if (tolower(dataFormat) == "maxquant") {
    removeContaminants <- if ("removeContaminants" %in% names(formatOpts)) {
      isTRUE(as.logical(formatOpts$removeContaminants))
    } else { FALSE }

    removeDecoys <- if ("removeDecoys" %in% names(formatOpts)) {
      isTRUE(as.logical(formatOpts$removeDecoys))
    } else { FALSE }
    removeOnlyBySite <- if ("removeOnlyBySite" %in% names(formatOpts)) {
      isTRUE(as.logical(formatOpts$removeOnlyBySite))
    } else { TRUE }

    minPeptides <- if ("minPeptides" %in% names(formatOpts)) {
      as.integer(formatOpts$minPeptides)
    } else { 1 }

    result <- ApplyMaxQuantFiltering(data, removeContaminants, removeDecoys, removeOnlyBySite, minPeptides)
    data <- result$data
    stats <- result$stats

  } else if (tolower(dataFormat) == "diann") {
    qvalueFilter <- if ("qvalueFilter" %in% names(formatOpts)) {
      isTRUE(as.logical(formatOpts$qvalueFilter))
    } else { FALSE }
    pepFilter <- if ("pepFilter" %in% names(formatOpts)) {
      isTRUE(as.logical(formatOpts$pepFilter))
    } else { FALSE }
    pepThreshold <- if ("pepThreshold" %in% names(formatOpts)) {
      as.numeric(formatOpts$pepThreshold)
    } else { 0.01 }
    minPeptides <- if ("minPeptides" %in% names(formatOpts)) {
      as.integer(formatOpts$minPeptides)
    } else { 0 }

    result <- ApplyDiannFiltering(data, qvalueFilter, pepFilter, pepThreshold, minPeptides)
    data <- result$data
    stats <- result$stats

  } else if (tolower(dataFormat) == "fragpipe") {
    removeContaminants <- if ("removeContaminants" %in% names(formatOpts)) {
      isTRUE(as.logical(formatOpts$removeContaminants))
    } else { FALSE }

    minProb <- if ("minProb" %in% names(formatOpts)) {
      as.numeric(formatOpts$minProb)
    } else { 0.0 }

    minPeptides <- if ("minPeptides" %in% names(formatOpts)) {
      as.integer(formatOpts$minPeptides)
    } else { 1 }

    result <- ApplyFragpipeFiltering(data, removeContaminants, minProb, minPeptides)
    data <- result$data
    stats <- result$stats

  } else if (tolower(dataFormat) == "spectronaut") {
    removeContaminants <- if ("removeContaminants" %in% names(formatOpts)) {
      isTRUE(as.logical(formatOpts$removeContaminants))
    } else { FALSE }
    qvalueFilter <- if ("spectronautQValueFilter" %in% names(formatOpts)) {
      isTRUE(as.logical(formatOpts$spectronautQValueFilter))
    } else { TRUE }
    qvalueThreshold <- if ("spectronautQValueThreshold" %in% names(formatOpts)) {
      as.numeric(formatOpts$spectronautQValueThreshold)
    } else { 0.01 }
    minPeptides <- if ("spectronautMinPeptides" %in% names(formatOpts)) {
      as.integer(formatOpts$spectronautMinPeptides)
    } else { 0 }

    result <- ApplySpectronautFiltering(data, removeContaminants, qvalueFilter, qvalueThreshold, minPeptides)
    data <- result$data
    stats <- result$stats
  }

  # Save the current filtered input separately from the immutable baseline.
  filtered_count <- nrow(data)

  # Create detailed pre-filtering message similar to phospho module
  if (!is.null(stats)) {
    msg_parts <- c(sprintf("Format-specific filtering (%s): Starting with %d proteins.", dataFormat, stats$n_total))

    # Add format-specific filter messages
    if (tolower(dataFormat) == "maxquant") {
      if (stats$n_contaminants > 0) {
        msg_parts <- c(msg_parts, sprintf("Removed %d contaminants.", stats$n_contaminants))
      }
      if (stats$n_decoys > 0) {
        msg_parts <- c(msg_parts, sprintf("Removed %d decoys/reverse hits.", stats$n_decoys))
      }
      if (stats$n_only_by_site > 0) {
        msg_parts <- c(msg_parts, sprintf("Removed %d proteins identified only by site.", stats$n_only_by_site))
      }
      if (stats$n_min_peptides > 0) {
        msg_parts <- c(msg_parts, sprintf("Removed %d proteins with insufficient peptides.", stats$n_min_peptides))
      }
    } else if (tolower(dataFormat) == "diann") {
      if (stats$n_qvalue > 0) {
        msg_parts <- c(msg_parts, sprintf("Removed %d proteins with Q-value >= 0.01.", stats$n_qvalue))
      }
      if (stats$n_pep > 0) {
        msg_parts <- c(msg_parts, sprintf("Removed %d proteins with PEP >= %.2f.", stats$n_pep, stats$pep_threshold))
      }
      if (stats$n_min_peptides > 0) {
        msg_parts <- c(msg_parts, sprintf("Removed %d proteins with insufficient peptide evidence.", stats$n_min_peptides))
      }
    } else if (tolower(dataFormat) == "fragpipe") {
      if (stats$n_contaminants > 0) {
        msg_parts <- c(msg_parts, sprintf("Removed %d contaminants.", stats$n_contaminants))
      }
      if (stats$n_min_prob > 0) {
        msg_parts <- c(msg_parts, sprintf("Removed %d proteins below probability threshold (%.2f).", stats$n_min_prob, stats$min_prob))
      }
      if (stats$n_min_peptides > 0) {
        msg_parts <- c(msg_parts, sprintf("Removed %d proteins with insufficient peptides.", stats$n_min_peptides))
      }
    } else if (tolower(dataFormat) == "spectronaut") {
      if (stats$n_contaminants > 0) {
        msg_parts <- c(msg_parts, sprintf("Removed %d contaminants.", stats$n_contaminants))
      }
      if (stats$n_qvalue > 0) {
        msg_parts <- c(msg_parts, sprintf("Removed %d proteins with q-value >= %.2f.", stats$n_qvalue, stats$qvalue_threshold))
      }
      if (stats$n_min_peptides > 0) {
        msg_parts <- c(msg_parts, sprintf("Removed %d proteins with insufficient peptide evidence.", stats$n_min_peptides))
      }
    }

    msg_parts <- c(msg_parts, sprintf("Remaining: %d proteins.", stats$n_final))
    prefilter_msg <- paste(msg_parts, collapse = " ")
    msgSet$prefilter.msg <- prefilter_msg
  } else {
    # Fallback to simple message if stats unavailable
    if (filtered_count < original_count) {
      removed_count <- original_count - filtered_count
      prefilter_msg <- paste0("Pre-filtering (", dataFormat, "): Removed ", removed_count,
                             " protein", ifelse(removed_count > 1, "s", ""),
                             " (", original_count, " → ", filtered_count, ")")
      msgSet$prefilter.msg <- prefilter_msg
    } else {
      prefilter_msg <- paste0("Pre-filtering (", dataFormat, "): No proteins removed (",
                             original_count, " proteins passed all filters)")
      msgSet$prefilter.msg <- prefilter_msg
    }
  }

  qs::qsave(data, "norm.input.anot.qs")
  saveSet(msgSet, "msgSet")

  return(1L)
}

PerformFiltering <- function(dataSet, var.thresh, count.thresh, filterUnmapped, countOpt,
                             removeMissing = "false", missingPercent = 50){
  msg <- "";
  remove.missing.flag <- isTRUE(removeMissing) ||
    identical(removeMissing, "true") ||
    identical(removeMissing, "TRUE")
  msgSet <- readSet(msgSet, "msgSet")

  if (isTRUE(dataSet$annotated)) {
    raw.data.anot <- dataSet$data.norm
    msg <- "[Filter] Using annotated matrix from dataset; skipping remap via annotation.qs."
  } else if(filterUnmapped == "false"){
    # need to update those with annotations
    data1 <- qs::qread("data.raw.qs");
    colnames(data1) <- colnames(dataSet$data.norm)
    anot.id <- qs::qread("annotation.qs");
    # Map by names when available to avoid positional mismatches
    if (!is.null(names(anot.id)) && length(intersect(rownames(data1), names(anot.id))) > 0) {
      mapped <- anot.id[rownames(data1)]
      hit.inx <- !is.na(mapped)
      rownames(data1)[hit.inx] <- mapped[hit.inx]
    } else {
      hit.inx <- !is.na(anot.id)
      rownames(data1)[hit.inx] <- anot.id[hit.inx]
    }
    res <- RemoveDuplicates(data1, "mean", quiet=T, paramSet, msgSet);
    data1 <- res[[1]];
    msgSet <- res[[2]];
    raw.data.anot <- data1;
    msg <- "Only features with annotations are kept for further analysis.";
  }else{
    # Load annotated data
    # IMPORTANT: Always read from orig.data.anot.qs to ensure re-running normalization
    # starts from the original data, not already-filtered data
    if(dataSet$type=="prot"){
      # For proteomics/phospho: prefer the current normalization input,
      # then fall back to the immutable annotated baseline.
      if (file.exists("norm.input.anot.qs")) {
        raw.data.anot <- qs::qread("norm.input.anot.qs")
        #msg("[Filter] Loading norm.input.anot.qs for proteomics")
      } else if (file.exists("orig.data.anot.qs")) {
        raw.data.anot <- qs::qread("orig.data.anot.qs")
        #msg("[Filter] Loading orig.data.anot.qs for proteomics baseline")
      } else if (file.exists("data.annotated.qs")) {
        raw.data.anot <- qs::qread("data.annotated.qs")
        #msg("[Filter] Loading data.annotated.qs for proteomics (fallback)")
      } else if (file.exists("data.anot.qs")) {
        raw.data.anot <- qs::qread("data.anot.qs")
        #msg("[Filter] Loading data.anot.qs for proteomics (fallback)")
      } else if (file.exists("data.missed.qs")) {
        # Fallback to data.missed.qs if it exists (from previous workflow)
        raw.data.anot <- qs::qread("data.missed.qs")
        #msg("[Filter] Loading data.missed.qs for proteomics (fallback)")
      } else {
        stop("No annotated data file found for proteomics filtering")
      }
    }else{
     if (file.exists("norm.input.anot.qs")) {
       raw.data.anot <- qs::qread("norm.input.anot.qs");
     } else {
       raw.data.anot <- qs::qread("orig.data.anot.qs");
     }
    }
    # Ensure annotated IDs applied if annotation vector available
    if (file.exists("annotation.qs")) {
      anot.id <- qs::qread("annotation.qs")
      if (!is.null(names(anot.id)) && length(intersect(rownames(raw.data.anot), names(anot.id))) > 0) {
        mapped <- anot.id[rownames(raw.data.anot)]
        hit.inx <- !is.na(mapped)
        rownames(raw.data.anot)[hit.inx] <- mapped[hit.inx]
      }
    }
   # Set column names from dataSet$data.norm if it has columns
   if (!is.null(dataSet$data.norm) && ncol(dataSet$data.norm) > 0) {
     #msg("[Filter] Setting colnames from dataSet$data.norm: ", paste(head(colnames(dataSet$data.norm), 5), collapse=", "))
     colnames(raw.data.anot) <- colnames(dataSet$data.norm)
   } else {
     #msg("[Filter] WARNING: dataSet$data.norm is empty or NULL, keeping original colnames")
   }
  }

  data <- raw.data.anot;

  common.cols <- which(colnames(data) %in% rownames(dataSet$meta.info))

  data <- data[, common.cols, drop = FALSE]

  #msg("[Filter] Starting filtering: ", nrow(data), " features x ", ncol(data), " samples")
  #msg("[Filter] Data type: ", dataSet$type, ", count.thresh=", count.thresh, ", var.thresh=", var.thresh)

  # PERFORMANCE FIX (Issue #1): Use vectorized row operations instead of apply()
  # rowSums/rowMeans are 60-100x faster than apply(data, 1, sum/mean)
  # Critical for normalization - affects 100% of users
  if (dataSet$type == "count"){
    if (countOpt == "sum") {
        # Sum approach: sum counts across samples for each gene
        sum.counts <- rowSums(data, na.rm = TRUE)
        rm.inx <- is.na(sum.counts) | sum.counts < count.thresh
        msg <- paste(msg, "Filtered ", sum(rm.inx), " features with low counts using sum method.", collapse = " ")
    } else if (countOpt == "average") {
        # Average approach: calculate average counts across samples for each gene
        avg.counts <- rowMeans(data, na.rm = TRUE)
        rm.inx <- is.na(avg.counts) | avg.counts < count.thresh
        msg <- paste(msg, "Filtered ", sum(rm.inx), " features with low counts using average method.", collapse = " ")
    }
  }else{
    avg.signal <- rowMeans(data, na.rm=TRUE)
    abundance.pct <- count.thresh/100;
    p05 <- quantile(avg.signal, abundance.pct, na.rm = TRUE)
    all.rows <- nrow(data)

    #msg("[Filter] Abundance filtering: p05 threshold=", round(p05, 4),
    #        ", avg.signal range=[", round(min(avg.signal, na.rm=T), 4), ", ", round(max(avg.signal, na.rm=T), 4), "]",
    #        ", NAs=", sum(is.na(avg.signal)))

    # Handle NA/NaN values in avg.signal: treat them as failing the filter
    rm.inx <- is.na(avg.signal) | avg.signal < p05;
    #msg("[Filter] Will filter ", sum(rm.inx), " / ", length(rm.inx), " features (",
    #        round(100*sum(rm.inx)/length(rm.inx), 1), "%)")
    msg <- paste(msg, "Filtered ", sum(rm.inx), " features with low relative abundance (average expression signal).", collapse=" ");
  }

  # Apply abundance/count filtering (always apply, not just when var.thresh > 0)
  data <- data[!rm.inx,];

  # Drop features with missing values above the user-defined threshold
  # Only apply when the UI checkbox is enabled.
  rm.miss.msg <- NULL
  miss.pct.num <- suppressWarnings(as.numeric(missingPercent))

  # Apply filtering only when enabled and threshold is valid (>0)
  if (remove.missing.flag && !is.na(miss.pct.num) && miss.pct.num > 0) {
    # Count features before filtering
    features_before <- nrow(data)

    miss.ind <- is.na(data)
    # Handle character-coded missing values without spamming warnings
    if (is.character(data)) {
      miss.ind <- miss.ind | data == "NA"
    }
    zero.ind <- !is.na(data) & data == 0
    miss.ind <- miss.ind | zero.ind
    miss.frac <- rowMeans(miss.ind)

    rm.miss <- miss.frac > (miss.pct.num/100)

    if (any(rm.miss)) {
      data <- data[!rm.miss, , drop = FALSE]
      rm.miss.msg <- paste("Removed", sum(rm.miss), "features with >", miss.pct.num, "% missing values.")
      msg <- paste(msg, rm.miss.msg)
    }
  }

  # Check if we have data left after filtering
  if (nrow(data) == 0) {
    #msg("[Filter] WARNING: All features filtered out! Check abundance/count thresholds.")
    msgSet$current.msg <- c(msgSet$current.msg, "Warning: All features were filtered out by abundance/count filtering. Relaxing filters to keep at least some data.")
    saveSet(msgSet, "msgSet")
    # Restore original data to prevent empty matrix
    data <- raw.data.anot
    data <- data[,which(colnames(data)%in% rownames(dataSet$meta.info))]
    filt.msg <<- "WARNING: Abundance filtering removed all features - filters were skipped."
  } else if(var.thresh > 0){
    data.before.var <- data
    filter.val <- apply(data, 1, IQR, na.rm=T);
    nm <- "Interquantile Range";
    filter.val <- -filter.val
    rk <- rank(filter.val, ties.method='random');
    # remove constant values
    good.inx <- -filter.val > 0;
    kp.pct <- (100 - var.thresh)/100;
    remain <- rk < nrow(data)*kp.pct;
    initial_gene_count <- nrow(data)
    data <- data[remain&good.inx,];
    # Calculate number of features filtered by IQR
    filtered_by_iqr <- initial_gene_count - nrow(data)

    # If IQR filtering removes everything (common in sparse phospho data), keep pre-IQR data.
    if (nrow(data) == 0 && initial_gene_count > 0) {
      data <- data.before.var
      filt.msg <<- paste(msg, "WARNING: IQR variance filtering removed all features - skipped variance filtering for this dataset.")
    } else {
      # Update message with correct number of filtered features
      filt.msg <<- paste(msg, "Filtered", filtered_by_iqr, "low variance features based on IQR.")
    }
  }else{
    filt.msg <<- paste(msg, paste("Filtered 0 low variance features based on IQR"), collapse=" ");
  }
  
  return(data);
}

NormalizeDataMetaMode <-function (nm, opt, colNorm="NA", scaleNorm="NA"){
  if(nm == "NA"){
    paramSet <- readSet(paramSet, "paramSet");
    mdata.all <- paramSet$mdata.all;
    sel.nms <- names(mdata.all);
    for(i in 1:length(sel.nms)){
      dataName <- sel.nms[i];
      dataSet = readDataset(dataName);
      data.filtered <- readDataQs("data.filtered.qs", paramSet$anal.type, dataName);
      data <- NormalizeData(data.filtered,opt, colNorm, scaleNorm);
      if(length(data) == 1){
        return(0);
      }
    dataSet$data.norm <- data;
    dataSet$norm.opt <- opt;
    fast.write(sanitizeSmallNumbers(dataSet$data.norm), file = "data_normalized.csv")
    RegisterData(dataSet);
    }
    return(1)
  }else{
    dataSet <- readDataset(nm);
    data.filtered <- readDataQs("data.filtered.qs", paramSet$anal.type, nm);
    data <- NormalizeData(data.filtered,opt, colNorm, scaleNorm);
    if(length(data) == 1){
      return(0);
    }
    dataSet$data.norm <- data;
    dataSet$norm.opt <- opt;
    qs::qsave(data, file="data.stat.qs");
    return(RegisterData(dataSet));
    
  }
}

NormalizeData <-function (data, norm.opt, colNorm="NA", scaleNorm="NA"){
  msg <- ""
  row.nms <- rownames(data);
  col.nms <- colnames(data);
  msgSet <- readSet(msgSet, "msgSet");
  
  data <- sanitizeSmallNumbers(data)

  # column(sample)-wise normalization
  if(colNorm=="SumNorm"){
    data<-t(apply(data, 2, SumNorm));
    rownm<-"Normalization to constant sum";
  }else if(colNorm=="MedianNorm"){
    data<-t(apply(data, 2, MedianNorm));
    rownm<-"Normalization to sample median";
  }else{
    # nothing to do
    rownm<-"N/A";
  }
  # norm.opt
  if(norm.opt=="log"){
    positiveVals <- data[data > 0]
    if(length(positiveVals) == 0){
      AddErrMsg("All values are non-positive; log normalization cannot proceed.")
      return(0)
    }
    min.pos <- max(min(positiveVals, na.rm=T)/10, 1e-6)
    numberOfNeg = sum(data<0, na.rm = TRUE) + 1; 
    totalNumber = length(data)
    if((numberOfNeg/totalNumber)>0.2){
      msg <- paste(msg, "Can't perform log2 normalization, over 20% of data are negative. Try a different method or maybe the data already normalized?", collapse=" ");
      msgSet$norm.msg <- msgSet$current.msg <- msg;
      saveSet(msgSet, "msgSet");
      return(0);
    }
    data[data<=0] <- min.pos;
    data <- log2(data);
    msg <- paste(msg, "Log2 transformation.", collapse=" ");
  }else if(norm.opt=="logMedian"){
    positiveVals <- data[data > 0]
    if(length(positiveVals) == 0){
      AddErrMsg("All values are non-positive; log normalization cannot proceed.")
      return(0)
    }
    min.pos <- max(min(positiveVals, na.rm=T)/10, 1e-6)
    numberOfNeg = sum(data<0, na.rm = TRUE) + 1;
    totalNumber = length(data)
    if((numberOfNeg/totalNumber)>0.2){
      msg <- paste(msg, "Can't perform log2 normalization, over 20% of data are negative. Try a different method or maybe the data already normalized?", collapse=" ");
      msgSet$norm.msg <- msgSet$current.msg <- msg;
      saveSet(msgSet, "msgSet");
      return(0);
    }
    data[data<=0] <- min.pos;
    data <- log2(data);
    col.medians <- apply(data, 2, median, na.rm=TRUE);
    global.median <- median(col.medians, na.rm=TRUE);
    data <- sweep(data, 2, col.medians - global.median);
    msg <- paste(msg, "Log2 transformation + median centering.", collapse=" ");
  }else if(norm.opt=="vsn"){
    require(limma);
    data <- tryCatch(
      normalizeVSN(data),
      error = function(e) {
        err.msg <- paste0(
          "VSN normalization failed: ",
          e$message,
          " Current matrix size: ",
          nrow(data), " features x ", ncol(data), " samples."
        )
        if (nrow(data) < 100) {
          err.msg <- paste0(
            err.msg,
            " VSN usually needs a larger number of features after filtering. ",
            "Try relaxing filtering thresholds or use log/quantile normalization instead."
          )
        }
        msgSet$norm.msg <- msgSet$current.msg <- err.msg;
        saveSet(msgSet, "msgSet");
        return(0)
      }
    );
    if (length(data) == 1 && identical(as.numeric(data), 0)) {
      return(0)
    }
    msg <- paste(msg, "VSN normalization.", collapse=" ");
  }else if(norm.opt=="quantile"){
    require('preprocessCore');
    data <- normalize.quantiles(as.matrix(data), copy=TRUE);
    msg <- paste(msg, "Quantile normalization.", collapse=" ");
  }else if(norm.opt=="combined"){
    require(limma);
    data <- tryCatch(
      normalizeVSN(data),
      error = function(e) {
        err.msg <- paste0(
          "Combined normalization failed at the VSN step: ",
          e$message,
          " Current matrix size: ",
          nrow(data), " features x ", ncol(data), " samples."
        )
        if (nrow(data) < 100) {
          err.msg <- paste0(
            err.msg,
            " VSN usually needs a larger number of features after filtering. ",
            "Try relaxing filtering thresholds or switch to quantile/log normalization."
          )
        }
        msgSet$norm.msg <- msgSet$current.msg <- err.msg;
        saveSet(msgSet, "msgSet");
        return(0)
      }
    );
    if (length(data) == 1 && identical(as.numeric(data), 0)) {
      return(0)
    }
    require('preprocessCore');
    data <- normalize.quantiles(as.matrix(data), copy=TRUE);
    msg <- paste(msg, "VSN followed by quantile normalization.", collapse=" ");
  }else if(norm.opt=="logcount"){ # for count data, do it in DE analysis, as it is dependent on design matrix
    require(edgeR);
    nf <- calcNormFactors(data, method = "none");
    y <- voom(data,plot=F,lib.size=colSums(data)*nf);
    data <- y$E; # copy per million
    msg <- paste(msg, "Limma based on log2-counts per million transformation.", collapse=" ");
  } else if(norm.opt=="RLE"){
    suppressMessages(require(edgeR))
    nf <- calcNormFactors(data,method="RLE");
    y <- voom(data,plot=F,lib.size=colSums(data)*nf);
    data <- y$E; # copy per million
    msg <- c(msg, paste("Performed RLE Normalization"));
  }else if(norm.opt=="TMM"){
    suppressMessages(require(edgeR))
    nf <- calcNormFactors(data,method="TMM");
    y <- voom(data,plot=F,lib.size=colSums(data)*nf);
    data <- y$E; # copy per million
    msg <- c(msg, paste("Performed TMM Normalization"));
  }else if(norm.opt=="clr"){
    data <- apply(data, 2, clr_transform);
    msg <- "Performed centered-log-ratio normalization.";
  }else if(norm.opt=='LogNorm'){
    min.val <- min(abs(data[data!=0]))/10;
    data<-apply(data, 2, LogNorm, min.val);
  }else if(norm.opt=='CrNorm'){
    norm.data <- abs(data)^(1/3);
    norm.data[data<0] <- - norm.data[data<0];
    data <- norm.data;
  }else if(norm.opt=='Rlr'){
    norm.data <- RLRNorm(data)
    msg <- paste(msg, "Performed Linear Regression Normalization.", collapse=" ");
  }else if(norm.opt=='Loess'){
    norm.data <- LoessNorm(data)
    msg <- paste(msg, "Performed Local Regression Normalization.", collapse=" ");
  }else if(norm.opt=='EigenMS'){
     msg <- paste(msg, "Performed EigenMS Normalization.", collapse=" ");
  }else if(norm.opt=='median'){
    data<- apply(data, 2, MedianNorm);
    msg <- paste(msg, "Normalization to sample median.", collapse=" ");
  }else if(norm.opt=='abundance_corr'){
    # Protein abundance correction for phosphoproteomics
    # Normalizes phosphosite intensity by dividing by protein abundance (log-ratio subtraction)
    paramSet <- readSet(paramSet, "paramSet");

    if (!is.null(paramSet$data.type) && paramSet$data.type == "phospho" &&
        !is.null(paramSet$has.protein.ref) && paramSet$has.protein.ref) {

      #msg("[NormalizeData] Applying protein abundance correction for phosphoproteomics...")
      protein_ref <- paramSet$protein.ref

      # Perform abundance correction
      corrected_data <- .correctPhosphoByProteinAbundance(data, protein_ref)

      if (!is.null(corrected_data)) {
        data <- corrected_data
        msg <- paste(msg, "Protein abundance correction applied (log-ratio subtraction).", collapse=" ");
        #msg("[NormalizeData] Protein abundance correction complete.")
      } else {
        msg <- paste(msg, "Warning: Protein abundance correction failed; keeping original data.", collapse=" ");
        #msg("[NormalizeData] Warning: Protein abundance correction failed.")
      }
    } else {
      msg <- paste(msg, "Warning: Protein abundance correction requires phospho data with protein reference.", collapse=" ");
      #msg("[NormalizeData] Warning: abundance_corr selected but requirements not met.")
    }
  }


  # scaling
  if(scaleNorm=='MeanCenter'){
    data<-apply(data, 1, MeanCenter);
    scalenm<-"Mean Centering";
  }else if(scaleNorm=='AutoNorm'){
    data<-apply(data, 1, AutoNorm);
    scalenm<-"Autoscaling";
  }else if(scaleNorm=='ParetoNorm'){
    data<-apply(data, 1, ParetoNorm);
    scalenm<-"Pareto Scaling";
  }else if(scaleNorm=='RangeNorm'){
    data<-apply(data, 1, RangeNorm);
    scalenm<-"Range Scaling";
  }else if(scaleNorm=="colsum"){
    data <- sweep(data, 2, colSums(data), FUN="/")
    data <- data*10000000;
    msg <- c(msg, paste("Performed total sum normalization."));
  }else if(scaleNorm=="upperquartile" || norm.opt == "upperquartile"){
    suppressMessages(require(edgeR))
    nf <- calcNormFactors(data,method="upperquartile");
    y <- voom(data,plot=F,lib.size=colSums(data)*nf);
    data <- y$E; # copy per million
    msg <- c(msg, paste("Performed upper quartile normalization"));
  }else if(scaleNorm=="CSS"){
    suppressMessages(require(metagenomeSeq))
    #biom and mothur data also has to be in class(matrix only not in phyloseq:otu_table)
    data1 <- as(data,"matrix");
    dataMR <- newMRexperiment(data1);
    data <- cumNorm(dataMR,p=cumNormStat(dataMR));
    data <- MRcounts(data,norm = T);
    msg <- c(msg, paste("Performed cumulative sum scaling normalization"));
  }else{
    scalenm<-"N/A";
  }
  
  if(scaleNorm %in% c('MeanCenter', 'AutoNorm', 'ParetoNorm', 'RangeNorm')){
    data <- t(data)
  }
  
  norm.msg <<- msg;

  rownames(data) <- row.nms;
  colnames(data) <- col.nms;

  msgSet$current.msg <- msg;
  saveSet(msgSet, "msgSet");
  return(data)
}

# Wrapper: peptide -> protein summarization using summarize_peptides()
# Expects a peptide-to-protein map saved as peptide_to_protein_map.qs (two columns: Peptide, Protein)
# Input: data.stat.qs (peptide-level matrix after normalization)
# Output: int.mat.qs + data.norm updated to protein-level for downstream annotation
SummarizeProteomicsData <- function(dataName = "",
                                    method = "median_polish",
                                    top_n = 3,
                                    min_peptides = 1) {
  dataSet <- readDataset(dataName)
  msgSet  <- readSet(msgSet, "msgSet")
  paramSet <- readSet(paramSet, "paramSet")
  msgSet$current.msg <- character(0)

  if (!(dataSet$type %in% c("prot", "peptide"))) {
    msgSet$current.msg <- "Peptide summarization skipped: dataset is not proteomics/peptide."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  if (!file.exists("data.stat.qs")) {
    msgSet$current.msg <- "Peptide summarization failed: normalized peptide matrix (data.stat.qs) not found."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Load peptide matrix
  pep.mat <- qs::qread("data.stat.qs")
  msgSet$current.msg <- paste0("Starting peptide summarization using method '", method,
                               "' with top_n=", top_n, ", min_peptides=", min_peptides, ".")
  saveSet(msgSet, "msgSet")

  # Load peptide-to-protein map
  pep.map <- NULL
  if (file.exists("peptide_to_protein_map.qs")) {
    pep.map <- qs::qread("peptide_to_protein_map.qs")
  } else if (!is.null(dataSet$prot.map)) {
    pep.map <- dataSet$prot.map
  } else if (!is.null(dataSet$pep.map)) {
    pep.map <- dataSet$pep.map
  }

  if (is.null(pep.map)) {
    msgSet$current.msg <- c(msgSet$current.msg, "Peptide summarization failed: peptide-to-protein map not found.")
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Run summarization
  prot.mat <- try(summarize_peptides(pep.mat,
                                     pep.map,
                                     method = method,
                                     top_n_count = top_n,
                                     min_peptides = min_peptides),
                  silent = TRUE)

  if (inherits(prot.mat, "try-error")) {
    msgSet$current.msg <- paste0("Peptide summarization error: ", prot.mat)
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Clean protein IDs to UniProt accessions (e.g., sp|Q96EU7| -> Q96EU7, drop _HUMAN suffix)
  clean_uniprot <- function(x) {
    if (is.na(x) || !nzchar(x)) return(NA_character_)
    first <- strsplit(as.character(x), ";")[[1]][1]
    first <- trimws(first)
    if (grepl("\\|", first)) {
      parts <- strsplit(first, "\\|")[[1]]
      if (length(parts) >= 2 && nzchar(parts[2])) {
        return(parts[2])
      }
    }
    # Fallback: strip species suffix/isoform
    first <- sub("_.*", "", first)
    first <- sub("-\\d+$", "", first)
    return(first)
  }
  prot.ids <- vapply(rownames(prot.mat), clean_uniprot, character(1))
  rownames(prot.mat) <- prot.ids

  # Persist protein-level matrix for downstream annotation
  qs::qsave(prot.mat, "int.mat.qs")
  qs::qsave(prot.mat, "orig.data.anot.qs")
  qs::qsave(prot.mat, "norm.input.anot.qs")
  qs::qsave(prot.mat, "data.missed.qs")
  qs::qsave(prot.mat, "data.raw.qs")

  # Save peptide-level data for peptide-level DE analysis (shadow save for Arrow)
  shadow_save(pep.mat, "peptide_level_data.qs")
  msgSet$current.msg <- paste0("Saved peptide-level data: ", nrow(pep.mat), " peptides x ", ncol(pep.mat), " samples.")

  dataSet$data.norm <- prot.mat
  msgSet$current.msg <- paste0("Peptide summarization (", method, ") completed: ",
                               nrow(prot.mat), " proteins x ", ncol(prot.mat), " samples.")

  # Perform protein-level annotation using UniProt IDs derived above
  dataSet <- PerformDataAnnotInternal(dataSet,
                                      dataName,
                                     org = paramSet$data.org,
                                      dataType = "prot",
                                      idType = "uniprot",
                                    lvlOpt = "mean")
  msgSet$current.msg <- c(msgSet$current.msg,
                          "Protein-level annotation completed (input IDs parsed to UniProt).")

  saveSet(msgSet, "msgSet")

  return(RegisterData(dataSet));
}

# Build MSstats long-format input directly from an expression matrix and metadata
.buildMSstatsInputFromMatrix <- function(int.mat, metadata) {
  if (is.null(int.mat) || is.null(metadata)) {
    return(NULL);
  }
  runs <- colnames(int.mat);
  meta <- as.data.frame(metadata, stringsAsFactors = FALSE);
  shared <- intersect(runs, rownames(meta));
  if (length(shared) == 0) {
    return(NULL);
  }
  runs <- shared;
  meta <- meta[runs, , drop = FALSE];
  if (!"Condition" %in% colnames(meta) && ncol(meta) >= 1) {
    meta$Condition <- meta[, 1];
  }
  if (!"BioReplicate" %in% colnames(meta)) {
    meta$BioReplicate <- rownames(meta);
  }
  if (is.null(meta$Condition) || is.null(meta$BioReplicate)) {
    return(NULL);
  }
  data.frame(
    ProteinName      = rep(rownames(int.mat), times = length(runs)),
    PeptideSequence  = rep(rownames(int.mat), times = length(runs)),
    PrecursorCharge  = NA_integer_,
    FragmentIon      = NA_character_,
    ProductCharge    = NA_integer_,
    IsotopeLabelType = "L",
    Condition        = rep(meta$Condition, each = nrow(int.mat)),
    BioReplicate     = rep(meta$BioReplicate, each = nrow(int.mat)),
    Run              = rep(runs, each = nrow(int.mat)),
    Intensity        = as.numeric(int.mat[, runs, drop = FALSE]),
    stringsAsFactors = FALSE
  )
}

# Normalize proteomics data using MSstats::dataProcess if available
.normalizeWithMSstats <- function(msstats.path, meta.info) {
  msgSet <- readSet(msgSet, "msgSet");
  if (!file.exists(msstats.path)) {
    return(list(mat = NULL, proc = NULL));
  }
  if (!requireNamespace("MSstats", quietly = TRUE)) {
    msgSet$current.msg <- c(msgSet$current.msg, "MSstats package not available for proteomics normalization.");
    saveSet(msgSet, "msgSet");
    return(list(mat = NULL, proc = NULL));
  }
  msin <- try(qs::qread(msstats.path), silent = TRUE);
  if (inherits(msin, "try-error") || is.null(msin)) {
    msgSet$current.msg <- c(msgSet$current.msg, "Failed to read msstats_input.qs.");
    saveSet(msgSet, "msgSet");
    return(list(mat = NULL, proc = NULL));
  }
  required.cols <- c("ProteinName", "Run", "BioReplicate", "Condition", "Intensity");
  if (!all(required.cols %in% colnames(msin))) {
    msgSet$current.msg <- c(msgSet$current.msg, "MSstats input missing required columns; skipping MSstats normalization.");
    saveSet(msgSet, "msgSet");
    return(list(mat = NULL, proc = NULL));
  }
  proc <- try(MSstats::dataProcess(msin,
                                  normalization = "equalizeMedians",
                                  summaryMethod = "TMP",
                                  censoredInt = "NA",
                                  MBimpute = TRUE,
                                  featureSubset = "top3"), silent = TRUE);
  if (inherits(proc, "try-error")) {
    msgSet$current.msg <- c(msgSet$current.msg, "MSstats dataProcess error; skipping MSstats normalization.");
    saveSet(msgSet, "msgSet");
    return(list(mat = NULL, proc = NULL));
  }
  prof <- proc$ProteinLevelData;
  if (is.null(prof)) {
    msgSet$current.msg <- c(msgSet$current.msg, "MSstats returned empty ProteinLevelData.");
    saveSet(msgSet, "msgSet");
    return(list(mat = NULL, proc = proc));
  }
  value.col <- if ("LogIntensities" %in% colnames(prof)) "LogIntensities" else if ("ABUNDANCE" %in% colnames(prof)) "ABUNDANCE" else NULL;
  if (is.null(value.col)) {
    msgSet$current.msg <- c(msgSet$current.msg, "MSstats output missing LogIntensities/ABUNDANCE.");
    saveSet(msgSet, "msgSet");
    return(list(mat = NULL, proc = proc));
  }
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    msgSet$current.msg <- c(msgSet$current.msg, "reshape2 package not available to reshape MSstats output.");
    saveSet(msgSet, "msgSet");
    return(list(mat = NULL, proc = proc));
  }
  wide <- reshape2::dcast(prof, Protein ~ RUN, value.var = value.col);
  if (!"Protein" %in% colnames(wide)) {
    msgSet$current.msg <- c(msgSet$current.msg, "MSstats output missing Protein column.");
    saveSet(msgSet, "msgSet");
    return(list(mat = NULL, proc = proc));
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

# Normalize with an already-built MSstats long data frame
.normalizeWithMSstatsInput <- function(msin) {
  msgSet <- readSet(msgSet, "msgSet");
  if (is.null(msin)) {
    return(list(mat = NULL, proc = NULL));
  }
  if (!requireNamespace("MSstats", quietly = TRUE) || !requireNamespace("reshape2", quietly = TRUE)) {
    msgSet$current.msg <- c(msgSet$current.msg, "MSstats/reshape2 not available for proteomics normalization.");
    saveSet(msgSet, "msgSet");
    return(list(mat = NULL, proc = NULL));
  }
  required.cols <- c("ProteinName", "Run", "BioReplicate", "Condition", "Intensity");
  if (!all(required.cols %in% colnames(msin))) {
    msgSet$current.msg <- c(msgSet$current.msg, "MSstats input missing required columns; skipping MSstats normalization.");
    saveSet(msgSet, "msgSet");
    return(list(mat = NULL, proc = NULL));
  }
  proc <- try(MSstats::dataProcess(msin,
                                  normalization = "equalizeMedians",
                                  summaryMethod = "TMP",
                                  censoredInt = "NA",
                                  MBimpute = TRUE,
                                  featureSubset = "top3"), silent = TRUE);
  if (inherits(proc, "try-error")) {
    msgSet$current.msg <- c(msgSet$current.msg, "MSstats dataProcess error; skipping MSstats normalization.");
    saveSet(msgSet, "msgSet");
    return(list(mat = NULL, proc = NULL));
  }
  prof <- proc$ProteinLevelData;
  if (is.null(prof)) {
    msgSet$current.msg <- c(msgSet$current.msg, "MSstats returned empty ProteinLevelData.");
    saveSet(msgSet, "msgSet");
    return(list(mat = NULL, proc = proc));
  }
  value.col <- if ("LogIntensities" %in% colnames(prof)) "LogIntensities" else if ("ABUNDANCE" %in% colnames(prof)) "ABUNDANCE" else NULL;
  if (is.null(value.col)) {
    msgSet$current.msg <- c(msgSet$current.msg, "MSstats output missing LogIntensities/ABUNDANCE.");
    saveSet(msgSet, "msgSet");
    return(list(mat = NULL, proc = proc));
  }
  wide <- reshape2::dcast(prof, Protein ~ RUN, value.var = value.col);
  if (!"Protein" %in% colnames(wide)) {
    msgSet$current.msg <- c(msgSet$current.msg, "MSstats output missing Protein column.");
    saveSet(msgSet, "msgSet");
    return(list(mat = NULL, proc = proc));
  }
  mat <- as.matrix(wide[, -1, drop = FALSE]);
  rownames(mat) <- wide$Protein;
  storage.mode(mat) <- "numeric";
  return(list(mat = mat, proc = proc));
}


########
#
#Normalization internal methods
#
########

# based on phyloseq post: https://github.com/joey711/shiny-phyloseq/blob/master/panels/paneldoc/Transform.md
clr_transform <- function(x, base=2){
  x <- log((x / gm_mean(x)), base)
  x[!is.finite(x) | is.na(x)] <- 0.0
  return(x)
}


# generalize log, tolerant to 0 and negative values
LogNorm<-function(x, min.val){
  log10((x + sqrt(x^2 + min.val^2))/2)
}


SumNorm<-function(x){
  1000*x/sum(x, na.rm=T);
}

# normalize by median
MedianNorm<-function(x){
  x/median(x, na.rm=T);
}


# normalize to zero mean and unit variance
AutoNorm<-function(x){
  (x - mean(x))/sd(x, na.rm=T);
}

# normalize to zero mean but variance/SE
ParetoNorm<-function(x){
  (x - mean(x))/sqrt(sd(x, na.rm=T));
}

# normalize to zero mean but variance/SE
MeanCenter<-function(x){
  x - mean(x);
}

# normalize to zero mean but variance/SE
RangeNorm<-function(x){
  if(max(x) == min(x)){
    x;
  }else{
    (x - mean(x))/(max(x)-min(x));
  }
}
########### adapted from NormalyzerDE (https://github.com/ComputationalProteomics/NormalyzerDE)

RLRNorm <- function(data) {

  # NOTE: For further optimization, consider using matrixStats::rowMedians() for 10-20x speedup
  # Current: apply(data, 1, median) works but is slower
  # Optimized alternative (requires matrixStats): matrixStats::rowMedians(data, na.rm=T)
  sampleLog2Median <- apply(data, 1, median,na.rm=T)
  
  calculateRLMForCol <- function(colIndex, sampleLog2Median, data) {
    
    lrFit <- MASS::rlm(as.matrix(data[, colIndex])~sampleLog2Median, na.action=stats::na.exclude)
    coeffs <- lrFit$coefficients
    coefIntercept <- coeffs[1]
    coefSlope <- coeffs[2]
    globalFittedRLRCol <- (data[, colIndex] - coefIntercept) / coefSlope
    globalFittedRLRCol
  }
  
  globalFittedRLR <- vapply(
    seq_len(ncol(data)),
    calculateRLMForCol,
    rep(0, nrow(data)),
    sampleLog2Median=sampleLog2Median,
    data=data
  )
  
  colnames(globalFittedRLR) <- colnames(data)
  
  return(globalFittedRLR)
}

LoessNorm <- function(x, weights = NULL, span=0.7, iterations = 3){
  x <- as.matrix(x)
  n <- ncol(x)
    for (k in 1:iterations) {
      a <- rowMeans(x,na.rm=TRUE)
      for (i in 1:n){
        m <- x[,i] - a
        f <- limma::loessFit(m, a, weights=weights, span=span)$fitted
        x[,i] <- x[,i] - f
      }
    }
  
  return(x)
}


# Prepare MORlog: save expression matrix into dat.in.qs
.prepare.morlog <- function(expr) {
  di <- list(expr = expr)
  qs::qsave(di, "dat.in.qs")
  return(1L)
}


# Apply MORlog back to the active dataSet
.apply.morlog <- function(dataName) {
  dataSet <- readDataset(dataName)
  di <- qs::qread("dat.in.qs")

  dataSet$expr <- di$expr
  dataSet$norm <- di$norm

  dataSet$data.norm <- dataSet$norm
  fast.write(dataSet$data.norm, file = "data_normalized.csv")
  qs::qsave(dataSet$data.norm, file = "data.stat.qs")

  msgSet <- readSet(msgSet, "msgSet")
  sample_norm_msg <- msgSet$sample.norm.msg
  if (is.null(sample_norm_msg)) {
    sample_norm_msg <- character(0)
  }
  msgSet$current.msg <- c(sample_norm_msg, "[MORlog] Applied DESeq2 size-factor normalization and log2(x+1).")
  saveSet(msgSet, "msgSet")

  return(RegisterData(dataSet))
}


# Microservice entrypoint: expects dat.in.qs in working dir,
# runs DESeq2 size factor normalization + log2(x+1)
morlog_micro_run <- function(expr_field = "expr", norm_field = "norm") {
  requireNamespace("DESeq2", quietly = TRUE)

  di <- qs::qread("dat.in.qs")
  m  <- as.matrix(di[[expr_field]])

  # basic checks
  if (any(m < 0, na.rm = TRUE)) stop("MORlog expects non-negative counts")
  if (!is.integer(m)) m <- round(m)

  # minimal DESeq2 object
  cd  <- S4Vectors::DataFrame(row.names = colnames(m))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = m, colData = cd, design = ~ 1)
  dds <- DESeq2::estimateSizeFactors(dds)

  norm_counts <- DESeq2::counts(dds, normalized = TRUE)
  di[[norm_field]] <- log2(norm_counts + 1)

  qs::qsave(di, "dat.in.qs")
  return(1L)
}


#' Peptide to Protein Summarization
#'
#' @param peptide_matrix A numeric matrix or data frame of peptide intensities. 
#'                       Rows = Peptides, Cols = Samples.
#'                       ASSUMPTION: Data is already Log2 transformed (typical for normalization steps).
#' @param peptide_to_protein_map A data frame with two columns: 'Peptide' and 'Protein'.
#' @param method The summarization method: "median_polish", "top_n", "mean", "median", "sum".
#' @param top_n_count Integer. Number of peptides to use for "top_n" method.
#' @param min_peptides Integer. Minimum number of peptides required to keep a protein.
#'
#' @return A matrix of protein abundances (Rows = Proteins, Cols = Samples).
summarize_peptides <- function(peptide_matrix, 
                               peptide_to_protein_map, 
                               method = "median_polish", 
                               top_n_count = 3, 
                               min_peptides = 1) {
  
library(dplyr)
library(tidyr)
library(tibble)
  # 1. Prepare Data
  # Ensure peptide_matrix has row names as Peptide sequences
  if(is.null(rownames(peptide_matrix))) {
    stop("Input matrix must have Peptide Sequences as row names.")
  }
  
  # Convert to long format for easier manipulation with dplyr
  long_data <- as.data.frame(peptide_matrix) %>%
    rownames_to_column(var = "Peptide") %>%
    pivot_longer(cols = -Peptide, names_to = "Sample", values_to = "Intensity") %>%
    mutate(Intensity = suppressWarnings(as.numeric(as.character(Intensity))))
  
  # Merge with Protein Map
  # Clean map to ensure column names are standard
  colnames(peptide_to_protein_map) <- c("Peptide", "Protein")
  
  combined_data <- left_join(long_data, peptide_to_protein_map, by = "Peptide") %>%
    filter(!is.na(Protein)) # Remove peptides not mapped to a protein
  
  #msg("[summarize_peptides] peptides: ", nrow(peptide_to_protein_map), 
  #        " mapped rows: ", nrow(combined_data), 
  #        " unique proteins: ", length(unique(combined_data$Protein)))

  # 2. Filter by Minimum Peptides per Protein
  # Count unique peptides per protein
  prot_counts <- combined_data %>%
    group_by(Protein) %>%
    summarise(n_peptides = n_distinct(Peptide)) %>%
    filter(n_peptides >= min_peptides)
  
  filtered_data <- combined_data %>%
    filter(Protein %in% prot_counts$Protein)
  
  if(nrow(filtered_data) == 0) {
    stop("No proteins remained after filtering for minimum peptide count.")
  }
  
  # 3. Apply Summarization Method
  
  #msg(paste("Performing summarization using method:", method))
  
  if (method == "median_polish") {    # FASTER ALTERNATIVE: Median Sweep (Approximation of Tukey's)
    
    # 1. Normalize peptides by their own median across samples
    normalized_peptides <- filtered_data %>%
      group_by(Peptide) %>%
      mutate(PepMedian = median(Intensity, na.rm = TRUE)) %>%
      mutate(NormIntensity = Intensity - PepMedian) %>%
      ungroup()
      
    # 2. Summarize proteins by taking median of normalized peptides
    final_df <- normalized_peptides %>%
      group_by(Protein, Sample) %>%
      summarise(Value = median(NormIntensity, na.rm = TRUE) + median(PepMedian, na.rm=TRUE), .groups="drop") %>%
      pivot_wider(names_from = Sample, values_from = Value)
    
  } else if (method == "mean") {
    # --- Simple Mean (Averaging Log intensities = Geometric Mean of Linear) ---
    final_df <- filtered_data %>%
      group_by(Protein, Sample) %>%
      summarise(Value = mean(Intensity, na.rm = TRUE), .groups="drop") %>%
      pivot_wider(names_from = Sample, values_from = Value)
      
  } else if (method == "median") {
    # --- Simple Median ---
    final_df <- filtered_data %>%
      group_by(Protein, Sample) %>%
      summarise(Value = median(Intensity, na.rm = TRUE), .groups="drop") %>%
      pivot_wider(names_from = Sample, values_from = Value)
      
  } else if (method == "sum") {
    # --- Sum all peptides (Linear scale sum) ---
    final_df <- filtered_data %>%
      mutate(LinearIntensity = 2^Intensity) %>%
      group_by(Protein, Sample) %>%
      summarise(SumLinear = sum(LinearIntensity, na.rm = TRUE), .groups="drop") %>%
      mutate(Log2Intensity = log2(SumLinear)) %>%
      select(Protein, Sample, Log2Intensity) %>%
      pivot_wider(names_from = Sample, values_from = Log2Intensity)
  }
  
  # 4. Final Formatting
  # Convert back to matrix
  final_matrix <- final_df %>%
    column_to_rownames("Protein") %>%
    as.matrix()
  #msg("[summarize_peptides] final matrix dim: ", paste(dim(final_matrix), collapse="x"))
  
  # Replace -Inf (from log2(0)) with NA if necessary
  final_matrix[is.infinite(final_matrix)] <- NA
  
  return(final_matrix)
}

PlotProteinProfile <- function(dataName = "", protein_id = "", imgName = "prot_profile_plot.png", dpi = 96, type = "png") {
  #msg("[PlotProteinProfile] start dataName=", dataName,
  #        " protein_id=", protein_id, " imgName=", imgName,
  #        " dpi=", dpi, " type=", type, " wd=", getwd())
  if (protein_id == "") {
    #msg("[PlotProteinProfile] abort: protein_id is empty")
    return("No protein_id provided")
  }
  if (!file.exists("int.mat.qs")) {
    #msg("[PlotProteinProfile] abort: int.mat.qs missing")
    return("Protein matrix (int.mat.qs) not found")
  }
  if (!file.exists("data.stat.qs")) {
    #msg("[PlotProteinProfile] abort: data.stat.qs missing")
    return("Peptide matrix (data.stat.qs) not found")
  }
  peptide_matrix <- qs::qread("data.stat.qs")
  protein_matrix <- qs::qread("int.mat.qs")
  #msg("[PlotProteinProfile] peptide_matrix dim: ", paste(dim(peptide_matrix), collapse = "x"),
  #        " protein_matrix dim: ", paste(dim(protein_matrix), collapse = "x"))
  ds <- readDataset(dataName)
  if (file.exists("peptide_to_protein_map.qs")) {
    #msg("[PlotProteinProfile] using peptide_to_protein_map.qs")
    peptide_to_protein_map <- qs::qread("peptide_to_protein_map.qs")
  } else if (!is.null(ds$prot.map)) {
    #msg("[PlotProteinProfile] using ds$prot.map")
    peptide_to_protein_map <- ds$prot.map
  } else if (!is.null(ds$pep.map)) {
    #msg("[PlotProteinProfile] using ds$pep.map")
    peptide_to_protein_map <- ds$pep.map
  } else {
    #msg("[PlotProteinProfile] abort: no peptide/protein map available")
    return("Peptide-to-protein map not found")
  }

  # normalize IDs to UniProt accession (same as summarization)
  clean_uniprot <- function(x) {
    if (is.na(x) || !nzchar(x)) return(NA_character_)
    first <- strsplit(as.character(x), ";")[[1]][1]
    first <- trimws(first)
    if (grepl("\\|", first)) {
      parts <- strsplit(first, "\\|")[[1]]
      if (length(parts) >= 2 && nzchar(parts[2])) return(parts[2])
    }
    first <- sub("_.*", "", first)
    first <- sub("-\\d+$", "", first)
    return(first)
  }
  map_ids <- vapply(peptide_to_protein_map[[2]], clean_uniprot, character(1))
  peptide_to_protein_map[[2]] <- map_ids
  #msg("[PlotProteinProfile] mapping rows: ", nrow(peptide_to_protein_map),
  #        " unique proteins: ", length(unique(peptide_to_protein_map[[2]])),
  #        " NA proteins: ", sum(is.na(peptide_to_protein_map[[2]])))

  # 1. protein profile
  prot_vals <- protein_matrix[protein_id, , drop = FALSE]
  if (nrow(prot_vals) == 0) {
    #msg("[PlotProteinProfile] abort: protein_id not found in protein matrix. available examples: ",
    #        paste(utils::head(rownames(protein_matrix), 5), collapse = ", "))
    return("Protein not found in protein matrix")
  }
  prot_df <- data.frame(
    Sample = factor(colnames(prot_vals), levels = colnames(prot_vals)),
    Intensity = as.numeric(prot_vals),
    Type = "Protein (Summarized)",
    Group = "Protein"
  )

  # 2. peptide profiles
  peps <- peptide_to_protein_map$Peptide[peptide_to_protein_map[[2]] == protein_id]
  if (length(peps) > 0) {
    #msg("[PlotProteinProfile] peptides mapped to protein: ", length(peps),
    #        " example: ", paste(utils::head(peps, 5), collapse = ", "))
    pep_mat <- peptide_matrix[rownames(peptide_matrix) %in% peps, , drop = FALSE]
    pep_df <- reshape2::melt(pep_mat)
    colnames(pep_df) <- c("Peptide", "Sample", "Intensity")
    pep_df$Sample <- factor(pep_df$Sample, levels = colnames(peptide_matrix))
    pep_df$Type <- "Peptide"
    pep_df$Group <- pep_df$Peptide
  } else {
    #msg("[PlotProteinProfile] no peptides mapped to protein_id")
    pep_df <- data.frame()
  }

  # 3. plot
  fileNm <- paste(imgName, "dpi", dpi, ".", sep = "")
  finalImg <- paste0(fileNm, type, sep = "")
  #msg("[PlotProteinProfile] writing plot to ", finalImg)
  if (tolower(type) == "pdf") {
    grDevices::pdf(file = finalImg, width = 12, height = 6)
  } else {
    grDevices::png(filename = finalImg, width = 900, height = 450, res = dpi)
  }
  p <- ggplot()
  if (nrow(pep_df) > 0) {
    p <- p + geom_line(data = pep_df, aes(x = Sample, y = Intensity, group = Group),
                       color = "grey80", size = 0.8, alpha = 0.7)
  }
  p <- p + geom_line(data = prot_df, aes(x = Sample, y = Intensity, group = Group),
                     color = "#2196F3", size = 1.5) +
    geom_point(data = prot_df, aes(x = Sample, y = Intensity),
               color = "#1976D2", size = 3) +
    theme_minimal() +
    labs(title = paste("Profile:", protein_id),
         subtitle = paste("Based on", length(peps), "peptides"),
         y = "Log2 Intensity", x = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
  grDevices::dev.off()
  #msg("[PlotProteinProfile] finished plot generation: ", finalImg)
  return(finalImg)
}

# Build protein dropdown options as "id||label" (label = gene symbol if annotated, else id)
BuildProteinOptions <- function() {
  save.image("build.RData");
  if (!file.exists("int.mat.qs")) {
    #msg("[BuildProteinOptions] int.mat.qs not found")
    return(character(0))
  }
  ids <- rownames(qs::qread("int.mat.qs"))
  labels <- ids
  mapped_entrez <- rep(NA_character_, length(ids))

  if (file.exists("annotation.qs")) {
    anot <- qs::qread("annotation.qs")
    if (!is.null(anot)) {
      if (!is.null(names(anot))) {
        mapped_entrez <- unname(anot[ids])
      } else if (length(anot) == length(ids)) {
        names(anot) <- ids
        mapped_entrez <- unname(anot[ids])
      }
    }
  }

  mapped_symbols <- rep(NA_character_, length(ids))
  if (file.exists("symbol.map.qs")) {
    paramSet <- readSet(paramSet, "paramSet");
    sm <- tryCatch(readDataQs("symbol.map.qs", paramSet$anal.type, paramSet$dataName),
                   error = function(e) NULL)
    # print(head(sm));
    if (is.data.frame(sm)) {
      sym.map <- NULL
      if ("gene_id" %in% colnames(sm)) {
        if ("symbol" %in% colnames(sm)) {
          sym.map <- setNames(as.character(sm$symbol), sm$gene_id)
        } else if ("Symbol" %in% colnames(sm)) {
          sym.map <- setNames(as.character(sm$Symbol), sm$gene_id)
        } else if ("name" %in% colnames(sm)) {
          sym.map <- setNames(as.character(sm$name), sm$gene_id)
        } else if ("Name" %in% colnames(sm)) {
          sym.map <- setNames(as.character(sm$Name), sm$gene_id)
        } else if ("accession" %in% colnames(sm)) {
          # Handle table with accession = gene symbol, gene_id = entrez
          sym.map <- setNames(as.character(sm$accession), sm$gene_id)
        }
      } else if (all(c("geneID", "Symbol") %in% colnames(sm))) {
        sym.map <- setNames(as.character(sm$Symbol), sm$geneID)
      }

      if (!is.null(sym.map)) {
        #msg("[BuildProteinOptions] symbol.map entries: ", length(sym.map))
        mapped_symbols <- sym.map[mapped_entrez]
      } else {
        #msg("[BuildProteinOptions] symbol.map.qs missing expected columns; using IDs only")
      }
    } else {
      #msg("[BuildProteinOptions] symbol.map.qs not a data.frame; using IDs only")
    }
  }

  # Prefer symbol; fall back to Entrez; else ID
  labels <- ifelse(!is.na(mapped_symbols) & mapped_symbols != "",
                   paste0(mapped_symbols, " (", ids, ")"),
                   labels)
  labels <- ifelse(labels == ids & !is.na(mapped_entrez) & mapped_entrez != "",
                   paste0(mapped_entrez, " (", ids, ")"),
                   labels)

  res <- paste(ids, labels, sep = "||")
  #msg("[BuildProteinOptions] returning ", length(res), " options")
  return(res)
}

#' Correct phosphosite intensities by protein abundance
#'
#' @param phospho_data Matrix of phosphosite intensities (log2-transformed), sites × samples
#' @param protein_ref Matrix of protein abundances (log2-transformed), proteins × samples
#' @return Corrected phosphosite matrix with protein abundance subtracted
.correctPhosphoByProteinAbundance <- function(phospho_data, protein_ref) {

  if (is.null(phospho_data) || is.null(protein_ref)) {
    #msg("[.correctPhosphoByProteinAbundance] Error: NULL input data")
    return(NULL)
  }

  if (nrow(phospho_data) == 0 || ncol(phospho_data) == 0) {
    #msg("[.correctPhosphoByProteinAbundance] Error: Empty phospho data")
    return(NULL)
  }

  # Extract protein IDs from phosphosite row names
  # Phosphosite format: "P12345_S_123" or "ProteinID_Residue_Position"
  site_ids <- rownames(phospho_data)

  # Extract protein ID (everything before the first underscore)
  # Handle cases where protein ID may contain multiple components (e.g., "sp|P12345|PROT_HUMAN")
  protein_ids <- sapply(strsplit(site_ids, "_"), function(x) x[1])

  #msg("[.correctPhosphoByProteinAbundance] Processing ", length(unique(protein_ids)), " unique proteins from ",
  #        nrow(phospho_data), " phosphosites")

  # Match samples between phospho and protein data
  common_samples <- intersect(colnames(phospho_data), colnames(protein_ref))
  if (length(common_samples) == 0) {
    #msg("[.correctPhosphoByProteinAbundance] Error: No common samples between phospho and protein data")
    return(NULL)
  }

  if (length(common_samples) < ncol(phospho_data)) {
    #msg("[.correctPhosphoByProteinAbundance] Warning: Only ", length(common_samples), " of ",
    #        ncol(phospho_data), " samples have protein reference data")
  }

  # Subset to common samples
  phospho_data_sub <- phospho_data[, common_samples, drop = FALSE]
  protein_ref_sub <- protein_ref[, common_samples, drop = FALSE]

  # Create corrected matrix
  corrected_data <- matrix(NA, nrow = nrow(phospho_data_sub), ncol = ncol(phospho_data_sub))
  rownames(corrected_data) <- rownames(phospho_data_sub)
  colnames(corrected_data) <- colnames(phospho_data_sub)

  # Match phosphosites to proteins and perform correction
  protein_ref_ids <- rownames(protein_ref_sub)
  sites_corrected <- 0
  sites_no_protein <- 0
  values_corrected <- 0
  values_kept_original <- 0

  for (i in 1:nrow(phospho_data_sub)) {
    site_protein_id <- protein_ids[i]

    # Try exact match first
    protein_match_idx <- which(protein_ref_ids == site_protein_id)

    # If no exact match, try partial match (for UniProt IDs with isoforms)
    if (length(protein_match_idx) == 0) {
      # Try matching the first part (main protein ID without isoform)
      main_protein_id <- strsplit(site_protein_id, "-")[[1]][1]
      protein_match_idx <- which(grepl(paste0("^", main_protein_id), protein_ref_ids))
    }

    if (length(protein_match_idx) > 0) {
      # Use first match if multiple matches found
      protein_abundances <- protein_ref_sub[protein_match_idx[1], ]

      # Perform log-ratio subtraction ONLY when BOTH values are present
      # If protein is NA, keep original phosphosite value (cannot correct)
      for (j in 1:ncol(phospho_data_sub)) {
        if (!is.na(phospho_data_sub[i, j]) && !is.na(protein_abundances[j])) {
          # Both present: perform correction
          corrected_data[i, j] <- phospho_data_sub[i, j] - protein_abundances[j]
          values_corrected <- values_corrected + 1
        } else {
          # Either phospho or protein is NA: keep original phospho value
          corrected_data[i, j] <- phospho_data_sub[i, j]
          if (!is.na(phospho_data_sub[i, j])) {
            values_kept_original <- values_kept_original + 1
          }
        }
      }

      sites_corrected <- sites_corrected + 1
    } else {
      # No protein reference found - keep original phosphosite intensity
      corrected_data[i, ] <- phospho_data_sub[i, ]
      sites_no_protein <- sites_no_protein + 1
    }
  }

  total_values <- nrow(phospho_data_sub) * ncol(phospho_data_sub)

  #msg("[.correctPhosphoByProteinAbundance] Correction summary:")
  #msg("  - Sites with protein match: ", sites_corrected)
  #msg("  - Sites without protein match: ", sites_no_protein)
  #msg("  - Values corrected (both phospho & protein present): ", values_corrected)
  #msg("  - Values kept original (protein missing): ", values_kept_original)
  #msg("  - Total values processed: ", total_values)
  #msg("  - Correction rate: ", round(100 * values_corrected / total_values, 1), "%")

  # For samples not in protein reference, keep original phospho values
  if (ncol(phospho_data) > length(common_samples)) {
    missing_samples <- setdiff(colnames(phospho_data), common_samples)
    #msg("[.correctPhosphoByProteinAbundance] Warning: ", length(missing_samples),
    #        " samples not corrected (no protein reference): ", paste(missing_samples, collapse = ", "))

    # Create full matrix with all original samples
    full_corrected <- phospho_data
    full_corrected[, common_samples] <- corrected_data
    return(full_corrected)
  }

  return(corrected_data)
}

##################################################
## Format-Specific Filtering Functions
## Applied during normalization step
##################################################

#' Apply MaxQuant-specific filtering
#' @param data Data matrix to filter
#' @param removeContaminants Remove contaminant proteins marked with "+"
#' @param removeDecoys Remove decoy/reverse proteins marked with "+"
#' @param minPeptides Minimum number of peptides required per protein
#' @return Filtered data matrix
ApplyMaxQuantFiltering <- function(data, removeContaminants = TRUE, removeDecoys = TRUE, removeOnlyBySite = TRUE, minPeptides = 2) {
  msgSet <- readSet(msgSet, "msgSet")
  paramSet <- readSet(paramSet, "paramSet")

  # Initialize statistics
  stats <- list(
    n_total = nrow(data),
    n_contaminants = 0,
    n_decoys = 0,
    n_only_by_site = 0,
    n_min_peptides = 0,
    n_final = nrow(data)
  )

  # Check if we have the original MaxQuant metadata
  if (!file.exists("maxquant_metadata.qs")) {
    msg <- "[MaxQuant Filter] No MaxQuant metadata found - skipping format-specific filtering"
    msgSet$current.msg <- c(msgSet$current.msg, msg)
    saveSet(msgSet, "msgSet")
    return(list(data = data, stats = stats))
  }

  mq_meta <- qs::qread("maxquant_metadata.qs")

  # Filter contaminants
  if (removeContaminants && "Potential.contaminant" %in% names(mq_meta)) {
    contaminant_ids <- mq_meta$Protein.IDs[mq_meta$Potential.contaminant == "+"]
    keep_rows <- !(rownames(data) %in% contaminant_ids)
    prev_count <- nrow(data)
    data <- data[keep_rows, , drop = FALSE]
    stats$n_contaminants <- prev_count - nrow(data)
    if (stats$n_contaminants > 0) {
      msg <- paste0("[MaxQuant Filter] Removed ", stats$n_contaminants, " contaminant proteins")
      msgSet$current.msg <- c(msgSet$current.msg, msg)
    }
  }

  # Filter decoys/reverse
  if (removeDecoys && "Reverse" %in% names(mq_meta)) {
    decoy_ids <- mq_meta$Protein.IDs[mq_meta$Reverse == "+"]
    keep_rows <- !(rownames(data) %in% decoy_ids)
    prev_count <- nrow(data)
    data <- data[keep_rows, , drop = FALSE]
    stats$n_decoys <- prev_count - nrow(data)
    if (stats$n_decoys > 0) {
      msg <- paste0("[MaxQuant Filter] Removed ", stats$n_decoys, " decoy/reverse proteins")
      msgSet$current.msg <- c(msgSet$current.msg, msg)
    }
  }

  # Filter only-by-site proteins
  if (removeOnlyBySite && "Only.identified.by.site" %in% names(mq_meta)) {
    site_ids <- mq_meta$Protein.IDs[mq_meta$Only.identified.by.site == "+"]
    keep_rows <- !(rownames(data) %in% site_ids)
    prev_count <- nrow(data)
    data <- data[keep_rows, , drop = FALSE]
    stats$n_only_by_site <- prev_count - nrow(data)
    if (stats$n_only_by_site > 0) {
      msg <- paste0("[MaxQuant Filter] Removed ", stats$n_only_by_site, " proteins identified only by site")
      msgSet$current.msg <- c(msgSet$current.msg, msg)
    }
  }

  # Filter by minimum peptides
  if (!is.na(minPeptides) && minPeptides > 1) {
    pep_col <- NULL
    for (col in c("Peptides", "Razor...unique.peptides", "Razor.unique.peptides")) {
      if (col %in% names(mq_meta)) {
        pep_col <- col
        break
      }
    }

    if (!is.null(pep_col)) {
      pep_counts <- mq_meta[[pep_col]]
      names(pep_counts) <- mq_meta$Protein.IDs

      # Match peptide counts to current data rows
      matched_counts <- pep_counts[rownames(data)]
      keep_rows <- !is.na(matched_counts) & matched_counts >= minPeptides
      prev_count <- nrow(data)
      data <- data[keep_rows, , drop = FALSE]
      stats$n_min_peptides <- prev_count - nrow(data)
      if (stats$n_min_peptides > 0) {
        msg <- paste0("[MaxQuant Filter] Removed ", stats$n_min_peptides, " proteins with < ", minPeptides, " peptides")
        msgSet$current.msg <- c(msgSet$current.msg, msg)
      }
    }
  }

  stats$n_final <- nrow(data)
  saveSet(msgSet, "msgSet")
  return(list(data = data, stats = stats))
}

#' Apply DIA-NN-specific filtering
#' @param data Data matrix to filter
#' @param qvalueFilter Filter proteins by Q-value < 0.01
#' @param pepFilter Filter proteins by PEP threshold
#' @param pepThreshold PEP threshold
#' @param minPeptides Minimum peptide/precursor evidence count
#' @return Filtered data matrix
ApplyDiannFiltering <- function(data, qvalueFilter = TRUE, pepFilter = FALSE, pepThreshold = 0.01, minPeptides = 0) {
  msgSet <- readSet(msgSet, "msgSet")

  # Initialize statistics
  stats <- list(
    n_total = nrow(data),
    n_qvalue = 0,
    n_pep = 0,
    n_min_peptides = 0,
    n_final = nrow(data),
    pep_threshold = pepThreshold
  )

  if (!qvalueFilter && !pepFilter && minPeptides <= 0) {
    return(list(data = data, stats = stats))
  }

  # Check if we have Q-value information
  if (!file.exists("diann_metadata.qs")) {
    msg <- "[DIA-NN Filter] No DIA-NN metadata found - skipping Q-value filtering"
    msgSet$current.msg <- c(msgSet$current.msg, msg)
    saveSet(msgSet, "msgSet")
    return(list(data = data, stats = stats))
  }

  diann_meta <- qs::qread("diann_metadata.qs")

  # Filter by Q-value
  if (qvalueFilter && ("Q.Value" %in% names(diann_meta) || "PG.Qvalue" %in% names(diann_meta))) {
    qval_col <- if ("Q.Value" %in% names(diann_meta)) "Q.Value" else "PG.Qvalue"
    qvalues <- diann_meta[[qval_col]]
    names(qvalues) <- diann_meta$Protein.IDs

    matched_qvals <- qvalues[rownames(data)]
    keep_rows <- !is.na(matched_qvals) & matched_qvals < 0.01
    prev_count <- nrow(data)
    data <- data[keep_rows, , drop = FALSE]
    stats$n_qvalue <- prev_count - nrow(data)

    if (stats$n_qvalue > 0) {
      msg <- paste0("[DIA-NN Filter] Removed ", stats$n_qvalue, " proteins with Q-value >= 0.01")
      msgSet$current.msg <- c(msgSet$current.msg, msg)
    }
  }

  if (pepFilter && "PEP" %in% names(diann_meta)) {
    peps <- as.numeric(diann_meta$PEP)
    names(peps) <- diann_meta$Protein.IDs

    matched_pep <- peps[rownames(data)]
    keep_rows <- !is.na(matched_pep) & matched_pep < pepThreshold
    prev_count <- nrow(data)
    data <- data[keep_rows, , drop = FALSE]
    stats$n_pep <- prev_count - nrow(data)

    if (stats$n_pep > 0) {
      msgSet$current.msg <- c(msgSet$current.msg,
                              paste0("[DIA-NN Filter] Removed ", stats$n_pep, " proteins with PEP >= ", pepThreshold))
    }
  }

  if (minPeptides > 0 && "Peptide.Count" %in% names(diann_meta)) {
    pep_counts <- as.numeric(diann_meta$Peptide.Count)
    names(pep_counts) <- diann_meta$Protein.IDs

    matched_counts <- pep_counts[rownames(data)]
    keep_rows <- !is.na(matched_counts) & matched_counts >= minPeptides
    prev_count <- nrow(data)
    data <- data[keep_rows, , drop = FALSE]
    stats$n_min_peptides <- prev_count - nrow(data)

    if (stats$n_min_peptides > 0) {
      msgSet$current.msg <- c(msgSet$current.msg,
                              paste0("[DIA-NN Filter] Removed ", stats$n_min_peptides, " proteins with peptide/precursor count < ", minPeptides))
    }
  }

  stats$n_final <- nrow(data)
  saveSet(msgSet, "msgSet")
  return(list(data = data, stats = stats))
}

#' Apply FragPipe-specific filtering
#' @param data Data matrix to filter
#' @param removeContaminants Remove contaminant proteins
#' @param minProb Minimum protein probability threshold
#' @param minPeptides Minimum combined total peptides
#' @return Filtered data matrix
ApplyFragpipeFiltering <- function(data, removeContaminants = TRUE, minProb = 0.99, minPeptides = 2) {
  msgSet <- readSet(msgSet, "msgSet")

  # Initialize statistics
  stats <- list(
    n_total = nrow(data),
    n_contaminants = 0,
    n_min_prob = 0,
    n_min_peptides = 0,
    n_final = nrow(data),
    min_prob = minProb
  )

  # Check if we have FragPipe metadata
  if (!file.exists("fragpipe_metadata.qs")) {
    msg <- "[FragPipe Filter] No FragPipe metadata found - skipping format-specific filtering"
    msgSet$current.msg <- c(msgSet$current.msg, msg)
    saveSet(msgSet, "msgSet")
    return(list(data = data, stats = stats))
  }

  fp_meta <- qs::qread("fragpipe_metadata.qs")

  # Filter contaminants
  if (removeContaminants && "Is.Contaminant" %in% names(fp_meta)) {
    contaminant_ids <- fp_meta$Protein[fp_meta$Is.Contaminant == "+"]
    keep_rows <- !(rownames(data) %in% contaminant_ids)
    prev_count <- nrow(data)
    data <- data[keep_rows, , drop = FALSE]
    stats$n_contaminants <- prev_count - nrow(data)
    if (stats$n_contaminants > 0) {
      msg <- paste0("[FragPipe Filter] Removed ", stats$n_contaminants, " contaminant proteins")
      msgSet$current.msg <- c(msgSet$current.msg, msg)
    }
  }

  # Filter by protein probability
  if (!is.na(minProb) && minProb > 0 && "Protein.Probability" %in% names(fp_meta)) {
    probs <- fp_meta$Protein.Probability
    names(probs) <- fp_meta$Protein

    matched_probs <- probs[rownames(data)]
    keep_rows <- !is.na(matched_probs) & matched_probs >= minProb

    prev_count <- nrow(data)
    data <- data[keep_rows, , drop = FALSE]
    stats$n_min_prob <- prev_count - nrow(data)

    if (stats$n_min_prob > 0) {
      msg <- paste0("[FragPipe Filter] Removed ", stats$n_min_prob, " proteins with probability < ", minProb)
      msgSet$current.msg <- c(msgSet$current.msg, msg)
    }
  }

  # Filter by minimum peptides
  if (!is.na(minPeptides) && minPeptides > 1 && "Combined.Total.Peptides" %in% names(fp_meta)) {
    pep_counts <- fp_meta$Combined.Total.Peptides
    names(pep_counts) <- fp_meta$Protein

    matched_counts <- pep_counts[rownames(data)]
    keep_rows <- !is.na(matched_counts) & matched_counts >= minPeptides
    prev_count <- nrow(data)
    data <- data[keep_rows, , drop = FALSE]
    stats$n_min_peptides <- prev_count - nrow(data)
    if (stats$n_min_peptides > 0) {
      msg <- paste0("[FragPipe Filter] Removed ", stats$n_min_peptides, " proteins with < ", minPeptides, " peptides")
      msgSet$current.msg <- c(msgSet$current.msg, msg)
    }
  }

  stats$n_final <- nrow(data)
  saveSet(msgSet, "msgSet")
  return(list(data = data, stats = stats))
}

#' Apply Spectronaut-specific filtering
#' @param data Data matrix to filter
#' @param removeContaminants Remove contaminant proteins
#' @return Filtered data matrix
ApplySpectronautFiltering <- function(data, removeContaminants = TRUE, qvalueFilter = TRUE, qvalueThreshold = 0.01, minPeptides = 0) {
  msgSet <- readSet(msgSet, "msgSet")

  # Initialize statistics
  stats <- list(
    n_total = nrow(data),
    n_contaminants = 0,
    n_qvalue = 0,
    n_min_peptides = 0,
    n_final = nrow(data),
    qvalue_threshold = qvalueThreshold
  )

  if (!removeContaminants && !qvalueFilter && minPeptides <= 0) {
    return(list(data = data, stats = stats))
  }

  # Check if we have Spectronaut metadata
  if (!file.exists("spectronaut_metadata.qs")) {
    msg <- "[Spectronaut Filter] No Spectronaut metadata found - skipping format-specific filtering"
    msgSet$current.msg <- c(msgSet$current.msg, msg)
    saveSet(msgSet, "msgSet")
    return(list(data = data, stats = stats))
  }

  spec_meta <- qs::qread("spectronaut_metadata.qs")

  # Filter contaminants (look for common contaminant markers)
  if (removeContaminants && "PG.IsContaminant" %in% names(spec_meta)) {
    contaminant_ids <- spec_meta$PG.ProteinGroups[spec_meta$PG.IsContaminant == TRUE]
    keep_rows <- !(rownames(data) %in% contaminant_ids)
    prev_count <- nrow(data)
    data <- data[keep_rows, , drop = FALSE]
    stats$n_contaminants <- prev_count - nrow(data)
    if (stats$n_contaminants > 0) {
      msg <- paste0("[Spectronaut Filter] Removed ", stats$n_contaminants, " contaminant proteins")
      msgSet$current.msg <- c(msgSet$current.msg, msg)
    }
  } else if (removeContaminants) {
    msgSet$current.msg <- c(msgSet$current.msg, "[Spectronaut Filter] Contaminant status not available in this Spectronaut file")
  }

  if (qvalueFilter && "PG.Qvalue" %in% names(spec_meta)) {
    qvals <- as.numeric(spec_meta$PG.Qvalue)
    names(qvals) <- spec_meta$PG.ProteinGroups
    matched_qvals <- qvals[rownames(data)]
    if (all(is.na(matched_qvals))) {
      msgSet$current.msg <- c(msgSet$current.msg, "[Spectronaut Filter] Protein q-values are unavailable for the current proteins")
    } else {
      keep_rows <- is.na(matched_qvals) | matched_qvals < qvalueThreshold
      prev_count <- nrow(data)
      data <- data[keep_rows, , drop = FALSE]
      stats$n_qvalue <- prev_count - nrow(data)
      if (stats$n_qvalue > 0) {
        msg <- paste0("[Spectronaut Filter] Removed ", stats$n_qvalue, " proteins with q-value >= ", qvalueThreshold)
        msgSet$current.msg <- c(msgSet$current.msg, msg)
      }
    }
  } else if (qvalueFilter) {
    msgSet$current.msg <- c(msgSet$current.msg, "[Spectronaut Filter] Protein q-values not available in this Spectronaut file")
  }

  if (minPeptides > 0 && "Peptide.Count" %in% names(spec_meta)) {
    pep_counts <- as.numeric(spec_meta$Peptide.Count)
    names(pep_counts) <- spec_meta$PG.ProteinGroups
    matched_counts <- pep_counts[rownames(data)]
    if (all(is.na(matched_counts))) {
      msgSet$current.msg <- c(msgSet$current.msg, "[Spectronaut Filter] Peptide evidence counts are unavailable for the current proteins")
    } else {
      keep_rows <- is.na(matched_counts) | matched_counts >= minPeptides
      prev_count <- nrow(data)
      data <- data[keep_rows, , drop = FALSE]
      stats$n_min_peptides <- prev_count - nrow(data)
      if (stats$n_min_peptides > 0) {
        msg <- paste0("[Spectronaut Filter] Removed ", stats$n_min_peptides, " proteins with peptide count < ", minPeptides)
        msgSet$current.msg <- c(msgSet$current.msg, msg)
      }
    }
  } else if (minPeptides > 0) {
    msgSet$current.msg <- c(msgSet$current.msg, "[Spectronaut Filter] Peptide evidence counts not available in this Spectronaut file")
  }

  stats$n_final <- nrow(data)
  saveSet(msgSet, "msgSet")
  return(list(data = data, stats = stats))
}

##############################################################
## Sample Normalization Functions (Row-wise Normalization) ##
## Based on MetaboAnalyst's sample normalization methods   ##
##############################################################

#' Apply sample normalization (row-wise normalization)
#' @param data Matrix with samples as columns and features as rows
#' @param method Sample normalization method: "none", "SumNorm", "MedianNorm", "ProbNorm", "SpecNorm"
#' @return Normalized data matrix
ApplySampleNormalization <- function(data, method = "none", sampleNormParam = "__manual__") {
  msgSet <- readSet(msgSet, "msgSet")
  paramSet <- readSet(paramSet, "paramSet")
  sample_norm_msg <- msgSet$sample.norm.msg
  if (is.null(sample_norm_msg)) {
    sample_norm_msg <- character(0)
  }

  if (method == "none" || is.null(method) || method == "NA") {
    sample_norm_msg <- c(sample_norm_msg, "[Sample Norm] No sample normalization applied")
    msgSet$sample.norm.msg <- sample_norm_msg
    msgSet$current.msg <- c(msgSet$current.msg, "[Sample Norm] No sample normalization applied")
    saveSet(msgSet, "msgSet")
    return(data)
  }

  # Apply the selected method
  if (method == "SumNorm") {
    data <- SampleNorm_Sum(data)
    sample_norm_msg <- c(sample_norm_msg, "[Sample Norm] Normalization to constant sum applied")
  } else if (method == "MedianNorm") {
    data <- SampleNorm_Median(data)
    sample_norm_msg <- c(sample_norm_msg, "[Sample Norm] Normalization to sample median applied")
  } else if (method == "ProbNorm") {
    data <- SampleNorm_PQN(data)
    sample_norm_msg <- c(sample_norm_msg, "[Sample Norm] Probabilistic Quotient Normalization (PQN) applied")
  } else if (method == "SpecNorm") {
    result <- SampleNorm_Specific(data, paramSet, sampleNormParam)
    data <- result$data
    sample_norm_msg <- c(sample_norm_msg, result$message)
  } else {
    sample_norm_msg <- c(sample_norm_msg, paste0("[Sample Norm] Unknown method: ", method, " - skipping"))
  }

  msgSet$sample.norm.msg <- sample_norm_msg
  msgSet$current.msg <- c(msgSet$current.msg, sample_norm_msg)
  saveSet(msgSet, "msgSet")
  return(data)
}

#' Normalization to constant sum
#' Each sample is divided by its sum and multiplied by a constant (1000)
#' @param data Matrix with samples as columns
#' @return Normalized matrix
SampleNorm_Sum <- function(data) {
  colSums <- apply(data, 2, sum, na.rm = TRUE)
  # Avoid division by zero
  colSums[!is.finite(colSums) | colSums == 0] <- 1
  # Normalize each sample to sum = 1000
  data_norm <- sweep(data, 2, colSums, FUN = "/") * 1000
  return(data_norm)
}

#' Normalization to sample median
#' Each sample is divided by its median value
#' @param data Matrix with samples as columns
#' @return Normalized matrix
SampleNorm_Median <- function(data) {
  colMedians <- apply(data, 2, median, na.rm = TRUE)
  # Avoid division by zero
  colMedians[!is.finite(colMedians) | colMedians == 0] <- 1
  # Normalize each sample by its median
  data_norm <- sweep(data, 2, colMedians, FUN = "/")
  return(data_norm)
}

#' Probabilistic Quotient Normalization (PQN)
#' 1. Perform an integral normalization (sum normalization)
#' 2. Choose/calculate a reference spectrum (median across samples)
#' 3. Calculate the quotients of all variables of interest relative to the reference
#' 4. Calculate the median of these quotients for each sample
#' 5. Divide all variables of each sample by this median quotient
#' @param data Matrix with samples as columns
#' @return Normalized matrix
SampleNorm_PQN <- function(data) {
  # Step 1: Integral normalization (sum normalization to same total)
  colSums <- apply(data, 2, sum, na.rm = TRUE)
  colSums[!is.finite(colSums) | colSums == 0] <- 1
  data_int <- sweep(data, 2, colSums, FUN = "/")

  # Step 2: Calculate reference spectrum (median across samples)
  ref_spectrum <- apply(data_int, 1, median, na.rm = TRUE)
  ref_spectrum[!is.finite(ref_spectrum) | ref_spectrum == 0] <- NA_real_

  # Step 3 & 4: Calculate quotients and their median for each sample
  quotients <- sweep(data_int, 1, ref_spectrum, FUN = "/")
  quotients[!is.finite(quotients)] <- NA_real_
  median_quotients <- apply(quotients, 2, median, na.rm = TRUE)

  # Avoid division by zero
  median_quotients[!is.finite(median_quotients) | median_quotients == 0] <- 1

  # Step 5: Normalize by median quotients
  data_norm <- sweep(data_int, 2, median_quotients, FUN = "/")

  return(data_norm)
}

#' Sample-specific normalization
#' Normalize each sample by a user-provided factor source
#' Factors come either from a manual vector or an explicitly selected metadata column
#' @param data Matrix with samples as columns
#' @param paramSet Parameter set containing metadata
#' @return List with normalized data and message
SampleNorm_Specific <- function(data, paramSet, sampleNormParam = "__manual__") {
  norm_vec <- NULL
  norm_col_name <- NULL

  if (!is.null(sampleNormParam) && sampleNormParam == "__manual__") {
    if (exists("norm.vec") && length(norm.vec) == ncol(data)) {
      norm_vec <- as.numeric(norm.vec)
      norm_col_name <- "manual factors"
    } else {
      message <- "[Sample Norm] Sample-specific normalization: Manual factors were not provided or did not match the sample count. Using default (no normalization)."
      return(list(data = data, message = message))
    }
  }
 
  if (is.null(norm_vec) && !is.null(paramSet$meta.info)) {
    meta <- paramSet$meta.info
    sample_names <- colnames(data)
    
    meta_sample_col <- if ("Sample" %in% colnames(meta)) "Sample" else NULL
    meta_sample_names <- if (!is.null(meta_sample_col)) meta[[meta_sample_col]] else rownames(meta)

    if (!is.null(sampleNormParam) && sampleNormParam != "" && sampleNormParam != "__manual__" && sampleNormParam %in% colnames(meta)) {
      norm_col_name <- sampleNormParam
      raw_norm_vals <- meta[[norm_col_name]][match(sample_names, meta_sample_names)]
      norm_vec <- suppressWarnings(as.numeric(raw_norm_vals))
      if (all(is.na(norm_vec)) && any(!is.na(raw_norm_vals))) {
        message <- paste0("[Sample Norm] Sample-specific normalization: Metadata column '", norm_col_name,
                          "' is not numeric and cannot be used as scaling factors. Using default (no normalization).")
        return(list(data = data, message = message))
      }
    }
  }
  
  if (is.null(norm_vec)) {
    message <- "[Sample Norm] Sample-specific normalization: No manual factors or selected metadata column were available. Using default (no normalization)."
    return(list(data = data, message = message))
  }
  
  # Check for missing or invalid values
  if (any(is.na(norm_vec))) {
    na_count <- sum(is.na(norm_vec))
    message <- paste0("[Sample Norm] Sample-specific normalization: ", na_count, 
                     " samples have missing normalization factors. Using default (1.0) for these samples.")
    norm_vec[is.na(norm_vec)] <- 1
  }
  
  # Check for zero values
  if (any(norm_vec == 0, na.rm = TRUE)) {
    zero_count <- sum(norm_vec == 0, na.rm = TRUE)
    message <- paste0("[Sample Norm] Sample-specific normalization: ", zero_count,
                     " samples have zero normalization factors. Using default (1.0) for these samples.")
    norm_vec[norm_vec == 0] <- 1
  }
  
  # Ensure positive values
  if (any(norm_vec < 0, na.rm = TRUE)) {
    message <- "[Sample Norm] Sample-specific normalization: Negative normalization factors detected. Using absolute values."
    norm_vec <- abs(norm_vec)
  }
  
  # Apply normalization: divide each sample (column) by its factor
  data_norm <- sweep(data, 2, norm_vec, FUN = "/")
  
  # Create success message with statistics
  message <- paste0("[Sample Norm] Sample-specific normalization applied using '", norm_col_name, 
                   "' column. Factors range: ", round(min(norm_vec, na.rm=TRUE), 2), 
                   " to ", round(max(norm_vec, na.rm=TRUE), 2))
  
  return(list(data = data_norm, message = message))
}

#' Apply phosphoproteomics format-specific filtering
#' Entry point for phospho-specific filtering - routes to appropriate format function
#' @param dataName Dataset name
#' @param dataFormat Format type (maxquant, diann, fragpipe)
#' @param formatOpts List of format-specific options
#' @return 1 on success, 0 on failure
ApplyPhosphoFormatSpecificFiltering <- function(dataName, dataFormat, formatOpts = list()) {
  msgSet <- readSet(msgSet, "msgSet")
  paramSet <- readSet(paramSet, "paramSet")
  dataSet <- readDataset(dataName)

  # Load current data
  if (file.exists("norm.input.anot.qs")) {
    data <- qs::qread("norm.input.anot.qs")
  } else if (file.exists("orig.data.anot.qs")) {
    data <- qs::qread("orig.data.anot.qs")
  } else {
    msgSet$current.msg <- c(msgSet$current.msg, "[Phospho Filter] No annotated data found - skipping format-specific filtering")
    saveSet(msgSet, "msgSet")
    return(1)
  }

  msg("[Phospho Filter] Applying ", dataFormat, " format-specific filtering to phospho data")

  # Route to appropriate format-specific function
  result <- NULL
  if (dataFormat == "maxquant") {
    removeContaminants <- if (!is.null(formatOpts$removeContaminants)) formatOpts$removeContaminants else TRUE
    removeDecoys <- if (!is.null(formatOpts$removeDecoys)) formatOpts$removeDecoys else TRUE
    minPeptides <- if (!is.null(formatOpts$minPeptides)) formatOpts$minPeptides else 0
    result <- ApplyMaxQuantPhosphoFiltering(data, removeContaminants, removeDecoys, minPeptides)
  } else if (dataFormat == "diann") {
    qvalueFilter <- if (!is.null(formatOpts$qvalueFilter)) formatOpts$qvalueFilter else TRUE
    minPeptides <- if (!is.null(formatOpts$minPeptides)) formatOpts$minPeptides else 0
    result <- ApplyDiannPhosphoFiltering(data, qvalueFilter, minPeptides)
  } else if (dataFormat == "fragpipe") {
    removeContaminants <- if (!is.null(formatOpts$removeContaminants)) formatOpts$removeContaminants else TRUE
    minProb <- if (!is.null(formatOpts$minProb)) formatOpts$minProb else 0.0
    minPeptides <- if (!is.null(formatOpts$minPeptides)) formatOpts$minPeptides else 0
    result <- ApplyFragpipePhosphoFiltering(data, removeContaminants, minProb, minPeptides)
  } else {
    msg <- paste0("[Phospho Filter] Unsupported format: ", dataFormat, " - skipping")
    msgSet$current.msg <- c(msgSet$current.msg, msg)
    saveSet(msgSet, "msgSet")
    return(1)
  }

  if (is.null(result) || is.null(result$data)) {
    msgSet$current.msg <- c(msgSet$current.msg, "[Phospho Filter] Format-specific filtering failed")
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Save filtered data back
  qs::qsave(result$data, "norm.input.anot.qs")

  # Generate summary message
  stats <- result$stats
  summary_msg <- paste0(
    "[Phospho Filter] Format-specific filtering complete: ",
    stats$n_total, " → ", stats$n_final, " phosphosites"
  )
  msgSet$current.msg <- c(msgSet$current.msg, summary_msg)
  saveSet(msgSet, "msgSet")

  return(1)
}

#' Apply MaxQuant phosphoproteomics filtering
#' @param data Data matrix to filter
#' @param removeContaminants Remove contaminant proteins
#' @param removeDecoys Remove decoy/reverse proteins
#' @param minPeptides Minimum number of peptides (not typically used for phospho)
#' @return List with filtered data and statistics
ApplyMaxQuantPhosphoFiltering <- function(data, removeContaminants = TRUE, removeDecoys = TRUE, minPeptides = 0) {
  msgSet <- readSet(msgSet, "msgSet")

  # Initialize statistics
  stats <- list(
    n_total = nrow(data),
    n_contaminants = 0,
    n_decoys = 0,
    n_min_peptides = 0,
    n_final = nrow(data)
  )

  # Check if we have MaxQuant phospho metadata
  if (!file.exists("maxquant_phospho_metadata.qs")) {
    msg <- "[MaxQuant Phospho Filter] No MaxQuant phospho metadata found - skipping format-specific filtering"
    msgSet$current.msg <- c(msgSet$current.msg, msg)
    saveSet(msgSet, "msgSet")
    return(list(data = data, stats = stats))
  }

  mq_meta <- qs::qread("maxquant_phospho_metadata.qs")

  # Filter contaminants
  if (removeContaminants && "Potential.contaminant" %in% names(mq_meta)) {
    contaminant_sites <- rownames(mq_meta)[which(mq_meta$Potential.contaminant == "+")]
    keep_rows <- !(rownames(data) %in% contaminant_sites)
    prev_count <- nrow(data)
    data <- data[keep_rows, , drop = FALSE]
    stats$n_contaminants <- prev_count - nrow(data)
    if (stats$n_contaminants > 0) {
      msg <- paste0("[MaxQuant Phospho Filter] Removed ", stats$n_contaminants, " contaminant phosphosites")
      msgSet$current.msg <- c(msgSet$current.msg, msg)
    }
  }

  # Filter decoys/reverse
  if (removeDecoys && "Reverse" %in% names(mq_meta)) {
    decoy_sites <- rownames(mq_meta)[which(mq_meta$Reverse == "+")]
    keep_rows <- !(rownames(data) %in% decoy_sites)
    prev_count <- nrow(data)
    data <- data[keep_rows, , drop = FALSE]
    stats$n_decoys <- prev_count - nrow(data)
    if (stats$n_decoys > 0) {
      msg <- paste0("[MaxQuant Phospho Filter] Removed ", stats$n_decoys, " decoy/reverse phosphosites")
      msgSet$current.msg <- c(msgSet$current.msg, msg)
    }
  }

  stats$n_final <- nrow(data)
  saveSet(msgSet, "msgSet")
  return(list(data = data, stats = stats))
}

#' Apply DIA-NN phosphoproteomics filtering
#' @param data Data matrix to filter
#' @param qvalueFilter Filter by Q-value < 0.01
#' @param minPeptides Minimum peptide/precursor evidence count
#' @return List with filtered data and statistics
ApplyDiannPhosphoFiltering <- function(data, qvalueFilter = TRUE, minPeptides = 0) {
  msgSet <- readSet(msgSet, "msgSet")

  # Initialize statistics
  stats <- list(
    n_total = nrow(data),
    n_qvalue = 0,
    n_min_peptides = 0,
    n_final = nrow(data)
  )

  # Check if we have DIA-NN phospho metadata
  if (!file.exists("diann_phospho_metadata.qs")) {
    msg <- "[DIA-NN Phospho Filter] No DIA-NN phospho metadata found - skipping format-specific filtering"
    msgSet$current.msg <- c(msgSet$current.msg, msg)
    saveSet(msgSet, "msgSet")
    return(list(data = data, stats = stats))
  }

  diann_meta <- qs::qread("diann_phospho_metadata.qs")

  # Filter by Q-value
  if (qvalueFilter && ("Q.Value" %in% names(diann_meta) || "PG.Qvalue" %in% names(diann_meta))) {
    qval_col <- if ("Q.Value" %in% names(diann_meta)) "Q.Value" else "PG.Qvalue"
    qvalues <- diann_meta[[qval_col]]
    names(qvalues) <- rownames(diann_meta)

    matched_qvals <- qvalues[rownames(data)]
    keep_rows <- !is.na(matched_qvals) & matched_qvals < 0.01
    prev_count <- nrow(data)
    data <- data[keep_rows, , drop = FALSE]
    stats$n_qvalue <- prev_count - nrow(data)
    if (stats$n_qvalue > 0) {
      msg <- paste0("[DIA-NN Phospho Filter] Removed ", stats$n_qvalue, " phosphosites with Q-value >= 0.01")
      msgSet$current.msg <- c(msgSet$current.msg, msg)
    }
  }

  stats$n_final <- nrow(data)
  saveSet(msgSet, "msgSet")
  return(list(data = data, stats = stats))
}

#' Apply FragPipe phosphoproteomics filtering
#' @param data Data matrix to filter
#' @param removeContaminants Remove contaminant proteins
#' @param minProb Minimum localization probability (0-1)
#' @param minPeptides Minimum number of peptides
#' @return List with filtered data and statistics
ApplyFragpipePhosphoFiltering <- function(data, removeContaminants = TRUE, minProb = 0.0, minPeptides = 0) {
  msgSet <- readSet(msgSet, "msgSet")

  # Initialize statistics
  stats <- list(
    n_total = nrow(data),
    n_contaminants = 0,
    n_min_prob = 0,
    n_min_peptides = 0,
    n_final = nrow(data),
    min_prob = minProb
  )

  # Check if we have FragPipe phospho metadata
  if (!file.exists("fragpipe_phospho_metadata.qs")) {
    msg <- "[FragPipe Phospho Filter] No FragPipe phospho metadata found - skipping format-specific filtering"
    msgSet$current.msg <- c(msgSet$current.msg, msg)
    saveSet(msgSet, "msgSet")
    return(list(data = data, stats = stats))
  }

  fp_meta <- qs::qread("fragpipe_phospho_metadata.qs")

  # Filter contaminants
  if (removeContaminants && "Is.Contaminant" %in% names(fp_meta)) {
    contaminant_sites <- rownames(fp_meta)[which(fp_meta$Is.Contaminant == "+")]
    keep_rows <- !(rownames(data) %in% contaminant_sites)
    prev_count <- nrow(data)
    data <- data[keep_rows, , drop = FALSE]
    stats$n_contaminants <- prev_count - nrow(data)
    if (stats$n_contaminants > 0) {
      msg <- paste0("[FragPipe Phospho Filter] Removed ", stats$n_contaminants, " contaminant phosphosites")
      msgSet$current.msg <- c(msgSet$current.msg, msg)
    }
  }

  stats$n_final <- nrow(data)
  saveSet(msgSet, "msgSet")
  return(list(data = data, stats = stats))
}
