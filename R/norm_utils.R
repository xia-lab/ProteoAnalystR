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

  # Read the peptide-level normalization input.  Priority:
  # 1. peptide.input.anot.qs — saved by SummarizeProteomicsData before it overwrites
  #    norm.input.anot.qs with protein-level data.  Ensures re-running normalization
  #    after summarization still starts from the original peptide matrix.
  # 2. norm.input.anot.qs   — the filtered peptide baseline (pre-summarization path).
  # 3. orig.data.anot.qs    — annotated baseline (fallback).
  # 4. data.anot.qs         — backward-compatibility fallback.
  if (file.exists("peptide.input.anot.qs")) {
    ds <- ov_qs_read("peptide.input.anot.qs")
  } else if(file.exists("norm.input.anot.qs")){
    ds <- ov_qs_read("norm.input.anot.qs");
  } else if(file.exists("orig.data.anot.qs")){
    ds <- ov_qs_read("orig.data.anot.qs");
  } else {
    ds <- ov_qs_read("data.anot.qs");
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
    AddErrMsg(error.msg);
    return(0);
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

    # DESeq2 quarantined — isolate in subprocess
    norm_counts <- rsclient_isolated_exec(
      func_body = function(input_data) {
        require(DESeq2)
        cd <- S4Vectors::DataFrame(row.names = colnames(input_data$m))
        dds <- DESeq2::DESeqDataSetFromMatrix(countData = input_data$m, colData = cd, design = ~ 1)
        dds <- DESeq2::estimateSizeFactors(dds)
        DESeq2::counts(dds, normalized = TRUE)
      },
      input_data = list(m = m),
      packages = c("DESeq2", "qs"),
      timeout = 120,
      output_type = "qs"
    )
    if (is.list(norm_counts) && isFALSE(norm_counts$success)) { AddErrMsg(norm_counts$message); return(0) }
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
    anot.id <- ov_qs_read("annotation.qs")
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
  ov_qs_save(data, file = "data.stat.qs")

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
    data <- ov_qs_read("orig.data.anot.qs")
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

  ov_qs_save(data, "norm.input.anot.qs")
  saveSet(msgSet, "msgSet")

  return(1L)
}

#' Require at least one proteotypic (unique) peptide per protein/peptide row.
#' Operates on norm.input.anot.qs, which must already exist.
#' Handles both protein-level and peptide-level input matrices.
ApplyProteotypicFilter <- function(dataName) {
  msgSet <- readSet(msgSet, "msgSet")

  skip <- function(reason) {
    msgSet$current.msg <<- c(msgSet$current.msg, paste0("[ProteotypicFilter] SKIP: ", reason))
    saveSet(msgSet, "msgSet")
    return(1L)
  }

  # After peptide summarization, norm.input.anot.qs holds protein-level data and must
  # not be filtered here.  peptide.input.anot.qs (written by SummarizeProteomicsData)
  # is the correct peptide-level target; PerformNormalization already prefers it.
  target_file <- if (file.exists("peptide.input.anot.qs")) "peptide.input.anot.qs" else "norm.input.anot.qs"

  if (!file.exists(target_file)) {
    return(skip("norm.input.anot.qs not found."))
  }

  data <- ov_qs_read(target_file)
  n_before <- nrow(data)
  rn <- rownames(data)
  message(sprintf("[ProteotypicFilter] Starting: %d rows, first row name: '%s'", n_before, rn[1]))
  message(sprintf("[ProteotypicFilter] Files present: diann_metadata=%s, maxquant_metadata=%s, peptide_map=%s",
      file.exists("diann_metadata.qs"), file.exists("maxquant_metadata.qs"), file.exists("peptide_to_protein_map.qs")))

  # --- DIA-NN pr_matrix: Proteotypic column is directly available ---
  if (file.exists("diann_metadata.qs")) {
    diann_meta <- ov_qs_read("diann_metadata.qs")
    message(sprintf("[ProteotypicFilter] diann_metadata cols: %s", paste(names(diann_meta), collapse = ", ")))
    if ("Proteotypic" %in% names(diann_meta)) {
      idx <- match(rn, diann_meta$Protein.IDs)
      matched <- diann_meta$Proteotypic[idx]
      n_matched <- sum(!is.na(matched))
      message(sprintf("[ProteotypicFilter] Proteotypic column: %d/%d row names matched metadata IDs", n_matched, n_before))
      if (n_matched > 0) {
        n_proteotypic <- sum(!is.na(matched) & matched == 1L)
        message(sprintf("[ProteotypicFilter] Of matched: %d proteotypic (==1), %d non-proteotypic", n_proteotypic, n_matched - n_proteotypic))
        keep <- !is.na(matched) & matched == 1L
        data <- data[keep, , drop = FALSE]
        n_removed <- n_before - nrow(data)
        if (nrow(data) == 0) {
          return(skip(sprintf("Filter would remove all %d rows using Proteotypic column.", n_before)))
        }
        msg <- sprintf("Proteotypic filter: Removed %d non-proteotypic precursors (%d → %d retained).", n_removed, n_before, nrow(data))
        msgSet$current.msg <- c(msgSet$current.msg, msg)
        msgSet$prefilter.msg <- paste(c(msgSet$prefilter.msg, msg), collapse = " ")
        ov_qs_save(data, target_file)
        saveSet(msgSet, "msgSet")
        return(1L)
      }
      message("[ProteotypicFilter] Proteotypic column present but 0 row names matched — falling through to other methods")
    } else {
      message("[ProteotypicFilter] diann_metadata.qs has no Proteotypic column")
    }
  }

  # --- MaxQuant protein-level: use Unique.peptides column directly ---
  if (file.exists("maxquant_metadata.qs")) {
    mq_meta <- ov_qs_read("maxquant_metadata.qs")
    message(sprintf("[ProteotypicFilter] maxquant_metadata cols: %s", paste(names(mq_meta), collapse = ", ")))
    if (!("Unique.peptides" %in% names(mq_meta))) {
      return(skip("Unique.peptides column not found in MaxQuant metadata."))
    }
    idx <- match(rn, mq_meta$Protein.IDs)
    matched <- mq_meta$Unique.peptides[idx]
    n_matched <- sum(!is.na(matched))
    message(sprintf("[ProteotypicFilter] MaxQuant: %d/%d row names matched Protein.IDs; first data row='%s', first meta ID='%s'",
        n_matched, n_before, rn[1], mq_meta$Protein.IDs[1]))
    if (n_matched == 0) {
      return(skip(sprintf("No row names matched MaxQuant Protein.IDs (first data row: '%s', first meta ID: '%s').", rn[1], mq_meta$Protein.IDs[1])))
    }
    n_with_unique <- sum(!is.na(matched) & matched >= 1)
    message(sprintf("[ProteotypicFilter] MaxQuant: %d proteins have Unique.peptides >= 1, %d would be removed", n_with_unique, n_matched - n_with_unique))
    keep <- !is.na(matched) & matched >= 1
    data <- data[keep, , drop = FALSE]

  # --- Peptide-to-protein map: handle both protein-level and peptide-level ---
  } else if (file.exists("peptide_to_protein_map.qs")) {
    prot_map <- ov_qs_read("peptide_to_protein_map.qs")
    message(sprintf("[ProteotypicFilter] prot_map: %d rows, cols: %s", nrow(prot_map), paste(names(prot_map), collapse = ", ")))
    message(sprintf("[ProteotypicFilter] prot_map first Peptide: '%s', first Protein: '%s'", prot_map$Peptide[1], prot_map$Protein[1]))
    is_dup <- duplicated(prot_map$Peptide) | duplicated(prot_map$Peptide, fromLast = TRUE)
    not_dup_idx       <- which(!is_dup)
    proteotypic_peps  <- prot_map$Peptide[not_dup_idx]
    proteotypic_prots <- unique(prot_map$Protein[not_dup_idx])
    message(sprintf("[ProteotypicFilter] prot_map: %d total peptides, %d proteotypic, %d proteotypic proteins",
        nrow(prot_map), length(proteotypic_peps), length(proteotypic_prots)))

    pep_overlap  <- mean(rn %in% prot_map$Peptide)
    prot_overlap <- mean(rn %in% prot_map$Protein)
    message(sprintf("[ProteotypicFilter] Row name overlap — peptide col: %.1f%%, protein col: %.1f%%",
        pep_overlap * 100, prot_overlap * 100))

    if (pep_overlap == 0 && prot_overlap == 0) {
      return(skip(sprintf("Row names match neither peptides nor proteins in map (first row: '%s').", rn[1])))
    }

    if (pep_overlap >= prot_overlap) {
      keep <- rn %in% proteotypic_peps
      message(sprintf("[ProteotypicFilter] Peptide-level mode: %d/%d rows are proteotypic peptides", sum(keep), n_before))
    } else {
      keep <- rn %in% proteotypic_prots
      message(sprintf("[ProteotypicFilter] Protein-level mode: %d/%d rows have a proteotypic peptide", sum(keep), n_before))
    }
    data <- data[keep, , drop = FALSE]

  } else {
    return(skip("No peptide-level information available (no maxquant_metadata.qs or peptide_to_protein_map.qs)."))
  }

  n_removed <- n_before - nrow(data)
  if (nrow(data) == 0) {
    return(skip(sprintf("Filter would remove all %d rows; something is wrong with ID matching.", n_before)))
  }

  message(sprintf("[ProteotypicFilter] Done: removed %d, retained %d", n_removed, nrow(data)))
  result_msg <- sprintf("Proteotypic filter: Removed %d entries without a unique peptide (%d → %d retained).", n_removed, n_before, nrow(data))
  msgSet$current.msg <- c(msgSet$current.msg, result_msg)
  msgSet$prefilter.msg <- paste(c(msgSet$prefilter.msg, result_msg), collapse = " ")

  ov_qs_save(data, target_file)
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
    data1 <- ov_qs_read("data.raw.qs");
    colnames(data1) <- colnames(dataSet$data.norm)
    anot.id <- ov_qs_read("annotation.qs");
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
        raw.data.anot <- ov_qs_read("norm.input.anot.qs")
        #msg("[Filter] Loading norm.input.anot.qs for proteomics")
      } else if (file.exists("orig.data.anot.qs")) {
        raw.data.anot <- ov_qs_read("orig.data.anot.qs")
        #msg("[Filter] Loading orig.data.anot.qs for proteomics baseline")
      } else if (file.exists("data.annotated.qs")) {
        raw.data.anot <- ov_qs_read("data.annotated.qs")
        #msg("[Filter] Loading data.annotated.qs for proteomics (fallback)")
      } else if (file.exists("data.anot.qs")) {
        raw.data.anot <- ov_qs_read("data.anot.qs")
        #msg("[Filter] Loading data.anot.qs for proteomics (fallback)")
      } else if (file.exists("data.missed.qs")) {
        # Fallback to data.missed.qs if it exists (from previous workflow)
        raw.data.anot <- ov_qs_read("data.missed.qs")
        #msg("[Filter] Loading data.missed.qs for proteomics (fallback)")
      } else {
        AddErrMsg("No annotated data file found for proteomics filtering"); return(0);
      }
    }else{
     if (file.exists("norm.input.anot.qs")) {
       raw.data.anot <- ov_qs_read("norm.input.anot.qs");
     } else {
       raw.data.anot <- ov_qs_read("orig.data.anot.qs");
     }
    }
    # Ensure annotated IDs applied if annotation vector available
    if (file.exists("annotation.qs")) {
      anot.id <- ov_qs_read("annotation.qs")
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
    ov_qs_save(data, file="data.stat.qs");
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
  }else if(norm.opt %in% c("logcount", "RLE", "TMM")){
    # calcNormFactors is from edgeR (quarantined) — isolate in subprocess; voom is limma (OK in Master)
    edger_method <- switch(norm.opt, logcount = "none", RLE = "RLE", TMM = "TMM")
    nf <- rsclient_isolated_exec(
      func_body = function(input_data) {
        require(edgeR)
        edgeR::calcNormFactors(input_data$data, method = input_data$method)
      },
      input_data = list(data = data, method = edger_method),
      packages = c("edgeR", "qs"),
      timeout = 120,
      output_type = "qs"
    )
    if (is.list(nf) && isFALSE(nf$success)) { AddErrMsg(nf$message); return(0) }
    y <- voom(data, plot=F, lib.size=colSums(data)*nf);
    data <- y$E; # copy per million
    if (norm.opt == "logcount") {
      msg <- paste(msg, "Limma based on log2-counts per million transformation.", collapse=" ")
    } else {
      msg <- c(msg, paste("Performed", norm.opt, "Normalization"))
    }
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
    # Normalizes phosphosite intensity by dividing by protein abundance (log-ratio subtraction).
    # The protein reference is already log2-transformed; ensure phospho data is also log2 first.
    paramSet <- readSet(paramSet, "paramSet");

    if (!is.null(paramSet$data.type) && paramSet$data.type == "phospho" &&
        !is.null(paramSet$has.protein.ref) && paramSet$has.protein.ref) {

      # Log2-transform phospho intensities if still on linear scale (values > 100 indicate linear)
      data_min <- min(data, na.rm = TRUE)
      data_max <- max(data, na.rm = TRUE)
      looks_linear <- data_min >= 0 && data_max > 100
      if (looks_linear) {
        data[data == 0] <- NA
        data <- log2(data)
        msg <- paste(msg, "Log2 transformation applied before protein abundance correction.", collapse=" ");
      }

      protein_ref <- paramSet$protein.ref
      corrected_data <- .correctPhosphoByProteinAbundance(data, protein_ref)

      if (!is.null(corrected_data)) {
        data <- corrected_data
        msg <- paste(msg, "Protein abundance correction applied (log-ratio subtraction).", collapse=" ");
      } else {
        msg <- paste(msg, "Warning: Protein abundance correction failed; keeping original data.", collapse=" ");
      }
    } else {
      msg <- paste(msg, "Warning: Protein abundance correction requires phospho data with protein reference.", collapse=" ");
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
    nf <- rsclient_isolated_exec(
      func_body = function(input_data) {
        require(edgeR)
        edgeR::calcNormFactors(input_data$data, method = "upperquartile")
      },
      input_data = list(data = data),
      packages = c("edgeR", "qs"),
      timeout = 120,
      output_type = "qs"
    )
    if (is.list(nf) && isFALSE(nf$success)) { AddErrMsg(nf$message); return(0) }
    y <- voom(data, plot=F, lib.size=colSums(data)*nf);
    data <- y$E; # copy per million
    msg <- c(msg, paste("Performed upper quartile normalization"));
  }else if(scaleNorm=="CSS"){
    # metagenomeSeq quarantined — isolate in subprocess
    data <- rsclient_isolated_exec(
      func_body = function(input_data) {
        require(metagenomeSeq)
        data1 <- as(input_data$data, "matrix")
        dataMR <- metagenomeSeq::newMRexperiment(data1)
        dataMR <- metagenomeSeq::cumNorm(dataMR, p = metagenomeSeq::cumNormStat(dataMR))
        metagenomeSeq::MRcounts(dataMR, norm = TRUE)
      },
      input_data = list(data = data),
      packages = c("metagenomeSeq", "qs"),
      timeout = 120,
      output_type = "qs"
    )
    if (is.list(data) && isFALSE(data$success)) { AddErrMsg(data$message); return(0) }
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
                                    min_peptides = 1,
                                    filter.unmapped = FALSE) {
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
  pep.mat <- ov_qs_read("data.stat.qs")
  msgSet$current.msg <- paste0("Starting peptide summarization using method '", method,
                               "' with top_n=", top_n, ", min_peptides=", min_peptides, ".")
  saveSet(msgSet, "msgSet")

  # Load peptide-to-protein map
  pep.map <- NULL
  if (file.exists("peptide_to_protein_map.qs")) {
    pep.map <- ov_qs_read("peptide_to_protein_map.qs")
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

  if (is.null(prot.mat) || !is.matrix(prot.mat) || nrow(prot.mat) == 0) {
    msgSet$current.msg <- paste0(
      "Peptide summarization failed: no proteins remained after applying the minimum peptide ",
      "count filter (", min_peptides, "). Try reducing 'Minimum Peptides per Protein'.")
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

  # Snapshot the pre-normalization peptide input before overwriting it with protein-level
  # data.  PerformNormalization reads norm.input.anot.qs as its starting point, so without
  # this copy re-running normalization after summarization would re-normalize protein-level
  # data instead of the original peptides.
  if (file.exists("norm.input.anot.qs")) {
    file.copy("norm.input.anot.qs", "peptide.input.anot.qs", overwrite = TRUE)
  }

  # Persist protein-level matrix for downstream annotation
  ov_qs_save(prot.mat, "int.mat.qs")
  ov_qs_save(prot.mat, "orig.data.anot.qs")
  ov_qs_save(prot.mat, "norm.input.anot.qs")
  ov_qs_save(prot.mat, "data.missed.qs")
  ov_qs_save(prot.mat, "data.raw.qs")

  # Save peptide-level data aligned to prot.mat column order.
  # pivot_wider inside summarize_peptides may reorder samples; the design will be built
  # from prot.mat's column order, so pep.mat must match or lmFit assigns wrong groups.
  pep.mat.save <- pep.mat[, colnames(prot.mat)[colnames(prot.mat) %in% colnames(pep.mat)], drop = FALSE]
  shadow_save(pep.mat.save, "peptide_level_data.qs")
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

  if (isTRUE(filter.unmapped)) {
    sym.map <- tryCatch(
      readDataQs("symbol.map.qs", paramSet$anal.type, dataName),
      error = function(e) NULL
    )
    if (!is.null(sym.map) && is.data.frame(sym.map) && nrow(sym.map) > 0 &&
        "uniprot" %in% colnames(sym.map) && "symbol" %in% colnames(sym.map)) {
      # symbol == uniprot means no real gene symbol was found (UniProt ID used as fallback)
      real.sym.inx <- !is.na(sym.map$symbol) & nzchar(sym.map$symbol) &
                      sym.map$symbol != sym.map$uniprot
      ids.with.sym <- unique(sym.map$uniprot[real.sym.inx])
      mat <- dataSet$data.norm
      keep <- intersect(rownames(mat), ids.with.sym)
      n.before <- nrow(mat)
      if (length(keep) > 0 && length(keep) < n.before) {
        mat <- mat[keep, , drop = FALSE]
        dataSet$data.norm <- mat
        ov_qs_save(mat, "int.mat.qs")
        ov_qs_save(mat, "orig.data.anot.qs")
        ov_qs_save(mat, "norm.input.anot.qs")
        ov_qs_save(mat, "data.missed.qs")
        ov_qs_save(mat, "data.raw.qs")
        n.removed <- n.before - nrow(mat)
        msgSet$current.msg <- c(msgSet$current.msg,
          paste0("Removed ", n.removed, " unannotated protein(s) without a gene symbol; ",
                 nrow(mat), " proteins retained."))
      } else if (length(keep) == 0) {
        msgSet$current.msg <- c(msgSet$current.msg,
          "Warning: all proteins lack a gene symbol; skipping unannotated filter.")
      }
    }
  }

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
  msin <- try(ov_qs_read(msstats.path), silent = TRUE);
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
  ov_qs_save(di, "dat.in.qs")
  return(1L)
}


# Apply MORlog back to the active dataSet
.apply.morlog <- function(dataName) {
  dataSet <- readDataset(dataName)
  di <- ov_qs_read("dat.in.qs")

  dataSet$expr <- di$expr
  dataSet$norm <- di$norm

  dataSet$data.norm <- dataSet$norm
  fast.write(dataSet$data.norm, file = "data_normalized.csv")
  ov_qs_save(dataSet$data.norm, file = "data.stat.qs")

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

  di <- ov_qs_read("dat.in.qs")
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

  ov_qs_save(di, "dat.in.qs")
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
    AddErrMsg("Input matrix must have Peptide Sequences as row names.");
    return(0);
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
    AddErrMsg("No proteins remained after filtering for minimum peptide count.");
    return(0);
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
  peptide_matrix <- ov_qs_read("data.stat.qs")
  protein_matrix <- ov_qs_read("int.mat.qs")
  #msg("[PlotProteinProfile] peptide_matrix dim: ", paste(dim(peptide_matrix), collapse = "x"),
  #        " protein_matrix dim: ", paste(dim(protein_matrix), collapse = "x"))
  ds <- readDataset(dataName)
  if (file.exists("peptide_to_protein_map.qs")) {
    #msg("[PlotProteinProfile] using peptide_to_protein_map.qs")
    peptide_to_protein_map <- ov_qs_read("peptide_to_protein_map.qs")
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
    grDevices::png(filename = finalImg, width = 900, height = 450, res = dpi, type = "cairo")
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
  #save.image("build.RData");
  if (!file.exists("int.mat.qs")) {
    #msg("[BuildProteinOptions] int.mat.qs not found")
    return(character(0))
  }
  ids <- rownames(ov_qs_read("int.mat.qs"))
  labels <- ids
  mapped_entrez <- rep(NA_character_, length(ids))

  if (file.exists("annotation.qs")) {
    anot <- ov_qs_read("annotation.qs")
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

  mq_meta <- ov_qs_read("maxquant_metadata.qs")

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

  diann_meta <- ov_qs_read("diann_metadata.qs")

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

  fp_meta <- ov_qs_read("fragpipe_metadata.qs")

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

  spec_meta <- ov_qs_read("spectronaut_metadata.qs")

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
    data <- ov_qs_read("norm.input.anot.qs")
  } else if (file.exists("orig.data.anot.qs")) {
    data <- ov_qs_read("orig.data.anot.qs")
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
  ov_qs_save(result$data, "norm.input.anot.qs")

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

  mq_meta <- ov_qs_read("maxquant_phospho_metadata.qs")

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

  diann_meta <- ov_qs_read("diann_phospho_metadata.qs")

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

  fp_meta <- ov_qs_read("fragpipe_phospho_metadata.qs")

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

# Proteoform analysis: genes where the ratio of proteoform-specific peptide intensities
# shifts between conditions while aggregate abundance remains stable.
# Requires: diann_metadata.qs (Proteotypic + Gene columns), peptide_to_protein_map.qs,
#           peptide_level_data.qs (post-summarization) or data.stat.qs (pre-summarization).
.proteoform_meta_info <- function(dataSet) {
  if (is.null(dataSet$meta.info)) return(NULL)
  if (is.data.frame(dataSet$meta.info)) return(dataSet$meta.info)
  if (!is.null(dataSet$meta.info$meta.info)) return(dataSet$meta.info$meta.info)
  NULL
}

.proteoform_disc_inx <- function(dataSet) {
  if (!is.null(dataSet$disc.inx)) return(dataSet$disc.inx)
  if (!is.null(dataSet$meta.info$disc.inx)) return(dataSet$meta.info$disc.inx)
  NULL
}

.proteoform_prepare_precursors <- function(meta, prot.map, pmat.ids) {
  if (!"Protein.IDs" %in% names(meta))
    return(list(error = "DIA-NN metadata does not contain Protein.IDs."))
  if (!"Proteotypic" %in% names(meta))
    return(list(error = "Proteotypic column missing from metadata. Re-upload the dataset."))
  if (!all(c("Peptide", "Protein") %in% names(prot.map)))
    return(list(error = "Peptide-to-protein map must contain Peptide and Protein columns."))

  prec.df <- meta[meta$Protein.IDs %in% pmat.ids, , drop = FALSE]
  if (nrow(prec.df) == 0) {
    return(list(error = paste0("No precursors in diann_metadata.qs match the normalized matrix row names. ",
                               "Re-upload the dataset to regenerate metadata.")))
  }

  prec.df$clean.prot <- prot.map$Protein[match(prec.df$Protein.IDs, prot.map$Peptide)]
  prec.df <- prec.df[!is.na(prec.df$Proteotypic) & prec.df$Proteotypic == 1L, , drop = FALSE]
  prec.df <- prec.df[!is.na(prec.df$clean.prot) & prec.df$clean.prot != "", , drop = FALSE]
  if (nrow(prec.df) == 0)
    return(list(error = "No proteotypic precursors with valid protein assignments found in the normalized matrix."))

  # Auto-detect grouping strategy:
  #   UniProt isoform FASTA: >=1% of protein IDs carry a -N suffix (Q09028-2) -> group by canonical base
  #   Gene-name mode (OpenProt / standard without isoform FASTA): group by gene name.
  frac.iso.suffix <- mean(grepl("-[0-9]+$", prec.df$clean.prot))
  if (isTRUE(frac.iso.suffix >= 0.01)) {
    prec.df$group.key <- sub("-[0-9]+$", "", prec.df$clean.prot)
    group.mode <- "uniprot"
  } else {
    gene.col <- NA_character_
    for (cand in c("GN", "Gene", "Genes")) {
      if (cand %in% names(prec.df)) {
        vals <- as.character(prec.df[[cand]])
        if (mean(!is.na(vals) & vals != "") > 0.3) {
          gene.col <- cand
          break
        }
      }
    }
    if (is.na(gene.col))
      return(list(error = "No UniProt isoform suffixes detected and no usable Gene/GN column. Cannot group proteoforms."))
    prec.df$group.key <- as.character(prec.df[[gene.col]])
    prec.df <- prec.df[!is.na(prec.df$group.key) & prec.df$group.key != "", , drop = FALSE]
    if (nrow(prec.df) == 0)
      return(list(error = "Gene column is empty for all proteotypic precursors. Cannot group by gene name."))
    group.mode <- "gene"
  }

  list(prec.df = prec.df, group.mode = group.mode)
}

.proteoform_plot_group <- function(dataSet, meta.info, samples) {
  if (!is.null(meta.info) && !is.null(rownames(meta.info))) {
    analysis.var <- dataSet$analysisVar
    if (!is.null(analysis.var) && analysis.var %in% colnames(meta.info)) {
      vals <- as.character(meta.info[samples, analysis.var])
      if (any(!is.na(vals) & vals != "")) return(factor(vals))
    }
    disc.inx <- .proteoform_disc_inx(dataSet)
    if (!is.null(disc.inx) && any(disc.inx)) {
      disc.cols <- names(which(disc.inx))
      disc.cols <- disc.cols[disc.cols %in% colnames(meta.info)]
      if (length(disc.cols) > 0) {
        vals <- as.character(meta.info[samples, disc.cols[1]])
        if (any(!is.na(vals) & vals != "")) return(factor(vals))
      }
    }
  }
  if (!is.null(dataSet$cls) && length(dataSet$cls) == length(samples)) return(factor(as.character(dataSet$cls)))
  factor(samples, levels = samples)
}

.proteoform_build_context <- function(dataSet, pep.mat) {
  meta.info <- .proteoform_meta_info(dataSet)
  samples.all <- colnames(pep.mat)

  design <- dataSet$design
  contrast.matrix <- dataSet$contrast.matrix
  if (!is.null(design)) {
    design <- as.matrix(design)
    if (!is.null(contrast.matrix)) {
      contrast.matrix <- as.matrix(contrast.matrix)
    } else {
      coef.cands <- unique(as.character(unlist(c(dataSet$active.comp.nm,
                                                 dataSet$par1,
                                                 dataSet$analysis.var))))
      coef.cands <- coef.cands[!is.na(coef.cands) & coef.cands != ""]
      coef.name <- coef.cands[coef.cands %in% colnames(design)][1]
      if (!is.na(coef.name)) {
        contrast.matrix <- matrix(0, nrow = ncol(design), ncol = 1,
                                  dimnames = list(colnames(design), coef.name))
        contrast.matrix[coef.name, 1] <- 1
      }
    }
  }
  has.de.design <- !is.null(design) && !is.null(contrast.matrix) &&
                   is.matrix(design) && is.matrix(contrast.matrix) &&
                   nrow(design) > 1 && ncol(contrast.matrix) >= 1

  if (has.de.design) {
    samples <- NULL

    if (!is.null(rownames(design)) && all(rownames(design) %in% samples.all)) {
      samples <- rownames(design)
    }
    if (is.null(samples) && !is.null(dataSet$rmidx) && length(dataSet$rmidx) > 0 &&
        length(samples.all) - length(dataSet$rmidx) == nrow(design)) {
      samples <- samples.all[-dataSet$rmidx]
    }
    if (is.null(samples) && !is.null(meta.info) && !is.null(rownames(meta.info))) {
      meta.samples <- intersect(rownames(meta.info), samples.all)
      if (length(meta.samples) == nrow(design)) samples <- meta.samples
    }
    if (is.null(samples) && nrow(design) == length(samples.all)) {
      samples <- samples.all
    }

    if (!is.null(samples) && length(samples) == nrow(design)) {
      if (!is.null(rownames(design)) && all(samples %in% rownames(design))) {
        design <- design[samples, , drop = FALSE]
      } else {
        rownames(design) <- samples
      }

      active <- dataSet$active.comp.nm
      if (is.null(active) || length(active) == 0 || !(active %in% colnames(contrast.matrix))) {
        active <- colnames(contrast.matrix)[1]
      }
      contrast <- contrast.matrix[, active, drop = FALSE]

      block <- dataSet$block
      if (!is.null(block) && length(block) > 0) {
        if (length(block) == nrow(design)) {
          names(block) <- rownames(design)
        } else {
          block <- NULL
        }
      }

      plot.group <- .proteoform_plot_group(dataSet, meta.info, samples)
      names(plot.group) <- samples

      coef.vec <- as.numeric(contrast[, 1])
      names(coef.vec) <- rownames(contrast)
      pos.cols <- names(coef.vec)[coef.vec > 0]
      neg.cols <- names(coef.vec)[coef.vec < 0]
      summary.group <- NULL
      cond1.label <- NA_character_
      cond2.label <- NA_character_
      if (length(pos.cols) == 1 && length(neg.cols) == 1 &&
          all(c(pos.cols, neg.cols) %in% colnames(design))) {
        summary.group <- rep(NA_character_, nrow(design))
        summary.group[design[, neg.cols] > 0] <- neg.cols
        summary.group[design[, pos.cols] > 0] <- pos.cols
        names(summary.group) <- rownames(design)
        cond1.label <- neg.cols
        cond2.label <- pos.cols
      }

      return(list(
        samples = samples,
        design = design,
        contrast = contrast,
        contrast.name = active,
        block = block,
        plot.group = plot.group,
        summary.group = summary.group,
        cond1.label = cond1.label,
        cond2.label = cond2.label,
        mode = "differential-analysis design"
      ))
    }
  }

  disc.inx <- .proteoform_disc_inx(dataSet)
  if (is.null(meta.info) || is.null(disc.inx) || sum(disc.inx) == 0)
    return(list(error = "No differential-analysis design or categorical condition found for proteoform analysis."))

  cond.col <- names(which(disc.inx))[1]
  groups <- as.character(meta.info[[cond.col]])
  names(groups) <- rownames(meta.info)
  samples <- intersect(samples.all, names(groups))
  groups <- groups[samples]
  keep <- !is.na(groups) & groups != ""
  samples <- samples[keep]
  groups <- groups[keep]
  unique.groups <- unique(groups)
  if (length(unique.groups) < 2)
    return(list(error = "At least two conditions are required for proteoform analysis."))

  selected.groups <- unique.groups[1:2]
  keep <- groups %in% selected.groups
  samples <- samples[keep]
  groups <- factor(groups[keep], levels = selected.groups)
  if (any(table(groups) < 2))
    return(list(error = "At least 2 replicates per condition are required."))

  design <- model.matrix(~ 0 + groups)
  colnames(design) <- make.names(levels(groups), unique = TRUE)
  rownames(design) <- samples
  contrast <- matrix(c(-1, 1), ncol = 1)
  rownames(contrast) <- colnames(design)
  colnames(contrast) <- paste0(colnames(design)[2], "-", colnames(design)[1])
  summary.group <- setNames(ifelse(groups == selected.groups[1], colnames(design)[1], colnames(design)[2]), samples)
  plot.group <- setNames(groups, samples)

  list(
    samples = samples,
    design = design,
    contrast = contrast,
    contrast.name = colnames(contrast)[1],
    block = NULL,
    plot.group = plot.group,
    summary.group = summary.group,
    cond1.label = colnames(design)[1],
    cond2.label = colnames(design)[2],
    mode = paste0("fallback first-two-group design: ", selected.groups[2], " vs ", selected.groups[1])
  )
}

.proteoform_limma <- function(mat, ctx, robustTrend = FALSE) {
  require(limma)
  mat <- mat[, ctx$samples, drop = FALSE]
  if (nrow(mat) == 0 || ncol(mat) == 0) stop("Empty proteoform matrix.")
  if (!is.fullrank(ctx$design)) stop("Proteoform design matrix is not full rank.")

  fit <- if (is.null(ctx$block)) {
    lmFit(mat, ctx$design)
  } else {
    corfit <- duplicateCorrelation(mat, ctx$design, block = ctx$block)
    lmFit(mat, ctx$design, block = ctx$block, correlation = corfit$consensus)
  }
  if (all(fit$df.residual == 0))
    stop("No residual degrees of freedom in the selected proteoform design.")

  fit2 <- contrasts.fit(fit, ctx$contrast)
  fit2 <- eBayes(fit2, trend = robustTrend, robust = robustTrend)
  tbl <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH", sort.by = "none")
  if (!is.null(tbl$ID)) {
    rownames(tbl) <- tbl$ID
    tbl$ID <- NULL
  }
  tbl
}

# Sample absorbs gene-level abundance; the interaction tests condition-specific isoform composition.
.proteoform_omnibus_anova <- function(iso.mat, ctx) {
  empty <- list(p.value = NA_real_, f.stat = NA_real_)
  if (is.null(iso.mat) || nrow(iso.mat) < 2 || ncol(iso.mat) < 4) return(empty)
  iso.mat <- iso.mat[, ctx$samples, drop = FALSE]

  groups <- ctx$summary.group
  if (is.null(groups) || !any(!is.na(groups) & groups != "")) {
    groups <- ctx$plot.group
  }
  if (is.null(groups)) return(empty)
  if (is.null(names(groups)) && length(groups) == length(ctx$samples)) names(groups) <- ctx$samples

  groups <- as.character(groups[colnames(iso.mat)])
  keep.samples <- !is.na(groups) & groups != ""
  if (sum(keep.samples) < 4) return(empty)
  iso.mat <- iso.mat[, keep.samples, drop = FALSE]
  groups <- factor(groups[keep.samples])
  if (nlevels(groups) < 2 || any(table(groups) < 2)) return(empty)

  df <- data.frame(
    Sample = rep(colnames(iso.mat), each = nrow(iso.mat)),
    Condition = rep(as.character(groups), each = nrow(iso.mat)),
    Isoform = rep(rownames(iso.mat), times = ncol(iso.mat)),
    Value = as.vector(iso.mat),
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$Value) & !is.na(df$Condition) & df$Condition != "", , drop = FALSE]
  if (nrow(df) == 0 || length(unique(df$Isoform)) < 2) return(empty)
  sample.counts <- tapply(df$Sample, df$Condition, function(x) length(unique(x)))
  if (length(sample.counts) < 2 || any(sample.counts < 2)) return(empty)

  df$Sample <- factor(df$Sample)
  df$Condition <- factor(df$Condition)
  df$Isoform <- factor(df$Isoform)

  tryCatch({
    reduced <- lm(Value ~ Sample + Isoform, data = df)
    full <- lm(Value ~ Sample + Isoform + Condition:Isoform, data = df)
    tab <- anova(reduced, full)
    if (nrow(tab) < 2) return(empty)
    p.val <- as.numeric(tab[["Pr(>F)"]][2])
    f.val <- as.numeric(tab[["F"]][2])
    if (!is.finite(p.val)) p.val <- NA_real_
    if (!is.finite(f.val)) f.val <- NA_real_
    list(p.value = p.val, f.stat = f.val)
  }, error = function(e) empty)
}

DetectProteoformAnalysis <- function(dataName) {
  msgSet  <- readSet(msgSet, "msgSet")
  dataSet <- readDataset(dataName)

  fail <- function(msg) {
    msgSet$current.msg <- msg
    saveSet(msgSet, "msgSet")
    return(0L)
  }

  pep.qs <- if (ov_qs_exists("peptide_level_data.qs")) "peptide_level_data.qs" else "data.stat.qs"
  if (!file.exists(pep.qs))
    return(fail("Peptide matrix not found. Run normalization (and summarization if applicable) first."))
  if (!ov_qs_exists("diann_metadata.qs"))
    return(fail("DIA-NN metadata not available. Re-upload the dataset to generate it."))
  if (!ov_qs_exists("peptide_to_protein_map.qs"))
    return(fail("Peptide-to-protein map not found."))

  meta <- ov_qs_read("diann_metadata.qs")
  pep.mat <- ov_qs_read(pep.qs)
  prot.map <- ov_qs_read("peptide_to_protein_map.qs")
  pmat.ids <- rownames(pep.mat)

  prep <- .proteoform_prepare_precursors(meta, prot.map, pmat.ids)
  if (!is.null(prep$error)) return(fail(prep$error))
  prec.df <- prep$prec.df
  group.mode <- prep$group.mode

  ctx <- .proteoform_build_context(dataSet, pep.mat)
  if (!is.null(ctx$error)) return(fail(ctx$error))

  group.isoforms <- tapply(prec.df$clean.prot, prec.df$group.key, function(x) unique(x))
  multi.groups <- names(group.isoforms)[sapply(group.isoforms, length) >= 2]
  if (length(multi.groups) == 0) {
    return(fail(paste0(
      "No ", if (group.mode == "uniprot") "canonical proteins" else "genes",
      " with >=2 distinct proteoforms having proteotypic peptides detected (",
      length(group.isoforms), " group(s) with any proteotypic coverage)."
    )))
  }

  get.gene <- function(grp) {
    if (group.mode == "gene") return(grp)
    if ("Gene" %in% names(prec.df)) {
      idx <- which(prec.df$group.key == grp & !is.na(prec.df$Gene) & prec.df$Gene != "")
      if (length(idx) > 0) return(as.character(prec.df$Gene[idx[1]]))
    }
    grp
  }

  get.ta <- function(grp, iso) {
    if (!"TA" %in% names(prec.df)) return(NA_character_)
    idx <- which(prec.df$group.key == grp & prec.df$clean.prot == iso & !is.na(prec.df$TA) & prec.df$TA != "")
    if (length(idx) > 0) prec.df$TA[idx[1]] else NA_character_
  }

  get.abd <- function(peps) {
    peps <- intersect(peps, pmat.ids)
    if (length(peps) == 0) return(NULL)
    mat <- pep.mat[peps, ctx$samples, drop = FALSE]
    if (length(peps) == 1) as.numeric(mat) else colMeans(mat, na.rm = TRUE)
  }

  mean_for <- function(x, label) {
    if (is.null(ctx$summary.group) || is.na(label)) return(NA_real_)
    idx <- ctx$summary.group == label
    if (!any(idx, na.rm = TRUE)) return(NA_real_)
    val <- mean(x[idx], na.rm = TRUE)
    if (is.nan(val)) NA_real_ else val
  }

  records <- list()
  ratio.rows <- list()
  total.rows <- list()
  min.peptides.per.isoform <- 2L

  for (grp in multi.groups) {
    isoforms <- group.isoforms[[grp]]
    pep.counts <- sapply(isoforms, function(iso)
      sum(prec.df$group.key == grp & prec.df$clean.prot == iso))
    isoforms <- isoforms[order(pep.counts, decreasing = TRUE)]
    eligible.isoforms <- isoforms[pep.counts[isoforms] >= min.peptides.per.isoform]
    if (length(eligible.isoforms) < 2) next

    iso.abd <- lapply(eligible.isoforms, function(iso) {
      peps <- prec.df$Protein.IDs[prec.df$group.key == grp & prec.df$clean.prot == iso]
      get.abd(peps)
    })
    names(iso.abd) <- eligible.isoforms
    valid.isoforms <- names(iso.abd)[vapply(iso.abd, function(x) {
      !is.null(x) && sum(is.finite(x)) >= 3
    }, logical(1))]
    if (length(valid.isoforms) < 2) next
    eligible.isoforms <- valid.isoforms

    iso.mat <- do.call(rbind, iso.abd[eligible.isoforms])
    rownames(iso.mat) <- eligible.isoforms
    colnames(iso.mat) <- ctx$samples
    omnibus <- .proteoform_omnibus_anova(iso.mat, ctx)

    all.abd <- get.abd(prec.df$Protein.IDs[prec.df$group.key == grp])
    if (is.null(all.abd)) next

    for (i in seq_len(length(eligible.isoforms) - 1L)) {
      for (j in seq.int(i + 1L, length(eligible.isoforms))) {
        isoA <- eligible.isoforms[i]
        isoB <- eligible.isoforms[j]

        abd.A <- iso.mat[isoA, ]
        abd.B <- iso.mat[isoB, ]

        log.ratio <- abd.A - abd.B
        if (sum(is.finite(log.ratio)) < 3) next

        pair.id <- make.unique(c(names(records), paste(grp, isoA, isoB, sep = "||")))
        pair.id <- pair.id[length(pair.id)]
        records[[pair.id]] <- list(
          Gene = get.gene(grp),
          Isoform.Group = grp,
          Isoform.Count = length(eligible.isoforms),
          Omnibus.F = omnibus$f.stat,
          Omnibus.Pvalue = omnibus$p.value,
          Proteoform.A = isoA,
          Transcript.A = get.ta(grp, isoA),
          Proteoform.B = isoB,
          Transcript.B = get.ta(grp, isoB),
          Peptides.A = unname(pep.counts[isoA]),
          Peptides.B = unname(pep.counts[isoB]),
          LogRatio.Cond1 = mean_for(log.ratio, ctx$cond1.label),
          LogRatio.Cond2 = mean_for(log.ratio, ctx$cond2.label)
        )
        ratio.rows[[pair.id]] <- log.ratio
        total.rows[[pair.id]] <- all.abd
      }
    }
  }

  if (length(ratio.rows) == 0)
    return(fail(paste0(
      "Insufficient data for proteoform analysis after requiring >= ",
      min.peptides.per.isoform, " proteotypic peptides per isoform."
    )))

  ratio.mat <- do.call(rbind, ratio.rows)
  total.mat <- do.call(rbind, total.rows)
  colnames(ratio.mat) <- ctx$samples
  colnames(total.mat) <- ctx$samples

  ratio.tbl <- tryCatch(.proteoform_limma(ratio.mat, ctx, robustTrend = isTRUE(dataSet$robustTrend)),
                        error = function(e) e)
  if (inherits(ratio.tbl, "error"))
    return(fail(paste("Proteoform ratio model failed:", ratio.tbl$message)))

  total.tbl <- tryCatch(.proteoform_limma(total.mat, ctx, robustTrend = isTRUE(dataSet$robustTrend)),
                        error = function(e) NULL)

  ids <- rownames(ratio.tbl)
  results <- do.call(rbind, lapply(ids, function(id) {
    rec <- records[[id]]
    if (is.null(rec)) return(NULL)
    total.p <- if (!is.null(total.tbl) && id %in% rownames(total.tbl)) total.tbl[id, "P.Value"] else NA_real_
    total.fc <- if (!is.null(total.tbl) && id %in% rownames(total.tbl) && "logFC" %in% colnames(total.tbl)) total.tbl[id, "logFC"] else NA_real_
    data.frame(
      Gene = rec$Gene,
      Comparison = ctx$contrast.name,
      Isoform.Group = rec$Isoform.Group,
      Isoform.Count = rec$Isoform.Count,
      Omnibus.F = rec$Omnibus.F,
      Omnibus.Pvalue = rec$Omnibus.Pvalue,
      Proteoform.A = rec$Proteoform.A,
      Transcript.A = rec$Transcript.A,
      Proteoform.B = rec$Proteoform.B,
      Transcript.B = rec$Transcript.B,
      Peptides.A = rec$Peptides.A,
      Peptides.B = rec$Peptides.B,
      LogRatio.Cond1 = rec$LogRatio.Cond1,
      LogRatio.Cond2 = rec$LogRatio.Cond2,
      Delta.LogRatio = ratio.tbl[id, "logFC"],
      Ratio.Pvalue = ratio.tbl[id, "P.Value"],
      Ratio.FDR = ratio.tbl[id, "adj.P.Val"],
      Total.LogFC = total.fc,
      Total.Pvalue = total.p,
      stringsAsFactors = FALSE
    )
  }))

  results <- results[!is.na(results$Ratio.Pvalue), , drop = FALSE]
  if (nrow(results) == 0)
    return(fail("Proteoform ratio model produced no valid p-values."))

  omnibus.p <- tapply(results$Omnibus.Pvalue, results$Isoform.Group, function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) NA_real_ else x[1]
  })
  omnibus.fdr <- setNames(rep(NA_real_, length(omnibus.p)), names(omnibus.p))
  valid.omnibus <- is.finite(omnibus.p)
  if (any(valid.omnibus)) {
    omnibus.fdr[valid.omnibus] <- p.adjust(omnibus.p[valid.omnibus], method = "BH")
  }
  results$Omnibus.FDR <- unname(omnibus.fdr[results$Isoform.Group])

  results$Is.Switch <- !is.na(results$Omnibus.FDR) &
                       results$Omnibus.FDR < 0.05 &
                       !is.na(results$Ratio.FDR) &
                       results$Ratio.FDR < 0.05 &
                       results$Peptides.A >= min.peptides.per.isoform &
                       results$Peptides.B >= min.peptides.per.isoform &
                       (is.na(results$Total.Pvalue) | results$Total.Pvalue > 0.1)
  results <- results[order(results$Omnibus.Pvalue, results$Ratio.Pvalue, na.last = TRUE), ]

  results$Omnibus.F <- round(results$Omnibus.F, 3)
  results$Omnibus.Pvalue <- signif(results$Omnibus.Pvalue, 3)
  results$Omnibus.FDR <- signif(results$Omnibus.FDR, 3)
  results$LogRatio.Cond1 <- round(results$LogRatio.Cond1, 3)
  results$LogRatio.Cond2 <- round(results$LogRatio.Cond2, 3)
  results$Delta.LogRatio <- round(results$Delta.LogRatio, 3)
  results$Ratio.Pvalue <- signif(results$Ratio.Pvalue, 3)
  results$Ratio.FDR <- signif(results$Ratio.FDR, 3)
  results$Total.LogFC <- round(results$Total.LogFC, 3)
  results$Total.Pvalue <- signif(results$Total.Pvalue, 3)

  result.cols <- c("Gene", "Comparison", "Isoform.Group", "Isoform.Count",
                   "Omnibus.F", "Omnibus.Pvalue", "Omnibus.FDR",
                   "Proteoform.A", "Transcript.A", "Proteoform.B", "Transcript.B",
                   "Peptides.A", "Peptides.B", "LogRatio.Cond1", "LogRatio.Cond2",
                   "Delta.LogRatio", "Ratio.Pvalue", "Ratio.FDR",
                   "Total.LogFC", "Total.Pvalue", "Is.Switch")
  results <- results[, result.cols[result.cols %in% names(results)], drop = FALSE]

  ov_qs_save(results, "proteoform_analysis_results.qs")
  ov_qs_save(list(mode = ctx$mode, contrast = ctx$contrast.name, samples = ctx$samples),
             "proteoform_analysis_context.qs")
  fast.write(results, "proteoform_analysis_results.csv", row.names = FALSE)

  msgSet$current.msg <- paste0(
    "Detected ", nrow(results), " isoform pair(s) across ",
    length(unique(results$Isoform.Group)), " gene/protein group(s) with >= ",
    min.peptides.per.isoform, " proteotypic peptides per isoform. ANOVA tests gene-level isoform usage; pairwise contrasts are post-hoc. Design: ",
    ctx$mode, " (", ctx$contrast.name, ")."
  )
  saveSet(msgSet, "msgSet")
  return(as.integer(nrow(results)))
}

PlotProteoformIntensityProfile <- function(dataName, imageName, isoA, isoB,
                                        format = "png", dpi = 96, paletteOpt = "default",
                                        plotType = "boxplot") {
  require(ggplot2)
  require(Cairo)
  require(grid)

  dataSet  <- readDataset(dataName)
  prot.map <- if (ov_qs_exists("peptide_to_protein_map.qs")) ov_qs_read("peptide_to_protein_map.qs") else NULL
  pep.qs   <- if (ov_qs_exists("peptide_level_data.qs")) "peptide_level_data.qs" else "data.stat.qs"
  if (!file.exists(pep.qs) || is.null(prot.map)) return(invisible(0))
  if (!ov_qs_exists("diann_metadata.qs")) return(invisible(0))

  pep.mat <- ov_qs_read(pep.qs)
  pmat.ids <- rownames(pep.mat)
  meta <- ov_qs_read("diann_metadata.qs")

  prep <- .proteoform_prepare_precursors(meta, prot.map, pmat.ids)
  if (!is.null(prep$error)) return(invisible(0))
  prec.df <- prep$prec.df

  ctx <- .proteoform_build_context(dataSet, pep.mat)
  if (!is.null(ctx$error)) return(invisible(0))

  grp.idx <- which(prec.df$clean.prot %in% c(isoA, isoB))
  if (length(grp.idx) == 0) return(invisible(0))
  grp <- prec.df$group.key[grp.idx[1]]

  peps.A <- prec.df$Protein.IDs[prec.df$group.key == grp & prec.df$clean.prot == isoA]
  peps.B <- prec.df$Protein.IDs[prec.df$group.key == grp & prec.df$clean.prot == isoB]
  peps.all <- prec.df$Protein.IDs[prec.df$group.key == grp]
  if (length(peps.A) == 0 && length(peps.B) == 0) return(invisible(0))

  samples <- ctx$samples
  plot.group <- ctx$plot.group[samples]
  if (all(is.na(plot.group))) plot.group <- factor(samples, levels = samples)
  plot.group <- factor(as.character(plot.group), levels = unique(as.character(plot.group)))

  avg_signal <- function(peps) {
    peps <- intersect(peps, pmat.ids)
    if (length(peps) == 0) return(NULL)
    mat <- pep.mat[peps, samples, drop = FALSE]
    if (length(peps) == 1) as.numeric(mat) else colMeans(mat, na.rm = TRUE)
  }

  abd.A <- avg_signal(peps.A)
  abd.B <- avg_signal(peps.B)
  abd.total <- avg_signal(peps.all)
  if (is.null(abd.A) || is.null(abd.B) || is.null(abd.total)) return(invisible(0))
  log.ratio <- abd.A - abd.B

  trunc_label <- function(x, n = 26) ifelse(nchar(x) > n, paste0(substr(x, 1, n), "..."), x)

  make_rows <- function(abd, label, panel) {
    if (is.null(abd)) return(NULL)
    abd <- as.numeric(abd)
    data.frame(
      Sample    = samples,
      Condition = plot.group,
      Feature   = label,
      Panel     = panel,
      Value     = abd,
      stringsAsFactors = FALSE
    )
  }

  df.intensity <- rbind(make_rows(abd.A, paste0("A: ", trunc_label(isoA, 20)), "Proteoform A/B intensity"),
                        make_rows(abd.B, paste0("B: ", trunc_label(isoB, 20)), "Proteoform A/B intensity"))
  df.ratio <- make_rows(log.ratio, "A - B score", "Proteoform A - B")
  df.total <- make_rows(abd.total, "Aggregate signal", "Aggregate signal")
  if (is.null(df.intensity) || is.null(df.ratio) || is.null(df.total)) return(invisible(0))

  clean_plot_df <- function(df) {
    df <- df[is.finite(df$Value) & !is.na(df$Condition), , drop = FALSE]
    df$Condition <- factor(as.character(df$Condition), levels = levels(plot.group))
    df
  }
  df.intensity <- clean_plot_df(df.intensity)
  df.ratio <- clean_plot_df(df.ratio)
  df.total <- clean_plot_df(df.total)
  if (nrow(df.intensity) == 0 || nrow(df.ratio) == 0 || nrow(df.total) == 0) return(invisible(0))
  df.intensity$Feature <- factor(df.intensity$Feature, levels = unique(df.intensity$Feature))

  abundance.ylim <- range(c(df.intensity$Value, df.total$Value), na.rm = TRUE)
  if (!all(is.finite(abundance.ylim))) abundance.ylim <- NULL
  if (!is.null(abundance.ylim)) {
    pad <- diff(abundance.ylim) * 0.06
    if (!is.finite(pad) || pad == 0) pad <- max(0.5, abs(abundance.ylim[1]) * 0.05)
    abundance.ylim <- abundance.ylim + c(-pad, pad)
  }

  group.cols <- GetGroupPalette(plot.group, paletteOpt)

  plotType <- tolower(trimws(plotType))
  if (!(plotType %in% c("violin", "boxplot"))) plotType <- "boxplot"
  use.violin <- plotType == "violin"

  p1.geom <- if (use.violin) {
    geom_violin(trim = FALSE, position = position_dodge(width = 0.72),
                aes(color = Condition), show.legend = TRUE)
  } else {
    geom_boxplot(position = position_dodge(width = 0.72), outlier.shape = NA,
                 width = 0.58, aes(color = Condition), show.legend = TRUE)
  }

  n.feat.p1 <- length(levels(df.intensity$Feature))
  band.df.p1 <- data.frame(
    xmin = c(-Inf, seq_len(n.feat.p1 - 1) + 0.5),
    xmax = c(seq_len(n.feat.p1 - 1) + 0.5, Inf),
    ymin = -Inf, ymax = Inf,
    Fill = rep(c("#d9d9d9", "#ececec"), length.out = n.feat.p1),
    stringsAsFactors = FALSE
  )

  p1 <- ggplot(df.intensity, aes(x = Feature, y = Value, fill = Condition)) +
    geom_rect(data = band.df.p1,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = band.df.p1$Fill, alpha = 0.7) +
    p1.geom +
    geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.72),
               color = "black", size = 1.25, alpha = 0.65, na.rm = TRUE, show.legend = FALSE) +
    stat_summary(fun = mean, colour = "yellow", geom = "point",
                 shape = 18, size = 2.4, na.rm = TRUE, show.legend = FALSE,
                 position = position_dodge(width = 0.72)) +
    scale_fill_manual(values = group.cols) +
    scale_color_manual(values = group.cols) +
    theme_gray(base_size = 10) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
          axis.ticks.x = element_line(colour = "#555555"),
          legend.position = "right",
          legend.title = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(color = "white", linewidth = 0.35),
          panel.background = element_rect(fill = "#e5e5e5", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          plot.title = element_blank()) +
    ylab("Normalized abundance")
  if (!is.null(abundance.ylim)) {
    p1 <- p1 + coord_cartesian(ylim = abundance.ylim)
  }

  make_group_plot <- function(df, ylab, title, ylims = NULL) {
    plot.geom <- if (use.violin) {
      geom_violin(trim = FALSE, aes(color = Condition), show.legend = FALSE)
    } else {
      geom_boxplot(aes(color = Condition), outlier.shape = NA, width = 0.55, show.legend = FALSE)
    }
    plot <- ggplot(df, aes(x = Condition, y = Value, fill = Condition)) +
      plot.geom +
      geom_jitter(width = 0.08, height = 0, color = "black", size = 1.25,
                  alpha = 0.65, na.rm = TRUE, show.legend = FALSE) +
      stat_summary(fun = mean, colour = "yellow", geom = "point",
                   shape = 18, size = 2.4, na.rm = TRUE, show.legend = FALSE) +
      scale_fill_manual(values = group.cols) +
      scale_color_manual(values = group.cols) +
      theme_gray(base_size = 10) +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_text(size = 8),
            axis.ticks.x = element_line(colour = "#555555"),
            legend.position = "none",
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_line(color = "white", linewidth = 0.35),
            panel.background = element_rect(fill = "#e5e5e5", color = NA),
            plot.background = element_rect(fill = "white", color = NA),
            plot.title = element_text(size = 10, hjust = 0.5, color = "#444444")) +
      ylab(ylab) +
      ggtitle(title)
    if (!is.null(ylims)) plot <- plot + coord_cartesian(ylim = ylims)
    plot
  }

  extract_legend <- function(plot) {
    grob <- ggplotGrob(plot + theme(legend.position = "right"))
    idx <- which(vapply(grob$grobs, function(x) x$name, character(1)) == "guide-box")
    if (length(idx) == 0) return(NULL)
    grob$grobs[[idx[1]]]
  }

  legend.grob <- extract_legend(p1)
  p1 <- p1 + theme(legend.position = "none")
  p2 <- make_group_plot(df.ratio, "A - B abundance score", "Proteoform A - B")
  p3 <- make_group_plot(df.total, "Normalized abundance", "Aggregate signal", abundance.ylim)

  imgName <- paste0(imageName, "dpi", dpi, ".", format)
  Cairo(file = imgName, width = 12.8, height = 5.2, unit = "in",
        dpi = dpi, bg = "white", type = format)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 4, widths = unit(c(1.35, 1, 1, 0.45), "null"))))
  print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(p3, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
  if (!is.null(legend.grob)) {
    grid.draw(editGrob(legend.grob, vp = viewport(layout.pos.row = 1, layout.pos.col = 4)))
  }
  invisible(dev.off())
  return(invisible(1))
}

PlotProteoformOverview <- function(dataName, isoformGroup, imageName,
                                   format = "png", dpi = 96, paletteOpt = "default",
                                   plotType = "boxplot") {
  require(ggplot2)
  require(Cairo)
  require(grid)

  dataSet  <- readDataset(dataName)
  prot.map <- if (ov_qs_exists("peptide_to_protein_map.qs")) ov_qs_read("peptide_to_protein_map.qs") else NULL
  pep.qs   <- if (ov_qs_exists("peptide_level_data.qs")) "peptide_level_data.qs" else "data.stat.qs"
  if (!file.exists(pep.qs) || is.null(prot.map)) return(invisible(0))
  if (!ov_qs_exists("diann_metadata.qs")) return(invisible(0))

  pep.mat  <- ov_qs_read(pep.qs)
  pmat.ids <- rownames(pep.mat)
  meta     <- ov_qs_read("diann_metadata.qs")

  prep <- .proteoform_prepare_precursors(meta, prot.map, pmat.ids)
  if (!is.null(prep$error)) return(invisible(0))
  prec.df <- prep$prec.df

  ctx <- .proteoform_build_context(dataSet, pep.mat)
  if (!is.null(ctx$error)) return(invisible(0))

  grp.inx <- prec.df$group.key == isoformGroup
  if (!any(grp.inx)) return(invisible(0))
  isoforms <- unique(prec.df$clean.prot[grp.inx])

  samples    <- ctx$samples
  plot.group <- ctx$plot.group[samples]
  if (all(is.na(plot.group))) plot.group <- factor(samples, levels = samples)
  plot.group <- factor(as.character(plot.group), levels = unique(as.character(plot.group)))

  avg_signal <- function(peps) {
    peps <- intersect(peps, pmat.ids)
    if (length(peps) == 0) return(NULL)
    mat <- pep.mat[peps, samples, drop = FALSE]
    if (length(peps) == 1) as.numeric(mat) else colMeans(mat, na.rm = TRUE)
  }

  trunc_label <- function(x, n = 24) ifelse(nchar(x) > n, paste0(substr(x, 1, n), "..."), x)

  df.list <- lapply(isoforms, function(iso) {
    peps <- prec.df$Protein.IDs[prec.df$group.key == isoformGroup & prec.df$clean.prot == iso]
    abd  <- avg_signal(peps)
    if (is.null(abd)) return(NULL)
    data.frame(Sample = samples, Condition = plot.group, Isoform = trunc_label(iso),
               Value = as.numeric(abd), stringsAsFactors = FALSE)
  })
  df.intensity <- do.call(rbind, Filter(Negate(is.null), df.list))
  if (is.null(df.intensity) || nrow(df.intensity) == 0) return(invisible(0))
  df.intensity <- df.intensity[is.finite(df.intensity$Value) & !is.na(df.intensity$Condition), ]
  if (nrow(df.intensity) == 0) return(invisible(0))
  df.intensity$Condition <- factor(as.character(df.intensity$Condition), levels = levels(plot.group))
  df.intensity$Isoform   <- factor(df.intensity$Isoform, levels = unique(df.intensity$Isoform))

  peps.all  <- prec.df$Protein.IDs[prec.df$group.key == isoformGroup]
  abd.total <- avg_signal(peps.all)
  df.total  <- if (!is.null(abd.total)) {
    df <- data.frame(Sample = samples, Condition = plot.group,
                     Value = as.numeric(abd.total), stringsAsFactors = FALSE)
    df <- df[is.finite(df$Value) & !is.na(df$Condition), ]
    df$Condition <- factor(as.character(df$Condition), levels = levels(plot.group))
    df
  } else NULL

  ylim.all <- range(c(df.intensity$Value, if (!is.null(df.total)) df.total$Value else NULL), na.rm = TRUE)
  if (all(is.finite(ylim.all))) {
    pad <- diff(ylim.all) * 0.06
    if (!is.finite(pad) || pad == 0) pad <- max(0.5, abs(ylim.all[1]) * 0.05)
    ylim.all <- ylim.all + c(-pad, pad)
  } else ylim.all <- NULL

  group.cols <- GetGroupPalette(plot.group, paletteOpt)
  plotType   <- tolower(trimws(plotType))
  if (!(plotType %in% c("violin", "boxplot"))) plotType <- "boxplot"
  use.violin <- plotType == "violin"
  n.iso      <- nlevels(df.intensity$Isoform)

  make_multi_geom <- function(show.legend = TRUE) {
    if (use.violin)
      geom_violin(trim = FALSE, position = position_dodge(width = 0.72),
                  aes(color = Condition), show.legend = show.legend)
    else
      geom_boxplot(position = position_dodge(width = 0.72), outlier.shape = NA,
                   width = 0.58, aes(color = Condition), show.legend = show.legend)
  }

  n.feat.ov <- nlevels(df.intensity$Isoform)
  band.df.ov <- data.frame(
    xmin = c(-Inf, seq_len(n.feat.ov - 1) + 0.5),
    xmax = c(seq_len(n.feat.ov - 1) + 0.5, Inf),
    ymin = -Inf, ymax = Inf,
    Fill = rep(c("#d9d9d9", "#ececec"), length.out = n.feat.ov),
    stringsAsFactors = FALSE
  )

  p1 <- ggplot(df.intensity, aes(x = Isoform, y = Value, fill = Condition)) +
    geom_rect(data = band.df.ov,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = band.df.ov$Fill, alpha = 0.7) +
    make_multi_geom(show.legend = FALSE) +
    geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.72),
               color = "black", size = 1.25, alpha = 0.65, na.rm = TRUE, show.legend = FALSE) +
    stat_summary(fun = mean, colour = "yellow", geom = "point",
                 shape = 18, size = 2.4, na.rm = TRUE, show.legend = FALSE,
                 position = position_dodge(width = 0.72)) +
    scale_fill_manual(values = group.cols) +
    scale_color_manual(values = group.cols) +
    theme_gray(base_size = 10) +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
          axis.ticks.x = element_line(colour = "#555555"),
          legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(color = "white", linewidth = 0.35),
          panel.background = element_rect(fill = "#e5e5e5", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          plot.title = element_blank()) +
    ylab("Normalized abundance")
  if (!is.null(ylim.all)) p1 <- p1 + coord_cartesian(ylim = ylim.all)

  extract_legend <- function(plot) {
    grob <- ggplotGrob(plot + theme(legend.position = "right", legend.title = element_blank()))
    idx  <- which(vapply(grob$grobs, function(x) x$name, character(1)) == "guide-box")
    if (length(idx) == 0) return(NULL)
    grob$grobs[[idx[1]]]
  }
  p1_leg <- ggplot(df.intensity, aes(x = Isoform, y = Value, fill = Condition)) +
    make_multi_geom(show.legend = TRUE) +
    scale_fill_manual(values = group.cols) + scale_color_manual(values = group.cols) +
    theme_gray(base_size = 10) + theme(legend.title = element_blank())
  legend.grob <- extract_legend(p1_leg)

  panels <- list(p1)
  widths  <- c(max(1, n.iso))

  if (!is.null(df.total) && nrow(df.total) > 0) {
    make_single_geom <- function() {
      if (use.violin) geom_violin(trim = FALSE, aes(color = Condition), show.legend = FALSE)
      else geom_boxplot(aes(color = Condition), outlier.shape = NA, width = 0.55, show.legend = FALSE)
    }
    p2 <- ggplot(df.total, aes(x = Condition, y = Value, fill = Condition)) +
      make_single_geom() +
      geom_jitter(width = 0.08, height = 0, color = "black", size = 1.25, alpha = 0.65,
                  na.rm = TRUE, show.legend = FALSE) +
      stat_summary(fun = mean, colour = "yellow", geom = "point",
                   shape = 18, size = 2.4, na.rm = TRUE, show.legend = FALSE) +
      scale_fill_manual(values = group.cols) + scale_color_manual(values = group.cols) +
      theme_gray(base_size = 10) +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_text(size = 8),
            axis.ticks.x = element_line(colour = "#555555"),
            legend.position = "none",
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_line(color = "white", linewidth = 0.35),
            panel.background = element_rect(fill = "#e5e5e5", color = NA),
            plot.background = element_rect(fill = "white", color = NA),
            plot.title = element_text(size = 10, hjust = 0.5, color = "#444444")) +
      ylab("Normalized abundance") +
      ggtitle("Aggregate signal")
    if (!is.null(ylim.all)) p2 <- p2 + coord_cartesian(ylim = ylim.all)
    panels <- c(panels, list(p2))
    widths  <- c(widths, 1)
  }
  widths <- c(widths, 0.45)

  img.width <- max(8, 2.5 * n.iso + 2.5)
  imgName   <- paste0(imageName, "dpi", dpi, ".", format)
  Cairo(file = imgName, width = img.width, height = 5.2, unit = "in",
        dpi = dpi, bg = "white", type = format)
  grid.newpage()
  n.cols <- length(panels) + 1L
  pushViewport(viewport(layout = grid.layout(1, n.cols, widths = unit(widths, "null"))))
  for (i in seq_along(panels)) {
    print(panels[[i]], vp = viewport(layout.pos.row = 1, layout.pos.col = i))
  }
  if (!is.null(legend.grob)) {
    grid.draw(editGrob(legend.grob, vp = viewport(layout.pos.row = 1, layout.pos.col = n.cols)))
  }
  invisible(dev.off())
  return(invisible(1))
}

# PTM Occupancy Analysis
# UniMod:4 (Carbamidomethyl) is a fixed sample-prep modification — excluded from occupancy pairing.
PTM_FIXED_MODS <- "UniMod:4"

PTM_UNIMOD_NAMES <- c(
  "UniMod:1"   = "Acetyl",
  "UniMod:21"  = "Phospho",
  "UniMod:35"  = "Oxidation",
  "UniMod:36"  = "Dimethyl",
  "UniMod:121" = "GlyGly",
  "UniMod:737" = "TMT6plex"
)

ptm.format.mod.name <- function(unimods) {
  nms <- PTM_UNIMOD_NAMES[unimods]
  nms[is.na(nms)] <- unimods[is.na(nms)]
  paste(sort(nms), collapse = "+")
}

ptm.extract.bio.mods <- function(seq) {
  all.mods <- regmatches(seq, gregexpr("UniMod:[0-9]+", seq, perl = TRUE))[[1]]
  setdiff(all.mods, PTM_FIXED_MODS)
}

.ptmAnnotateOccupancyCompartments <- function(results, org, lib.path) {
  if (is.null(results) || nrow(results) == 0 || !"Gene" %in% names(results)) return(results)

  results$Parent.Compartment <- "Unknown"
  results$Parent.All.Compartments <- "Unknown"
  results$Landscape.Group <- "Unknown"

  entrez.db <- try(queryGeneDB("entrez", org), silent = TRUE)
  if (inherits(entrez.db, "try-error") || is.null(entrez.db) ||
      !all(c("gene_id", "symbol") %in% colnames(entrez.db))) {
    return(results)
  }
  entrez.db[] <- lapply(entrez.db, as.character)
  sym.to.entrez <- setNames(as.character(entrez.db$gene_id), toupper(as.character(entrez.db$symbol)))

  genes <- as.character(results$Gene)
  gene.entrez <- vector("list", length(genes))
  for (i in seq_along(genes)) {
    gene.tokens <- unlist(strsplit(genes[i], "[;, ]+", perl = TRUE), use.names = FALSE)
    gene.tokens <- gene.tokens[nzchar(gene.tokens)]
    gene.entrez[[i]] <- unique(na.omit(sym.to.entrez[toupper(gene.tokens)]))
  }

  loc.path <- paste0(lib.path, org, "/", org, "_localization.qs")
  if (!file.exists(loc.path)) loc.path <- paste0(lib.path, org, "/localization.qs")
  if (!file.exists(loc.path)) {
    results$Landscape.Group <- results$Parent.Compartment
    return(results)
  }

  loc <- try(ov_qs_read(loc.path), silent = TRUE)
  if (inherits(loc, "try-error") || is.null(loc) ||
      !all(c("EntrezID", "Broad.category") %in% colnames(loc))) {
    results$Landscape.Group <- results$Parent.Compartment
    return(results)
  }
  if (!"Main.location" %in% colnames(loc)) loc$Main.location <- "Unknown"
  primary <- rep("Unknown", length(genes))
  all.comp <- rep("Unknown", length(genes))
  loc.entrez <- as.character(loc$EntrezID)

  for (i in seq_along(genes)) {
    entrez.ids <- gene.entrez[[i]]
    if (length(entrez.ids) == 0) next

    loc.rows <- loc[loc.entrez %in% entrez.ids, , drop = FALSE]
    if (nrow(loc.rows) == 0) next

    broad <- paste(unique(as.character(loc.rows$Broad.category)), collapse = "; ")
    main <- paste(unique(as.character(loc.rows$Main.location)), collapse = "; ")
    resolved <- .paPrimaryCompartment(broad, main)
    primary[i] <- resolved$primary
    all.comp[i] <- resolved$all_categories
  }

  results$Parent.Compartment <- primary
  results$Parent.All.Compartments <- all.comp
  results$Landscape.Group <- primary
  results
}

# For each gene with peptides detected in both modified and unmodified forms,
# computes per-sample occupancy (fraction modified) and tests whether it shifts
# between conditions independently of total peptide abundance.
DetectPTMOccupancy <- function(dataName) {
  message("[PTMOccupancy] Starting DetectPTMOccupancy for dataset: ", dataName)
  msgSet  <- readSet(msgSet, "msgSet")
  dataSet <- readDataset(dataName)

  fail <- function(msg) {
    message("[PTMOccupancy] FAIL: ", msg)
    msgSet$current.msg <- msg
    saveSet(msgSet, "msgSet")
    return(0L)
  }

  if (!ov_qs_exists("diann_metadata.qs"))
    return(fail("DIA-NN metadata not available. Re-upload the pr_matrix file."))

  meta <- ov_qs_read("diann_metadata.qs")
  message("[PTMOccupancy] diann_metadata.qs cols: ", paste(names(meta), collapse = ", "))
  message("[PTMOccupancy] diann_metadata rows: ", nrow(meta))

  if (!all(c("Stripped.Sequence", "Modified.Sequence") %in% names(meta)))
    return(fail("Sequence columns missing from metadata. Re-upload the pr_matrix file to regenerate it."))

  pep.qs <- if (ov_qs_exists("peptide_level_data.qs")) "peptide_level_data.qs" else "data.stat.qs"
  message("[PTMOccupancy] Using peptide matrix: ", pep.qs, " (exists=", file.exists(pep.qs), ")")
  if (!file.exists(pep.qs))
    return(fail("Peptide matrix not found. Run normalization first."))

  pep.mat <- ov_qs_read(pep.qs)
  message("[PTMOccupancy] pep.mat dim: ", nrow(pep.mat), " x ", ncol(pep.mat))

  if (!is.null(dataSet$disc.inx)) {
    disc.inx  <- dataSet$disc.inx
    meta.info <- if (is.data.frame(dataSet$meta.info)) dataSet$meta.info else dataSet$meta.info$meta.info
  } else {
    meta.info <- dataSet$meta.info$meta.info
    disc.inx  <- dataSet$meta.info$disc.inx
  }
  if (is.null(meta.info) || is.null(disc.inx) || sum(disc.inx) == 0)
    return(fail("No categorical condition found in dataset metadata."))

  cond.col <- names(which(disc.inx))[1]
  groups   <- as.character(meta.info[[cond.col]])
  names(groups) <- rownames(meta.info)
  message("[PTMOccupancy] cond.col=", cond.col, " groups=", paste(unique(groups), collapse=","))

  samples <- intersect(colnames(pep.mat), names(groups))
  pep.mat <- pep.mat[, samples, drop = FALSE]
  groups  <- groups[samples]
  unique.groups <- unique(groups)
  message("[PTMOccupancy] samples matched: ", length(samples))

  if (length(unique.groups) < 2) return(fail("At least two conditions required."))
  g1.idx <- groups == unique.groups[1]
  g2.idx <- groups == unique.groups[2]
  message("[PTMOccupancy] g1=", sum(g1.idx), " g2=", sum(g2.idx))
  if (sum(g1.idx) < 2 || sum(g2.idx) < 2)
    return(fail("At least 2 replicates per condition required."))

  pmat.ids <- rownames(pep.mat)
  meta <- meta[meta$Protein.IDs %in% pmat.ids, , drop = FALSE]
  message("[PTMOccupancy] meta rows after matching to pep.mat: ", nrow(meta))
  if (nrow(meta) == 0)
    return(fail("No precursors match the peptide matrix. Re-upload the dataset."))

  # Vectorized mod extraction: one gregexpr call over all sequences instead of lapply
  all.matches  <- regmatches(meta$Modified.Sequence,
                              gregexpr("UniMod:[0-9]+", meta$Modified.Sequence, perl = TRUE))
  meta$bio.mods   <- lapply(all.matches, function(m) setdiff(m, PTM_FIXED_MODS))
  meta$is.bio.mod <- lengths(meta$bio.mods) > 0L
  meta$mod.sig    <- vapply(meta$bio.mods, function(m) {
    if (length(m) == 0L) "unmodified" else paste(sort(unique(m)), collapse = "+")
  }, character(1L))
  message("[PTMOccupancy] bio.modified rows: ", sum(meta$is.bio.mod),
          " unmodified rows: ", sum(!meta$is.bio.mod))
  message("[PTMOccupancy] mod.sig table: ", paste(names(table(meta$mod.sig)), table(meta$mod.sig), sep="=", collapse=", "))

  # Vectorized intensity aggregation with rowsum() — avoids per-sequence matrix operations.
  # 1. Exponentiate the whole matrix once to linear scale.
  # 2. rowsum() group-sums by (Stripped.Sequence + mod.sig) in a single C pass.
  # 3. Per-sample NA tracking: if every precursor in a group is NA for a sample → NA result.
  grp.key  <- paste0(meta$Stripped.Sequence, "___", meta$mod.sig)
  lin.mat  <- 2^pep.mat[meta$Protein.IDs, samples, drop = FALSE]
  lin.mat[is.nan(lin.mat)] <- NA

  lin.mat0 <- lin.mat; lin.mat0[is.na(lin.mat0)] <- 0  # NA→0 for rowsum
  det.mat  <- (!is.na(lin.mat)) + 0L                    # 1 if detected, 0 if NA

  grp.sum  <- rowsum(lin.mat0, grp.key)
  grp.det  <- rowsum(det.mat,  grp.key)
  grp.sum[grp.det == 0] <- NA_real_   # restore NA where nothing was detected

  grp.rn  <- rownames(grp.sum)
  grp.seq <- sub("___.*$",  "", grp.rn)
  grp.sig <- sub("^[^_]*___", "", grp.rn)

  prec.count <- table(grp.key)

  # Gene annotation: build lookup table in one pass via match() instead of per-seq scan
  uniq.seqs <- unique(grp.seq)
  gene.for.seq <- local({
    gene.col <- NA_character_
    for (gc in c("Gene", "GN")) {
      if (gc %in% names(meta)) { gene.col <- gc; break }
    }
    if (!is.na(gene.col)) {
      valid    <- !is.na(meta[[gene.col]]) & meta[[gene.col]] != ""
      meta.sub <- meta[valid, ]
      idx      <- match(uniq.seqs, meta.sub$Stripped.Sequence)
      setNames(meta.sub[[gene.col]][idx], uniq.seqs)
    } else {
      setNames(rep(NA_character_, length(uniq.seqs)), uniq.seqs)
    }
  })

  seqs.with.both <- intersect(
    unique(grp.seq[grp.sig == "unmodified"]),
    unique(grp.seq[grp.sig != "unmodified"])
  )
  message("[PTMOccupancy] peptide sequences with both forms: ", length(seqs.with.both))

  if (length(seqs.with.both) == 0)
    return(fail(paste0(
      "No peptide sequences found with both modified and unmodified precursors. ",
      "PTM occupancy requires both forms to be quantified in the same run."
    )))

  mod.sigs.by.seq <- split(
    grp.sig[grp.sig != "unmodified"],
    grp.seq[grp.sig != "unmodified"]
  )

  results.list <- lapply(seqs.with.both, function(sq) {
    int.unmod <- grp.sum[paste0(sq, "___unmodified"), ]
    mod.sigs  <- unique(mod.sigs.by.seq[[sq]])

    lapply(mod.sigs, function(sig) {
      int.mod   <- grp.sum[paste0(sq, "___", sig), ]
      total.int <- int.unmod + int.mod
      occ <- ifelse(is.na(int.unmod) | is.na(int.mod) | total.int == 0,
                    NA_real_, int.mod / total.int)

      valid.g1 <- occ[g1.idx][!is.na(occ[g1.idx])]
      valid.g2 <- occ[g2.idx][!is.na(occ[g2.idx])]
      if (length(valid.g1) < 2 || length(valid.g2) < 2) return(NULL)

      logit <- function(p) log((p + 1e-6) / (1 - p + 1e-6))
      tt.occ  <- tryCatch(t.test(logit(valid.g1), logit(valid.g2)), error = function(e) NULL)
      if (is.null(tt.occ)) return(NULL)

      total.log <- log2(total.int + 1)
      tt.total  <- tryCatch(t.test(total.log[g1.idx], total.log[g2.idx]), error = function(e) NULL)

      data.frame(
        Gene             = gene.for.seq[sq],
        Peptide          = sq,
        Modification     = ptm.format.mod.name(strsplit(sig, "+", fixed = TRUE)[[1]]),
        Mod.Sig          = sig,
        Precursors.Unmod = as.integer(prec.count[paste0(sq, "___unmodified")]),
        Precursors.Mod   = as.integer(prec.count[paste0(sq, "___", sig)]),
        Occupancy.Cond1  = round(mean(valid.g1), 3),
        Occupancy.Cond2  = round(mean(valid.g2), 3),
        Delta.Occupancy  = round(mean(valid.g2) - mean(valid.g1), 3),
        Occ.Pvalue       = signif(tt.occ$p.value, 3),
        Total.Pvalue     = if (!is.null(tt.total)) signif(tt.total$p.value, 3) else NA_real_,
        Total.LogFC      = round(mean(total.log[g2.idx]) - mean(total.log[g1.idx]), 3),
        stringsAsFactors = FALSE
      )
    })
  })

  results <- do.call(rbind, Filter(Negate(is.null),
                                   unlist(results.list, recursive = FALSE)))
  if (is.null(results) || nrow(results) == 0)
    return(fail(paste0(
      "Insufficient data for PTM occupancy analysis. ",
      "Need ≥2 samples with both modified and unmodified forms detected per condition."
    )))

  results$Occ.FDR <- signif(p.adjust(results$Occ.Pvalue, method = "BH"), 3)
  results <- results[order(results$Occ.Pvalue), ]
  results$Cond1.Label <- unique.groups[1]
  results$Cond2.Label <- unique.groups[2]
  paramSet <- readSet(paramSet, "paramSet")
  results <- .ptmAnnotateOccupancyCompartments(results, paramSet$data.org, paramSet$lib.path)

  message("[PTMOccupancy] Final results: ", nrow(results), " pairs")
  ov_qs_save(results, "ptm_occupancy_results.qs")
  fast.write(results, "ptm_occupancy_results.csv")

  msgSet$current.msg <- paste0(
    "Analyzed PTM occupancy for ", nrow(results), " peptide-modification pair(s)."
  )
  saveSet(msgSet, "msgSet")
  message("[PTMOccupancy] Done. Returning ", nrow(results))
  return(as.integer(nrow(results)))
}

PlotPTMOccupancyProfile <- function(dataName, imageName, peptide, modSig,
                                    format = "png", dpi = 96, paletteOpt = "default") {
  require(ggplot2)

  dataSet <- readDataset(dataName)

  if (!ov_qs_exists("diann_metadata.qs")) return(invisible(0))
  meta <- ov_qs_read("diann_metadata.qs")
  if (!all(c("Stripped.Sequence", "Modified.Sequence") %in% names(meta))) return(invisible(0))

  pep.qs <- if (ov_qs_exists("peptide_level_data.qs")) "peptide_level_data.qs" else "data.stat.qs"
  if (!file.exists(pep.qs)) return(invisible(0))
  pep.mat <- ov_qs_read(pep.qs)

  if (!is.null(dataSet$disc.inx)) {
    disc.inx  <- dataSet$disc.inx
    meta.info <- if (is.data.frame(dataSet$meta.info)) dataSet$meta.info else dataSet$meta.info$meta.info
  } else {
    meta.info <- dataSet$meta.info$meta.info
    disc.inx  <- dataSet$meta.info$disc.inx
  }
  cond.col <- names(which(disc.inx))[1]
  groups   <- as.character(meta.info[[cond.col]])
  names(groups) <- rownames(meta.info)
  samples  <- intersect(colnames(pep.mat), names(groups))
  pep.mat  <- pep.mat[, samples, drop = FALSE]
  groups   <- groups[samples]

  # Filter to the target peptide first — avoids running mod extraction on all 100k+ rows
  seq.rows <- meta[meta$Stripped.Sequence == peptide & meta$Protein.IDs %in% rownames(pep.mat), , drop = FALSE]
  if (nrow(seq.rows) == 0) return(invisible(0))

  all.matches        <- regmatches(seq.rows$Modified.Sequence,
                                    gregexpr("UniMod:[0-9]+", seq.rows$Modified.Sequence, perl = TRUE))
  seq.rows$bio.mods   <- lapply(all.matches, function(m) setdiff(m, PTM_FIXED_MODS))
  seq.rows$is.bio.mod <- lengths(seq.rows$bio.mods) > 0L
  seq.rows$mod.sig    <- vapply(seq.rows$bio.mods, function(m) {
    if (length(m) == 0L) "unmodified" else paste(sort(unique(m)), collapse = "+")
  }, character(1L))
  unmod.ids <- seq.rows$Protein.IDs[!seq.rows$is.bio.mod]
  mod.ids   <- seq.rows$Protein.IDs[seq.rows$mod.sig == modSig]

  sum.lin <- function(ids) {
    ids <- intersect(ids, rownames(pep.mat))
    if (length(ids) == 0) return(rep(NA_real_, length(samples)))
    mat.lin <- 2^pep.mat[ids, samples, drop = FALSE]
    if (nrow(mat.lin) == 1) return(as.numeric(mat.lin))
    apply(mat.lin, 2, function(col) if (all(is.na(col))) NA_real_ else sum(col, na.rm = TRUE))
  }

  int.unmod <- sum.lin(unmod.ids)
  int.mod   <- sum.lin(mod.ids)
  total     <- int.unmod + int.mod
  occ       <- ifelse(is.na(int.unmod) | is.na(int.mod) | total == 0,
                      NA_real_, int.mod / total)

  df <- data.frame(
    Sample    = samples,
    Condition = groups[samples],
    Occupancy = occ,
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$Occupancy), ]
  if (nrow(df) == 0) return(invisible(0))

  col <- GetGroupPalette(df$Condition, paletteOpt)

  mod.label <- ptm.format.mod.name(strsplit(modSig, "+", fixed = TRUE)[[1]])
  trunc.seq <- if (nchar(peptide) > 28) paste0(substr(peptide, 1, 28), "…") else peptide

  df$Condition <- factor(df$Condition, levels = names(col))

  myplot <- ggplot(df, aes(x = Condition, y = Occupancy, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.7, color = "#555555") +
    scale_fill_manual(values = col, guide = "none") +
    scale_y_continuous(
      limits = c(0, 1), expand = expansion(mult = c(0.02, 0.05)),
      labels = function(x) paste0(round(x * 100), "%")
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.title.x     = element_blank(),
      axis.text.x      = element_text(size = 9),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.title       = element_text(size = 11, hjust = 0.5),
      plot.subtitle    = element_text(size = 9,  hjust = 0.5, color = "#666666")
    ) +
    ylab("Modification Occupancy") +
    ggtitle(trunc.seq, subtitle = mod.label)

  imgName <- paste0(imageName, "dpi", dpi, ".", format)
  Cairo(file = imgName, width = 5, height = 5, unit = "in", dpi = dpi, bg = "white", type = format)
  invisible(print(myplot))
  invisible(dev.off())
  return(invisible(1))
}

PlotPTMOccupancyLandscape <- function(imageName, format = "png", dpi = 96) {
  require(ggplot2)

  res <- ov_qs_read("ptm_occupancy_results.qs")
  if (is.null(res) || nrow(res) == 0) return(invisible(0))

  cond1.lbl <- if ("Cond1.Label" %in% names(res)) res$Cond1.Label[1] else "Cond1"
  cond2.lbl <- if ("Cond2.Label" %in% names(res)) res$Cond2.Label[1] else "Cond2"

  # Per-modification summary for ordering and significance
  mods <- unique(res$Modification)
  summary_df <- do.call(rbind, lapply(mods, function(mod) {
    sub <- res[res$Modification == mod, , drop = FALSE]
    data.frame(
      Modification = mod,
      Tested       = nrow(sub),
      Median.Occ   = median(c(as.numeric(sub$Occupancy.Cond1),
                               as.numeric(sub$Occupancy.Cond2)), na.rm = TRUE),
      Best.FDR     = min(as.numeric(sub$Occ.FDR), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  summary_df <- summary_df[order(summary_df$Median.Occ), ]
  summary_df$Label <- paste0(summary_df$Modification, "  (n=", summary_df$Tested, ")")
  summary_df$Asterisk <- ifelse(summary_df$Best.FDR <= 0.01, "**",
                         ifelse(summary_df$Best.FDR <= 0.05, "*", ""))

  # Long form: one row per feature per condition
  df_long <- rbind(
    data.frame(Modification = res$Modification,
               Occupancy    = as.numeric(res$Occupancy.Cond1),
               Condition    = cond1.lbl,
               stringsAsFactors = FALSE),
    data.frame(Modification = res$Modification,
               Occupancy    = as.numeric(res$Occupancy.Cond2),
               Condition    = cond2.lbl,
               stringsAsFactors = FALSE)
  )
  df_long <- df_long[!is.na(df_long$Occupancy), ]
  df_long$Label     <- summary_df$Label[match(df_long$Modification, summary_df$Modification)]
  df_long$Label     <- factor(df_long$Label, levels = summary_df$Label)
  df_long$Condition <- factor(df_long$Condition, levels = c(cond1.lbl, cond2.lbl))

  cond_colors <- c(setNames("#3182bd", cond1.lbl), setNames("#e6550d", cond2.lbl))

  # Significance annotation: place at top of y scale per modification
  sig_df <- summary_df[summary_df$Asterisk != "", , drop = FALSE]
  sig_df$Label <- factor(sig_df$Label, levels = summary_df$Label)

  p <- ggplot(df_long, aes(x = Label, y = Occupancy, fill = Condition)) +
    geom_violin(position = position_dodge(0.8), width = 0.7,
                scale = "width", alpha = 0.75, trim = TRUE) +
    stat_summary(aes(group = Condition),
                 fun = median, geom = "point", shape = 18, size = 2.5,
                 position = position_dodge(0.8), color = "white") +
    geom_text(data = sig_df,
              aes(x = Label, y = 1.02, label = Asterisk),
              inherit.aes = FALSE,
              size = 4, color = "#c0392b", fontface = "bold", vjust = 0) +
    scale_fill_manual(values = cond_colors, name = "Condition") +
    scale_y_continuous(labels = function(x) sprintf("%.0f%%", x * 100),
                       limits = c(0, 1.08), expand = c(0, 0)) +
    coord_flip() +
    labs(x = NULL, y = "Occupancy",
         caption = "◆ median   * FDR ≤ 0.05   ** FDR ≤ 0.01") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.y        = element_text(size = 9),
      plot.caption       = element_text(size = 8, color = "#666666", hjust = 0),
      legend.position    = "right",
      legend.title       = element_text(size = 9),
      legend.text        = element_text(size = 8)
    )

  imgNm <- paste0(imageName, "dpi", dpi, ".", format)
  h <- max(3, nrow(summary_df) * 1.1 + 1.5)
  require(Cairo)
  Cairo(file = imgNm, width = 6.5, height = h, unit = "in", dpi = dpi, bg = "white",
        type = if (format == "pdf") "pdf" else "png")
  print(p)
  dev.off()
  return(imgNm)
}

PrepPTMOccupancyVolcano <- function(dataName, fileNm) {
  res <- ov_qs_read("ptm_occupancy_results.qs")
  if (is.null(res) || nrow(res) == 0) {
    message("[PTMVolcano] No PTM occupancy results found")
    return(invisible(0))
  }

  res$label <- make.unique(paste0(res$Gene, " [", res$Modification, "]"))

  delta  <- setNames(as.list(res$Delta.Occupancy), res$label)
  p.log  <- setNames(as.list(-log10(pmax(res$Occ.FDR, 1e-300))), res$label)
  inx.up   <- !is.na(res$Delta.Occupancy) & !is.na(res$Occ.FDR) &
               res$Delta.Occupancy > 0.1 & res$Occ.FDR < 0.05
  inx.down <- !is.na(res$Delta.Occupancy) & !is.na(res$Occ.FDR) &
               res$Delta.Occupancy < -0.1 & res$Occ.FDR < 0.05
  inx.p    <- inx.up | inx.down

  sig.up.ids   <- res$label[inx.up]
  sig.down.ids <- res$label[inx.down]
  non.sig.ids  <- res$label[!inx.p]

  conv <- data.frame(anot.id = res$Gene, symbol = res$Gene, stringsAsFactors = FALSE)

  json.obj <- list(
    fc.log     = delta,
    fc.log.uniq = delta,
    p.log      = p.log,
    p.raw      = setNames(as.list(res$Occ.FDR), res$label),
    inx.up     = setNames(as.list(inx.up),   res$label),
    inx.down   = setNames(as.list(inx.down), res$label),
    inx.p      = setNames(as.list(inx.p),    res$label),
    sigUpIds   = sig.up.ids,
    sigDownIds = sig.down.ids,
    nonSigIds  = non.sig.ids,
    conv       = conv,
    raw.threshx = 0.1,
    raw.threshy = 0.05,
    thresh.y    = -log10(0.05),
    fc.symb     = res$label,
    analType    = "ptm_occupancy",
    org         = paramSet$data.org,
    naviString  = "PTM Occupancy Volcano"
  )

  json.path <- paste0(fileNm, ".json")
  write(RJSONIO::toJSON(json.obj), json.path)
  return(invisible(1))
}
