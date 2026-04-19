
ReadPhosphoData <- function(fileName, metafileName, phosphoLocProb = 0, dataFormat="maxquant", proteinRefFile="", peptideFdr = 0.01, proteinFdr = 0.01, removeContaminants = TRUE, removeDecoys = TRUE, residues = c("S", "T", "Y")) {

  #msg("[ReadPhosphoData] Starting phosphoproteomics data processing.")

  # Set options for the specialized reader
  opts <- list(
    loc.prob = phosphoLocProb,
    removeContaminants = removeContaminants,
    removeDecoys = removeDecoys,
    peptide.fdr = peptideFdr,
    protein.fdr = proteinFdr,
    residues = residues
  )

  # Support multiple phospho data formats
  if (dataFormat == "fragpipe") {
      msg("[ReadPhosphoData] Detected FragPipe format. Using FragPipe reader.")
      dataSet <- .readFragPipePhospho(fileName, opts)
  } else if (dataFormat == "diann") {
      msg("[ReadPhosphoData] Detected DIA-NN format. Using DIA-NN reader.")
      dataSet <- .readDIANNPhospho(fileName, opts)
  } else if (dataFormat %in% c("text", "table", "plain")) {
    msg("[ReadPhosphoData] Detected plain text format. Using tab reader.")
    dataSet <- .readTabData(fileName)
  } else if (dataFormat == "maxquant" || grepl("Phospho \\(STY\\)Sites\\.txt$", fileName, ignore.case = TRUE)) {
    #msg("[ReadPhosphoData] Using MaxQuant phosphosite reader.")
    dataSet <- .readMaxQuantPhospho(fileName, opts)
  } else {
    msgSet$current.msg <- "Unsupported data format for phosphoproteomics analysis."
    saveSet(msgSet, "msgSet")
    return(0)
  }
  
  if (is.null(dataSet)) {
    msgSet$current.msg <- "Failed to read and process the phosphoproteomics data file."
    saveSet(msgSet, "msgSet")
    return(0)
  }
  
  # --- Metadata Processing ---
  # For MaxQuant phospho data, the metadata file may have "Intensity " prefix in sample names
  # We need to preprocess the metadata to strip this prefix to match our column names

  # Check if metadata file is empty or contains only BOM/control characters
  metadata_is_valid <- FALSE
  if (metafileName != "" && file.exists(metafileName)) {
    file_size <- file.info(metafileName)$size
    if (!is.na(file_size) && file_size > 10) {  # At least 10 bytes for minimal content
      metadata_is_valid <- TRUE
    }
  }

  if (!metadata_is_valid) {
    # If metadata file is empty or invalid, set to empty string to trigger synthesis
    metafileName <- ""
  } else if (dataFormat == "maxquant" && metafileName != "") {
    # Read metadata file and check if it has "Intensity " prefix
    if (file.exists(metafileName)) {
      meta_raw <- try(data.table::fread(metafileName, header = TRUE, data.table = FALSE), silent = TRUE)
      if (!inherits(meta_raw, "try-error") && ncol(meta_raw) > 0) {
        # Check if first column (usually "Sample") has "Intensity " prefix
        first_col <- meta_raw[, 1]
        if (any(grepl("^Intensity ", first_col))) {
          #msg("[ReadPhosphoData] Removing 'Intensity ' prefix from metadata sample names.")
          meta_raw[, 1] <- sub("^Intensity ", "", meta_raw[, 1])
          # Save the preprocessed metadata to a temp file
          temp_meta_file <- paste0(tempfile(), "_meta.txt")
          data.table::fwrite(meta_raw, temp_meta_file, sep = "\t", quote = FALSE)
          metafileName <- temp_meta_file
        }
      } else {
        # If reading failed, set to empty to trigger synthesis
        metafileName <- ""
      }
    }
  }

  # This part is similar to ReadTabExpressData
  msg("[ReadPhosphoData] ========================================")
  msg("[ReadPhosphoData] About to load metadata...")
  msg("[ReadPhosphoData] metafileName: '", metafileName, "'")
  msg("[ReadPhosphoData] File exists: ", file.exists(metafileName))
  if (file.exists(metafileName)) {
    first_line <- readLines(metafileName, n = 1)
    msg("[ReadPhosphoData] Metadata file first line: ", first_line)
  }
  msg("[ReadPhosphoData] dataSet$data_orig class: ", class(dataSet$data_orig))
  msg("[ReadPhosphoData] dataSet$data_orig dimensions: ", paste(dim(dataSet$data_orig), collapse=" rows x "), " cols")
  msg("[ReadPhosphoData] dataSet$data_orig column names (ALL): ", paste(colnames(dataSet$data_orig), collapse=", "))
  msg("[ReadPhosphoData] dataSet$data column names (first 5): ", paste(head(colnames(dataSet$data), 5), collapse=", "))
  msg("[ReadPhosphoData] ========================================")

  meta.info <- if (metafileName != "") {
    msg("[ReadPhosphoData] Calling .readMetaData()...")
    result <- .readMetaData(metafileName, dataSet$data_orig, metaContain = "false")
    msg("[ReadPhosphoData] .readMetaData() returned: ", if(is.null(result)) "NULL" else "data")
    if (!is.null(result) && !is.null(result$meta.info)) {
      msg("[ReadPhosphoData] Metadata loaded successfully: ", paste(dim(result$meta.info), collapse=" rows x "), " columns")
      msg("[ReadPhosphoData] Metadata columns: ", paste(colnames(result$meta.info), collapse=", "))
    }
    result
  } else {
    msg("[ReadPhosphoData] metafileName is empty string - skipping metadata load")
    NULL
  }
  # If a plain data.frame slipped through, wrap it into the expected structure
  if (!is.null(meta.info) && is.data.frame(meta.info) && is.null(meta.info$meta.info)) {
    meta.df <- meta.info
    # Use GetDiscreteInx to properly detect discrete vs continuous columns
    disc.inx <- GetDiscreteInx(meta.df, min.rep = 2)
    cont.inx <- setNames(!disc.inx, names(disc.inx))

    # DEBUG: Print metadata structure
    msg(paste0("[ReadPhosphoData] Metadata columns: ", paste(colnames(meta.df), collapse=", ")))
    msg(paste0("[ReadPhosphoData] Discrete columns: ", paste(names(disc.inx)[disc.inx], collapse=", ")))
    msg(paste0("[ReadPhosphoData] First 3 rows of metadata:"))
    print(head(meta.df, 3))

    meta.info <- list(meta.info = meta.df, disc.inx = disc.inx, cont.inx = cont.inx)
  }
  
  if (is.null(meta.info) || is.null(meta.info$meta.info)) {
    msg("[ReadPhosphoData][ERROR] ========================================")
    msg("[ReadPhosphoData][ERROR] METADATA FILE FAILED TO LOAD!")
    msg("[ReadPhosphoData][ERROR] metafileName provided: ", metafileName)
    msg("[ReadPhosphoData][ERROR] Falling back to synthesized metadata (sample names as groups)")
    msg("[ReadPhosphoData][ERROR] This will cause incorrect groupings in analysis!")
    msg("[ReadPhosphoData][ERROR] ========================================")
    meta.info <- .synthesizeMetaFromRuns(colnames(dataSet$data))
  }
  
  dataSet$meta.info <- meta.info$meta.info
  dataSet$disc.inx <- meta.info$disc.inx
  dataSet$cont.inx <- meta.info$cont.inx

  # DEBUG: Verify metadata was loaded correctly
  msg(paste0("[ReadPhosphoData] FINAL metadata structure:"))
  msg(paste0("[ReadPhosphoData] Columns: ", paste(colnames(dataSet$meta.info), collapse=", ")))
  msg(paste0("[ReadPhosphoData] Discrete: ", paste(names(dataSet$disc.inx)[dataSet$disc.inx], collapse=", ")))
  msg(paste0("[ReadPhosphoData] Continuous: ", paste(names(dataSet$cont.inx)[dataSet$cont.inx], collapse=", ")))
  # Convert factors to character to show labels, not numeric codes
  msg(paste0("[ReadPhosphoData] Sample 1 metadata: ", paste(as.character(unlist(dataSet$meta.info[1,])), collapse=", ")))
  if (nrow(dataSet$meta.info) >= 3) {
    msg(paste0("[ReadPhosphoData] Sample 3 metadata: ", paste(as.character(unlist(dataSet$meta.info[3,])), collapse=", ")))
  }
  # Also show the actual data types and factor levels for Group column
  if ("Group" %in% colnames(dataSet$meta.info)) {
    msg(paste0("[ReadPhosphoData] Group column class: ", class(dataSet$meta.info$Group)))
    msg(paste0("[ReadPhosphoData] Group factor levels: ", paste(levels(dataSet$meta.info$Group), collapse=", ")))
    msg(paste0("[ReadPhosphoData] Group unique values: ", paste(unique(as.character(dataSet$meta.info$Group)), collapse=", ")))
  }
  
  # --- Final Alignment and Processing ---
  int.mat <- dataSet$data
  # Treat strict zeros as missing (common in MaxQuant/Proteomics exports)
  zero.count <- sum(int.mat == 0, na.rm = TRUE)
  na.count <- sum(is.na(int.mat))
  int.mat[int.mat == 0] <- NA
  #msg("[ReadPhosphoData] Missing values before zero->NA: NA=", na.count, " zeros=", zero.count,
  #        " total=", na.count + zero.count)
  
  # Align matrix columns with metadata rows
  common_samples <- intersect(colnames(int.mat), rownames(dataSet$meta.info))
  if (length(common_samples) == 0) {
    # Add diagnostic information to help debug sample name mismatches
    #msg(paste0("[ReadPhosphoData][ERROR] No matching samples between data and metadata!"))
    #msg(paste0("[ReadPhosphoData] Data columns (first 5): ", paste(head(colnames(int.mat), 5), collapse=", ")))
    #msg(paste0("[ReadPhosphoData] Metadata rows (first 5): ", paste(head(rownames(dataSet$meta.info), 5), collapse=", ")))

    msgSet$current.msg <- paste0("No samples could be matched between data and metadata. ",
                                  "Data columns: ", paste(head(colnames(int.mat), 3), collapse=", "),
                                  ". Metadata rows: ", paste(head(rownames(dataSet$meta.info), 3), collapse=", "))
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Log successful alignment (commented out for production, uncomment for debugging)
  #msg(paste0("[ReadPhosphoData] Successfully matched ", length(common_samples), " / ",
  #           ncol(int.mat), " samples between data and metadata"))
  
  int.mat <- int.mat[, common_samples, drop = FALSE]
  dataSet$meta.info <- dataSet$meta.info[common_samples, , drop = FALSE]

  # Reorder matrix to match metadata order
  int.mat <- int.mat[, rownames(dataSet$meta.info), drop = FALSE]

  # Ensure metadata index vectors are present
  if (is.null(dataSet$disc.inx) || length(dataSet$disc.inx) == 0) {
    dataSet$disc.inx <- setNames(rep(TRUE, ncol(dataSet$meta.info)), colnames(dataSet$meta.info))
  }
  if (is.null(dataSet$cont.inx) || length(dataSet$cont.inx) == 0) {
    dataSet$cont.inx <- setNames(rep(FALSE, ncol(dataSet$meta.info)), colnames(dataSet$meta.info))
  }

  #msg("[ReadPhosphoData] Aligned data: ", nrow(int.mat), " sites and ", ncol(int.mat), " samples.")

  # Report missing values after alignment
  total_values <- nrow(int.mat) * ncol(int.mat)
  na_count_aligned <- sum(is.na(int.mat))
  missing_pct <- round(100 * na_count_aligned / total_values, 1)
  #msg("[ReadPhosphoData] Missing values after alignment: ", na_count_aligned, " / ", total_values,
  #        " (", missing_pct, "%)")

  # Report missing values per sample
  missing_per_sample <- colSums(is.na(int.mat))
  #msg("[ReadPhosphoData] Missing values per sample: min=", min(missing_per_sample),
  #        ", max=", max(missing_per_sample), ", mean=", round(mean(missing_per_sample), 1))

  # --- Load Global Proteome Reference (Optional) ---
  # For phosphoproteomics, users can optionally provide total protein abundance data
  # (e.g., proteinGroups.txt from MaxQuant) for protein-level normalization
  paramSet <- readSet(paramSet, "paramSet")

  if (proteinRefFile != "" && file.exists(proteinRefFile)) {
    #msg("[ReadPhosphoData] Loading global proteome reference file: ", proteinRefFile)
    protein_data <- .readMaxQuantProteome(proteinRefFile, colnames(int.mat))

    if (!is.null(protein_data) && nrow(protein_data) > 0) {
      paramSet$protein.ref <- protein_data
      paramSet$has.protein.ref <- TRUE
      #msg("[ReadPhosphoData] Global proteome file loaded: ", nrow(protein_data), " proteins × ", ncol(protein_data), " samples.")

      # Report missing values in protein reference data
      total_protein_values <- nrow(protein_data) * ncol(protein_data)
      na_protein <- sum(is.na(protein_data))
      protein_missing_pct <- round(100 * na_protein / total_protein_values, 1)
      #msg("[ReadPhosphoData] Protein reference missing values: ", na_protein, " / ", total_protein_values,
      #        " (", protein_missing_pct, "%)")

      missing_per_protein_sample <- colSums(is.na(protein_data))
      #msg("[ReadPhosphoData] Protein reference missing per sample: min=", min(missing_per_protein_sample),
      #        ", max=", max(missing_per_protein_sample), ", mean=", round(mean(missing_per_protein_sample), 1))

      msgSet$current.msg <- paste0(msgSet$current.msg, " Global proteome reference loaded (", nrow(protein_data), " proteins).")
    } else {
      paramSet$has.protein.ref <- FALSE
      #msg("[ReadPhosphoData] Warning: Failed to process global proteome file.")
    }
  } else {
    paramSet$has.protein.ref <- FALSE
    if (proteinRefFile != "" && !file.exists(proteinRefFile)) {
      #msg("[ReadPhosphoData] Warning: Protein reference file specified but not found: ", proteinRefFile)
    } else {
      #msg("[ReadPhosphoData] No global proteome reference provided.")
    }
  }

  # --- Save and Register ---
  paramSet$data.type <- "phospho"
  paramSet$anal.type <- "onedata"
  paramSet$data.format <- dataSet$format
  paramSet$oneDataAnalType <- "phospho";

  ov_qs_save(int.mat, "data.raw.qs")
  ov_qs_save(int.mat, "int.mat.qs");
  fast.write(sanitizeSmallNumbers(int.mat), file="data_original.csv")
  
  dataSet$data.norm <- int.mat
  dataSet$name <- fileName
  paramSet$dataName <- fileName;

  dataSet$cls <- dataSet$meta.info[, 1]
  
  # Clean up and save
  dataSet$data <- NULL
  dataSet$data_orig <- NULL

  # Final missing value summary
  final_na <- sum(is.na(int.mat))
  final_total <- nrow(int.mat) * ncol(int.mat)
  final_pct <- round(100 * final_na / final_total, 1)
  #msg("[ReadPhosphoData] ========================================")
  #msg("[ReadPhosphoData] FINAL DATA SUMMARY:")
  #msg("[ReadPhosphoData]   Sites: ", nrow(int.mat))
  #msg("[ReadPhosphoData]   Samples: ", ncol(int.mat))
  #msg("[ReadPhosphoData]   Total values: ", final_total)
  #msg("[ReadPhosphoData]   Missing values: ", final_na, " (", final_pct, "%)")
  #msg("[ReadPhosphoData]   Present values: ", final_total - final_na, " (", round(100 - final_pct, 1), "%)")
  if (paramSet$has.protein.ref) {
    #msg("[ReadPhosphoData]   Protein reference: AVAILABLE (", nrow(paramSet$protein.ref), " proteins)")
  } else {
    #msg("[ReadPhosphoData]   Protein reference: NOT AVAILABLE")
  }
  #msg("[ReadPhosphoData] ========================================")

  # CRITICAL: Phosphosite IDs (e.g., P12345_S_123) are NOT Entrez IDs
  # and should NOT be annotated as such. Set annotated flag to FALSE
  # to prevent doEntrezIDAnot() from being called on phosphosite IDs during DE analysis.
  dataSet$annotated <- FALSE
  #msg("[ReadPhosphoData] Set dataSet$annotated = FALSE (phosphosites are not Entrez IDs)")

  # Generate filtering message if statistics are available
  if (exists("filterStatsPhospho", envir = .GlobalEnv)) {
    stats <- get("filterStatsPhospho", envir = .GlobalEnv)
    if (!is.null(stats) && stats$n_total > 0) {
      msg_parts <- c(sprintf("Format-specific filtering (MaxQuant phospho): Starting with %d phosphosites.", stats$n_total))
      if (stats$n_residues > 0) {
        msg_parts <- c(msg_parts, sprintf("Removed %d phosphosites (excluded residue types).", stats$n_residues))
      }
      if (stats$n_contaminants > 0) {
        msg_parts <- c(msg_parts, sprintf("Removed %d contaminants.", stats$n_contaminants))
      }
      if (stats$n_decoys > 0) {
        msg_parts <- c(msg_parts, sprintf("Removed %d decoys/reverse hits.", stats$n_decoys))
      }
      if (stats$n_loc_prob > 0) {
        msg_parts <- c(msg_parts, sprintf("Removed %d phosphosites below localization probability threshold (%.2f).", stats$n_loc_prob, opts$loc.prob))
      }
      msg_parts <- c(msg_parts, sprintf("Remaining: %d phosphosites.", stats$n_final))
      filter_msg <- paste(msg_parts, collapse = " ")

      # Store in msgSet for Java to retrieve
      msgSet$current.msg <- filter_msg
      saveSet(msgSet, "msgSet")
    }
  }

  saveSet(paramSet, "paramSet")

  return(RegisterData(dataSet))
}

# Convert intensity-like columns to numeric safely (handles integer64 and character vectors).
.coerceNumericDF <- function(df) {
  if (is.null(df) || ncol(df) == 0) {
    return(matrix(nrow = nrow(df), ncol = 0))
  }
  to_numeric_col <- function(x) {
    if (inherits(x, "integer64")) {
      return(as.numeric(as.character(x)))
    }
    suppressWarnings(as.numeric(x))
  }
  mat <- vapply(df, to_numeric_col, numeric(nrow(df)))
  if (is.null(dim(mat))) {
    mat <- matrix(mat, ncol = 1)
  }
  rownames(mat) <- rownames(df)
  colnames(mat) <- colnames(df)
  return(mat)
}

.readMaxQuantPhospho <- function(filePath, opts) {

  if (!file.exists(filePath)) {
    #msg("[.readMaxQuantPhospho] Error: File not found at ", filePath)
    return(NULL)
  }

  dat <- try(data.table::fread(
    filePath,
    header = TRUE,
    check.names = FALSE,
    data.table = FALSE,
    integer64 = "double"
  ))
  if (inherits(dat, "try-error")) {
    #msg("[.readMaxQuantPhospho] Error: Failed to read the file.")
    return(NULL)
  }

  # Track filtering statistics
  n_total <- nrow(dat)
  n_contaminants <- 0
  n_decoys <- 0
  n_residues <- 0
  n_loc_prob <- 0

  # 1. Filter by residue type (S/T/Y)
  if (!is.null(opts$residues) && length(opts$residues) > 0 && "Amino acid" %in% colnames(dat)) {
    n_before <- nrow(dat)
    dat <- dat[dat$`Amino acid` %in% opts$residues, ]
    n_residues <- n_before - nrow(dat)
  }

  # 2. Filter contaminants and reverse hits
  if (isTRUE(opts$removeContaminants) && "Potential contaminant" %in% colnames(dat)) {
    n_before <- nrow(dat)
    dat <- dat[dat$`Potential contaminant` != "+", ]
    n_contaminants <- n_before - nrow(dat)
  }
  if (isTRUE(opts$removeDecoys) && "Reverse" %in% colnames(dat)) {
    n_before <- nrow(dat)
    dat <- dat[dat$Reverse != "+", ]
    n_decoys <- n_before - nrow(dat)
  }

  # Note: "Number of Phospho (STY)" column represents phosphorylation multiplicity (biological state),
  # not peptide count, so it should not be used for quality filtering

  # 4. Filter by localization probability
  if ("Localization prob" %in% colnames(dat)) {
    n_before <- nrow(dat)
    dat <- dat[dat$`Localization prob` >= opts$loc.prob, ]
    n_loc_prob <- n_before - nrow(dat)
  }

  # Store filtering statistics for reporting
  if (!exists("filterStatsPhospho", envir = .GlobalEnv)) {
    assign("filterStatsPhospho", list(), envir = .GlobalEnv)
  }
  filterStatsPhospho <- list(
    n_total = n_total,
    n_residues = n_residues,
    n_contaminants = n_contaminants,
    n_decoys = n_decoys,
    n_loc_prob = n_loc_prob,
    n_final = nrow(dat)
  )
  assign("filterStatsPhospho", filterStatsPhospho, envir = .GlobalEnv)

  if (nrow(dat) == 0) {
    #msg("[.readMaxQuantPhospho] Error: No data remains after filtering.")
    return(NULL)
  }
  
  # 3. Handle Multiplicity
  intensity_cols <- grep("^Intensity ", colnames(dat), value = TRUE)
  if (length(intensity_cols) == 0) {
    return(NULL)
  }

  intensity_df <- dat[, intensity_cols, drop = FALSE]
  intensity_num <- .coerceNumericDF(intensity_df)
  
  # Extract base sample names by removing multiplicity suffix (___1, ___2, etc.)
  base_sample_names <- unique(sub("___[0-9]+$", "", intensity_cols))
  
  # Aggregate intensities by summing up multiplicities for each sample
  agg_intens <- sapply(base_sample_names, function(sample_name) {
    multiplicity_cols <- grep(paste0("^", sample_name, "(___[0-9]+)?$"), intensity_cols, value = TRUE)
    rowSums(intensity_num[, multiplicity_cols, drop = FALSE], na.rm = TRUE)
  })
  
  # Clean column names
  colnames(agg_intens) <- sub("^Intensity ", "", colnames(agg_intens))
  
  # Convert 0 back to NA (as rowSums with na.rm=T produces 0 for all-NA rows)
  agg_intens[agg_intens == 0] <- NA
  
  # 4. Create Composite Key for phosphosites
  protein_col <- dat$Proteins
  # Take the first protein ID if multiple are listed
  protein_ids <- sapply(strsplit(protein_col, ";"), function(x) x[1])
  
  residue_col <- dat$`Amino acid`
  position_col <- dat$Position
  
  # Create a unique identifier for each site
  site_ids <- paste(protein_ids, residue_col, position_col, sep = "_")
  
  # Prepare feature info BEFORE aggregation
  # Only select columns that actually exist in the file
  possible_feature_cols <- c("id", "Proteins", "Gene names", "Amino acid", "Position",
                             "Sequence window", "Localization prob", "Fasta headers")
  feature_cols <- possible_feature_cols[possible_feature_cols %in% colnames(dat)]

  # Ensure we have at least the essential columns
  if (!"Proteins" %in% feature_cols || !"Amino acid" %in% feature_cols || !"Position" %in% feature_cols) {
    #msg("[.readMaxQuantPhospho] Error: Missing essential feature columns.")
    return(NULL)
  }

  # Add site_id to feature info for tracking
  feature_info <- dat[, feature_cols, drop = FALSE]
  feature_info$site_id <- site_ids

  # The same site can be on different peptides, leading to duplicate IDs.
  # We aggregate these by taking the mean intensity.
  if (any(duplicated(site_ids))) {
    #msg("[.readMaxQuantPhospho] Aggregating duplicate sites by mean.")
    agg_df <- as.data.frame(agg_intens)
    agg_df$site_id <- site_ids

    # Use split/lapply instead of aggregate(. ~ site_id, ...).
    final_intens <- do.call(rbind, lapply(
      split(agg_df[, colnames(agg_df) != "site_id", drop = FALSE],
            agg_df$site_id),
      function(g) colMeans(g, na.rm = TRUE)
    ))
    final_intens <- as.data.frame(final_intens)
    agg_intens <- as.matrix(final_intens)

    # Also aggregate feature info - take first occurrence of each unique site
    # For numeric columns (like Localization prob), take the mean
    feature_info_agg <- feature_info[!duplicated(feature_info$site_id), , drop = FALSE]
    rownames(feature_info_agg) <- feature_info_agg$site_id

    # Match the order of aggregated intensities
    feature_info_final <- feature_info_agg[rownames(agg_intens), , drop = FALSE]
    feature_info_final$site_id <- NULL  # Remove temporary column
  } else {
    rownames(agg_intens) <- site_ids
    feature_info_final <- feature_info
    feature_info_final$site_id <- NULL  # Remove temporary column
    rownames(feature_info_final) <- site_ids
  }

  # 5. Prepare the object to return
  # For .readMetaData to work, data_orig needs a feature column
  data_orig_df <- cbind(data.frame(SiteID = rownames(agg_intens)), as.data.frame(agg_intens))

  return(list(
    data = agg_intens,
    data_orig = data_orig_df,
    type = "phospho",
    format = "maxquant-phospho",
    meta.info = NULL, # To be filled by the caller
    feature.info = feature_info_final
  ))
}

#' Read MaxQuant protein abundance data for phosphoproteomics normalization
#'
#' @param proteinFile Path to MaxQuant proteinGroups.txt file
#' @param phospho_samples Character vector of sample names from phospho data
#' @return Matrix of protein abundances with proteins as rows, samples as columns
.readMaxQuantProteome <- function(proteinFile, phospho_samples) {

  if (!file.exists(proteinFile)) {
    #msg("[.readMaxQuantProteome] Error: File not found at ", proteinFile)
    return(NULL)
  }

  # Read proteinGroups.txt
  dat <- try(data.table::fread(
    proteinFile,
    header = TRUE,
    check.names = FALSE,
    data.table = FALSE,
    integer64 = "double"
  ))
  if (inherits(dat, "try-error")) {
    #msg("[.readMaxQuantProteome] Error: Failed to read the file.")
    return(NULL)
  }

  #msg("[.readMaxQuantProteome] Read proteinGroups.txt with ", nrow(dat), " rows and ", ncol(dat), " columns.")

  # Filter out contaminants, reverse hits, and proteins only identified by site
  if ("Potential contaminant" %in% colnames(dat)) {
    contam_filter <- is.na(dat$`Potential contaminant`) | dat$`Potential contaminant` == ""
    dat <- dat[contam_filter, ]
    #msg("[.readMaxQuantProteome] Removed contaminants: ", nrow(dat), " proteins remaining.")
  }

  if ("Reverse" %in% colnames(dat)) {
    reverse_filter <- is.na(dat$Reverse) | dat$Reverse == ""
    dat <- dat[reverse_filter, ]
    #msg("[.readMaxQuantProteome] Removed reverse hits: ", nrow(dat), " proteins remaining.")
  }

  if ("Only identified by site" %in% colnames(dat)) {
    site_only_filter <- is.na(dat$`Only identified by site`) | dat$`Only identified by site` == ""
    dat <- dat[site_only_filter, ]
    #msg("[.readMaxQuantProteome] Removed 'only by site' proteins: ", nrow(dat), " proteins remaining.")
  }

  # Extract "Intensity " columns
  intensity_cols <- grep("^Intensity ", colnames(dat), value = TRUE)
  if (length(intensity_cols) == 0) {
    #msg("[.readMaxQuantProteome] Error: No 'Intensity ' columns found.")
    return(NULL)
  }

  intensity_mat <- .coerceNumericDF(dat[, intensity_cols, drop = FALSE])

  # Strip "Intensity " prefix from column names to match phospho sample names
  colnames(intensity_mat) <- sub("^Intensity ", "", colnames(intensity_mat))

  # Match sample names with phospho data
  common_samples <- intersect(colnames(intensity_mat), phospho_samples)
  if (length(common_samples) == 0) {
    #msg("[.readMaxQuantProteome] Warning: No matching samples between protein and phospho data.")
    return(NULL)
  }

  if (length(common_samples) < length(phospho_samples)) {
    #msg("[.readMaxQuantProteome] Warning: Only ", length(common_samples), " of ",
    #        length(phospho_samples), " phospho samples found in protein data.")
  }

  # Subset to common samples and reorder to match phospho sample order
  intensity_mat <- intensity_mat[, common_samples, drop = FALSE]
  phospho_order <- phospho_samples[phospho_samples %in% common_samples]
  intensity_mat <- intensity_mat[, phospho_order, drop = FALSE]

  # Log2 transform (MaxQuant intensities are linear scale)
  intensity_mat[intensity_mat == 0] <- NA
  intensity_mat <- log2(intensity_mat)

  # Use Protein IDs or Gene names as row identifiers
  if ("Protein IDs" %in% colnames(dat)) {
    protein_ids <- dat$`Protein IDs`
  } else if ("Majority protein IDs" %in% colnames(dat)) {
    protein_ids <- dat$`Majority protein IDs`
  } else {
    protein_ids <- paste0("Protein_", 1:nrow(dat))
  }

  rownames(intensity_mat) <- protein_ids

  # Remove proteins with all NA values
  valid_rows <- rowSums(!is.na(intensity_mat)) > 0
  intensity_mat <- intensity_mat[valid_rows, , drop = FALSE]

  #msg("[.readMaxQuantProteome] Returning protein abundance matrix: ",
  #        nrow(intensity_mat), " proteins × ", ncol(intensity_mat), " samples.")

  return(intensity_mat)
}

#' Collapse phosphosite-level data to protein-level
#'
#' This function collapses multiple phosphosites to one value per protein/gene,
#' using a local collapse function. The site with the maximum absolute
#' statistic is chosen to represent each protein.
# Local fallback for PhosR::phosCollapse (max-abs stat per protein/gene).
phosphoCollapseLocal <- function(mat, id, stat, by = "max") {
  if (by != "max") {
    stop("Only by = 'max' is supported for local phospho collapse.")
  }

  if (is.null(mat) || nrow(mat) == 0) {
    return(mat)
  }

  if (is.null(id) || length(id) != nrow(mat)) {
    stop("id must be a vector with length equal to nrow(mat).")
  }

  if (is.null(stat) || length(stat) != nrow(mat)) {
    stop("stat must be a vector with length equal to nrow(mat).")
  }

  ids <- as.character(id)
  stats <- as.numeric(stat)
  stats[is.na(stats) | is.infinite(stats)] <- -Inf

  idx_by_id <- split(seq_len(nrow(mat)), ids)
  chosen_idx <- vapply(idx_by_id, function(idxs) {
    if (length(idxs) == 1) {
      return(idxs)
    }

    local_stats <- stats[idxs]
    if (all(is.infinite(local_stats))) {
      return(idxs[1])
    }
    idxs[which.max(local_stats)]
  }, integer(1))

  collapsed <- mat[chosen_idx, , drop = FALSE]
  rownames(collapsed) <- names(chosen_idx)

  return(collapsed)
}

#'
#' @return 1 if successful, 0 if failed
CollapsePhosphoToProtein <- function() {
  # Read paramSet and msgSet at the start (like NormalizeData)
  paramSet <- readSet(paramSet, "paramSet");
  msgSet <- readSet(msgSet, "msgSet");

  msg <- "";
  msg("[CollapsePhosphoToProtein] ======================================")
  msg("[CollapsePhosphoToProtein] Starting phosphosite to protein-level collapse")
  msg("[CollapsePhosphoToProtein] ======================================")

  # Load current data
  dataName <- paramSet$dataName;
  msg("[CollapsePhosphoToProtein] Loading dataset: ", dataName)

  dataSet <- readDataset(dataName);
  if (is.null(dataSet)) {
    msg <- "ERROR: Could not load dataset."
    msg("[CollapsePhosphoToProtein] ", msg)
    msgSet$current.msg <- msg;
    saveSet(msgSet, "msgSet");
    return(0)
  }
  msg("[CollapsePhosphoToProtein] Dataset loaded successfully")

  if (is.null(dataSet$data.norm)) {
    msg <- "ERROR: No normalized phosphosite data available for collapsing."
    msg("[CollapsePhosphoToProtein] ", msg)
    msgSet$current.msg <- msg;
    saveSet(msgSet, "msgSet");
    return(0)
  }

  phospho_mat <- dataSet$data.norm;
  row.nms <- rownames(phospho_mat);
  col.nms <- colnames(phospho_mat);

  msg("[CollapsePhosphoToProtein] Input matrix: ", nrow(phospho_mat), " phosphosites × ",
          ncol(phospho_mat), " samples")

  # Check feature info
  if (is.null(dataSet$feature.info)) {
    msg <- "ERROR: Feature information is missing. Cannot extract protein IDs from phosphosites."
    msg("[CollapsePhosphoToProtein] ", msg)
    msgSet$current.msg <- msg;
    saveSet(msgSet, "msgSet");
    return(0)
  }

  feature_info <- dataSet$feature.info;
  msg("[CollapsePhosphoToProtein] Feature info columns: ", paste(colnames(feature_info), collapse=", "))

  # Extract protein IDs/gene names
  msg("[CollapsePhosphoToProtein] Extracting protein IDs from phosphosites...")

  if ("Proteins" %in% colnames(feature_info)) {
    msg("[CollapsePhosphoToProtein] Using 'Proteins' column from feature info")
    protein_ids <- sapply(strsplit(as.character(feature_info$Proteins), ";"), function(x) x[1])
  } else {
    msg("[CollapsePhosphoToProtein] Extracting protein IDs from row names")
    site_ids <- rownames(phospho_mat)
    protein_ids <- sapply(strsplit(site_ids, "_"), function(x) x[1])
  }

  # Check for gene names
  use_gene_names <- FALSE;
  if ("Gene names" %in% colnames(feature_info) && !all(is.na(feature_info$`Gene names`))) {
    gene_names <- feature_info$`Gene names`;
    non_empty_genes <- sum(!is.na(gene_names) & gene_names != "");
    msg("[CollapsePhosphoToProtein] Found ", non_empty_genes, " gene names in feature info")

    if (non_empty_genes > 0) {
      use_gene_names <- TRUE;
      gene_ids <- ifelse(is.na(gene_names) | gene_names == "", protein_ids, gene_names);
      id_vector <- gene_ids;
      msg("[CollapsePhosphoToProtein] Using gene names where available, protein IDs as fallback")
    } else {
      id_vector <- protein_ids;
      msg("[CollapsePhosphoToProtein] No valid gene names found, using protein IDs only")
    }
  } else {
    id_vector <- protein_ids;
    msg("[CollapsePhosphoToProtein] Using protein IDs for collapsing (no gene names column)")
  }

  unique_proteins <- length(unique(id_vector));
  msg("[CollapsePhosphoToProtein] Identified ", unique_proteins, " unique proteins/genes")

  if (unique_proteins == 0) {
    msg <- "ERROR: Could not extract any protein IDs from phosphosite data."
    msg("[CollapsePhosphoToProtein] ", msg)
    msgSet$current.msg <- msg;
    saveSet(msgSet, "msgSet");
    return(0)
  }

  # Show distribution of sites per protein
  sites_per_protein <- table(id_vector);
  msg("[CollapsePhosphoToProtein] Sites per protein: min=", min(sites_per_protein),
          ", max=", max(sites_per_protein),
          ", median=", median(sites_per_protein))

  # Determine statistic for collapsing
  msg("[CollapsePhosphoToProtein] Determining statistic for site selection...")

  if (!is.null(paramSet$sig.mat) && "t" %in% colnames(paramSet$sig.mat)) {
    msg("[CollapsePhosphoToProtein] Found DE test statistics in paramSet$sig.mat")

    common_sites <- intersect(rownames(phospho_mat), rownames(paramSet$sig.mat));
    msg("[CollapsePhosphoToProtein] Matched ", length(common_sites), " sites with DE statistics")

    if (length(common_sites) > 0) {
      stat_vec <- abs(paramSet$sig.mat[common_sites, "t"]);
      names(stat_vec) <- common_sites;

      # Create stat_vec for all sites
      all_stat_vec <- rep(0, nrow(phospho_mat));
      names(all_stat_vec) <- rownames(phospho_mat);
      all_stat_vec[common_sites] <- stat_vec;

      # For sites not in sig.mat, use row variance
      missing_sites <- setdiff(rownames(phospho_mat), common_sites);
      if (length(missing_sites) > 0) {
        msg("[CollapsePhosphoToProtein] Computing variance for ", length(missing_sites),
                " sites without DE statistics")
        missing_stats <- apply(phospho_mat[missing_sites, , drop = FALSE], 1, var, na.rm = TRUE);
        all_stat_vec[missing_sites] <- missing_stats;
      }

      stat_vec <- all_stat_vec;
      msg("[CollapsePhosphoToProtein] Using DE test statistics (t-statistic) for collapsing")
    } else {
      msg("[CollapsePhosphoToProtein] No overlap with DE statistics, using row variance")
      stat_vec <- apply(phospho_mat, 1, var, na.rm = TRUE);
    }
  } else {
    msg("[CollapsePhosphoToProtein] No DE statistics found, using row variance")
    stat_vec <- apply(phospho_mat, 1, var, na.rm = TRUE);
  }

  # Handle NA/Inf in statistics
  na_count <- sum(is.na(stat_vec) | is.infinite(stat_vec));
  if (na_count > 0) {
    msg("[CollapsePhosphoToProtein] Replacing ", na_count, " NA/Inf values with 0")
    stat_vec[is.na(stat_vec) | is.infinite(stat_vec)] <- 0;
  }

  msg("[CollapsePhosphoToProtein] Statistic range: min=", round(min(stat_vec), 4),
          ", max=", round(max(stat_vec), 4))

  # Perform collapse using local implementation (no PhosR dependency)
  msg("[CollapsePhosphoToProtein] ======================================")
  msg("[CollapsePhosphoToProtein] Collapsing phosphosites to protein level...")
  msg("[CollapsePhosphoToProtein] ======================================")

  protein_mat <- tryCatch({
    result <- phosphoCollapseLocal(
      mat = phospho_mat,
      id = id_vector,
      stat = stat_vec,
      by = "max"
    )
    msg("[CollapsePhosphoToProtein] Collapse completed successfully")
    result
  }, error = function(e) {
    msg <- paste0("ERROR during phosphosite collapse: ", e$message);
    msg("[CollapsePhosphoToProtein] ", msg)
    msgSet$current.msg <- msg;
    saveSet(msgSet, "msgSet");
    return(NULL)
  });

  if (is.null(protein_mat)) {
    return(0)
  }

  msg("[CollapsePhosphoToProtein] Collapsed result: ", nrow(protein_mat), " proteins × ",
          ncol(protein_mat), " samples")
  msg("[CollapsePhosphoToProtein] Reduction: ", nrow(phospho_mat), " sites → ",
          nrow(protein_mat), " proteins (",
          round(100 * nrow(protein_mat) / nrow(phospho_mat), 1), "%)")

  # Update dataSet
  msg("[CollapsePhosphoToProtein] Updating dataSet with protein-level data...")
  dataSet$data.norm <- protein_mat;

  # Create new feature info
  msg("[CollapsePhosphoToProtein] Creating protein-level feature info...")
  protein_feature_info <- data.frame(
    Protein = rownames(protein_mat),
    row.names = rownames(protein_mat)
  );

  if (use_gene_names) {
    protein_feature_info$Gene <- rownames(protein_mat);
  }

  dataSet$feature.info <- protein_feature_info;

  # Mark as protein-level
  paramSet$phospho.data.level <- "protein";
  paramSet$phospho.original.sites <- nrow(phospho_mat);
  paramSet$phospho.collapsed.proteins <- nrow(protein_mat);
  msg("[CollapsePhosphoToProtein] Marked data as protein-level in paramSet")

  # Save files (needed for annotation to work)
  msg("[CollapsePhosphoToProtein] Saving protein-level data to disk...")
  ov_qs_save(protein_mat, "data.raw.qs");
  ov_qs_save(protein_mat, "int.mat.qs");
  ov_qs_save(protein_mat, "data.stat.qs");
  fast.write(sanitizeSmallNumbers(protein_mat), file="data_collapsed_protein.csv");
  msg("[CollapsePhosphoToProtein] Saved to data.raw.qs, int.mat.qs, data.stat.qs, data_collapsed_protein.csv")

  # Update dataSet with protein-level data
  dataSet$data.norm <- protein_mat;
  dataSet$type <- "prot";
  dataSet$annotated <- FALSE;

  # Save paramSet before annotation
  saveSet(paramSet, "paramSet");
  msg("[CollapsePhosphoToProtein] Saved paramSet")

  # Perform annotation: convert protein IDs/gene names to Entrez IDs
  msg("[CollapsePhosphoToProtein] ======================================")
  msg("[CollapsePhosphoToProtein] Starting annotation: Protein/Gene → Entrez")
  msg("[CollapsePhosphoToProtein] ======================================")

  # Determine ID type for annotation
  idType <- "uniprot";  # Default to UniProt
  #if (use_gene_names) {
    # If we're using gene names, try symbol first
  #  idType <- "symbol";
  #  msg("[CollapsePhosphoToProtein] Using 'symbol' as ID type for annotation")
  #} else {
    msg("[CollapsePhosphoToProtein] Using 'uniprot' as ID type for annotation")
  #}

  # Get organism from paramSet
  org <- paramSet$data.org;
  if (is.null(org) || org == "NA") {
    org <- "hsa";  # Default to human
    msg("[CollapsePhosphoToProtein] No organism set, defaulting to human (hsa)")
  }
  msg("[CollapsePhosphoToProtein] Organism: ", org)

  # Call PerformDataAnnotInternal to map to Entrez IDs
  tryCatch({
    dataSet <- PerformDataAnnotInternal(
      dataSet = dataSet,
      dataName = dataName,
      org = org,
      dataType = "prot",
      idType = idType,
      lvlOpt = "mean"
    );
    msg("[CollapsePhosphoToProtein] Annotation completed successfully")
    msg("[CollapsePhosphoToProtein] Annotated data dimensions: ",
            nrow(dataSet$data.norm), " features × ", ncol(dataSet$data.norm), " samples")
  }, error = function(e) {
    msg("[CollapsePhosphoToProtein] WARNING: Annotation failed: ", e$message)
    msg("[CollapsePhosphoToProtein] Continuing without annotation (protein IDs will be used as-is)")
  });

  # Build success message
  if (dataSet$annotated) {
    msg <- paste0(
      "Successfully collapsed ", nrow(phospho_mat), " phosphosites to ",
      nrow(protein_mat), " proteins using local collapse, then annotated to ",
      nrow(dataSet$data.norm), " Entrez IDs. Data is ready for enrichment analysis."
    );
  } else {
    msg <- paste0(
      "Successfully collapsed ", nrow(phospho_mat), " phosphosites to ",
      nrow(protein_mat), " proteins using local collapse. ",
      "Note: Annotation to Entrez IDs was not performed. Enrichment may be limited."
    );
  }

  msgSet$current.msg <- msg;
  saveSet(msgSet, "msgSet");
  msg("[CollapsePhosphoToProtein] ", msg)

  msg("[CollapsePhosphoToProtein] ======================================")
  msg("[CollapsePhosphoToProtein] COLLAPSE COMPLETE - Registering dataset")
  msg("[CollapsePhosphoToProtein] ======================================")

  return(RegisterData(dataSet))
}

#' Map phosphosite UniProt IDs to gene symbols for downstream visualization
#'
#' This function keeps phosphosite-level data intact while creating a symbol
#' mapping that can be used in downstream visualizations (volcano, enrichment,
#' heatmap, ridgeline, upset, GSEA). The phosphosite IDs remain as-is, but
#' symbols are extracted and stored for display purposes.
#'
#' @return 1 if successful, 0 if failed
MapPhosphositeUniprotToSymbol <- function() {
  paramSet <- readSet(paramSet, "paramSet")
  msgSet <- readSet(msgSet, "msgSet")

  msg <- ""
  msg("[MapPhosphositeUniprotToSymbol] ======================================")
  msg("[MapPhosphositeUniprotToSymbol] Starting UniProt to Symbol mapping")
  msg("[MapPhosphositeUniprotToSymbol] ======================================")

  # Load current data
  dataName <- paramSet$dataName
  msg("[MapPhosphositeUniprotToSymbol] Loading dataset: ", dataName)

  dataSet <- readDataset(dataName)
  if (is.null(dataSet)) {
    msg <- "ERROR: Could not load dataset."
    msg("[MapPhosphositeUniprotToSymbol] ", msg)
    msgSet$current.msg <- msg
    saveSet(msgSet, "msgSet")
    return(0)
  }
  msg("[MapPhosphositeUniprotToSymbol] Dataset loaded successfully")

  if (is.null(dataSet$data.norm)) {
    msg <- "ERROR: No normalized phosphosite data available."
    msg("[MapPhosphositeUniprotToSymbol] ", msg)
    msgSet$current.msg <- msg
    saveSet(msgSet, "msgSet")
    return(0)
  }

  phospho_mat <- dataSet$data.norm
  phosphosite_ids <- rownames(phospho_mat)

  msg("[MapPhosphositeUniprotToSymbol] Data matrix: ", nrow(phospho_mat),
          " phosphosites × ", ncol(phospho_mat), " samples")
  msg("[MapPhosphositeUniprotToSymbol] Sample phosphosite IDs: ",
          paste(head(phosphosite_ids, 3), collapse=", "))

  # Extract UniProt IDs from phosphosite IDs (e.g., "P12345_S123" -> "P12345")
  msg("[MapPhosphositeUniprotToSymbol] Extracting UniProt IDs from phosphosites...")

  # Store original UniProt IDs WITH isoforms for later reconstruction
  uniprot_ids_with_isoform <- sapply(strsplit(as.character(phosphosite_ids), "_"), function(x) x[1])

  # Extract isoform suffix if present (e.g., "A2AJT9-3" -> "-3")
  isoform_suffix <- ifelse(grepl("-\\d+$", uniprot_ids_with_isoform),
                           sub("^[^-]+(-\\d+)$", "\\1", uniprot_ids_with_isoform),
                           "")

  # Extract site suffix (everything after first underscore, e.g., "_S_123")
  site_suffix <- ifelse(grepl("_", phosphosite_ids),
                        sub("^[^_]+(_.*)", "\\1", phosphosite_ids),
                        "")

  # Strip isoform suffixes for mapping (e.g., "A2AJT9-3" -> "A2AJT9")
  # This is essential for matching against annotation databases
  uniprot_ids <- sub("-\\d+$", "", uniprot_ids_with_isoform)

  msg("[MapPhosphositeUniprotToSymbol] Sample UniProt IDs: ",
          paste(head(uniprot_ids, 3), collapse=", "))

  # Get organism
  org <- paramSet$data.org
  if (is.null(org) || org == "NA") {
    org <- "hsa"  # Default to human
    msg("[MapPhosphositeUniprotToSymbol] No organism set, defaulting to human (hsa)")
  }
  msg("[MapPhosphositeUniprotToSymbol] Organism: ", org)

  # Map UniProt IDs to gene symbols using existing mapping infrastructure
  msg("[MapPhosphositeUniprotToSymbol] Mapping UniProt IDs to gene symbols...")

  tryCatch({
    # Convert UniProt to Entrez first
    entrez_ids <- .doGeneIDMapping(uniprot_ids, "uniprot", paramSet, "vec", keepNA = FALSE)

    # Then convert Entrez to Symbols
    base_symbols <- doEntrez2SymbolMapping(entrez_ids, org, "entrez")

    # Reconstruct display names: SYMBOL-isoform_site_position
    # e.g., "A2AJT9-3_S_78" -> "CXORF23-3_S_78"
    symbols <- base_symbols
    has_symbol <- !is.na(symbols) & symbols != "" & symbols != "NA"

    if (any(has_symbol)) {
      # For successful mappings, reconstruct with isoform and site suffix
      symbols[has_symbol] <- paste0(base_symbols[has_symbol],
                                     isoform_suffix[has_symbol],
                                     site_suffix[has_symbol])
    }

    # For sites without symbol mapping, keep the original phosphosite ID
    no_symbol <- !has_symbol
    symbols[no_symbol] <- phosphosite_ids[no_symbol]

    msg("[MapPhosphositeUniprotToSymbol] Mapped ", sum(!no_symbol), " phosphosites to symbols")
    msg("[MapPhosphositeUniprotToSymbol] ", sum(no_symbol), " phosphosites kept original IDs")
    msg("[MapPhosphositeUniprotToSymbol] Sample symbols: ",
            paste(head(symbols, 5), collapse=", "))

    # Store the mapping in paramSet for downstream use
    phospho_symbol_map <- data.frame(
      phosphosite_id = phosphosite_ids,
      uniprot_id = uniprot_ids,
      entrez_id = entrez_ids,
      symbol = symbols,
      stringsAsFactors = FALSE
    )
    rownames(phospho_symbol_map) <- phosphosite_ids

    # Save to file for downstream visualizations
    saveDataQs(phospho_symbol_map, "phospho_symbol_map.qs", paramSet$anal.type, dataName)
    msg("[MapPhosphositeUniprotToSymbol] Saved mapping to phospho_symbol_map.qs")

    # Also update feature.info if it exists
    if (!is.null(dataSet$feature.info)) {
      # Check if feature.info dimensions match current data
      if (nrow(dataSet$feature.info) == length(phosphosite_ids)) {
        # Dimensions match - safe to update
        dataSet$feature.info$Symbol <- symbols
        dataSet$feature.info$UniProt <- uniprot_ids
        msg("[MapPhosphositeUniprotToSymbol] Updated feature.info with symbols")
      } else {
        # Dimensions don't match - subset feature.info to current data or create new one
        msg("[MapPhosphositeUniprotToSymbol] feature.info has ", nrow(dataSet$feature.info),
                " rows but data has ", length(phosphosite_ids), " rows - subsetting to match")

        # Keep only rows that are in current data
        matching_rows <- rownames(dataSet$feature.info) %in% phosphosite_ids
        if (sum(matching_rows) > 0) {
          dataSet$feature.info <- dataSet$feature.info[matching_rows, , drop = FALSE]
          # Now add symbols in the correct order
          reorder_idx <- match(phosphosite_ids, rownames(dataSet$feature.info))
          dataSet$feature.info <- dataSet$feature.info[reorder_idx, , drop = FALSE]
          dataSet$feature.info$Symbol <- symbols
          dataSet$feature.info$UniProt <- uniprot_ids
          msg("[MapPhosphositeUniprotToSymbol] Updated feature.info after subsetting")
        } else {
          # No matching rows - create new feature.info
          msg("[MapPhosphositeUniprotToSymbol] Creating new feature.info")
          dataSet$feature.info <- data.frame(
            Symbol = symbols,
            UniProt = uniprot_ids,
            row.names = phosphosite_ids,
            stringsAsFactors = FALSE
          )
        }
      }
    }

    # Save updated dataset
    RegisterData(dataSet)

    msg <- paste0(
      "Successfully mapped ", nrow(phospho_mat), " phosphosites to gene symbols. ",
      sum(!no_symbol), " sites have gene symbols, ", sum(no_symbol), " sites use UniProt IDs. ",
      "Phosphosite-level data is preserved for all downstream analyses."
    )

    msgSet$current.msg <- msg
    saveSet(msgSet, "msgSet")
    msg("[MapPhosphositeUniprotToSymbol] ", msg)

    msg("[MapPhosphositeUniprotToSymbol] ======================================")
    msg("[MapPhosphositeUniprotToSymbol] MAPPING COMPLETE")
    msg("[MapPhosphositeUniprotToSymbol] ======================================")

    return(1)

  }, error = function(e) {
    msg <- paste0("ERROR during symbol mapping: ", e$message)
    msg("[MapPhosphositeUniprotToSymbol] ", msg)
    msgSet$current.msg <- msg
    saveSet(msgSet, "msgSet")
    return(0)
  })
}

.readFragPipePhospho <- function(filePath, opts) {
  
  if (!file.exists(filePath)) {
    return(NULL)
  }
  
  msg("[.readFragPipePhospho] Reading FragPipe file: ", filePath)
  dat <- try(data.table::fread(filePath, header = TRUE, check.names = FALSE, data.table = FALSE))
  
  if (inherits(dat, "try-error")) {
    return(NULL)
  }
  
  # --- 1. Basic Filtering ---
  # FragPipe usually handles decoys/contaminants upstream, but we check standard columns if they exist
  if (isTRUE(opts$removeContaminants)) {
     # FragPipe often puts contaminants in "Protein" column with "CONT_" prefix
     if ("Protein" %in% colnames(dat)) {
         dat <- dat[!grepl("^CONT_", dat$Protein), ]
     }
  }
  
  # Filter by Localization Probability
  # FragPipe column is often "MaxPepProb" or "Probability" depending on version
  prob_col <- grep("Probability", colnames(dat), value = TRUE, ignore.case = TRUE)[1]
  if (!is.na(prob_col)) {
      msg("[.readFragPipePhospho] Filtering by localization probability column: ", prob_col)
      dat <- dat[dat[[prob_col]] >= opts$loc.prob, ]
  }

  if (nrow(dat) == 0) return(NULL)

  # --- 2. Extract Intensity Data ---
  # FragPipe intensity columns usually look like "SampleName Intensity" or "SampleName"
  # We look for columns that end in " Intensity" or " Spectral Count" depending on settings.
  # Assuming standard TMT/LFQ intensity columns:
  
  # Strategy: Find columns that are numeric and not metadata
  meta_cols <- c("Index", "Gene", "Protein", "Peptide", "Modified Peptide", "Probability", 
                 "MaxPepProb", "Reference", "Organism", "Entry Name", "Amino Acid", "Position")
  
  # Regex to capture intensity columns. FragPipe often names them "Sample 1 Intensity", "Sample 2 Intensity"
  intensity_cols <- grep(" Intensity$", colnames(dat), value = TRUE)
  
  # If no " Intensity" suffix found, this might be a custom table. 
  # Fallback: Exclude known metadata columns
  if (length(intensity_cols) == 0) {
     intensity_cols <- setdiff(colnames(dat), meta_cols)
     # Keep only numeric columns
     intensity_cols <- intensity_cols[sapply(dat[, intensity_cols], is.numeric)]
  }
  
  if (length(intensity_cols) == 0) {
      msg("[.readFragPipePhospho] Error: Could not identify intensity columns.")
      return(NULL)
  }

  msg("[.readFragPipePhospho] Found ", length(intensity_cols), " intensity columns.")
  
  # Extract Data Matrix
  int.mat <- as.matrix(dat[, intensity_cols, drop = FALSE])
  
  # Clean column names: Remove " Intensity" suffix to match Metadata sample names
  colnames(int.mat) <- sub(" Intensity$", "", colnames(int.mat))
  colnames(int.mat) <- sub(" MaxLFQ Intensity$", "", colnames(int.mat)) # Handle MaxLFQ specific suffix

  # Zeros to NA
  int.mat[int.mat == 0] <- NA

  # --- 3. Construct Site IDs ---
  # Your downstream code needs "Protein_Residue_Position" (e.g., P12345_S_10)
  
  # FragPipe standard columns:
  # Protein: "P00533"
  # Amino Acid: "S" (sometimes called "Code")
  # Position: "1056"
  
  # Handle variations in header names
  prot_col <- if("Protein" %in% colnames(dat)) "Protein" else grep("Protein", colnames(dat), value=T)[1]
  aa_col <- if("Amino Acid" %in% colnames(dat)) "Amino Acid" else if("Code" %in% colnames(dat)) "Code" else NA
  pos_col <- if("Position" %in% colnames(dat)) "Position" else NA
  
  # If Position/AA columns are missing (common in some aggregated reports), try to parse "Index"
  # Example Index: P00533_1056_S
  if (is.na(aa_col) || is.na(pos_col)) {
      if ("Index" %in% colnames(dat)) {
          msg("[.readFragPipePhospho] Using 'Index' column to parse site info.")
          # Split FragPipe Index: Protein_Position_AA (Note: FragPipe often uses Prot_Pos_AA format)
          split_index <- strsplit(as.character(dat$Index), "_")
          
          # We need to robustly map this to Protein, AA, Pos
          # Assuming standard FragPipe Index
          protein_ids <- sapply(split_index, `[`, 1) 
          positions <- sapply(split_index, `[`, 2)
          residues <- sapply(split_index, `[`, 3)
      } else {
          return(NULL) # Cannot build IDs
      }
  } else {
      protein_ids <- dat[[prot_col]]
      residues <- dat[[aa_col]]
      positions <- dat[[pos_col]]
  }
  
  # Handle multiple proteins (P123;P456) -> Take first
  protein_ids <- sapply(strsplit(as.character(protein_ids), ";"), `[`, 1)
  
  # Create the ID
  site_ids <- paste(protein_ids, residues, positions, sep = "_")
  
  # --- 4. Handle Duplicates ---
  # Just like your MaxQuant reader, if site_ids are duplicated, aggregate by mean
  if (any(duplicated(site_ids))) {
      msg("[.readFragPipePhospho] Aggregating duplicate sites.")
      agg_df <- as.data.frame(int.mat)
      agg_df$site_id <- site_ids
      
      # Aggregate intensities
      final_intens <- aggregate(. ~ site_id, data = agg_df, FUN = mean, na.rm = TRUE)
      
      rownames(final_intens) <- final_intens$site_id
      final_intens$site_id <- NULL
      int.mat <- as.matrix(final_intens)
      
      # Re-align site_ids after aggregation
      site_ids <- rownames(int.mat)
  } else {
      rownames(int.mat) <- site_ids
  }
  
  # --- 5. Feature Information ---
  # Create the feature info table required by your `CollapsePhosphoToProtein` function
  # It needs: "Proteins", "Gene names", "Amino acid", "Position"
  
  feature_info <- data.frame(
      id = site_ids,
      Proteins = sapply(strsplit(site_ids, "_"), `[`, 1),
      "Amino acid" = sapply(strsplit(site_ids, "_"), `[`, 2),
      Position = sapply(strsplit(site_ids, "_"), `[`, 3),
      check.names = FALSE
  )
  
  # Add Gene Names if available
  if ("Gene" %in% colnames(dat)) {
      # We need to match the gene names to the aggregated sites. 
      # Simplest way: Create a lookup from original data
      gene_map <- setNames(dat$Gene, paste(protein_ids, residues, positions, sep="_"))
      feature_info$`Gene names` <- gene_map[site_ids]
  } else {
      feature_info$`Gene names` <- NA
  }
  
  # Add Localization Prob if available (useful for downstream filtering)
  if (!is.na(prob_col)) {
      prob_map <- setNames(dat[[prob_col]], paste(protein_ids, residues, positions, sep="_"))
      feature_info$`Localization prob` <- prob_map[site_ids]
  }
  
  rownames(feature_info) <- site_ids
  
  # --- 6. Return Data Object ---
  # Create data_orig for metadata matching
  data_orig_df <- cbind(data.frame(SiteID = rownames(int.mat)), as.data.frame(int.mat))

  return(list(
    data = int.mat,
    data_orig = data_orig_df,
    type = "phospho",
    format = "fragpipe",
    meta.info = NULL,
    feature.info = feature_info
  ))
}


#' .readDIANNPhospho
#'
#' Parse DIA-NN phospho report.tsv file
#'
#' Expected columns:
#'   - Run (sample name)
#'   - Protein.Group (protein ID)
#'   - Modified.Sequence (peptide with modifications, e.g. PEPT(Phospho)IDE)
#'   - Precursor.Id (unique peptide+charge ID)
#'   - Stripped.Sequence (unmodified peptide)
#'   - Genes (gene names)
#'   - Precursor.Quantity or Precursor.Normalised (intensity)
#'   - Q.Value (peptide FDR)
#'   - Global.PG.Q.Value (protein FDR)
#'   - PTM.Site.Confidence (localization probability, 0-1)
#'
#' @param filePath Path to DIA-NN report file
#' @param opts List with filtering options:
#'   - loc.prob: Minimum PTM.Site.Confidence (default 0.75)
#'   - peptide.fdr: Maximum Q.Value (default 0.01)
#'   - protein.fdr: Maximum Global.PG.Q.Value (default 0.01)
#'   - removeContaminants: Remove contaminants (default TRUE)
#'   - removeDecoys: Remove decoys (default TRUE)
#'
#' @return List with data, data_orig, type, format, meta.info, feature.info
#'
.readDIANNPhospho <- function(filePath, opts) {

  if (!file.exists(filePath)) {
    msg("[.readDIANNPhospho] Error: File not found at ", filePath)
    return(NULL)
  }

  msg("[.readDIANNPhospho] Reading DIA-NN file: ", filePath)

  # Read DIA-NN report (long format: one row per precursor per run)
  dat <- try(data.table::fread(filePath, header = TRUE, check.names = FALSE, data.table = FALSE))

  if (inherits(dat, "try-error")) {
    msg("[.readDIANNPhospho] Error: Failed to read the file.")
    return(NULL)
  }

  # --- 1. Quality Filtering ---

  # Filter by protein FDR (Global.PG.Q.Value)
  protein_fdr_threshold <- if (!is.null(opts$protein.fdr)) opts$protein.fdr else 0.01
  if ("Global.PG.Q.Value" %in% colnames(dat)) {
    msg("[.readDIANNPhospho] Filtering by Global.PG.Q.Value < ", protein_fdr_threshold)
    dat <- dat[dat$Global.PG.Q.Value < protein_fdr_threshold, ]
  }

  # Filter by peptide FDR (Q.Value)
  peptide_fdr_threshold <- if (!is.null(opts$peptide.fdr)) opts$peptide.fdr else 0.01
  if ("Q.Value" %in% colnames(dat)) {
    msg("[.readDIANNPhospho] Filtering by Q.Value < ", peptide_fdr_threshold)
    dat <- dat[dat$Q.Value < peptide_fdr_threshold, ]
  }

  # Filter by PTM localization confidence
  loc_prob_threshold <- opts$loc.prob
  if ("PTM.Site.Confidence" %in% colnames(dat)) {
    msg("[.readDIANNPhospho] Filtering by PTM.Site.Confidence >= ", loc_prob_threshold)
    dat <- dat[!is.na(dat$PTM.Site.Confidence) & dat$PTM.Site.Confidence >= loc_prob_threshold, ]
  }

  # Remove contaminants (DIA-NN marks them as "cont_" in Protein.Group)
  if (isTRUE(opts$removeContaminants) && "Protein.Group" %in% colnames(dat)) {
    dat <- dat[!grepl("^cont_|^CONT_", dat$Protein.Group, ignore.case = TRUE), ]
  }

  # Remove decoys (typically marked with "DECOY_" or "REV_")
  if (isTRUE(opts$removeDecoys) && "Protein.Group" %in% colnames(dat)) {
    dat <- dat[!grepl("^DECOY_|^REV_", dat$Protein.Group, ignore.case = TRUE), ]
  }

  if (nrow(dat) == 0) {
    msg("[.readDIANNPhospho] Error: No data remains after filtering.")
    return(NULL)
  }

  msg("[.readDIANNPhospho] After filtering: ", nrow(dat), " precursors")

  # --- 2. Extract Phosphosites from Modified.Sequence ---
  # Modified.Sequence example: "PEPT(Phospho (STY))IDE" or "PEPTIDE(UniMod:21)"
  # We need to parse this to get: Protein_Residue_Position format

  if (!"Modified.Sequence" %in% colnames(dat)) {
    msg("[.readDIANNPhospho] Error: Missing Modified.Sequence column.")
    return(NULL)
  }

  # Parse modified sequences to extract phosphosites
  # This is the key function that extracts P12345_S_123 from modified peptide
  # NOTE: This function can return MORE rows than input if peptides have multiple phosphosites
  # The parser returns a row_idx column indicating which input row each phosphosite came from
  phospho_sites <- .parseDIANNModifiedSequence(
    dat$Modified.Sequence,
    dat$Protein.Group,
    dat$Stripped.Sequence
  )

  if (is.null(phospho_sites) || nrow(phospho_sites) == 0) {
    msg("[.readDIANNPhospho] Error: Could not parse phosphosites from Modified.Sequence.")
    return(NULL)
  }

  # Expand the original data to match the phosphosites
  # We need to replicate rows when a peptide has multiple phosphosites
  # Use the row_idx from the parser to know which precursor each phosphosite came from
  dat_expanded <- dat[phospho_sites$row_idx, ]

  # Now add site info to the expanded data
  dat_expanded$site_id <- phospho_sites$site_id
  dat_expanded$protein_id <- phospho_sites$protein_id
  dat_expanded$residue <- phospho_sites$residue
  dat_expanded$position <- phospho_sites$position

  # Remove rows without valid phosphosites
  dat_expanded <- dat_expanded[!is.na(dat_expanded$site_id), ]

  if (nrow(dat_expanded) == 0) {
    msg("[.readDIANNPhospho] Error: No phosphosites found in Modified.Sequence column.")
    return(NULL)
  }

  msg("[.readDIANNPhospho] Found ", length(unique(dat_expanded$site_id)), " unique phosphosites")

  # --- 3. Pivot to Wide Format (Sites × Samples) ---
  # DIA-NN is long format (one row per precursor per run)
  # We need to aggregate to wide format for downstream analysis

  # Use Precursor.Normalised if available, otherwise Precursor.Quantity
  quant_col <- if ("Precursor.Normalised" %in% colnames(dat_expanded)) {
    "Precursor.Normalised"
  } else if ("Precursor.Quantity" %in% colnames(dat_expanded)) {
    "Precursor.Quantity"
  } else {
    msg("[.readDIANNPhospho] Error: Neither Precursor.Normalised nor Precursor.Quantity found.")
    return(NULL)
  }

  msg("[.readDIANNPhospho] Using quantification column: ", quant_col)

  # Aggregate: For each site_id + Run combination, sum intensities across precursors
  # (multiple precursors can have the same phosphosite)
  agg_data <- aggregate(
    as.formula(paste(quant_col, "~ site_id + Run")),
    data = dat_expanded,
    FUN = sum,
    na.rm = TRUE
  )

  # Pivot to wide format
  library(tidyr)
  wide_data <- tidyr::pivot_wider(
    agg_data,
    names_from = "Run",
    values_from = quant_col,
    values_fill = 0  # DIA-NN: missing = not detected
  )

  # Convert to matrix
  site_ids <- wide_data$site_id
  int.mat <- as.matrix(wide_data[, -1])
  rownames(int.mat) <- site_ids

  # Treat zeros as missing (proteomics convention)
  int.mat[int.mat == 0] <- NA

  msg("[.readDIANNPhospho] Final matrix: ", nrow(int.mat), " sites × ", ncol(int.mat), " samples")

  # --- 4. Build Feature Info Table ---
  # Extract metadata for each phosphosite

  # Get unique site info (take first occurrence of each site)
  site_info <- dat_expanded[!duplicated(dat_expanded$site_id), c("site_id", "protein_id", "residue", "position")]

  # Add gene names
  if ("Genes" %in% colnames(dat_expanded)) {
    gene_map <- setNames(dat_expanded$Genes, dat_expanded$site_id)
    site_info$gene_name <- gene_map[site_info$site_id]
  } else {
    site_info$gene_name <- NA
  }

  # Add PTM confidence (average across precursors if multiple)
  if ("PTM.Site.Confidence" %in% colnames(dat_expanded)) {
    conf_map <- aggregate(
      PTM.Site.Confidence ~ site_id,
      data = dat_expanded,
      FUN = mean,
      na.rm = TRUE
    )
    site_info$localization_prob <- conf_map$PTM.Site.Confidence[match(site_info$site_id, conf_map$site_id)]
  }

  # Format as expected by downstream code
  feature_info <- data.frame(
    id = site_info$site_id,
    Proteins = site_info$protein_id,
    "Amino acid" = site_info$residue,
    Position = site_info$position,
    "Gene names" = site_info$gene_name,
    "Localization prob" = if ("localization_prob" %in% colnames(site_info)) site_info$localization_prob else NA,
    check.names = FALSE
  )

  rownames(feature_info) <- feature_info$id

  # Align feature_info with int.mat rownames
  feature_info <- feature_info[rownames(int.mat), ]

  # --- 5. Create data_orig for Metadata Matching ---
  data_orig_df <- cbind(data.frame(SiteID = rownames(int.mat)), as.data.frame(int.mat))

  # --- 6. Return Data Object ---
  return(list(
    data = int.mat,
    data_orig = data_orig_df,
    type = "phospho",
    format = "diann",
    meta.info = NULL,
    feature.info = feature_info
  ))
}


#' .parseDIANNModifiedSequence
#'
#' Parse DIA-NN Modified.Sequence to extract phosphosite positions
#'
#' DIA-NN format examples:
#'   "PEPT(Phospho (STY))IDE"  -> Phospho on T at position 4
#'   "PEPS(UniMod:21)EQRES(UniMod:21)" -> Two phosphosites (UniMod:21 = Phospho)
#'   "M(Oxidation (M))PEPTIDE" -> Oxidation (not phospho, skip)
#'   "(UniMod:1)PEPS(UniMod:21)TIDE" -> N-term acetyl + phospho
#'
#' @param mod_seq Vector of modified sequences
#' @param protein_ids Vector of protein IDs
#' @param stripped_seq Vector of stripped sequences (no modifications)
#'
#' @return Data frame with columns: site_id, protein_id, residue, position
#'
.parseDIANNModifiedSequence <- function(mod_seq, protein_ids, stripped_seq) {

  results <- list()

  for (i in seq_along(mod_seq)) {
    seq <- mod_seq[i]
    prot <- as.character(protein_ids[i])

    # Take first protein if semicolon-separated list
    prot <- strsplit(prot, ";")[[1]][1]

    # Find all phosphorylation sites
    # Support two notation styles:
    # 1. Named format: S(Phospho (STY)), T(Phospho (STY)), Y(Phospho (STY))
    # 2. UniMod format: S(UniMod:21), T(UniMod:21), Y(UniMod:21)
    #    UniMod:21 = Phosphorylation

    # Pattern 1: Named phospho notation
    phospho_pattern1 <- "([STY])\\(Phospho[^)]*\\)"
    # Pattern 2: UniMod:21 notation (phosphorylation)
    phospho_pattern2 <- "([STY])\\(UniMod:21\\)"

    # Try both patterns
    matches1 <- gregexpr(phospho_pattern1, seq, perl = TRUE)[[1]]
    matches2 <- gregexpr(phospho_pattern2, seq, perl = TRUE)[[1]]

    # Combine matches from both patterns
    all_matches <- c()
    if (matches1[1] != -1) {
      all_matches <- c(all_matches, matches1)
    }
    if (matches2[1] != -1) {
      all_matches <- c(all_matches, matches2)
    }

    if (length(all_matches) == 0) {
      # No phosphosites found
      next
    }

    # Remove duplicates and sort
    all_matches <- sort(unique(all_matches))

    # Extract each phosphosite
    for (match_start in all_matches) {
      # Extract the residue (S, T, or Y)
      residue <- substr(seq, match_start, match_start)

      # Calculate position in stripped sequence
      # Remove modification annotations to count position
      seq_before_mod <- substr(seq, 1, match_start)

      # Remove all modification annotations from sequence before this point
      clean_seq_before <- gsub("\\([^)]+\\)", "", seq_before_mod)

      position <- nchar(clean_seq_before)

      # Create site ID
      site_id <- paste(prot, residue, position, sep = "_")

      results[[length(results) + 1]] <- data.frame(
        row_idx = i,  # Track which input row this came from
        site_id = site_id,
        protein_id = prot,
        residue = residue,
        position = position,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(results) == 0) {
    return(NULL)
  }

  # Combine all results
  do.call(rbind, results)
}
