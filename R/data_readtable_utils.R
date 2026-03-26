##################################################
## R script for ProteoAnalyst
## Description: functions for reading data table 
##
## Authors: 
## Jeff Xia, jeff.xia@mcgill.ca
## Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

.safe_numeric_column <- function(x) {
  if (inherits(x, "integer64") || is.numeric(x)) {
    return(as.numeric(x))
  }
  suppressWarnings(as.numeric(as.character(x)))
}

.safe_numeric_matrix <- function(df) {
  df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
  out <- lapply(df, .safe_numeric_column)
  out <- as.data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
  as.matrix(out)
}

#' Read Tabular Expression Data and Metadata
#'
#' This function reads tabular expression data along with metadata and processes the data.
#'
#' @param fileName A character string specifying the name of the expression data file.
#' @param metafileName A character string specifying the name of the metadata file, ignore if metaContain = true
#' @param metaContain A logical value indicating whether metadata is contained in the data file if false, metadata is in separate file
#' @param oneDataAnalType The type of analysis to perform on the one-data setup.
#' @param path The path to the files if they are located in a different directory.
#'
#' @return A processed dataset object containing expression data and metadata.
#'
#' @author Guangyan Zhou \email{guangyan.zhou@mail.mcgill.ca}
#' @details Additional details about the function, if needed.
#'
#' @examples
#' \dontrun{
#' ReadTabExpressData(fileName = "expression_data.csv", metafileName = "metadata.csv",
#'                    metaContain = TRUE, oneDataAnalType = "default", path = "")
#' }
#'
#' @export
#' @license MIT License

ReadTabExpressData <- function(fileName, metafileName="",metaContain="true",oneDataAnalType="default", path="", dataFormat="default") {
  #msg("[ReadTabExpressData] start file=", fileName, " metaContain=", metaContain, " path=", path, " dataFormat=", dataFormat)
  is.maxquant <- dataFormat == "maxquant" || grepl("proteinGroups", fileName, ignore.case = TRUE)
  is.fragpipe <- dataFormat == "fragpipe" || grepl("fragpipe", fileName, ignore.case = TRUE) || grepl("fragpipe", metafileName, ignore.case = TRUE)
  is.diann <- dataFormat == "diann" || grepl("diann", fileName, ignore.case = TRUE) || grepl("diann", metafileName, ignore.case = TRUE)
  is.spectronaut <- dataFormat == "spectronaut" || grepl("spectronaut", fileName, ignore.case = TRUE)

  if (is.maxquant) {
    #msg("[ReadTabExpressData] MaxQuant format detected")
    mx.opts <- getOption("pa.maxquant.opts", list(quantType = "lfq", removeContaminants = TRUE, removeDecoys = TRUE))
    # Route to appropriate reader based on quantType
    qtype <- if (!is.null(mx.opts$quantType)) tolower(mx.opts$quantType) else "lfq"
    if (qtype == "peptide") {
      # Use peptide-level reader for peptides.txt
      #msg("[ReadTabExpressData] Using peptide-level reader for MaxQuant (quantType=", qtype, ")")
      dataSet <- .readMaxQuantPeptides(paste0(path, fileName), mx.opts)
    } else {
      # Use protein-level reader for proteinGroups.txt (lfq, ibaq, intensity)
      #msg("[ReadTabExpressData] Using protein-level reader for MaxQuant (quantType=", qtype, ")")
      dataSet <- .readMaxQuantProteinGroups(paste0(path, fileName), mx.opts)
    }
    if (is.null(dataSet)) {
      #msg("[ReadTabExpressData] MaxQuant reader failed; falling back to standard reader")
      dataSet <- .readTabData(paste0(path, fileName))
    } else {
      #msg("[ReadTabExpressData] MaxQuant data loaded with meta.info rows=", nrow(dataSet$meta.info$meta.info))
      metaContain <- "false"  # metadata loaded by MaxQuant reader
    }
  } else if (is.fragpipe) {
  fp.opts <- getOption("pa.fragpipe.opts", list(quantType = "protein_maxlfq", removeContaminants = TRUE, minProb = NA, minPeptides = NA))
  # Route to appropriate reader based on quantType
  qtype <- if (!is.null(fp.opts$quantType)) fp.opts$quantType else "protein_maxlfq"
  if (qtype == "peptide") {
    # Use peptide-level reader for combined_peptide.tsv or combined_ion.tsv
    #msg("[ReadTabExpressData] Using peptide-level reader for FragPipe (quantType=", qtype, ")")
    dataSet <- .readFragpipePeptide(paste0(path, fileName), paste0(path, metafileName), fp.opts)
  } else {
    # Use protein-level reader for combined_protein.tsv
    #msg("[ReadTabExpressData] Using protein-level reader for FragPipe (quantType=", qtype, ")")
    dataSet <- .readFragpipeLfq(paste0(path, fileName), paste0(path, metafileName), fp.opts)
  }
  metaContain <- "false"  # metadata loaded separately
  } else if (is.diann) {
    di.opts <- getOption("pa.diann.opts", list(fileType = "protein_matrix", qvalueFilter = TRUE))
    dataSet <- .readDiannReport(paste0(path, fileName), paste0(path, metafileName), di.opts)
    metaContain <- "false"  # metadata loaded separately
  } else if (is.spectronaut) {
    spec.opts <- getOption("pa.spectronaut.opts", list(inputType = "protein"))
    dataSet <- .readSpectronaut(paste0(path, fileName), spec.opts)
    metaContain <- "false"  # metadata inferred from sample names
  } else {
    dataSet <- .readTabData(paste0(path, fileName));
  }
  if(is.null(dataSet)){
    return(0);
  }  
  datOrig <- dataSet$data_orig;

  row.num <- nrow(datOrig);
  col.num <- ncol(datOrig);
   if(row.num > 100){
       row.num <- 100;
   }
   if(col.num > 10){
        col.num <- 10;
   }

   write.csv(sanitizeSmallNumbers(datOrig[1:row.num, 1:col.num]), file="raw_dataview.csv");  

  use.external.meta <- (!is.null(metafileName) && nchar(metafileName) > 0)

  if (use.external.meta) {
    #msg("[ReadTabExpressData] external metadata provided; loading meta file (metaContain=", metaContain, ")")
    meta.info <- .readMetaData(metafileName, dataSet$data_orig, metaContain)
  } else if (!is.null(dataSet$meta.info)) {
    #msg("[ReadTabExpressData] using meta.info injected by reader; metaContain set to false")
    meta.info <- dataSet$meta.info;
    metaContain <- "false";
  } else if (metaContain == "true") {
    #msg("[ReadTabExpressData] calling .readMetaData with metaContain=", metaContain, " metaFile=", metafileName)
    meta.info <- .readMetaData(metafileName,dataSet$data_orig,metaContain);
  } else {
    #msg("[ReadTabExpressData] no meta provided; synthesizing metadata from column names.")
    meta.info <- .synthesizeMetaFromRuns(colnames(dataSet$data));
  }
  if (is.null(meta.info) || is.null(meta.info$meta.info)) {
    if (!is.null(meta.info) && is.data.frame(meta.info)) {
      # Wrap plain data.frame metadata into expected list structure
      meta.df <- meta.info
      disc.inx <- setNames(rep(TRUE, ncol(meta.df)), colnames(meta.df))
      cont.inx <- setNames(rep(FALSE, ncol(meta.df)), colnames(meta.df))
      meta.info <- list(meta.info = meta.df, disc.inx = disc.inx, cont.inx = cont.inx)
    }
  }
  if (is.null(meta.info) || is.null(meta.info$meta.info)) {
    #msg("[ReadTabExpressData][WARN] meta.info missing; synthesizing from run/sample names.")
    meta.info <- .synthesizeMetaFromRuns(colnames(dataSet$data));
  }

  # Normalize metadata slots so downstream Java/JSF always see data frames and index vectors
  meta.df <- as.data.frame(meta.info$meta.info, stringsAsFactors = FALSE);
  meta.df[] <- lapply(meta.df, function(col) { factor(col) });
  rownames(meta.df) <- rownames(meta.info$meta.info);
  meta.info$meta.info <- meta.df;
  dataSet$meta.info <- meta.df;
  dataSet$cont.inx <- meta.info$cont.inx;
  dataSet$disc.inx <- meta.info$disc.inx;
  
  # Persist a legacy-friendly matrix for in-house processing (avoid MSstats normalization by default)
  if (dataSet$type == "prot") {
    prot.mat <- dataSet$data;
    if (!is.null(prot.mat)) {
      qs::qsave(prot.mat, "int.mat.qs");
      qs::qsave(prot.mat, "data.raw.qs");
      fast.write(sanitizeSmallNumbers(prot.mat), file="data_original.csv");
    }
    # NOTE: pepcount will be saved later after data filtering/alignment
  } else {
    qs::qsave(dataSet$data, "data.raw.qs");
    fast.write(sanitizeSmallNumbers(dataSet$data), file="data_original.csv");
  }
  
  # Check metadata to determine if it resembles a "dose" pattern with at least 3 replicates per group
  metadata <- meta.info$meta.info
  if(is.null(oneDataAnalType) || oneDataAnalType == "") {
    if(is.numeric(metadata[, 1]) && length(unique(metadata[, 1])) > 1) {
      # Check if values consistently increase without requiring a fixed increment size
      sorted_values <- sort(unique(metadata[, 1]))
      consistent_increase <- all(diff(sorted_values) > 0)
      
      # Check if each unique value has at least 3 replicates
      replicate_counts <- table(metadata[, 1])
      sufficient_replicates <- all(replicate_counts >= 3)
      
      # If values consistently increase and each group has at least 3 replicates, assume dose
      if(consistent_increase && sufficient_replicates) {
        oneDataAnalType <- "dose"
      } else {
        oneDataAnalType <- "default"
      }
    } else {
      oneDataAnalType <- "default"
    }
  }
  
  msgSet <- readSet(msgSet, "msgSet");
  paramSet <- readSet(paramSet, "paramSet");
  paramSet$data.format <- if (!is.null(dataSet$format)) dataSet$format else "text";
  paramSet$isMetaContain <- metaContain
  paramSet$oneDataAnalType <- oneDataAnalType;
  
  # clean memory
  dataSet$data_orig <- NULL; 

  # rename data to data.orig
  int.mat <- dataSet$data;
  matched.cols <- colnames(int.mat)[colnames(int.mat) %in% rownames(meta.info$meta.info)]
  int.mat <- int.mat[, matched.cols, drop = FALSE]
  if (ncol(int.mat) == 0) {
    #msg("[ReadTabExpressData][ERROR] No columns matched between expression matrix and metadata runs.")
    msgSet$current.msg <- "No samples matched between expression matrix and metadata.";
    saveSet(msgSet, "msgSet");
    return(0);
  }
  # align to metadata order
  int.mat <- int.mat[, match(rownames(meta.info$meta.info), colnames(int.mat)), drop = FALSE]
  #msg("[ReadTabExpressData] aligned intensity dims=", paste(dim(int.mat), collapse = "x"),
  #        " meta rows=", nrow(meta.info$meta.info), " meta cols=", ncol(meta.info$meta.info))
  dataSet$data <- NULL;
  dataSet$name <- fileName;
  if (nrow(int.mat) == 0 || ncol(int.mat) == 0) {
    msgSet$current.msg <- "Expression matrix is empty after aligning to metadata.";
    saveSet(msgSet, "msgSet");
    return(0);
  }
  qs::qsave(int.mat, "int.mat.qs");
  msg <- paste("a total of ", ncol(int.mat), " samples and ", nrow(int.mat), " features were found");
  # remove NA, null
  # OPTIMIZED: Use rowSums instead of apply for 60-100x speedup
  row.nas <- rowSums(is.na(int.mat)|is.null(int.mat));
  #good.inx<- row.nas/ncol(int.mat) < 0.5;
  #if(sum(!good.inx) > 0){
  #  int.mat <- int.mat[good.inx,];
  #  msg <- c(msg, paste("removed ", sum(!good.inx), " features with over 50% missing values"));
  #}
  
  finite_vals <- int.mat[is.finite(int.mat)];
  if (length(finite_vals) == 0) {
    msgSet$current.msg <- "All intensity values are missing or non-finite after filtering.";
    saveSet(msgSet, "msgSet");
    return(0);
  }
  minVal <- min(finite_vals, na.rm=TRUE);
  na.inx <- is.na(int.mat);

  # Store original missing value count BEFORE imputation for later reporting
  original.missing.count <- sum(na.inx);
  paramSet$original.missing.count <- original.missing.count;
  saveSet(paramSet, "paramSet");

  # Save original data WITH missing values for QC plots (before imputation)
  qs::qsave(int.mat, "data.with.missing.qs");

  if(sum(na.inx) > 0){
    int.mat[na.inx] <- minVal/2;
    # msg <- c(msg, "the remaining", sum(na.inx), "missing variables were replaced with data min");
  }
  if (nrow(int.mat) == 0 || ncol(int.mat) == 0) {
    msgSet$current.msg <- "No features or samples remain after filtering missing values.";
    saveSet(msgSet, "msgSet");
    return(0);
  }
  #msg("[ReadTabExpressData] dims after pruning missing values: ", paste(dim(int.mat), collapse = "x"))
  msgSet$current.msg <- paste(msg, collapse="; ");
  #res <- RemoveDuplicates(int.mat, "mean", quiet=T, paramSet, msgSet);Gene-level summarization
  data.proc <- int.mat #res[[1]];
  #msgSet <- res[[2]];
  paramSet$smpl.num <- ncol(data.proc);
  
  metadata <- meta.info$meta.info;
  dataSet$meta.info <- metadata;
  if(oneDataAnalType == "dose"){
    # re-order everything numerically by dose
    dose <- as.numeric(gsub(".*_", "", as.character(metadata[,1])))
    int.mat <- int.mat[ ,order(dose)]
    meta.reorder <- as.data.frame(metadata[order(dose),])
    colnames(meta.reorder) <- colnames(metadata)
    rownames(meta.reorder) <- rownames(meta.info$meta.info)[order(dose)]
    dataSet$meta.info <- meta.reorder
    
    # re-level the factor to be numeric instead of alphabetic
    dataSet$meta.info[,1] <- format(as.numeric(as.character(dataSet$meta.info[,1])), scientific = FALSE) # remove scientific notation
    dataSet$meta.info[,1] <- gsub(" ", "", dataSet$meta.info[,1])
    dataSet$meta.info[,1] <- factor(dataSet$meta.info[,1], levels = unique(dataSet$meta.info[,1]))
    
    # rename data to data.orig
    data.proc <- int.mat;
    paramSet$dataSet$meta.info <- dataSet$meta.info;
    dataSet$cls <- dataSet$meta.info[,1];
    dataSet$data <- NULL;
    dataSet$listData <- FALSE;
    
    dataSet$imgSet <- list();
    dataSet$reportSummary <- list();
    
  }
  
  # save processed data for download user option
  data.proc <- sanitizeSmallNumbers(data.proc);
  fast.write(sanitizeSmallNumbers(data.proc), file="data_original.csv");
  qs::qsave(data.proc, "data.raw.qs");
  dataSet$data.norm  <- data.proc;
  metaInx = which(rownames(dataSet$meta.info) %in% colnames(data.proc))
  
  paramSet$dataSet <- list();
  meta.types <- rep("disc", ncol(dataSet$meta.info));
  meta.types[meta.info$cont.inx] <- "cont";
  names(meta.types) <- colnames(dataSet$meta.info);

  paramSet$dataSet$meta.types <- meta.types;
  paramSet$dataSet$meta.info <- dataSet$metaOrig <- dataSet$meta.info[metaInx,,drop=F]
  paramSet$dataSet$disc.inx <- dataSet$disc.inx <-dataSet$disc.inx.orig <- meta.info$disc.inx
  paramSet$dataSet$cont.inx <- dataSet$cont.inx <-dataSet$cont.inx.orig  <- meta.info$cont.inx
  
  meta.types <- rep("disc", ncol(dataSet$meta.info));
  meta.types[meta.info$cont.inx] <- "cont";
  names(meta.types) <- colnames(dataSet$meta.info);
  dataSet$meta.types <-meta.types;
  paramSet$anal.type <- "onedata";
  if (is.null(paramSet$data.idType) || paramSet$data.idType == "") {
    paramSet$data.idType <- "NA";
  }
  paramSet$partialToBeSaved <- c(paramSet$partialToBeSaved, fileName);
  paramSet$jsonNms$dataName <- fileName;
  paramSet$dataName <- fileName;

  # Build an MSstats-ready long-format table for downstream proteomics processing
  msstats.input <- .buildMSstatsInput(int.mat = data.proc,
                                      metadata = dataSet$meta.info,
                                      fileName = fileName);
  if (!is.null(msstats.input)) {
    qs::qsave(msstats.input, "msstats_input.qs");
    fast.write(msstats.input, file = "msstats_input.csv");
    paramSet$dataSet$msstats.input <- "msstats_input.qs";
  }

  # Persist peptide-to-protein map for peptide-level data summarization (shadow save for Arrow)
  if (!is.null(dataSet$prot.map)) {
    shadow_save(dataSet$prot.map, "peptide_to_protein_map.qs");
    #msg("[ReadTabExpressData] Saved peptide-to-protein map: ", nrow(dataSet$prot.map), " peptides mapped")
  }

  # Save peptide/PSM counts aligned with final processed data (for DEqMS)
  if (!is.null(dataSet$pepcount) && dataSet$type == "prot") {
    #msg("[ReadTabExpressData] Aligning and saving pepcount for DEqMS")
    # Align pepcount with final data.proc rownames
    pepcount_aligned <- dataSet$pepcount[rownames(data.proc)]
    # Handle any missing values (proteins that were filtered out)
    pepcount_aligned[is.na(pepcount_aligned)] <- 1
    #msg("[ReadTabExpressData] DEBUG: Original pepcount length=", length(dataSet$pepcount))
    #msg("[ReadTabExpressData] DEBUG: Aligned pepcount length=", length(pepcount_aligned))
    #msg("[ReadTabExpressData] DEBUG: data.proc rows=", nrow(data.proc))
    #msg("[ReadTabExpressData] DEBUG: pepcount range: ", min(pepcount_aligned, na.rm=TRUE), " - ", max(pepcount_aligned, na.rm=TRUE))
    #msg("[ReadTabExpressData] DEBUG: pepcount unique values: ", length(unique(pepcount_aligned)))
    #msg("[ReadTabExpressData] DEBUG: First 10 values: ", paste(head(pepcount_aligned, 10), collapse=", "))
    qs::qsave(pepcount_aligned, "pepcount.qs");
  } else if (is.null(dataSet$pepcount)) {
    #msg("[ReadTabExpressData] INFO: dataSet$pepcount is NULL - no peptide counts available")
  }

  saveSet(paramSet, "paramSet");
  saveSet(msgSet, "msgSet");
  return(RegisterData(dataSet));
}

.readMSstatsLong <- function(filePath) {
  msgSet <- readSet(msgSet, "msgSet");
  if (!file.exists(filePath)) {
    return(NULL);
  }
  # Peek header to see if it looks like MSstats long
  header <- try(readLines(filePath, n = 1), silent = TRUE);
  if (inherits(header, "try-error")) {
    return(NULL);
  }
  required <- c("ProteinName", "Run", "Condition", "BioReplicate", "Intensity");
  cols0 <- strsplit(header, "[,\t]")[[1]];
  if (!all(required %in% cols0)) {
    return(NULL);
  }
  #msg("[MSstats-long] reading ", filePath)
  dat <- try(data.table::fread(filePath, header = TRUE, check.names = FALSE, data.table = FALSE));
  if (inherits(dat, "try-error")) {
    msgSet$current.msg <- "Failed to read MSstats long-format file.";
    saveSet(msgSet, "msgSet");
    return(NULL);
  }
  if (!all(required %in% colnames(dat))) {
    msgSet$current.msg <- "MSstats long-format file missing required columns.";
    saveSet(msgSet, "msgSet");
    return(NULL);
  }
  value.col <- if ("LogIntensities" %in% colnames(dat)) "LogIntensities" else "Intensity";
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    msgSet$current.msg <- "reshape2 package not available to reshape MSstats long-format.";
    saveSet(msgSet, "msgSet");
    return(NULL);
  }
  wide <- reshape2::dcast(dat, ProteinName ~ Run, value.var = value.col);
  int.mat <- as.matrix(wide[, -1, drop = FALSE]);
  storage.mode(int.mat) <- "numeric";
  rownames(int.mat) <- wide$ProteinName;
  # Build metadata from Condition/BioReplicate per Run
  meta.df <- unique(dat[, c("Run", "Condition", "BioReplicate")]);
  rownames(meta.df) <- meta.df$Run;
  meta.df$Run <- NULL;
  disc.inx <- rep(TRUE, ncol(meta.df)); names(disc.inx) <- colnames(meta.df);
  cont.inx <- rep(FALSE, ncol(meta.df)); names(cont.inx) <- colnames(meta.df);
  # Order matrix columns to match metadata row order
  run.order <- intersect(rownames(meta.df), colnames(int.mat));
  int.mat <- int.mat[, run.order, drop = FALSE];
  meta.df <- meta.df[run.order, , drop = FALSE];
  meta.list <- list(meta.info = meta.df, cont.inx = cont.inx, disc.inx = disc.inx);
  #msg("[MSstats-long] matrix dims: ", paste(dim(int.mat), collapse = "x"))
  list(
    name = basename(filePath),
    data_orig = dat,
    data = int.mat,
    type = "prot",
    meta.info = meta.list,
    format = "msstats-long"
  );
}

# Synthesize minimal metadata when missing
.synthesizeMetaFromRuns <- function(run.names) {
  if (is.null(run.names) || length(run.names) == 0) {
    return(list(meta.info = data.frame(Condition = factor(), row.names = character(0)),
                disc.inx = c(Condition = TRUE),
                cont.inx = c(Condition = FALSE)))
  }
  cond.vec <- run.names;
  meta.df <- data.frame(Condition = factor(cond.vec), row.names = run.names, stringsAsFactors = FALSE);
  disc.inx <- c(Condition = TRUE);
  cont.inx <- c(Condition = FALSE);
  list(meta.info = meta.df, disc.inx = disc.inx, cont.inx = cont.inx);
}

#read annotation table file when user selects custom annotdataSet$meta.info <-ation option
ReadAnnotationTable <- function(fileName) {
  anot.data <- .readDataTable(fileName);
  msgSet <- readSet(msgSet, "msgSet");
  
  if(length(colnames(anot.data)) != 3){
    msgSet$current.msg <- "Please make sure the annotation contains exactly 3 columns";
  }
  colnames(anot.data) = c("gene_id", "symbol", "name");
  qs::qsave(anot.data, "anot_table.qs");
  saveSet(msgSet, "msgSet");
  return(1);
}

# Internal helper: convert wide expression matrix + metadata to MSstats long format
.buildMSstatsInput <- function(int.mat, metadata, fileName) {
  if (is.null(int.mat) || is.null(metadata)) {
    return(NULL)
  }

  runs <- colnames(int.mat);
  if (is.null(runs)) {
    return(NULL)
  }

  # Ensure metadata rows match runs
  shared <- intersect(runs, rownames(metadata));
  if (length(shared) == 0) {
    return(NULL)
  }
  runs <- shared

  # Pick first metadata column as condition; fall back to "Condition"
  cond.col <- metadata[runs, 1, drop = TRUE];
  if (is.null(cond.col)) {
    cond.col <- rep("Condition", length(runs));
  }

  # BioReplicate defaults to sample name
  bio.rep <- rownames(metadata)[match(runs, rownames(metadata))];
  if (is.null(bio.rep)) {
    bio.rep <- runs;
  }

  # Expand to long format
  long.df <- data.frame(
    ProteinName = rep(rownames(int.mat), times = length(runs)),
    PeptideSequence = rep(rownames(int.mat), times = length(runs)),
    PrecursorCharge = NA_integer_,
    FragmentIon = NA_character_,
    ProductCharge = NA_integer_,
    IsotopeLabelType = "L",
    Condition = rep(cond.col, each = nrow(int.mat)),
    BioReplicate = rep(bio.rep, each = nrow(int.mat)),
    Run = rep(runs, each = nrow(int.mat)),
    Intensity = as.numeric(c(int.mat[, runs, drop = FALSE])),
    stringsAsFactors = FALSE
  );

  # Remove rows with missing intensity
  long.df <- long.df[!is.na(long.df$Intensity), , drop = FALSE];

  # Add source name for traceability
  attr(long.df, "source_file") <- fileName;

  return(long.df);
}


ReadCustomLib <- function(fileName) {
  # Load the msgSet and paramSet objects
  msgSet <- readSet(msgSet, "msgSet")
  paramSet <- readSet(paramSet, "paramSet")
  
  paramSet$init.lib <- "custom";
  org <- paramSet$data.org
  idType <- paramSet$data.idType
  
  # Sanity check: Ensure the file exists
  if (!file.exists(fileName)) {
    msgSet$current.msg <- paste("Error: Gene set file", fileName, "not found.")
    saveSet(msgSet, "msgSet")
    return(0)
  }
  
  # Read the file into a character vector, line by line
  gene_data <- readLines(fileName)
  
  # Sanity check: Ensure the file is not empty
  if (length(gene_data) == 0) {
    msgSet$current.msg <- "Error: The input gene set file is empty."
    saveSet(msgSet, "msgSet")
    return(0)
  }
  
  # Step 1: Get the universe of all features across all sets
  universe_features <- unlist(lapply(gene_data, function(line) {
    line_parts <- strsplit(line, "\t+")[[1]]
    return(line_parts[-1])  # Return only the features, exclude the set name
  }))
  
  # Step 2: Perform ID conversion if idType is not "NA"
  if (idType != "NA") {
    # Perform ID conversion on the entire universe of features
    converted_universe <- .doAnnotation(universe_features, idType, paramSet)
    converted_universe <- unname(converted_universe)
    
    # Create a mapping of original gene IDs to converted IDs (for subsetting later)
    universe_map <- setNames(converted_universe, universe_features)
    
    if(!file.exists("symbol.map.qs")){
      if(idType %in% c("s2f", "generic", "ko")){
        symbol.map <- .doGeneIDMapping(anot.id, idType, paramSet, "matrix");
      }else{
        symbol.map <- .doGeneIDMapping(anot.id, "entrez", paramSet, "matrix");
      }
      symbol.map <- symbol.map[which(symbol.map$gene_id %in% anot.id),];
      saveDataQs(symbol.map, "symbol.map.qs", paramSet$anal.type, dataName);
      
    }
    
  } else {
    # If no ID conversion, keep the original gene IDs as the universe map
    universe_map <- setNames(universe_features, universe_features)
    symbol.map <- data.frame(gene_id=universe_features,symbol=universe_features);
    saveDataQs(symbol.map, "symbol.map.qs", paramSet$anal.type, dataName);
  }
  
  
  # Step 3: Initialize an empty list to store gene sets by cell line
  gene_set_list <- list()
  
  # Step 4: Loop through each line and use the pre-processed universe map
  for (line in gene_data) {
    # Split the line by tabs
    line_parts <- strsplit(line, "\t+")[[1]]
    
    # Sanity check: Ensure the line has at least one gene
    if (length(line_parts) < 2) {
      msgSet$current.msg <- paste("Error: Gene set for", line_parts[1], "is missing gene entries.")
      saveSet(msgSet, "msgSet")
      return(0)
    }
    
    # The first part is the set (e.g., cell line), the rest are the features
    set_name <- line_parts[1]
    features <- line_parts[-1]  # Everything after the first element
    
    # Step 5: Subset the universe map to get the converted IDs (or original IDs if no conversion)
    converted_ids <- universe_map[features]
    
    # Filter out any NA values (unmatched IDs)
    hit_inx <- !is.na(converted_ids)
    
    # Store the successfully converted IDs (or original IDs) for the gene set
    gene_set_list[[set_name]] <- converted_ids[hit_inx]
  }
  
  # Step 6: Save the gene set list into a .qs file
  qs::qsave(gene_set_list, "custom_lib.qs")
  
  # Update paramSet with the custom library file name
  paramSet$custom.lib <- fileName
  
  # Save msgSet and paramSet after the process
  saveSet(msgSet, "msgSet")
  saveSet(paramSet, "paramSet")
  
  return(1)
}




GetCustomLibColNames <- function(){
  paramSet <- readSet(paramSet, "paramSet")
  return(paramSet$custom.lib);
}

#' Read Metadata for Meta-Analysis Mode
#'
#' This function reads metadata for meta-analysis mode and performs necessary checks.
#'
#' @param metafilename A character string specifying the name of the metadata file.
#'
#' @return An integer indicating the success of reading and processing metadata.
#'
#' @author Guangyan Zhou \email{guangyan.zhou@mail.mcgill.ca}
#' @details Additional details about the function, if needed.
#'
#' @examples
#' \dontrun{
#' ReadMetaData(metafilename = "metadata.csv")
#' }
#'
#' @export
#' @license MIT License
#'
ReadMetaData <- function(metafilename){
  paramSet <- readSet(paramSet, "paramSet");
  msgSet <- readSet(msgSet,"msgSet");
  metadata <- try(data.table::fread(metafilename, header=TRUE, check.names=FALSE, data.table=FALSE));
  metadata[is.na(metadata)] = "NA"
  if(class(metadata) == "try-error"){
    msgSet$current.msg = "Failed to read in the metadata file! Please make sure that the metadata file is in the right format and does not have empty cells or contains NA."
    saveSet(msgSet, "msgSet");
    return(NULL);
  }
  # look for #NAME, store in a list
  sam.inx <- grep("^#NAME", colnames(metadata)[1]);
  if(length(sam.inx) > 0){
    smpl_nms<-metadata[,1];
    smpl_var<-colnames(metadata[-1]);
  }else{
    msgSet$current.msg = "Please make sure you have the label #NAME in your sample data file!"
    saveSet(msgSet, "msgSet");
    return(NULL);
  }
  # converting to character matrix as duplicate row names not allowed in data frame.
  metadata <-data.frame(lapply(1:ncol(metadata),function(x){
    metadata[,x]=unlist(ClearFactorStrings(metadata[,x]))
  }))
  metadata <- metadata[,-1,drop=F];
  if(nrow(metadata)==1){
    msgSet$current.msg = "Only one sample in the dataset or the metadata file must be transposed!"
    saveSet(msgSet, "msgSet");
    return(NULL);
  }
  rownames(metadata) <- smpl_nms;
  colnames(metadata) <- smpl_var;
  
  na.msg <- ""
  disc.inx <- GetDiscreteInx(metadata);
  if(sum(disc.inx) == length(disc.inx)){
    msgSet$na.msg <- "All metadata columns are OK!"
  }else{
    bad.meta<- paste(names(disc.inx)[!disc.inx], collapse="; ");
    msgSet$na.msg <- paste0("<font style=\"color:red\">Detected presence of unique values in the following columns: <b>", bad.meta, "</b></font>","Please make sure the metadata is in right format! You can use meta editor to update the information !");
  }
  
  
  cont.inx <- GetNumbericalInx(metadata);
  cont.inx <- !disc.inx & cont.inx; # discrete is first
  
  if(sum(cont.inx)>0){
    # make sure the discrete data is on the left side
    metadata <- cbind(metadata[,disc.inx, drop=FALSE], metadata[,cont.inx, drop=FALSE]);
  }
  
  metadata$Dataset <- rep("NA", nrow(metadata));
  
  mdata.all <- paramSet$mdata.all;
  # need to add metadata sanity check
  # are sample names identical to data$orig
  # order samples in same way as in abundance table
  sel.nms <- names(mdata.all);
  
  for(i in 1:length(sel.nms)){
    dataSet <- readDataset(sel.nms[i]);
    data.smpl.nms <- colnames(dataSet$data.norm)
    nm.hits <- data.smpl.nms %in% smpl_nms;
    if(!all(nm.hits)){
      msgSet$current.msg = paste0("Some sample names including ",paste(mis.nms, collapse="; ") ," in your data are not in the metadata file!")
      saveSet(msgSet, "msgSet");
      return(NULL);
    }
    
    # now remove extra meta if present, and order them
    nm.hits2 <- which(smpl_nms %in% data.smpl.nms);
    metadata$Dataset[nm.hits2] <- sel.nms[i];
    metadata1 <- metadata[nm.hits2,,drop=F];
    metadata1[] <- lapply( metadata1, factor)
    dataSet$meta.info <- dataSet$metaOrig <- metadata1
    dataSet$disc.inx <-dataSet$disc.inx.orig <- disc.inx[colnames(metadata1)]
    dataSet$cont.inx <-dataSet$cont.inx.orig  <- cont.inx[colnames(metadata1)]
    dataSet$cls <- dataSet$meta.info[,1];
    meta.types <- rep("disc", ncol(dataSet$meta.info));
    meta.types[cont.inx[colnames(metadata1)]] <- "cont";
    names(meta.types) <- colnames(dataSet$meta.info);
    dataSet$meta.types <-meta.types;
    
    RegisterData(dataSet);
  }
  
  paramSet$dataSet <- list();
  meta.types <- rep("disc", ncol(metadata));
  meta.types[cont.inx] <- "cont";
  names(meta.types) <- colnaymes(metadata);
  
  paramSet$dataSet$meta.types <- meta.types;
  paramSet$dataSet$meta.status <- rep("OK", ncol(metadata));
  paramSet$dataSet$cont.inx <- cont.inx;
  paramSet$dataSet$disc.inx <- disc.inx;
  paramSet$dataSet$meta.info <- metadata;
  paramSet$dataSet$metaOrig <- metadata;
  saveSet(msgSet, "msgSet");
  saveSet(paramSet, "paramSet");
  return(1);
}

# read tab delimited file
# can have many classes, stored in meta.info (starts with #)
# return a list (data.name, data.frame, meta.data)
.readTabData <- function(dataName) {

  msgSet <- readSet(msgSet, "msgSet");
  #msg("[readTabData] entry: ", dataName)
  # MaxQuant handling moved to ReadTabExpressData

  # Special handling for Proteome Discoverer allProteins exports (example and uploads)
  pd.res <- .readProteomeDiscoverer(dataName);
  if (!is.null(pd.res)) {
    #msg("[readTabData] returning Proteome Discoverer object with meta.info rows=", nrow(pd.res$meta.info$meta.info))
    return(pd.res);
  }

  # Special handling for MSstats long-format inputs
  msstats.res <- .readMSstatsLong(dataName);
  if (!is.null(msstats.res)) {
    #msg("[readTabData] returning MSstats-long object with meta.info rows=", nrow(msstats.res$meta.info$meta.info))
    return(msstats.res);
  }
  if(length(grep('\\.zip$',dataName,perl=TRUE))>0){
    dataName <- unzip(dataName);
    if(length(dataName) > 1){
      # test if "__MACOSX" or ".DS_Store"
      osInx <- grep('MACOSX',dataName,perl=TRUE);
      if(length(osInx) > 0){
        dataName <- dataName[-osInx];
      }
      dsInx <- grep('DS_Store',dataName,perl=TRUE);
      if(length(dsInx) > 0){
        dataName <- dataName[-dsInx];
      }
      dat.inx <- grep(".[Tt][Xx][Tt]$", dataName);
      if(length(dat.inx) != 1){
        msgSet$current.msg <- "More than one text files (.txt) found in the zip file.";
        saveSet(msgSet, "msgSet");        
        return(NULL);
      }
    }
  }
  
  msg <- NULL;
  # using the powerful fread function, 10 times faster, note: default return data.table, turn off
  datOrig <- .readDataTable(dataName);
  if(is.null(datOrig)){
    return(NULL);
  }
  dat1 <- .to.numeric.mat(datOrig);
  list(
    name= basename(dataName),
    data_orig = datOrig,
    data=dat1,
    type="count", # to be updated later
    format = "text"
  );
}
.readMaxQuantProteinGroups <- function(filePath, opts = list(quantType = "lfq", removeContaminants = TRUE, removeDecoys = TRUE, removeOnlyBySite = FALSE, minPeptides = NA)) {
  
  #msg("[MaxQuant] reading proteinGroups via MSnbase at ", filePath)
  
  if (!file.exists(filePath)) {
    warning("MaxQuant file not found: ", filePath)
    return(NULL)
  }

  # --- Step 1: Detect Columns (Quantification & IDs) ---
  # We read just the header first to find the column indices for MSnbase
  header <- names(data.table::fread(filePath, nrows = 1, header = TRUE, check.names = FALSE))
  
  # Regex logic to find Quant columns (LFQ, iBAQ, or Intensity)
  lfq.cols <- grep("(?i)^lfq[ ._]intensity[ ._]", header, perl = TRUE)
  ibaq.cols <- grep("(?i)^ibaq[ ._]", header, perl = TRUE)
  int.cols <- grep("(?i)^intensity[ ._]", header, perl = TRUE)
  
  qtype <- if (!is.null(opts$quantType)) tolower(opts$quantType) else "lfq"
  
  if (qtype == "lfq" && length(lfq.cols) > 0) {
    ecols <- lfq.cols; source.label <- "LFQ.intensity"
  } else if (qtype == "ibaq" && length(ibaq.cols) > 0) {
    ecols <- ibaq.cols; source.label <- "iBAQ"
  } else if (qtype == "intensity" && length(int.cols) > 0) {
    ecols <- int.cols; source.label <- "Intensity"
  } else {
    # Fallback priority
    if (length(lfq.cols) > 0) { ecols <- lfq.cols; source.label <- "LFQ.intensity" }
    else if (length(ibaq.cols) > 0) { ecols <- ibaq.cols; source.label <- "iBAQ" }
    else if (length(int.cols) > 0) { ecols <- int.cols; source.label <- "Intensity" }
    else { warning("No quantification columns found."); return(NULL) }
  }

  # Regex logic to find ID column
  id.col.cands <- c("^Protein IDs$", "^Protein\\.IDs$", "^Majority protein IDs$", "^Majority\\.protein\\.IDs$")
  id.col.name <- NULL
  for(cand in id.col.cands) {
    hits <- grep(cand, header, perl=TRUE, value=TRUE)
    if(length(hits) > 0) { id.col.name <- hits[1]; break }
  }
  if(is.null(id.col.name)) {
    # Broad fallback
    broad <- grep("(?i)protein[ ._]IDs?$", header, perl = TRUE, value = TRUE)
    if(length(broad) > 0) id.col.name <- broad[1]
    else { warning("No Protein ID column found."); return(NULL) }
  }

  # --- Step 2: Ingest Data using fread ---
  dat <- tryCatch(
    data.table::fread(filePath, header = TRUE, check.names = FALSE, data.table = FALSE),
    error = function(e) {
      warning("Failed to read MaxQuant file: ", e$message)
      NULL
    }
  )

  if (is.null(dat) || nrow(dat) == 0) {
    warning("MaxQuant file is empty or failed to load.")
    return(NULL)
  }

  if (!id.col.name %in% colnames(dat)) {
    warning("Protein ID column not found in MaxQuant file.")
    return(NULL)
  }

  quant.names <- header[ecols]
  quant.names <- quant.names[quant.names %in% colnames(dat)]
  if (length(quant.names) == 0) {
    warning("No quantification columns found in MaxQuant file.")
    return(NULL)
  }

  # --- Step 3: Filtering (Contaminants & Decoys) ---
  fd <- dat
  
  # Helper to find column case-insensitive
  find_col <- function(pattern, cols) { grep(pattern, cols, ignore.case=TRUE, value=TRUE)[1] }
  
  # Filter Contaminants
  if (isTRUE(opts$removeContaminants)) {
    cont_col <- find_col("^Potential[ .]contaminant$", colnames(fd))
    if (!is.na(cont_col)) {
      # Keep rows where column is NOT "+" and NOT "TRUE"
      keep <- fd[[cont_col]] != "+" & fd[[cont_col]] != "TRUE"
      # Handle NAs (sometimes MaxQuant leaves blank instead of -)
      keep[is.na(keep)] <- TRUE 
      dat <- dat[keep, , drop = FALSE]
    }
  }

  # Filter Decoys (Reverse)
  if (isTRUE(opts$removeDecoys)) {
    rev_col <- find_col("^Reverse$", colnames(fd))
    if (!is.na(rev_col)) {
      keep <- fd[[rev_col]] != "+" & fd[[rev_col]] != "TRUE"
      keep[is.na(keep)] <- TRUE
      dat <- dat[keep, , drop = FALSE]
    }
  }
  
  # Keep only-by-site entries at import; let normalization-view filtering decide later.
  site_col <- find_col("^Only[ .]identified[ .]by[ .]site$", colnames(fd))
  if (isTRUE(opts$removeOnlyBySite) && !is.na(site_col)) {
    keep <- fd[[site_col]] != "+" & fd[[site_col]] != "TRUE"
    keep[is.na(keep)] <- TRUE
    dat <- dat[keep, , drop = FALSE]
  }
  # refresh feature data after filtering
  fd <- dat

  # Optional min peptides per protein (razor + unique)
  if (!is.null(opts$minPeptides) && !is.na(opts$minPeptides)) {
    pep_col <- find_col("^Razor[ .+]unique[ .+]peptides$", colnames(fd))
    if (is.na(pep_col)) {
      # try sanitized name
      pep_col <- find_col("^Razor\\.\\.\\.unique\\.peptides$", colnames(fd))
    }
    if (!is.na(pep_col)) {
      keep <- suppressWarnings(as.numeric(fd[[pep_col]]) >= as.numeric(opts$minPeptides))
      keep[is.na(keep)] <- FALSE
      dat <- dat[keep, , drop = FALSE]
      fd <- dat
    }
  }

  # --- Step 4: Clean IDs and Quantification Matrix ---
  
  # Clean IDs: MaxQuant returns "ID1;ID2", we want just "ID1"
  current_ids <- as.character(dat[[id.col.name]])
  clean_ids <- vapply(strsplit(current_ids, ";"), function(x) x[1], character(1))
  
  # Assign clean IDs back. Note: If duplicates exist after splitting, MSnbase might complain.
  # We make unique just in case, matching your original script logic
  if(any(duplicated(clean_ids))) clean_ids <- make.unique(clean_ids)
  rownames(dat) <- clean_ids
  
  # Clean Quant Data: Replace 0 with NA (Standard MaxQuant practice)
  int.mat <- as.matrix(dat[, quant.names, drop = FALSE])
  for (i in seq_len(ncol(int.mat))) {
    int.mat[, i] <- suppressWarnings(as.numeric(int.mat[, i]))
  }
  int.mat[int.mat == 0] <- NA

  # --- Step 5: Format Output to Match Original Function ---
  
  # 1. Intensity Matrix
  # Remove rows that are all NA after zero->NA
  keep.non.na <- rowSums(!is.na(int.mat)) > 0
  if (!any(keep.non.na)) {
    warning("No data left after filtering (all rows NA).")
    return(NULL)
  }
  dat <- dat[keep.non.na, , drop = FALSE]
  int.mat <- int.mat[keep.non.na, , drop = FALSE]
  clean_ids <- clean_ids[keep.non.na]
  rownames(int.mat) <- clean_ids
  fd_final <- dat[, setdiff(colnames(dat), quant.names), drop = FALSE]
  rownames(fd_final) <- clean_ids
  
  # Clean column names (Remove "LFQ intensity " prefix)
  # The column names in MSnbase correspond to the selected ecols
  raw_colnames <- colnames(int.mat)
  run_ids <- sub("(?i)^lfq[ ._]intensity[ ._]|(?i)^intensity[ ._]|(?i)^ibaq[ ._]", "", raw_colnames, perl = TRUE)
  
  # Handle duplicates in run IDs
  if(any(duplicated(run_ids))) run_ids <- make.unique(run_ids, sep="_")
  colnames(int.mat) <- run_ids

  # 2. Metadata (Conditions)
  cond.vec <- sub("\\.[^.]+$", "", run_ids)
  meta.df <- data.frame(Condition = cond.vec, row.names = run_ids, check.names = FALSE, stringsAsFactors = FALSE)
  disc.inx <- rep(TRUE, ncol(meta.df)); names(disc.inx) <- colnames(meta.df)
  cont.inx <- rep(FALSE, ncol(meta.df)); names(cont.inx) <- colnames(meta.df)
  meta.info <- list(meta.info = meta.df, cont.inx = cont.inx, disc.inx = disc.inx)

  # 3. Peptide Counts (for DEqMS)
  # We look inside fd_final for the count columns
  pepcount <- NULL

  # Priority list for peptide counts
  # Note: "Razor...unique.peptides" is how R usually sanitizes "Razor + unique peptides"
  # But since we used check.names=FALSE in readMSnSet2, we might find "Razor + unique peptides"
  count_candidates <- c("Peptides", "Razor + unique peptides", "Razor...unique.peptides", "MS.MS.count", "MS.MS.Count")

  #msg("[MaxQuant] DEBUG: Looking for peptide count columns")
  #msg("[MaxQuant] DEBUG: Available fData columns: ", paste(colnames(fd_final), collapse=", "))
  #msg("[MaxQuant] DEBUG: Checking candidates: ", paste(count_candidates, collapse=", "))

  for(col in count_candidates) {
    if(col %in% colnames(fd_final)) {
      pepcount <- as.numeric(fd_final[[col]])
      names(pepcount) <- clean_ids
      #msg("[MaxQuant] Extracted '", col, "' for DEqMS.")
      #msg("[MaxQuant] DEBUG: pepcount range: ", min(pepcount, na.rm=TRUE), " - ", max(pepcount, na.rm=TRUE))
      #msg("[MaxQuant] DEBUG: pepcount unique values: ", length(unique(pepcount)))
      #msg("[MaxQuant] DEBUG: First 10 values: ", paste(head(pepcount, 10), collapse=", "))
      break
    }
  }

  if(is.null(pepcount)) {
    #msg("[MaxQuant] WARNING: No peptide count column found. Available columns: ", paste(colnames(fd_final), collapse=", "))
  }

  # Save MaxQuant metadata for later filtering during normalization
  # Include Protein IDs, contaminant flags, reverse flags, and peptide counts
  mq_metadata <- data.frame(
    Protein.IDs = clean_ids,
    stringsAsFactors = FALSE
  )

  # Add contaminant column if present
  cont_col <- grep("(?i)^Potential[ .]contaminant$", colnames(fd_final), perl = TRUE, value = TRUE)
  if (length(cont_col) > 0) {
    mq_metadata$Potential.contaminant <- fd_final[[cont_col[1]]]
  }

  # Add reverse column if present
  rev_col <- grep("(?i)^Reverse$", colnames(fd_final), perl = TRUE, value = TRUE)
  if (length(rev_col) > 0) {
    mq_metadata$Reverse <- fd_final[[rev_col[1]]]
  }

  # Add peptide count if available
  if (!is.null(pepcount)) {
    pep_col_name <- NULL
    for (col in c("Peptides", "Razor...unique.peptides", "Razor.unique.peptides")) {
      if (col %in% colnames(fd_final)) {
        pep_col_name <- col
        break
      }
    }
    if (!is.null(pep_col_name)) {
      mq_metadata[[pep_col_name]] <- fd_final[[pep_col_name]]
    }
  }

  qs::qsave(mq_metadata, "maxquant_metadata.qs")
  #msg("[MaxQuant] Saved metadata with ", nrow(mq_metadata), " proteins for later filtering")

  # Return the exact list structure your workflow expects
  list(
    name = basename(filePath),
    # Combine feature data and expression data to mimic the "original full dataframe"
    data_orig = cbind(fd_final, int.mat),
    data = int.mat,
    type = "prot",
    meta.info = meta.info,
    format = "maxquant",
    pepcount = pepcount
  )
}

# Helper: read MaxQuant peptides.txt for peptide-level analysis
.readMaxQuantPeptides <- function(filePath, opts = list(quantType = "intensity", removeContaminants = TRUE, removeDecoys = TRUE)) {

  #msg("[MaxQuantPep] reading peptides.txt at ", filePath)

  if (!file.exists(filePath)) {
    warning("MaxQuantPep peptides file not found: ", filePath)
    return(NULL)
  }

  # Read the full table
  dat <- tryCatch(
    data.table::fread(filePath, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) {
      warning("Failed to read MaxQuantPep peptides.txt: ", e$message)
      return(NULL)
    }
  )

  if (is.null(dat) || nrow(dat) == 0) {
    warning("MaxQuantPep peptides.txt is empty or failed to load")
    return(NULL)
  }

  # Convert to data.frame for easier manipulation
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)

  # --- Step 1: Filter Contaminants and Decoys ---
  if (isTRUE(opts$removeContaminants)) {
    cont_col <- grep("(?i)^Potential[ .]contaminant$", colnames(dat), perl = TRUE, value = TRUE)
    if (length(cont_col) > 0) {
      keep <- dat[[cont_col[1]]] != "+" & dat[[cont_col[1]]] != "TRUE"
      keep[is.na(keep)] <- TRUE
      dat <- dat[keep, , drop = FALSE]
    }
    # Also check Proteins column for CON__ contaminants
    if ("Proteins" %in% colnames(dat)) {
      keep <- !grepl("CON__|^CON_", dat$Proteins)
      dat <- dat[keep, , drop = FALSE]
    }
  }

  if (isTRUE(opts$removeDecoys)) {
    rev_col <- grep("(?i)^Reverse$", colnames(dat), perl = TRUE, value = TRUE)
    if (length(rev_col) > 0) {
      keep <- dat[[rev_col[1]]] != "+" & dat[[rev_col[1]]] != "TRUE"
      keep[is.na(keep)] <- TRUE
      dat <- dat[keep, , drop = FALSE]
    }
    # Also check Proteins column for REV__ decoys
    if ("Proteins" %in% colnames(dat)) {
      keep <- !grepl("REV__|^REV_", dat$Proteins)
      dat <- dat[keep, , drop = FALSE]
    }
  }

  if (nrow(dat) == 0) {
    warning("No peptides remain after filtering contaminants/decoys")
    return(NULL)
  }

  # --- Step 2: Identify Intensity Columns ---
  # MaxQuant peptides.txt uses "Intensity {sample}" columns
  int_cols <- grep("^Intensity ", colnames(dat), value = TRUE)

  if (length(int_cols) == 0) {
    warning("No 'Intensity' columns found in peptides.txt")
    return(NULL)
  }

  # Extract sample names from column names (remove "Intensity " prefix)
  sample_names <- sub("^Intensity ", "", int_cols)

  # --- Step 3: Build Intensity Matrix ---
  intens <- dat[, int_cols, drop = FALSE]
  intens <- .safe_numeric_matrix(intens)

  # Replace 0 with NA (MaxQuant convention)
  intens[intens == 0] <- NA

  # Remove rows that are all NA
  keep_rows <- rowSums(!is.na(intens)) > 0
  if (!any(keep_rows)) {
    warning("All peptide intensities are NA/zero")
    return(NULL)
  }

  intens <- intens[keep_rows, , drop = FALSE]
  dat <- dat[keep_rows, , drop = FALSE]

  # --- Step 4: Set Row IDs ---
  # Prefer "Modified sequence" to distinguish PTMs, fallback to "Sequence"
  if ("Modified sequence" %in% colnames(dat)) {
    peptide_ids <- dat[["Modified sequence"]]
  } else if ("Sequence" %in% colnames(dat)) {
    peptide_ids <- dat[["Sequence"]]
  } else {
    warning("No 'Sequence' or 'Modified sequence' column found")
    return(NULL)
  }

  # Handle duplicates
  if (any(duplicated(peptide_ids))) {
    #msg("[MaxQuantPep] Warning: Duplicate peptide sequences found. Making IDs unique.")
    peptide_ids <- make.unique(as.character(peptide_ids))
  }

  rownames(intens) <- peptide_ids
  colnames(intens) <- sample_names

  # --- Step 5: Build Peptide-to-Protein Map ---
  prot_col <- NULL
  if ("Proteins" %in% colnames(dat)) {
    prot_col <- dat[["Proteins"]]
  } else if ("Leading razor protein" %in% colnames(dat)) {
    prot_col <- dat[["Leading razor protein"]]
  } else if ("Leading proteins" %in% colnames(dat)) {
    prot_col <- dat[["Leading proteins"]]
  }

  prot_map <- NULL
  if (!is.null(prot_col)) {
    # Clean protein IDs: MaxQuant uses semicolon-separated lists, take first
    clean_prots <- vapply(strsplit(as.character(prot_col), ";"), function(x) x[1], character(1))
    prot_map <- data.frame(
      Peptide = peptide_ids,
      Protein = clean_prots,
      stringsAsFactors = FALSE
    )
  }

  # --- Step 6: Infer Metadata from Sample Names ---
  # Use same logic as protein reader: split by last dot/underscore
  cond_vec <- sub("\\.[^.]+$", "", sample_names)
  meta_df <- data.frame(
    Condition = factor(cond_vec),
    stringsAsFactors = FALSE
  )
  rownames(meta_df) <- sample_names

  disc_inx <- c(Condition = TRUE)
  cont_inx <- c(Condition = FALSE)

  #msg("[MaxQuant] Peptide data loaded: ", nrow(intens), " peptides × ", ncol(intens), " samples")

  # --- Step 7: Return Structure ---
  list(
    data = intens,
    data_orig = intens,
    type = "peptide",  # Flag as peptide-level
    format = "maxquant-peptide",
    meta.info = list(meta.info = meta_df, disc.inx = disc_inx, cont.inx = cont_inx),
    prot.map = prot_map  # Peptide-to-protein mapping
  )
}

# Record MaxQuant options from UI (stored in R options)
SetMaxQuantOptions <- function(quantType = "lfq", removeContaminants = TRUE, removeDecoys = TRUE, minPeptides = NA) {
  paramSet <- readSet(paramSet, "paramSet")
  opts <- list(
    quantType = tolower(quantType),
    removeContaminants = isTRUE(as.logical(removeContaminants)),
    removeDecoys = isTRUE(as.logical(removeDecoys)),
    minPeptides = suppressWarnings(as.numeric(minPeptides))
  )
  options(pa.maxquant.opts = opts)
  paramSet$maxquant <- opts
  saveSet(paramSet, "paramSet")
  #msg("[MaxQuant] options set: quantType=", opts$quantType,
  #        " removeContaminants=", opts$removeContaminants,
  #        " removeDecoys=", opts$removeDecoys,
  #        " minPeptides=", opts$minPeptides)
  return(1L)
}

# Record DIA-NN options from UI (stored in R options + paramSet)
SetDiannOptions <- function(fileType = "protein_matrix", qvalueFilter = TRUE) {
  paramSet <- readSet(paramSet, "paramSet")
  opts <- list(
    fileType = fileType,
    qvalueFilter = isTRUE(as.logical(qvalueFilter))
  )
  options(pa.diann.opts = opts)
  paramSet$diann <- opts
  saveSet(paramSet, "paramSet")
  #msg("[DIA-NN] options set: fileType=", opts$fileType, " qvalueFilter=", opts$qvalueFilter)
  return(1L)
}

# Record FragPipe options from UI (stored in R options + paramSet)
SetFragpipeOptions <- function(quantType = "protein_maxlfq", removeContaminants = TRUE, minProb = NA, minPeptides = NA) {
  paramSet <- readSet(paramSet, "paramSet")
  opts <- list(
    quantType = quantType,
    removeContaminants = isTRUE(as.logical(removeContaminants)),
    minProb = suppressWarnings(as.numeric(minProb)),
    minPeptides = suppressWarnings(as.numeric(minPeptides))
  )
  options(pa.fragpipe.opts = opts)
  paramSet$fragpipe <- opts
  saveSet(paramSet, "paramSet")
  #msg("[FragPipe] options set: quantType=", opts$quantType, " removeContaminants=", opts$removeContaminants,
  #        " minProb=", opts$minProb, " minPeptides=", opts$minPeptides)
  return(1L)
}

SetSpectronautOptions <- function(inputType = "protein") {
  paramSet <- readSet(paramSet, "paramSet")
  opts <- list(
    inputType = if (is.null(inputType) || !nzchar(inputType)) "protein" else tolower(inputType)
  )
  options(pa.spectronaut.opts = opts)
  paramSet$spectronaut <- opts
  saveSet(paramSet, "paramSet")
  return(1L)
}

.readProteomeDiscoverer <- function(filePath) {
  msgSet <- readSet(msgSet, "msgSet");
  # peek header, support .gz
  header <- tryCatch({
    con <- if (grepl("\\.gz$", filePath, ignore.case = TRUE)) gzfile(filePath, "rt") else file(filePath, "rt")
    on.exit(close(con), add = TRUE)
    readLines(con, n = 1)
  }, error = function(e) NA_character_)
  if (length(header) == 0 || all(is.na(header)) || !grepl("Abundance\\.", header)) {
    return(NULL);
  }
  #msg("[ProteomeDiscoverer] reading allProteins at ", filePath)
  # only keep accession + abundance columns
  dat <- try({
    cols <- data.table::fread(filePath, header = TRUE, nrows = 0, data.table = FALSE, check.names = FALSE);
    ab.cols <- grep("^Abundance[[:space:].]", colnames(cols), value = TRUE);
    keep.cols <- c("Accession", ab.cols);
    if (length(ab.cols) == 0) stop("no Abundance columns");
    data.table::fread(filePath, header = TRUE, check.names = FALSE, data.table = FALSE, select = keep.cols)
  }, silent = TRUE);
  if (inherits(dat, "try-error")) {
    msgSet$current.msg <- c(msgSet$current.msg, "Failed to read Proteome Discoverer file.");
    saveSet(msgSet, "msgSet");
    return(NULL);
  }
  abund.idx <- grep("^Abundance", colnames(dat));
  if (length(abund.idx) == 0) {
    return(NULL);
  }
  ab.names <- colnames(dat)[abund.idx];
  ab.names <- sub("^Abundance[[:space:].]*", "", ab.names);
  samp.names <- ab.names;
  cond <- gsub("[._]?rep.*$", "", samp.names, ignore.case = TRUE);
  cond[cond == ""] <- "Condition";
  # Attempt ingestion via MSstats PD converter if available
  int.mat <- NULL; meta.list <- NULL; used_msstats <- FALSE;
  if (requireNamespace("MSstats", quietly = TRUE) && requireNamespace("reshape2", quietly = TRUE)) {
    ann <- data.frame(Run = samp.names,
                      Condition = cond,
                      BioReplicate = samp.names,
                      stringsAsFactors = FALSE);
    msin <- try(MSstats::PDtoMSstatsFormat(raw = dat,
                                           annotation = ann,
                                           useUniquePeptide = TRUE,
                                           removeFewMeasurements = FALSE,
                                           removeOxidationMpeptide = FALSE,
                                           removeProtein_with1Peptide = FALSE,
                                           summaryMethod = "TMP",
                                           MBimpute = TRUE,
                                           logTrans = 2),
                silent = TRUE);
    if (!inherits(msin, "try-error") &&
        all(c("ProteinName", "Run", "Condition", "BioReplicate", "Intensity") %in% colnames(msin))) {
      wide <- reshape2::dcast(msin, ProteinName ~ Run, value.var = "Intensity");
      int.mat <- as.matrix(wide[, -1, drop = FALSE]);
      rownames(int.mat) <- wide$ProteinName;
      storage.mode(int.mat) <- "numeric";
      meta.df <- unique(msin[, c("Run", "Condition", "BioReplicate")]);
      rownames(meta.df) <- meta.df$Run;
      meta.df$Run <- NULL;
      disc.inx <- c(Condition = TRUE, BioReplicate = TRUE);
      cont.inx <- c(Condition = FALSE, BioReplicate = FALSE);
      meta.list <- list(meta.info = meta.df, disc.inx = disc.inx, cont.inx = cont.inx);
      used_msstats <- TRUE;
    }
  }

  # Fallback to manual matrix build
  if (is.null(int.mat)) {
    int.mat <- as.matrix(dat[, abund.idx, drop = FALSE]);
    storage.mode(int.mat) <- "numeric";
    colnames(int.mat) <- samp.names;

    # feature IDs: prefer Accession, else rownames
    if ("Accession" %in% colnames(dat)) {
      row.ids <- dat$Accession;
    } else {
      row.ids <- seq_len(nrow(int.mat));
    }
    rownames(int.mat) <- row.ids;

    # metadata: Condition inferred from prefix before first digit/underscore; BioReplicate = sample name
    meta.df <- data.frame(
      Condition = factor(cond),
      BioReplicate = factor(samp.names),
      row.names = samp.names,
      stringsAsFactors = FALSE
    );
    disc.inx <- c(Condition = TRUE, BioReplicate = TRUE);
    cont.inx <- c(Condition = FALSE, BioReplicate = FALSE);
    meta.list <- list(meta.info = meta.df, disc.inx = disc.inx, cont.inx = cont.inx);
  }

  list(
    name = basename(filePath),
    data_orig = dat,
    data = int.mat,
    type = "prot",
    meta.info = meta.list,
    format = "proteome-discoverer"
  );
}


# note, try to use the fread, however, it has issues with 
# some windows 10 files "Line ending is \r\r\n. .... appears to add the extra \r in text mode on Windows"
# in such as, use the slower read.table method
.readDataTable <- function(fileName){
  msgSet <- readSet(msgSet, "msgSet");
  if(length(grep('\\.zip$',fileName,perl=TRUE))>0){
    fileName <- unzip(fileName);
    if(length(fileName) > 1){
      # test if "__MACOSX" or ".DS_Store"
      osInx <- grep('MACOSX',fileName,perl=TRUE);
      if(length(osInx) > 0){
        fileName <- fileName[-osInx];
      }
      dsInx <- grep('DS_Store',fileName,perl=TRUE);
      if(length(dsInx) > 0){
        fileName <- fileName[-dsInx];
      }
      dat.inx <- grep(".[Tt][Xx][Tt]$", fileName);
      if(length(dat.inx) != 1){
        msgSet$current.msg <- "More than one text files (.txt) found in the zip file.";
        return(NULL);
      }
    }
  }
  dat <- try(data.table::fread(fileName, header=TRUE, check.names=FALSE, data.table=FALSE));
  rm.inx <- apply(dat,2,function(x){all(is.na(x))});
  dat <- dat[,!rm.inx];
  if(class(dat) == "try-error"){
    #try to use "tr" to remove double return characters
    trFileName <- paste("tr -d \'\\r\' <", fileName);
    dat <- try(data.table::fread(trFileName, header=TRUE, check.names=FALSE, data.table=FALSE));
    if(class(dat) == "try-error"){
      print("Using slower file reader ...");
      formatStr <- substr(fileName, nchar(fileName)-2, nchar(fileName))
      if(formatStr == "txt"){
        dat <-try(read.table(fileName,header=TRUE,comment.char = "", check.names=F, as.is=T));
      }else{ # note, read.csv is more than read.table with sep=","
        dat <-try(read.csv(fileName,header=TRUE,comment.char = "", check.names=F, as.is=T));
      }  
    }
  }
  if(class(dat) == "try-error"){
    msgSet$current.msg <- "Failed to read the data table! Please check your data format.";
    saveSet(msgSet, "msgSet");
    return(NULL);
  }
  
  
  # need to remove potential empty columns
  dat <- dat[!sapply(dat, function(x) all(x == "" | is.na(x)))];
  
  return(dat);
}


####read meta file
#### return a list
.readMetaData <- function(metafileName, datOrig, metaContain) {
  #msg("[.readMetaData] entered function; metafileName=", metafileName)
  msgSet <- readSet(msgSet, "msgSet");
  na.msg <- ""
  if(is.null(msgSet$current.msg)){
    msg <-""
  }else{
    msg <- msgSet$current.msg
  }
  match.msg <- "";
  
  if(metaContain=="true"){
    meta.info <- list();
    # look for #CLASS, could have more than 1 class labels, store in a list
    cls.inx <- grep("^#CLASS", datOrig[,1]);
    if(length(cls.inx) > 0){ 
      for(i in 1:length(cls.inx)){
        inx <- cls.inx[i];
        cls.nm <- substring(datOrig[inx, 1],2); # discard the first char #
        if(nchar(cls.nm) > 6){
          cls.nm <- substring(cls.nm, 7); # remove class
        }
        if(grepl("[[:blank:]]", cls.nm)){
          cls.nm<- gsub("\\s+","_", cls.nm);
          msg <- c(msg, " Blank spaces in group names are replaced with underscore '_'! ");
        }
        cls.lbls <- setNames(as.character(datOrig[inx, -1]),colnames(datOrig)[-1]);
        # test NA
        na.inx <- is.na(cls.lbls);
        cls.lbls[na.inx] <- "NA";
        cls.lbls <- ClearFactorStrings(cls.lbls);
        meta.info[[cls.nm]] <- cls.lbls;
      }
    }else{
      msgSet$current.msg <- "No metadata labels #CLASS found in your data!";
      saveSet(msgSet, "msgSet");
      return(NULL);
    }
    
    meta.info <- data.frame(meta.info);
    rownames(meta.info) = colnames(datOrig)[-1]
  }else{ # metadata input as an individual table
    #msg("[.readMetaData] reading external metadata file.")
    mydata <- try(data.table::fread(metafileName, header=TRUE, check.names=FALSE, data.table=FALSE));
    if(class(mydata) == "try-error"){
      msgSet$current.msg <- "Failed to read the metadata table! Please check your data format.";
      saveSet(msgSet, "msgSet");
      return(NULL);
    }
    #msg("[.readMetaData] mydata head:\n", paste(capture.output(head(mydata)), collapse="\n"))

    # standardize column name for #NAME
    name_col_inx <- which(tolower(colnames(mydata)) == "#name")
    if (length(name_col_inx) > 0) {
        colnames(mydata)[name_col_inx[1]] <- "#NAME"
        #msg("[.readMetaData] Found and standardized #NAME column.")
    }

    # trim whitespace from sample names before matching
    if (!("#NAME" %in% colnames(mydata))) {
      if ("SampleID" %in% colnames(mydata)) {
        mydata$`#NAME` <- mydata$SampleID
      } else if ("Sample" %in% colnames(mydata)) {
        mydata$`#NAME` <- mydata$Sample
      } else {
        msgSet$current.msg <- "Metadata file must contain a '#NAME' or 'SampleID' column for sample identifiers.";
        saveSet(msgSet, "msgSet");
        return(NULL);
      }
    }
    # Ensure #NAME is the first column for downstream processing
    mydata <- mydata[, c("#NAME", setdiff(colnames(mydata), "#NAME")), drop = FALSE]
    # Drop duplicate sample identifier columns (e.g., Sample / SampleID) to avoid treating them as metadata
    dup_cols <- setdiff(which(vapply(mydata, function(col) {
      identical(trimws(as.character(col)), mydata$`#NAME`)
    }, logical(1))), 1)
    if (length(dup_cols) > 0) {
      mydata <- mydata[, -dup_cols, drop = FALSE]
    }
    mydata$`#NAME` <- trimws(as.character(mydata$`#NAME`))
    colnames(datOrig) <- trimws(colnames(datOrig))
    #msg("[.readMetaData] datOrig dim=", paste(dim(datOrig), collapse="x"), "; is.matrix=", is.matrix(datOrig))
    #msg("[.readMetaData] head(colnames(datOrig)): ", paste(head(colnames(datOrig)), collapse=", "))

    # Align to data columns (exclude feature column)
    smpl_nms <- colnames(datOrig)
    if (!is.matrix(datOrig)) {
      smpl_nms <- smpl_nms[-1]
    }
    #msg("[.readMetaData] smpl_nms to match: ", paste(head(smpl_nms), collapse=", "))
    
    hit.inx <- match(smpl_nms, mydata$`#NAME`)
    #msg("[.readMetaData] hit.inx: ", paste(head(hit.inx), collapse=", "))

    miss.meta <- is.na(hit.inx)
    if (any(miss.meta)) {
      drop.smpl <- smpl_nms[miss.meta]
      if (is.matrix(datOrig)) {
        datOrig <- datOrig[, !miss.meta, drop=FALSE]
      } else {
        datOrig <- datOrig[, c(TRUE, !miss.meta), drop=FALSE]
      }
      match.msg <- paste0(match.msg, length(drop.smpl), " samples not found in metadata dropped: ", paste(drop.smpl, collapse=", "), ". ")
      smpl_nms <- smpl_nms[!miss.meta]
      #msg("[.readMetaData] smpl_nms updated after dropping non-meta samples: ", paste(head(smpl_nms), collapse=", "))
      hit.inx <- hit.inx[!miss.meta]
    }
    miss.data <- is.na(match(mydata$`#NAME`, smpl_nms))
    if (any(miss.data)) {
      drop.meta <- mydata$`#NAME`[miss.data]
      mydata <- mydata[!miss.data, , drop=FALSE]
      match.msg <- paste0(match.msg, length(drop.meta), " metadata rows removed (no matching data): ", paste(head(drop.meta, 3), collapse=", "), ". ")
    }
    mydata <- mydata[match(smpl_nms, mydata$`#NAME`), , drop=FALSE]
    rownames(mydata) <- mydata$`#NAME`
    # Replace common Excel error tokens with NA so downstream checks keep the samples
    errTokens <- c("#VALUE!", "#DIV/0!", "#N/A", "#NULL!", "#NUM!", "#REF!")
    for(token in errTokens) {
      mydata[mydata == token] <- NA
    }
    mydata[mydata == ""] <- NA
    mydata[is.na(mydata)] <- "NA";
    # look for #NAME, store in a list
    sam.inx <- which(colnames(mydata) == "#NAME")
    if(length(sam.inx) > 0){
      name_col <- sam.inx[1]
      smpl_nm<-mydata[, name_col];
    smpl_var<-colnames(mydata[-name_col]);
    }else{
      msgSet$current.msg <- "Please make sure you have the label #NAME in your sample data file!";
      saveSet(msgSet, "msgSet");
      return(NULL);
    }

    empt <- remove_empty_cols(mydata)
    mydata <- empt$cleaned
    if (length(empt$droppedNames)) {
      match.msg <- paste0(
        match.msg, "Columns ",
        paste(empt$droppedNames, collapse = ", "),
        " were completely empty and removed from metadata file.   "
      )
    }
    
    # covert to factor
    mydata <-data.frame(lapply(1:ncol(mydata),function(x){
      mydata[,x]=unlist(ClearFactorStrings(mydata[,x]))
    }))
    mydata <- mydata[,-name_col,drop=F]; # converting to character matrix as duplicate row names not allowed in data frame.
    if(nrow(mydata)==1){
      msgSet$current.msg <- "Only one sample in the dataset or the metadata file must be transposed!";
      saveSet(msgSet, "msgSet");
      return(NULL);
    }
    rownames(mydata) <- smpl_nm;
    colnames(mydata) <- smpl_var;
    
    # empty cell or NA cannot be tolerated in metadata
    na.inx  <- is.na(mydata);
    na.msg <- na.msg1 <- NULL;
    if(sum(na.inx) > 0){
      na.msg1 <- paste("A total of", sum(na.inx), "empty or NA values were detected. Please update in using metadata editor");
    }
    
    #Check group label names for spaces and replace with underscore
    meta.info <- data.frame(mydata,check.names=FALSE);
    if(any(grepl("[[:blank:]]", names(meta.info)))){
      names(meta.info) <- gsub("\\s+","_", names(meta.info));
      na.msg1 <- c(na.msg1, "Blank spaces in group names are replaced with underscore '_'");
    }
    
  }
  
  # Ensure meta.info is a data.frame
  if (!is.data.frame(meta.info)) {
    meta.info <- as.data.frame(meta.info, stringsAsFactors = FALSE)
    if (is.null(rownames(meta.info)) || any(rownames(meta.info) == "")) {
      rownames(meta.info) <- colnames(datOrig)[-1][seq_len(nrow(meta.info))]
    }
  }

  # expose sanitized metadata for debugging
  fast.write(meta.info, file = "metadata_processed.csv", row.names = TRUE)
  disc.inx <- GetDiscreteInx(meta.info);
  if(sum(disc.inx) == length(disc.inx)){
    na.msg <- c(na.msg,"All metadata columns are OK!")
  }else{
    bad.meta<- paste(names(disc.inx)[!disc.inx], collapse="; ");
    na.msg <- c(na.msg, paste0("<font style=\"color:red\">Detected presence of unique values in the following columns: <b>", bad.meta, "</b></font>","Please make sure the metadata is in right format! You can use meta editor to update the information !"));
  }
  
  cont.inx <- GetNumbericalInx(meta.info);
  cont.inx <- !disc.inx & cont.inx; # discrete is first
  
  rmcol <- intersect(which(!disc.inx),which(!cont.inx ))
  
  if(length(rmcol)==1){
    match.msg <- paste0(match.msg, "Column ",names(meta.info)[rmcol]," is removed due to lack of replicates!   " )
  }else if(length(rmcol)>1){
    match.msg <- paste0(match.msg, "Columns ",paste(names(meta.info)[rmcol],collapse = ", ")," are removed due to lack of replicates!   " )
  }
  
  if(sum(cont.inx)>0){
    # make sure the discrete data is on the left side
    meta.info <- cbind(meta.info[,disc.inx, drop=FALSE], meta.info[,cont.inx, drop=FALSE]);
  }
  disc.inx <- disc.inx[colnames(meta.info)]
  cont.inx <- cont.inx[colnames(meta.info)]
  msgSet$match.msg <-match.msg
  msgSet$na.msg <- na.msg
  saveSet(msgSet, "msgSet");  
  return(list(meta.info=meta.info,disc.inx=disc.inx,cont.inx=cont.inx))
}

# Helper: read FragPipe LFQ protein table + experiment annotation
.readFragpipeLfq <- function(data.file, meta.file, opts = list(quantType = "protein_maxlfq", removeContaminants = TRUE, minProb = NA, minPeptides = NA)) {
  if (!file.exists(data.file)) {
    #msg("[FragPipe] data file not found: ", data.file)
    return(NULL)
  }
  if (is.null(meta.file) || meta.file == "" || !file.exists(meta.file)) {
    #msg("[FragPipe][ERROR] metadata file is required: ", meta.file)
    return(NULL)
  }
  dat <- tryCatch(read.delim(data.file, check.names = FALSE, stringsAsFactors = FALSE, comment.char = ""), error = function(e) NULL)
  if (is.null(dat)) {
    #msg("[FragPipe][ERROR] unable to read data file")
    return(NULL)
  }
  meta <- tryCatch(read.delim(meta.file, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(meta)) {
    #msg("[FragPipe][ERROR] unable to read metadata file")
    return(NULL)
  }

  # Identify run column
  run.col.cands <- c("Run", "Filename", "Raw.file", "Raw", "Experiment", "experiment", "BioReplicate", "BioRep",
                     "sample_name", "Sample.Name", "Sample", "sample", "#NAME", "NAME")
  run.col <- run.col.cands[run.col.cands %in% colnames(meta)]
  if (length(run.col) == 0) {
    run.col <- colnames(meta)[1]
  } else {
    run.col <- run.col[1]
  }
  runs <- as.character(meta[[run.col]])

  # Condition and BioReplicate
  cond.col.cands <- c("Condition", "Group", "BioGroup", "condition", "group")
  cond.col <- cond.col.cands[cond.col.cands %in% colnames(meta)]
  condition <- if (length(cond.col) > 0) meta[[cond.col[1]]] else rep("Group1", length(runs))
  biorep.col <- c("BioReplicate", "BioRep", "Replicate", "replicate")[c("BioReplicate", "BioRep", "Replicate", "replicate") %in% colnames(meta)]
  biorep <- if (length(biorep.col) > 0) meta[[biorep.col[1]]] else seq_along(runs)

  meta.df <- data.frame(Condition = factor(condition), BioReplicate = factor(biorep), stringsAsFactors = FALSE)
  rownames(meta.df) <- runs
  disc.inx <- c(Condition = TRUE, BioReplicate = TRUE)
  cont.inx <- c(Condition = FALSE, BioReplicate = FALSE)

  # optional contaminant filtering (marked by "is contaminant" or similar)
  if (isTRUE(opts$removeContaminants)) {
    contam.cols <- grep("contaminant", tolower(colnames(dat)), value = TRUE)
    if (length(contam.cols) > 0) {
      keep <- dat[[contam.cols[1]]] != "+" & dat[[contam.cols[1]]] != TRUE
      if (any(keep)) {
        dat <- dat[keep, , drop = FALSE]
      }
    }
  }

  # optional evidence thresholds
  if (!is.null(opts$minProb) && !is.na(opts$minProb)) {
    prob.cols <- grep("(?i)^protein[ ._]probability$", colnames(dat), perl = TRUE, value = TRUE)
    if (length(prob.cols) > 0) {
      dat <- dat[as.numeric(dat[[prob.cols[1]]]) >= as.numeric(opts$minProb), , drop = FALSE]
    }
  }
  if (!is.null(opts$minPeptides) && !is.na(opts$minPeptides)) {
    pep.cols <- grep("(?i)^combined[ ._]total[ ._]peptides$", colnames(dat), perl = TRUE, value = TRUE)
    if (length(pep.cols) > 0) {
      dat <- dat[as.numeric(dat[[pep.cols[1]]]) >= as.numeric(opts$minPeptides), , drop = FALSE]
    }
  }

  # Intensity columns
  norm <- function(x) gsub("[^A-Za-z0-9]", "", tolower(x))
  runs.norm <- norm(runs)
  cols <- colnames(dat)
  cols.norm <- norm(cols)

  pick_col_for_run <- function(rn) {
    idx <- which(grepl(rn, cols.norm, fixed = TRUE))
    if (length(idx) == 0) return(NA_character_)
    qtype <- if (!is.null(opts$quantType)) opts$quantType else "protein_maxlfq"
    if (qtype == "protein_maxlfq") {
      maxlfq <- idx[grepl("maxlfq", cols.norm[idx])]
      if (length(maxlfq) > 0) return(cols[maxlfq[1]])
      inten <- idx[grepl("intensity", cols.norm[idx])]
      if (length(inten) > 0) return(cols[inten[1]])
      return(cols[idx[1]])
    } else if (qtype == "protein_intensity") {
      inten <- idx[grepl("intensity", cols.norm[idx])]
      if (length(inten) > 0) return(cols[inten[1]])
      maxlfq <- idx[grepl("maxlfq", cols.norm[idx])]
      if (length(maxlfq) > 0) return(cols[maxlfq[1]])
      return(cols[idx[1]])
    } else if (qtype == "ion_msstats") {
      # allow precursor/ion level: prefer Ion or Precursor columns if present
      ion <- idx[grepl("ion|precursor", cols.norm[idx])]
      if (length(ion) > 0) return(cols[ion[1]])
      inten <- idx[grepl("intensity", cols.norm[idx])]
      if (length(inten) > 0) return(cols[inten[1]])
      return(cols[idx[1]])
    } else {
      cols[idx[1]]
    }
  }

  match.cols <- vapply(runs.norm, pick_col_for_run, character(1))
  missing <- is.na(match.cols) | match.cols == ""
  if (any(missing)) {
    #msg("[FragPipe][WARN] could not find intensity columns for runs: ", paste(runs[missing], collapse = ", "))
    match.cols <- match.cols[!missing]
    runs <- runs[!missing]
    runs.norm <- runs.norm[!missing]
  }
  if (length(match.cols) == 0) {
    # fallback: any column starting with Intensity or LFQ.intensity
    match.cols <- grep("^(Intensity|LFQ.intensity)", colnames(dat), value = TRUE, ignore.case = TRUE)
  }
  if (length(match.cols) == 0) {
    #msg("[FragPipe][ERROR] could not match intensity columns to runs (checked normalized run names and Intensity*/LFQ.intensity*)")
    return(NULL)
  }

  # keep in the same order as runs vector (mapping run -> selected column)
  run_cols <- setNames(match.cols, runs)
  intens <- dat[, unname(run_cols), drop = FALSE]
  intens <- .safe_numeric_matrix(intens)
  id.col <- if ("Protein ID" %in% colnames(dat)) {
    dat[["Protein ID"]]
  } else if ("Protein" %in% colnames(dat)) {
    dat$Protein
  } else if ("Protein.IDs" %in% colnames(dat)) {
    dat$Protein.IDs
  } else if ("Protein.ID" %in% colnames(dat)) {
    dat$Protein.ID
  } else {
    dat[[1]]
  }
  rownames(intens) <- id.col
  colnames(intens) <- names(run_cols)[seq_len(ncol(intens))]

  # Extract peptide/spectral count for DEqMS
  pepcount <- NULL
  count.col.candidates <- c("Combined Total Peptides", "Total Peptides", "Peptides",
                            "Combined Spectral Count", "Spectral Count", "Total Spectral Count")
  for (col.name in count.col.candidates) {
    if (col.name %in% colnames(dat)) {
      pepcount <- as.numeric(dat[[col.name]])
      names(pepcount) <- id.col
      #msg("[FragPipe] Extracted ", col.name, " for DEqMS (range: ", min(pepcount, na.rm=TRUE), "-", max(pepcount, na.rm=TRUE), ")")
      break
    }
  }

  # Save FragPipe metadata for later filtering during normalization
  fp_metadata <- data.frame(
    Protein = id.col,
    stringsAsFactors = FALSE
  )

  # Helper function to find column with flexible name matching (space vs dot)
  find_column <- function(dat, candidates) {
    for (name in candidates) {
      if (name %in% colnames(dat)) return(name)
    }
    return(NULL)
  }

  # Add contaminant flag if present (check both space and dot versions)
  contam_col <- find_column(dat, c("Is Contaminant", "Is.Contaminant"))
  if (!is.null(contam_col)) {
    fp_metadata$Is.Contaminant <- dat[[contam_col]]
  }

  # Add protein probability if present (check both space and dot versions)
  prob_col <- find_column(dat, c("Protein Probability", "Protein.Probability"))
  if (!is.null(prob_col)) {
    fp_metadata$Protein.Probability <- as.numeric(dat[[prob_col]])
  }

  # Add peptide counts if present
  if (!is.null(pepcount)) {
    for (col.name in count.col.candidates) {
      if (col.name %in% colnames(dat)) {
        fp_metadata$Combined.Total.Peptides <- as.numeric(dat[[col.name]])
        break
      }
    }
  }

  qs::qsave(fp_metadata, "fragpipe_metadata.qs")
  #msg("[FragPipe] Saved metadata with ", nrow(fp_metadata), " proteins for later filtering")

  list(
    data = intens,
    data_orig = intens,
    type = "prot",
    format = "fragpipe-lfq",
    meta.info = list(meta.info = meta.df, disc.inx = disc.inx, cont.inx = cont.inx),
    pepcount = pepcount
  )
}

  remove_empty_cols <- function(df) {
    keep <- vapply(df, function(col) {
      any(nzchar(trimws(col)) & !is.na(col))
    }, logical(1))
    list(
      cleaned = df[ , keep, drop = FALSE],
      droppedNames = names(keep)[!keep]
    )
  }

GetAnalysisType <- function(){
  paramSet <- readSet(paramSet, "paramSet");
  return(paramSet$oneDataAnalType);
}

 GetViewData <- function(dataname){
    dat <- .readDataTable(dataname);   
        row.num <- nrow(dat);
        col.num <- ncol(dat);
        if(row.num > 100){
            row.num <- 100;
        }
        if(col.num > 10){
            col.num <- 10;
        }
        write.csv(dat[1:row.num, 1:col.num], file="raw_dataview.csv");
 
}


 


# Helper: read DIA-NN report.tsv (toy example from diann-rpackage)
.readDiannReport <- function(data.file, meta.file = "", opts = list(fileType = "protein_matrix", qvalueFilter = TRUE)) {
  if (!file.exists(data.file)) {
    #msg('[DIA-NN] data file not found: ', data.file)
    return(NULL)
  }
  dat <- tryCatch(read.delim(data.file, check.names = FALSE, stringsAsFactors = FALSE, comment.char = ""), error = function(e) NULL)
  if (is.null(dat)) {
    #msg('[DIA-NN][ERROR] unable to read data file')
    return(NULL)
  }

  file.type <- if (!is.null(opts$fileType)) opts$fileType else "protein_matrix"
  q.filter <- isTRUE(opts$qvalueFilter)

  # Read metadata file if provided
  meta <- NULL
  has.external.meta <- !is.null(meta.file) && meta.file != "" && file.exists(meta.file)
  if (has.external.meta) {
    #msg('[DIA-NN] Reading external metadata file: ', meta.file)
    # Try tab-separated first (most common for .tsv), then comma-separated
    meta <- tryCatch(read.delim(meta.file, check.names = FALSE, stringsAsFactors = FALSE, sep = "\t", comment.char = ""),
                     error = function(e) NULL)

    # Check if tab-separated reading resulted in only 1 column - likely a CSV file
    if (!is.null(meta) && ncol(meta) == 1) {
      #msg('[DIA-NN] Tab-separated read yielded only 1 column; trying comma-separated')
      meta <- tryCatch(read.delim(meta.file, check.names = FALSE, stringsAsFactors = FALSE, sep = ",", comment.char = ""),
                       error = function(e) NULL)
    }

    # If still null or failed, try auto-detect
    if (is.null(meta)) {
      meta <- tryCatch(read.table(meta.file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, comment.char = ""),
                       error = function(e) NULL)
    }

    if (is.null(meta)) {
      #msg('[DIA-NN][WARN] unable to read metadata file; will use defaults')
      has.external.meta <- FALSE
    }
  }

  # Check if this is matrix format (has Protein.Ids column + multiple intensity columns)
  cols <- colnames(dat)
  is.matrix.format <- "Protein.Ids" %in% cols &&
                      ("Protein.Names" %in% cols || "features" %in% cols) &&
                      ncol(dat) > 3  # More than just annotation columns
  if (file.type == "protein_matrix" && !is.matrix.format) {
    #msg("[DIA-NN] WARNING: fileType=protein_matrix but could not detect matrix format; attempting anyway")
  }

  if (is.matrix.format || file.type == "protein_matrix") {
    #msg('[DIA-NN] Detected matrix format with Protein.Ids column')

    # Use Protein.Ids as rownames (keep full IDs, don't split)
    if (!"Protein.Ids" %in% cols) {
      #msg('[DIA-NN][ERROR] Protein.Ids column not found in matrix format')
      return(NULL)
    }

    # Identify annotation columns to remove (more comprehensive list)
    annot.cols <- c("Protein.Ids", "Protein.Names", "features", "Protein.Groups",
                    "First.Protein.Description", "Proteotypic", "Stripped.Sequence",
                    "Modified.Sequence", "Precursor.Charge", "Q.Value", "PEP",
                    "Global.Q.Value", "Protein.Q.Value", "PG.Q.Value", "GG.Q.Value",
                    "Translated.Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value",
                    "Protein.Group", "Protein.Names", "Protein.Descriptions",
                    "Gene.Names", "Gene.Symbols", "Gene", "Description",
                    "Organism", "Species", "UniProt", "Accession", "Entry",
                    "Length", "Mass", "Sequence.Coverage", "Coverage",
                    "Unique.Peptides", "Peptides", "PSMs", "Razor.Peptides",
                    "MW", "Calc.pI", "Score", "iBAQ", "Fasta.headers")
    annot.cols.present <- intersect(annot.cols, cols)

    # Intensity columns are columns that:
    # 1. Are NOT in the annotation list
    # 2. Contain file paths (have "/" or "\") OR are all numeric
    candidate.cols <- setdiff(cols, annot.cols.present)

    intensity.cols <- sapply(candidate.cols, function(col) {
      col.data <- dat[[col]]
      # Check if column name looks like a file path
      has.path <- grepl("[/\\\\]", col)
      # Check if values are numeric (try converting first few non-NA values)
      sample.vals <- head(col.data[!is.na(col.data)], 10)
      if (length(sample.vals) == 0) return(FALSE)
      is.numeric.col <- !any(is.na(suppressWarnings(as.numeric(sample.vals))))
      # Accept if it's a file path OR all numeric
      return(has.path || is.numeric.col)
    })

    intensity.cols <- names(intensity.cols)[intensity.cols]

    if (length(intensity.cols) == 0) {
      #msg('[DIA-NN][ERROR] No intensity columns found in matrix format')
      #msg('[DIA-NN] Candidate columns were: ', paste(candidate.cols, collapse=", "))
      return(NULL)
    }

    #msg('[DIA-NN] Found ', length(intensity.cols), ' intensity columns')

    # Extract intensity matrix
    intens <- as.matrix(dat[, intensity.cols, drop = FALSE])

    # Split Protein.Ids and take first protein (leading protein)
    # This is consistent with the workflow assumption of one protein per row
    protein.ids <- sapply(as.character(dat$Protein.Ids), function(x) {
      parts <- strsplit(x, "[;|,]", perl = TRUE)[[1]]
      if (length(parts) == 0) return(NA_character_)
      parts[1]  # Take first/leading protein
    })

    # Remove rows with NA protein IDs
    valid.rows <- !is.na(protein.ids)
    if (any(!valid.rows)) {
      #msg('[DIA-NN] Removing ', sum(!valid.rows), ' rows with invalid Protein.Ids')
      intens <- intens[valid.rows, , drop = FALSE]
      protein.ids <- protein.ids[valid.rows]
    }

    rownames(intens) <- protein.ids

    # Convert to numeric and handle missing values
    storage.mode(intens) <- "numeric"

    # If there are duplicate protein IDs after splitting, aggregate them (sum intensities)
    if (any(duplicated(protein.ids))) {
      n.dup <- sum(duplicated(protein.ids))
      #msg('[DIA-NN] Found ', n.dup, ' duplicate protein IDs after splitting; aggregating by sum')
      intens.df <- as.data.frame(intens)
      intens.df$Protein <- protein.ids
      intens.agg <- stats::aggregate(. ~ Protein, data = intens.df, FUN = sum, na.rm = TRUE)
      rownames(intens.agg) <- intens.agg$Protein
      intens.agg$Protein <- NULL
      intens <- as.matrix(intens.agg)
    }

    # Remove rows that are all NA
    all.na.rows <- apply(intens, 1, function(x) all(is.na(x)))
    if (any(all.na.rows)) {
      #msg('[DIA-NN] Removing ', sum(all.na.rows), ' rows with all NA values')
      intens <- intens[!all.na.rows, , drop = FALSE]
    }

    #msg('[DIA-NN] Final matrix dimensions: ', nrow(intens), ' x ', ncol(intens))

    # Create or use metadata
    runs <- colnames(intens)

    if (has.external.meta) {
      #msg('[DIA-NN] Using external metadata file')

      # Identify run column in metadata
      run.col.cands <- c("Run", "Filename", "Raw.file", "Raw", "File.Name", "Experiment", "experiment",
                         "BioReplicate", "BioRep", "sample_name", "Sample.Name", "Sample", "sample",
                         "Name", "#Name", "X.NAME", "NAME")
      run.col <- run.col.cands[run.col.cands %in% colnames(meta)]
      if (length(run.col) == 0) {
        run.col <- colnames(meta)[1]
        #msg('[DIA-NN] Using first column as run identifier: ', run.col)
      } else {
        run.col <- run.col[1]
        #msg('[DIA-NN] Using column as run identifier: ', run.col)
      }

      meta.runs <- as.character(meta[[run.col]])

      # Try to match metadata runs to intensity column names
      # Match by basename or full path
      run.map <- sapply(runs, function(r) {
        # Try exact match first
        idx <- which(meta.runs == r)
        if (length(idx) > 0) return(idx[1])
        # Try basename match
        r.base <- basename(r)
        idx <- which(meta.runs == r.base)
        if (length(idx) > 0) return(idx[1])
        # Try basename without extension
        r.base.noext <- gsub("\\.[^.]*$", "", r.base)
        idx <- which(meta.runs == r.base.noext)
        if (length(idx) > 0) return(idx[1])
        # Try matching just the core filename (before first dot)
        r.core <- gsub("\\..*$", "", r.base)
        for (i in seq_along(meta.runs)) {
          meta.core <- gsub("\\..*$", "", basename(as.character(meta.runs[i])))
          if (r.core == meta.core) return(i)
        }
        return(NA)
      })

      # If all matches failed, provide detailed error message
      if (all(is.na(run.map))) {
        msg("[DIA-NN][ERROR] Could not match ANY data columns to metadata rows!")
        msg("[DIA-NN] Data column examples: ", paste(head(runs, 3), collapse=", "))
        msg("[DIA-NN] Metadata row examples: ", paste(head(meta.runs, 3), collapse=", "))
        stop("[DIA-NN] No matches found between data columns and metadata rows. Check that sample names match between files.")
      }

      if (any(is.na(run.map))) {
        msg('[DIA-NN][WARN] Could not match all runs to metadata. Matched: ', sum(!is.na(run.map)), '/', length(run.map))
        msg('[DIA-NN][WARN] Unmatched samples: ', paste(runs[is.na(run.map)], collapse=", "))
      }

      # Build metadata data frame with matching rows
      meta.df <- meta[run.map[!is.na(run.map)], , drop = FALSE]
      # Remove the run column
      meta.df[[run.col]] <- NULL
      rownames(meta.df) <- runs[!is.na(run.map)]

      # Subset intensity matrix to matched runs only
      if (any(is.na(run.map))) {
        intens <- intens[, !is.na(run.map), drop = FALSE]
        msg('[DIA-NN] Retained ', ncol(intens), ' matched samples out of ', length(runs), ' total')
      }

      # Determine which columns are discrete vs continuous
      if (ncol(meta.df) == 0 || nrow(meta.df) == 0) {
        msg("[DIA-NN][ERROR] Metadata is empty after processing (", nrow(meta.df), " rows, ", ncol(meta.df), " cols)")
        msg("[DIA-NN][ERROR] Original metadata had ", nrow(meta), " rows and ", ncol(meta), " columns")
        msg("[DIA-NN][ERROR] Run column used: ", run.col)
        stop("[DIA-NN] Metadata is empty after processing. Check that metadata file has valid data and sample names match between data and metadata files.")
      }

      disc.inx <- sapply(meta.df, function(col) {
        is.character(col) || is.factor(col) || (is.numeric(col) && length(unique(col)) <= 10)
      })

      # Ensure disc.inx is a properly named logical vector
      if (!is.logical(disc.inx) || is.null(names(disc.inx))) {
        names(disc.inx) <- colnames(meta.df)
      }

      cont.inx <- !disc.inx
      names(cont.inx) <- colnames(meta.df)

      # Convert discrete to factors
      if (any(disc.inx)) {
        for (i in which(disc.inx)) {
          if (!is.factor(meta.df[[i]])) {
            meta.df[[i]] <- factor(meta.df[[i]])
          }
        }
      }

    } else {
      # No external metadata; create defaults
      #msg('[DIA-NN] No external metadata; creating default metadata')
      meta.df <- data.frame(Condition = factor(rep('Group1', length(runs))), stringsAsFactors = FALSE)
      rownames(meta.df) <- runs
      disc.inx <- setNames(rep(TRUE, ncol(meta.df)), colnames(meta.df))
      cont.inx <- setNames(rep(FALSE, ncol(meta.df)), colnames(meta.df))
    }

    # Extract peptide/precursor count for DEqMS (from original data before filtering)
    pepcount <- NULL
    count.col.candidates <- c("Unique.Peptides", "Peptides", "PSMs", "Precursors")
    for (col.name in count.col.candidates) {
      if (col.name %in% annot.cols.present) {
        # Get count for valid rows only
        count.vals <- as.numeric(dat[[col.name]][valid.rows])
        # After aggregation, we need to map back to aggregated protein IDs
        if (any(duplicated(protein.ids))) {
          # If there were duplicates, sum the counts for aggregated proteins
          count.df <- data.frame(Protein = rownames(intens), stringsAsFactors = FALSE)
          count.df$Count <- 0  # Initialize
          # Match back to original protein IDs and sum
          for (i in seq_along(protein.ids)) {
            idx <- which(rownames(intens) == protein.ids[i])
            if (length(idx) > 0) {
              count.df$Count[idx] <- count.df$Count[idx] + count.vals[i]
            }
          }
          pepcount <- count.df$Count
          names(pepcount) <- count.df$Protein
        } else {
          pepcount <- count.vals
          names(pepcount) <- protein.ids
        }
        #msg('[DIA-NN] Extracted ', col.name, ' for DEqMS (range: ', min(pepcount, na.rm=TRUE), '-', max(pepcount, na.rm=TRUE), ')')
        break
      }
    }

    # Save DIA-NN metadata for later filtering during normalization
    # Note: For matrix format, Q-values are typically pre-filtered, but we save what we have
    diann_metadata <- data.frame(
      Protein.IDs = rownames(intens),
      stringsAsFactors = FALSE
    )
    if (!is.null(pepcount)) {
      diann_metadata$Peptide.Count <- pepcount[rownames(intens)]
    }
    qs::qsave(diann_metadata, "diann_metadata.qs")
    #msg("[DIA-NN] Saved metadata with ", nrow(diann_metadata), " proteins for later filtering")

    return(list(
      data = intens,
      data_orig = intens,
      type = 'prot',
      format = 'diann-matrix',
      idType = 'uniprot',
      organism = 'human',
      meta.info = list(meta.info = meta.df, disc.inx = disc.inx, cont.inx = cont.inx),
      pepcount = pepcount
    ))
  }

  # Original report.tsv format processing
  # Prefer MSstatsConvert/MSstats importer when available
  if (q.filter && "Q.Value" %in% colnames(dat)) {
    before <- nrow(dat)
    dat <- dat[dat$Q.Value < 0.01, , drop = FALSE]
    #msg("[DIA-NN] Q.Value filter < 0.01 removed ", before - nrow(dat), " rows")
  }
  if (!"Run" %in% colnames(dat)) dat$Run <- dat$File.Name
  if (!"BioReplicate" %in% colnames(dat)) dat$BioReplicate <- dat$Run
  if (!"Condition" %in% colnames(dat)) dat$Condition <- factor(make.names(basename(as.character(dat$Run))))
  if (!"LibPGQValue" %in% colnames(dat)) dat$LibPGQValue <- 0
  annot <- unique(dat[, c("Run", "BioReplicate", "Condition")])

  fmt <- NULL
  if (requireNamespace("MSstatsConvert", quietly = TRUE)) {
    fmt <- tryCatch({
      MSstatsConvert::DIANNtoMSstatsFormat(dat,
                                           annotation = annot,
                                           useUniquePeptide = TRUE,
                                           removeFewMeasurements = TRUE,
                                           summary = FALSE)
    }, error = function(e) NULL)
  } else if (requireNamespace("MSstats", quietly = TRUE)) {
    fmt <- tryCatch({
      MSstats::DIANNtoMSstatsFormat(dat,
                                    useUniquePeptide = TRUE,
                                    removeFewMeasurements = TRUE,
                                    summary = FALSE,
                                    annotation = annot)
    }, error = function(e) NULL)
  }

  if (!is.null(fmt)) {
    diann.df <- as.data.frame(fmt, stringsAsFactors = FALSE)
    diann.df$Intensity <- as.numeric(diann.df$Intensity)
    diann.df <- diann.df[!is.na(diann.df$Intensity), ]
    agg <- aggregate(Intensity ~ ProteinName + Run, data = diann.df, FUN = mean)
    wide <- reshape(agg, idvar = "ProteinName", timevar = "Run", direction = "wide")
    if (nrow(wide) == 0) {
      #msg("[DIA-NN][ERROR] MSstats returned empty data after filtering")
    } else if (ncol(wide) >= 2) {
      rownames(wide) <- wide$ProteinName
      wide$ProteinName <- NULL
      colnames(wide) <- sub("^Intensity\\.", "", colnames(wide))
      intens <- as.matrix(wide)

      meta.df <- unique(diann.df[, c("Run", "BioReplicate", "Condition")])
      rownames(meta.df) <- meta.df$Run
      meta.df$Run <- NULL
      meta.df$BioReplicate <- factor(meta.df$BioReplicate)
      meta.df$Condition <- factor(meta.df$Condition)
      disc.inx <- setNames(rep(TRUE, ncol(meta.df)), colnames(meta.df))
      cont.inx <- setNames(rep(FALSE, ncol(meta.df)), colnames(meta.df))

      diann_meta <- data.frame(Protein.IDs = rownames(intens), stringsAsFactors = FALSE)
      if ("Q.Value" %in% colnames(dat)) {
        qv <- stats::aggregate(as.numeric(dat$Q.Value), by = list(Protein = sapply(as.character(dat$Protein.Ids), function(x) strsplit(x, "[;|,]", perl = TRUE)[[1]][1])), FUN = min, na.rm = TRUE)
        diann_meta$Q.Value <- qv$x[match(diann_meta$Protein.IDs, qv$Protein)]
      }
      if ("PEP" %in% colnames(dat)) {
        pepv <- stats::aggregate(as.numeric(dat$PEP), by = list(Protein = sapply(as.character(dat$Protein.Ids), function(x) strsplit(x, "[;|,]", perl = TRUE)[[1]][1])), FUN = min, na.rm = TRUE)
        diann_meta$PEP <- pepv$x[match(diann_meta$Protein.IDs, pepv$Protein)]
      }
      count_col <- NULL
      for (cand in c("Unique.Peptides", "Peptides", "PSMs", "Precursors")) {
        if (cand %in% colnames(dat)) {
          count_col <- cand
          break
        }
      }
      if (!is.null(count_col)) {
        cnt <- stats::aggregate(as.numeric(dat[[count_col]]), by = list(Protein = sapply(as.character(dat$Protein.Ids), function(x) strsplit(x, "[;|,]", perl = TRUE)[[1]][1])), FUN = max, na.rm = TRUE)
        diann_meta$Peptide.Count <- cnt$x[match(diann_meta$Protein.IDs, cnt$Protein)]
      }
      qs::qsave(diann_meta, "diann_metadata.qs")

      return(list(
        data = intens,
        data_orig = intens,
        type = 'prot',
        format = 'diann-report',
        idType = 'uniprot',
        organism = 'human',
        meta.info = list(meta.info = meta.df, disc.inx = disc.inx, cont.inx = cont.inx)
      ))
    } else {
      #msg("[DIA-NN][WARN] MSstats returned <2 runs; falling back to manual parser.")
    }
  }

  # Fallback manual parse
  cols <- colnames(dat)
  run.vec <- dat$File.Name
  if (is.null(run.vec)) {
    #msg("[DIA-NN][ERROR] File.Name column missing")
    return(NULL)
  }
  qty_col <- NULL
  for (cand in c("features.MaxLFQ", "PG.Quantity", "features.Quantity", "Precursor.Quantity")) {
    if (cand %in% cols) { qty_col <- cand; break }
  }
  if (is.null(qty_col)) {
    #msg("[DIA-NN][ERROR] no quantity column found")
    return(NULL)
  }
  df <- data.frame(Protein = dat$Protein.Ids, Run = run.vec, Qty = as.numeric(dat[[qty_col]]), stringsAsFactors = FALSE)
  df$Protein <- sapply(as.character(df$Protein), function(x) {
    parts <- strsplit(x, "[;|,]", perl = TRUE)[[1]]
    if (length(parts) == 0) return(NA_character_)
    parts[1]
  })
  df <- df[!is.na(df$Protein) & !is.na(df$Qty), ]
  agg <- stats::aggregate(Qty ~ Protein + Run, data = df, FUN = mean)
  wide <- reshape(agg, idvar = "Protein", timevar = "Run", direction = "wide")
  rownames(wide) <- wide$Protein
  wide$Protein <- NULL
  colnames(wide) <- sub("^Qty\\.", "", colnames(wide))
  intens <- as.matrix(wide)
  diann_meta <- data.frame(Protein.IDs = rownames(intens), stringsAsFactors = FALSE)
  if ("Q.Value" %in% colnames(dat)) {
    qv <- stats::aggregate(as.numeric(dat$Q.Value), by = list(Protein = df$Protein), FUN = min, na.rm = TRUE)
    diann_meta$Q.Value <- qv$x[match(diann_meta$Protein.IDs, qv$Protein)]
  }
  if ("PEP" %in% colnames(dat)) {
    pepv <- stats::aggregate(as.numeric(dat$PEP), by = list(Protein = df$Protein), FUN = min, na.rm = TRUE)
    diann_meta$PEP <- pepv$x[match(diann_meta$Protein.IDs, pepv$Protein)]
  }
  count_col <- NULL
  for (cand in c("Unique.Peptides", "Peptides", "PSMs", "Precursors")) {
    if (cand %in% colnames(dat)) {
      count_col <- cand
      break
    }
  }
  if (!is.null(count_col)) {
    cnt <- stats::aggregate(as.numeric(dat[[count_col]]), by = list(Protein = df$Protein), FUN = max, na.rm = TRUE)
    diann_meta$Peptide.Count <- cnt$x[match(diann_meta$Protein.IDs, cnt$Protein)]
  }
  qs::qsave(diann_meta, "diann_metadata.qs")
  runs <- colnames(intens)
  meta.df <- data.frame(Condition = factor(rep('Group1', length(runs))), stringsAsFactors = FALSE)
  rownames(meta.df) <- runs
  disc.inx <- setNames(rep(TRUE, ncol(meta.df)), colnames(meta.df))
  cont.inx <- setNames(rep(FALSE, ncol(meta.df)), colnames(meta.df))
  list(
    data = intens,
    data_orig = intens,
    type = 'prot',
    format = 'diann-report',
    idType = 'uniprot',
    organism = 'human',
    meta.info = list(meta.info = meta.df, disc.inx = disc.inx, cont.inx = cont.inx)
  )
}

# Helper: read Spectronaut report
.readSpectronaut <- function(data.file, opts = list(inputType = "protein")) {
  if (!file.exists(data.file)) {
    #msg('[Spectronaut] data file not found: ', data.file)
    return(NULL)
  }

  pepcount <- NULL

  dat <- try(data.table::fread(data.file, header = TRUE, check.names = FALSE, data.table = FALSE, quote = "\""))
  if (inherits(dat, "try-error") || is.null(dat)) {
    #msg('[Spectronaut] Failed to read data file: ', data.file)
    return(NULL)
  }

  protein.id.candidates <- c('PG.ProteinAccessions', 'PG.ProteinGroups', 'PG.features',
                             'Accession', 'Protein.Group', 'ProteinName')
  get.spectronaut.protein.ids <- function(x, protein.col) {
    proteins <- as.character(x[[protein.col]])
    proteins <- vapply(strsplit(proteins, "[;|,]", perl = TRUE), function(parts) {
      if (length(parts) == 0) return(NA_character_)
      trimws(parts[1])
    }, character(1))
    proteins
  }

  build.spectronaut.metadata <- function(raw.dat, protein.col, protein.ids, keep.ids, pepcount = NULL) {
    spec_metadata <- data.frame(
      PG.ProteinGroups = keep.ids,
      stringsAsFactors = FALSE
    )

    valid <- !is.na(protein.ids) & protein.ids != ""
    if (!any(valid)) {
      if (!is.null(pepcount)) {
        spec_metadata$Peptide.Count <- pepcount[keep.ids]
      }
      return(spec_metadata)
    }

    protein.ids <- protein.ids[valid]
    raw.sub <- raw.dat[valid, , drop = FALSE]

    if ("PG.IsContaminant" %in% colnames(raw.sub)) {
      cont.df <- stats::aggregate(as.logical(raw.sub$PG.IsContaminant), by = list(Protein = protein.ids), FUN = function(x) any(isTRUE(x), na.rm = TRUE))
      spec_metadata$PG.IsContaminant <- cont.df$x[match(keep.ids, cont.df$Protein)]
      spec_metadata$PG.IsContaminant[is.na(spec_metadata$PG.IsContaminant)] <- FALSE
    }

    qval.col <- if ("PG.Qvalue" %in% colnames(raw.sub)) {
      "PG.Qvalue"
    } else if ("Qvalue" %in% colnames(raw.sub)) {
      "Qvalue"
    } else if ("EG.Qvalue" %in% colnames(raw.sub)) {
      "EG.Qvalue"
    } else {
      NULL
    }
    if (!is.null(qval.col)) {
      qvals <- suppressWarnings(as.numeric(raw.sub[[qval.col]]))
      qval.df <- stats::aggregate(qvals, by = list(Protein = protein.ids), FUN = function(x) suppressWarnings(min(x, na.rm = TRUE)))
      spec_metadata$PG.Qvalue <- qval.df$x[match(keep.ids, qval.df$Protein)]
      spec_metadata$PG.Qvalue[!is.finite(spec_metadata$PG.Qvalue)] <- NA_real_
    }

    if (!is.null(pepcount)) {
      spec_metadata$Peptide.Count <- pepcount[keep.ids]
    } else {
      pep.col <- if ("PG.NrOfStrippedSequencesIdentified" %in% colnames(raw.sub)) {
        "PG.NrOfStrippedSequencesIdentified"
      } else if ("Stripped.Sequence" %in% colnames(raw.sub)) {
        "Stripped.Sequence"
      } else if ("Sequence" %in% colnames(raw.sub)) {
        "Sequence"
      } else if ("EG.ModifiedSequence" %in% colnames(raw.sub)) {
        "EG.ModifiedSequence"
      } else if ("ModifiedSequence" %in% colnames(raw.sub)) {
        "ModifiedSequence"
      } else if ("PG.Quantity" %in% colnames(raw.sub)) {
        NULL
      } else {
        NULL
      }
      if (!is.null(pep.col)) {
        if (pep.col %in% c("Stripped.Sequence", "Sequence", "EG.ModifiedSequence", "ModifiedSequence")) {
          pepvals <- as.character(raw.sub[[pep.col]])
          pep.df <- stats::aggregate(pepvals, by = list(Protein = protein.ids), FUN = function(x) length(unique(stats::na.omit(x))))
        } else {
          pepvals <- suppressWarnings(as.numeric(raw.sub[[pep.col]]))
          pep.df <- stats::aggregate(pepvals, by = list(Protein = protein.ids), FUN = function(x) suppressWarnings(max(x, na.rm = TRUE)))
        }
        spec_metadata$Peptide.Count <- pep.df$x[match(keep.ids, pep.df$Protein)]
        spec_metadata$Peptide.Count[!is.finite(spec_metadata$Peptide.Count)] <- NA_real_
      }
    }

    spec_metadata
  }

  read.spectronaut.peptide <- function(raw.dat) {
    peptide.col <- if ("ModifiedSequence" %in% colnames(raw.dat)) {
      "ModifiedSequence"
    } else if ("EG.ModifiedSequence" %in% colnames(raw.dat)) {
      "EG.ModifiedSequence"
    } else if ("Sequence" %in% colnames(raw.dat)) {
      "Sequence"
    } else if ("Stripped.Sequence" %in% colnames(raw.dat)) {
      "Stripped.Sequence"
    } else {
      NULL
    }

    protein.col <- protein.id.candidates[protein.id.candidates %in% colnames(raw.dat)][1]
    if (is.null(peptide.col) || is.na(protein.col) || !nzchar(protein.col)) {
      return(NULL)
    }

    run.col <- if ("R.FileName" %in% colnames(raw.dat)) {
      "R.FileName"
    } else if ("Run" %in% colnames(raw.dat)) {
      "Run"
    } else {
      NULL
    }

    quantity.col <- c("FG.Quantity", "F.PeakArea", "PeakArea", "EG.Quantity", "Quantity")
    quantity.col <- quantity.col[quantity.col %in% colnames(raw.dat)][1]

    if (!is.na(run.col) && !is.null(run.col) && !is.na(quantity.col) && !is.null(quantity.col)) {
      long.df <- data.frame(
        Peptide = as.character(raw.dat[[peptide.col]]),
        Protein = get.spectronaut.protein.ids(raw.dat, protein.col),
        Run = as.character(raw.dat[[run.col]]),
        Intensity = suppressWarnings(as.numeric(raw.dat[[quantity.col]])),
        stringsAsFactors = FALSE
      )
      long.df <- long.df[!is.na(long.df$Peptide) & nzchar(long.df$Peptide) &
                           !is.na(long.df$Protein) & nzchar(long.df$Protein) &
                           !is.na(long.df$Run) & nzchar(long.df$Run) &
                           is.finite(long.df$Intensity), , drop = FALSE]
      if (nrow(long.df) == 0) {
        return(NULL)
      }

      agg <- stats::aggregate(Intensity ~ Peptide + Run, data = long.df, FUN = median, na.rm = TRUE)
      if (!requireNamespace("reshape2", quietly = TRUE)) {
        return(NULL)
      }
      wide <- reshape2::dcast(agg, Peptide ~ Run, value.var = "Intensity")
      intens <- as.matrix(wide[, -1, drop = FALSE])
      rownames(intens) <- make.unique(as.character(wide$Peptide))
      storage.mode(intens) <- "numeric"

      prot.map <- unique(long.df[, c("Peptide", "Protein"), drop = FALSE])
      prot.map <- prot.map[match(rownames(intens), prot.map$Peptide), , drop = FALSE]
      runs <- colnames(intens)
      meta.df <- data.frame(Condition = factor(runs), stringsAsFactors = FALSE)
      rownames(meta.df) <- runs
      disc.inx <- setNames(rep(TRUE, ncol(meta.df)), colnames(meta.df))
      cont.inx <- setNames(rep(FALSE, ncol(meta.df)), colnames(meta.df))

      return(list(
        data = intens,
        data_orig = intens,
        type = "peptide",
        format = "spectronaut-peptide",
        meta.info = list(meta.info = meta.df, disc.inx = disc.inx, cont.inx = cont.inx),
        prot.map = prot.map
      ))
    }

    annot.cols <- c('PG.ProteinGroups', 'PG.ProteinAccessions', 'PG.features', 'PG.ProteinNames',
                    'PG.ProteinDescriptions', 'PG.Organisms', 'PG.FastaFiles',
                    'EG.ModifiedSequence', 'EG.Cscore', 'EG.Qvalue', 'FG.Quantity',
                    'PG.Quantity', 'PG.NrOfStrippedSequencesIdentified',
                    'Experiment', 'Run', 'Accession', 'ModifiedSequence', 'Sequence', 'Precursor',
                    'PrecursorCharge', 'Protein.Group', 'ProteinName')
    intensity.cols <- grep("^[^P].*\\.[HEW]", colnames(raw.dat), value = TRUE, ignore.case = FALSE)
    if (length(intensity.cols) == 0) {
      candidate.cols <- setdiff(colnames(raw.dat), annot.cols)
      candidate.cols <- grep("^RtoF|^Ratio|^FC_", candidate.cols, value = TRUE, invert = TRUE, ignore.case = TRUE)
      intensity.cols <- sapply(candidate.cols, function(col) {
        col.data <- raw.dat[[col]]
        sample.vals <- head(col.data[!is.na(col.data)], 10)
        !any(is.na(suppressWarnings(as.numeric(sample.vals))))
      })
      intensity.cols <- names(intensity.cols[intensity.cols])
    }
    if (length(intensity.cols) == 0) {
      return(NULL)
    }

    intens <- as.matrix(raw.dat[, intensity.cols, drop = FALSE])
    storage.mode(intens) <- "numeric"
    peptide.ids <- as.character(raw.dat[[peptide.col]])
    proteins <- get.spectronaut.protein.ids(raw.dat, protein.col)
    valid.rows <- !is.na(peptide.ids) & nzchar(peptide.ids) & !is.na(proteins) & nzchar(proteins)
    if (!any(valid.rows)) {
      return(NULL)
    }
    intens <- intens[valid.rows, , drop = FALSE]
    peptide.ids <- make.unique(peptide.ids[valid.rows])
    proteins <- proteins[valid.rows]
    rownames(intens) <- peptide.ids
    prot.map <- data.frame(Peptide = peptide.ids, Protein = proteins, stringsAsFactors = FALSE)
    runs <- colnames(intens)
    meta.df <- data.frame(Condition = factor(runs), stringsAsFactors = FALSE)
    rownames(meta.df) <- runs
    disc.inx <- setNames(rep(TRUE, ncol(meta.df)), colnames(meta.df))
    cont.inx <- setNames(rep(FALSE, ncol(meta.df)), colnames(meta.df))

    list(
      data = intens,
      data_orig = intens,
      type = "peptide",
      format = "spectronaut-peptide",
      meta.info = list(meta.info = meta.df, disc.inx = disc.inx, cont.inx = cont.inx),
      prot.map = prot.map
    )
  }

  # If no obvious protein column is present, Spectronaut header likely misaligned (e.g., missing a column)
  if (!any(protein.id.candidates %in% colnames(dat))) {
    dat <- tryCatch({
      raw.lines <- readLines(data.file, n = 2)
      header.tokens.orig <- if (length(raw.lines) >= 1) strsplit(raw.lines[1], "\t")[[1]] else character(0)
      header.tokens.orig <- gsub('^\"|\"$', '', header.tokens.orig)

      first.data.tokens <- if (length(raw.lines) >= 2) strsplit(raw.lines[2], "\t")[[1]] else header.tokens.orig
      data.len <- length(first.data.tokens)

      build.header <- function(tokens, insert.after, filler) {
        if (length(tokens) == 0) return(NULL)
        insert.after <- ifelse(length(insert.after) == 0, length(tokens), insert.after[1])
        if (data.len > length(tokens)) {
          tokens <- append(tokens, rep(filler, data.len - length(tokens)), after = insert.after)
        } else if (data.len < length(tokens)) {
          tokens <- tokens[seq_len(data.len)]
        }
        tokens
      }

      header.tries <- list(
        build.header(header.tokens.orig, which(header.tokens.orig %in% c('Precursor', 'Precursor.Id', 'ModifiedSequence', 'Sequence')), 'PrecursorCharge'),
        build.header(header.tokens.orig, 1, 'Run')
      )

      dat.repaired <- NULL
      for (hdr in header.tries) {
        if (is.null(hdr)) next
        dat.try <- try(data.table::fread(
          data.file,
          header = FALSE,
          skip = 1,
          col.names = hdr,
          check.names = FALSE,
          data.table = FALSE,
          quote = "\"",
          fill = TRUE
        ), silent = TRUE)

        protein.matches <- intersect(protein.id.candidates, colnames(dat.try))
        has.valid.protein.ids <- FALSE
        if (!inherits(dat.try, "try-error") && length(protein.matches) > 0) {
          uniq.counts <- sapply(protein.matches, function(col) length(unique(stats::na.omit(dat.try[[col]]))))
          has.valid.protein.ids <- any(uniq.counts > 1)
        }

        if (!inherits(dat.try, "try-error") && has.valid.protein.ids) {
          #msg('[Spectronaut] Repairing header mismatch: ', length(header.tokens.orig), ' -> ', length(hdr), ' columns')
          dat.repaired <- dat.try
          break
        }
      }

      dat.repaired
    }, error = function(e) {
      #msg('[Spectronaut] Failed to repair Spectronaut header: ', e$message)
      NULL
    })

    if (inherits(dat, "try-error") || is.null(dat)) {
      return(NULL)
    }
  }

  #msg('[Spectronaut] Read data table with ', nrow(dat), ' rows, ', ncol(dat), ' columns')

  input.type <- if (!is.null(opts$inputType)) tolower(opts$inputType) else "protein"
  if (identical(input.type, "peptide")) {
    peptide.res <- read.spectronaut.peptide(dat)
    if (!is.null(peptide.res)) {
      return(peptide.res)
    }
  }

  # Check if this is a pivot table format (peptides x samples) or long format
  is.pivot.format <- ("Accession" %in% colnames(dat) || "Protein.Group" %in% colnames(dat)) &&
                      ("Sequence" %in% colnames(dat) || "ModifiedSequence" %in% colnames(dat)) &&
                      !("R.FileName" %in% colnames(dat))

  if (is.pivot.format) {
    #msg('[Spectronaut] Pivot table format detected; using manual parser')
  }

  # Prefer MSstatsConvert/MSstats pipeline when available and data is in long format
  if (!is.pivot.format &&
      requireNamespace("MSstatsConvert", quietly = TRUE) &&
      requireNamespace("MSstats", quietly = TRUE)) {

    #msg('[Spectronaut] Attempting MSstats pipeline')

    # Prepare annotation if needed
    # Spectronaut reports typically have R.Condition, R.Replicate columns
    # If not present, infer from R.FileName
    annot <- NULL
    if ("R.Condition" %in% colnames(dat) && "R.Replicate" %in% colnames(dat)) {
      annot <- unique(dat[, c("R.FileName", "R.Condition", "R.Replicate")])
      colnames(annot) <- c("Run", "Condition", "BioReplicate")
    } else if ("R.FileName" %in% colnames(dat)) {
      runs <- unique(dat$R.FileName)
      conditions <- gsub("_[0-9]+$", "", runs)
      conditions <- gsub("_DIA.*$", "", conditions, ignore.case = TRUE)
      annot <- data.frame(Run = runs,
                          Condition = conditions,
                          BioReplicate = runs,
                          stringsAsFactors = FALSE)
    }

    msstats.fmt <- tryCatch({
      MSstatsConvert::SpectronauttoMSstatsFormat(
        input = dat,
        annotation = annot,
        intensity = 'PeakArea',  # or 'NormalizedArea' depending on export
        filter_with_Qvalue = FALSE,
        qvalue_cutoff = 0.01,
        useUniquePeptide = TRUE,
        removeFewMeasurements = TRUE,
        removeProtein_with1Feature = FALSE
      )
    }, error = function(e) {
      #msg('[Spectronaut] MSstatsConvert failed: ', e$message)
      NULL
    })

    if (!is.null(msstats.fmt)) {
      # Save msstats input for later use (normalization, imputation, summarization in ProteoAnalyst pipeline)
      qs::qsave(msstats.fmt, "msstats_input.qs")
      #msg('[Spectronaut] MSstats format saved with ', length(unique(msstats.fmt$ProteinName)), ' proteins')

      # Simple aggregation to protein level (median intensity per protein per run)
      # Full MSstats processing (normalization, TMP summarization) will be done later
      if (requireNamespace("reshape2", quietly = TRUE)) {
        # Aggregate peptide intensities to protein level using median
        agg <- stats::aggregate(Intensity ~ ProteinName + Run, data = msstats.fmt, FUN = median, na.rm = TRUE)
        wide <- reshape2::dcast(agg, ProteinName ~ Run, value.var = "Intensity")
        intens <- as.matrix(wide[, -1, drop = FALSE])
        rownames(intens) <- wide$ProteinName
        storage.mode(intens) <- "numeric"

        runs <- colnames(intens)
        # Extract condition info from MSstats annotation
        if (!is.null(annot)) {
          meta.map <- unique(msstats.fmt[, c("Run", "Condition")])
          rownames(meta.map) <- meta.map$Run
          conditions <- meta.map[runs, "Condition"]
        } else {
          conditions <- gsub("_[0-9]+$", "", runs)
          conditions <- gsub("_DIA.*$", "", conditions, ignore.case = TRUE)
        }

        meta.df <- data.frame(Condition = factor(conditions), stringsAsFactors = FALSE)
        rownames(meta.df) <- runs
        disc.inx <- setNames(rep(TRUE, ncol(meta.df)), colnames(meta.df))
        cont.inx <- setNames(rep(FALSE, ncol(meta.df)), colnames(meta.df))

        # Extract peptide count per protein for DEqMS
        pepcount <- NULL
        if ("PeptideSequence" %in% colnames(msstats.fmt) || "FullPeptideName" %in% colnames(msstats.fmt)) {
          pep.col <- if ("PeptideSequence" %in% colnames(msstats.fmt)) "PeptideSequence" else "FullPeptideName"
          pep.counts <- stats::aggregate(as.formula(paste(pep.col, "~ ProteinName")),
                                        data = msstats.fmt,
                                        FUN = function(x) length(unique(x)))
          colnames(pep.counts)[2] <- "Count"
          pepcount <- pep.counts$Count
          names(pepcount) <- pep.counts$ProteinName
          pepcount <- pepcount[rownames(intens)]  # Match to final protein order
          #msg('[Spectronaut] Extracted unique peptide counts for DEqMS (range: ', min(pepcount, na.rm=TRUE), '-', max(pepcount, na.rm=TRUE), ')')
        }

        protein.col <- protein.id.candidates[protein.id.candidates %in% colnames(dat)][1]
        if (!is.na(protein.col) && nzchar(protein.col)) {
          spec_metadata <- build.spectronaut.metadata(dat, protein.col, get.spectronaut.protein.ids(dat, protein.col), rownames(intens), pepcount)
          qs::qsave(spec_metadata, "spectronaut_metadata.qs")
        }

        #msg('[Spectronaut] MSstatsConvert successful: ', nrow(intens), ' proteins, ', ncol(intens), ' runs')
        return(list(
          data = intens,
          data_orig = intens,
          type = 'prot',
          format = 'spectronaut-msstats',
          idType = 'uniprot',
          pepcount = pepcount,
          organism = 'ecoli',
          meta.info = list(meta.info = meta.df, disc.inx = disc.inx, cont.inx = cont.inx)
        ))
      }
    }

    #msg('[Spectronaut] MSstats pipeline failed; falling back to manual parser')
  }

  # Fallback manual parser

  # Identify annotation columns (typically at the start)
  annot.cols <- c('PG.ProteinGroups', 'PG.ProteinAccessions', 'PG.features', 'PG.ProteinNames',
                  'PG.ProteinDescriptions', 'PG.Organisms', 'PG.FastaFiles',
                  'EG.ModifiedSequence', 'EG.Cscore', 'EG.Qvalue', 'FG.Quantity',
                  'PG.Quantity', 'PG.NrOfStrippedSequencesIdentified',
                  'Experiment', 'Run', 'Accession', 'ModifiedSequence', 'Sequence', 'Precursor',
                  'PrecursorCharge', 'Protein.Group', 'ProteinName')

  # Find intensity columns (typically end with sample names and contain numeric data)
  # Spectronaut directLFQ columns or regular intensity columns
  intensity.cols <- grep("^[^P].*\\.[HEW]", colnames(dat), value = TRUE, ignore.case = FALSE)

  # If no intensity columns found with that pattern, try broader search
  if (length(intensity.cols) == 0) {
    candidate.cols <- setdiff(colnames(dat), annot.cols)
    # Exclude ratio/calculated columns
    candidate.cols <- grep("^RtoF|^Ratio|^FC_", candidate.cols, value = TRUE, invert = TRUE, ignore.case = TRUE)
    intensity.cols <- sapply(candidate.cols, function(col) {
      col.data <- dat[[col]]
      sample.vals <- head(col.data[!is.na(col.data)], 10)
      is.numeric.col <- !any(is.na(suppressWarnings(as.numeric(sample.vals))))
      return(is.numeric.col)
    })
    intensity.cols <- names(intensity.cols[intensity.cols])
  }

  if (length(intensity.cols) == 0) {
    #msg('[Spectronaut][ERROR] No intensity columns found')
    return(NULL)
  }

  #msg('[Spectronaut] Found ', length(intensity.cols), ' intensity columns')

  # Use protein identifier column (try multiple common Spectronaut column names)
  if ('PG.ProteinAccessions' %in% colnames(dat)) {
    protein.col <- 'PG.ProteinAccessions'
  } else if ('PG.ProteinGroups' %in% colnames(dat)) {
    protein.col <- 'PG.ProteinGroups'
  } else if ('PG.features' %in% colnames(dat)) {
    protein.col <- 'PG.features'
  } else if ('Accession' %in% colnames(dat)) {
    protein.col <- 'Accession'
  } else if ('Protein.Group' %in% colnames(dat)) {
    protein.col <- 'Protein.Group'
  } else if ('ProteinName' %in% colnames(dat)) {
    protein.col <- 'ProteinName'
  } else {
    #msg('[Spectronaut][ERROR] No protein identifier column found. Available columns: ', paste(colnames(dat)[1:min(10, ncol(dat))], collapse=', '))
    return(NULL)
  }

  # Extract intensity matrix
  intens <- as.matrix(dat[, intensity.cols, drop = FALSE])
  storage.mode(intens) <- 'numeric'

  # Set row names from protein identifiers, handle semicolon-separated lists
  proteins <- get.spectronaut.protein.ids(dat, protein.col)
  raw.protein.ids <- proteins

  # Remove rows with missing protein IDs
  valid.rows <- !is.na(proteins) & proteins != ""
  if (sum(valid.rows) == 0) {
    #msg('[Spectronaut][ERROR] No valid protein IDs found')
    return(NULL)
  }

  intens <- intens[valid.rows, , drop = FALSE]
  proteins <- proteins[valid.rows]
  raw.protein.ids <- raw.protein.ids[valid.rows]

  # Handle duplicate protein IDs by taking the maximum intensity
  if (any(duplicated(proteins))) {
    #msg('[Spectronaut] Aggregating ', sum(duplicated(proteins)), ' duplicate protein IDs')
    intens <- as.data.frame(intens)
    intens$Protein <- proteins
    intens <- stats::aggregate(. ~ Protein, data = intens, FUN = function(x) max(x, na.rm = TRUE))
    proteins <- intens$Protein
    intens$Protein <- NULL
    intens <- as.matrix(intens)
  }

  rownames(intens) <- proteins

  # Infer metadata from sample names
  runs <- colnames(intens)
  #msg('[Spectronaut] Creating metadata for ', length(runs), ' samples')

  # Try to extract condition from sample names (e.g., HEof_n600_DIA_1 -> HEof_n600)
  conditions <- gsub("_[0-9]+$", "", runs)  # Remove trailing numbers
  conditions <- gsub("_DIA.*$", "", conditions, ignore.case = TRUE)  # Remove DIA suffix
  conditions <- factor(conditions)

  meta.df <- data.frame(Condition = conditions, stringsAsFactors = FALSE)
  rownames(meta.df) <- runs

  disc.inx <- setNames(rep(TRUE, ncol(meta.df)), colnames(meta.df))
  cont.inx <- setNames(rep(FALSE, ncol(meta.df)), colnames(meta.df))

  # Save Spectronaut metadata for later filtering during normalization
  spec_metadata <- build.spectronaut.metadata(dat[valid.rows, , drop = FALSE], protein.col, raw.protein.ids, rownames(intens), pepcount)
  qs::qsave(spec_metadata, "spectronaut_metadata.qs")
  #msg("[Spectronaut] Saved metadata with ", nrow(spec_metadata), " proteins for later filtering")

  return(list(
    data = intens,
    data_orig = intens,
    type = 'prot',
    format = 'spectronaut',
    idType = 'uniprot',
    organism = 'ecoli',
    pepcount = pepcount,
    meta.info = list(meta.info = meta.df, disc.inx = disc.inx, cont.inx = cont.inx)
  ))
}


# Record Phospho options from UI (stored in R options + paramSet)
SetPhosphoOptions <- function(format = "diann_report", ptm = "phospho", locProb = 0.75, residues = c("S","T","Y")) {
  paramSet <- readSet(paramSet, "paramSet")
  opts <- list(
    format = format,
    ptm = ptm,
    locProb = suppressWarnings(as.numeric(locProb)),
    residues = residues
  )
  options(pa.phospho.opts = opts)
  paramSet$phospho <- opts
  saveSet(paramSet, "paramSet")
  #msg("[Phospho] options set: format=", opts$format,
  #        " ptm=", opts$ptm,
  #        " locProb=", opts$locProb,
  #        " residues=", paste(residues, collapse=","))
  return(1L)
}

# Helper: read FragPipe combined_peptide.tsv + experiment annotation
.readFragpipePeptide <- function(data.file, meta.file, opts = list(quantType = "intensity", removeContaminants = TRUE)) {
  
  if (!file.exists(data.file)) {
    #msg("[FragPipe] data file not found: ", data.file)
    return(NULL)
  }
  
  # 1. Read Data
  # check.names=FALSE is crucial to keep "Sample Intensity" headers intact
  dat <- tryCatch(read.delim(data.file, check.names = FALSE, stringsAsFactors = FALSE, comment.char = ""), error = function(e) NULL)
  if (is.null(dat)) return(NULL)
  
  # 2. Read Metadata
  if (is.null(meta.file) || meta.file == "" || !file.exists(meta.file)) {
    #msg("[FragPipe][ERROR] metadata file is required.")
    return(NULL)
  }
  meta <- tryCatch(read.delim(meta.file, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(meta)) return(NULL)

  # 3. Handle Metadata Columns (Run, Condition, BioRep)
  # (Kept identical to your logic for consistency)
  run.col.cands <- c("Run", "Filename", "Raw.file", "Raw", "Experiment", "experiment", "sample_name", "Sample.Name", "Sample", "#NAME", "NAME")
  run.col <- run.col.cands[run.col.cands %in% colnames(meta)][1]
  runs <- as.character(meta[[run.col]])

  cond.col.cands <- c("Condition", "Group", "BioGroup", "condition", "group")
  cond.col <- cond.col.cands[cond.col.cands %in% colnames(meta)][1]
  condition <- if (!is.na(cond.col)) meta[[cond.col]] else rep("Group1", length(runs))
  
  biorep.col.cands <- c("BioReplicate", "BioRep", "Replicate", "replicate")
  biorep.col <- biorep.col.cands[biorep.col.cands %in% colnames(meta)][1]
  if (length(biorep.col) == 0 || is.na(biorep.col)) {
    biorep <- seq_along(runs)
  } else {
    biorep <- meta[[biorep.col]]
  }

  meta.df <- data.frame(Condition = factor(condition), BioReplicate = factor(biorep), stringsAsFactors = FALSE)
  rownames(meta.df) <- runs
  
  # 4. Filter Contaminants
  # FragPipe peptide files usually mark contaminants in "Protein" column (e.g. contam_sp|...) 
  # or "Mapped Proteins".
  if (isTRUE(opts$removeContaminants)) {
    # Check Protein column for "CONT_" or "contam"
    prot_col <- grep("Protein", colnames(dat), value = TRUE)[1]
    if (!is.na(prot_col)) {
        # Common FragPipe contaminant flags
        is_contam <- grepl("CONT_|REV_|contam", dat[[prot_col]])
        dat <- dat[!is_contam, , drop = FALSE]
    }
  }

  # 5. Map Run Names to Data Columns
  # FragPipe combined_peptide.tsv columns usually look like: "Experiment1 Intensity" or "Experiment1 MaxLFQ Intensity"
  
  # Normalize function for matching
  norm <- function(x) gsub("[^A-Za-z0-9]", "", tolower(x))
  cols <- colnames(dat)
  cols.norm <- norm(cols)
  runs.norm <- norm(runs)
  
  pick_col_for_run <- function(rn) {
    # Find columns containing the run name
    idx <- which(grepl(rn, cols.norm, fixed = TRUE))
    if (length(idx) == 0) return(NA_character_)
    
    # Filter for Intensity columns
    # We look for "intensity" string in the column name
    int_idx <- idx[grepl("intensity", cols.norm[idx])]
    
    if (length(int_idx) == 0) return(NA_character_)
    
    cand_cols <- cols[int_idx]
    
    # Specific Logic for LFQ vs Intensity
    # If user wants MaxLFQ, prioritize columns with "MaxLFQ"
    if (!is.null(opts$quantType) && opts$quantType == "maxlfq") {
        lfq_match <- grep("maxlfq", tolower(cand_cols), value = TRUE)
        if (length(lfq_match) > 0) return(lfq_match[1])
    }
    
    # Default: Return the shortest match containing "Intensity" (usually "Exp1 Intensity")
    # or just the first one found
    return(cand_cols[1])
  }

  match.cols <- vapply(runs.norm, pick_col_for_run, character(1))
  
  # Error handling for missing columns
  missing <- is.na(match.cols) | match.cols == ""
  if (any(missing)) {
    #msg("[FragPipe][WARN] Missing columns for runs: ", paste(runs[missing], collapse = ", "))
    # Filter metadata to match available data
    match.cols <- match.cols[!missing]
    runs <- runs[!missing]
    meta.df <- meta.df[!missing, , drop=FALSE]
  }
  if (length(match.cols) == 0) {
    #msg("[FragPipe][ERROR] No intensity columns matched any run names. Please verify run naming in the FragPipe report vs metadata.")
    return(NULL)
  }

  # 6. Extract Intensity Matrix
  run_cols <- setNames(match.cols, runs)
  intens <- dat[, unname(run_cols), drop = FALSE]
  
  # Ensure numeric
  intens <- .safe_numeric_matrix(intens)
  
  # 7. Set Row IDs (Peptide Sequence)
  # Prefer "Modified Sequence" if available to distinguish PTMs, otherwise "Peptide Sequence"
  if ("Modified Sequence" %in% colnames(dat)) {
    id.col <- dat[["Modified Sequence"]]
  } else if ("Peptide Sequence" %in% colnames(dat)) {
    id.col <- dat[["Peptide Sequence"]]
  } else {
    id.col <- rownames(dat)
  }
  
  # Handle duplicate peptides (if any remain, though combined_peptide should be unique)
  if (any(duplicated(id.col))) {
    #msg("[FragPipe] Warning: Duplicate peptide sequences found. Making IDs unique.")
    id.col <- make.unique(as.character(id.col))
  }
  
  rownames(intens) <- id.col
  colnames(intens) <- names(run_cols)

  # 8. Create Peptide-to-Protein Map (Vital for your new modular step)
  # We extract this so you can use it in the 'Summarization' step later
  prot.col.name <- grep("Protein", colnames(dat), value=TRUE)[1]
  prot_map <- NULL
  if (!is.na(prot.col.name)) {
      prot_map <- data.frame(
          Peptide = id.col,
          Protein = dat[[prot.col.name]],
          stringsAsFactors = FALSE
      )
  }

  # Return compatible list structure
  list(
    data = intens,             # Wide Matrix (Peptide x Sample)
    data_orig = intens,
    type = "peptide",          # Flag as peptide
    format = "fragpipe-peptide",
    meta.info = list(meta.info = meta.df, disc.inx = c(Condition=TRUE, BioReplicate=TRUE), cont.inx = c(Condition=FALSE, BioReplicate=FALSE)),
    prot.map = prot_map        # New field: Peptide-Protein Mapping
  )
}
