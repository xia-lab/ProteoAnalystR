##################################################
## R script for ProteoAnalyst
## Description: functions for id annotation onedata and metadata
##
## Authors: 
## Jeff Xia, jeff.xia@mcgill.ca
## Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

# read in the data and perform
# gene ID mapping using built in libraries
# matchMin is minimal matched probe (%)
# idType: INVEX supported ID types
# lvlOpt: "NA" to keep original, other values will merge original ID to entrez gene IDs

# return the total matched gene number
# note: unmapped IDs will be retained as 
# original label (i.e. intergenic regions) in further analysis
#'Perform data annotation
#'@description Read data and perform gene ID mapping using built in databases
#'@param dataSetObj Input the name of the created datasetObj (see Init.Data).
#'@param org three letters annotation of organism (i.e hsa, mmu)
#'@param dataType Either "array" (microarray) or "count" (rna-seq)
#'@param idType original id type of features
#'@param lvlOpt merging original ID to entrez gene IDs. "NA" to keep original IDs without merging
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: MIT
#'@export
#'
PerformDataAnnot <- function(dataName="", org="hsa", dataType="array", idType="entrez", lvlOpt="mean"){
  dataSet <- readDataset(dataName);
  dataSet <- PerformDataAnnotInternal(dataSet, dataName, org, dataType, idType, lvlOpt);
  return(RegisterData(dataSet));   

}
PerformDataAnnotInternal <- function(dataSet, dataName=NULL, org="hsa", dataType="array", idType="entrez", lvlOpt="mean"){

  paramSet <- readSet(paramSet, "paramSet");
  msgSet <- readSet(msgSet, "msgSet");
  current.msg <- "";
  
  if(org == "custom"){
    idType <- "custom";
  }
  paramSet$data.org <- org;
  paramSet$data.idType <- idType;
  
  dataSet$type <- dataType;
  dataSet$id.orig <- dataSet$id.current <- idType;
  dataSet$annotated <- F;
  # Use legacy-friendly matrix produced at ingestion
  if(dataType=="prot"){
    data.proc <- ov_qs_read("int.mat.qs");
  }else{
    data.proc <- ov_qs_read("data.raw.qs");
  }
  #msg("[Annot] Loaded data matrix: ", nrow(data.proc), " features x ", ncol(data.proc), " samples")
  data.anot <- data.proc;
  # print(head(data.anot));
  # Save intensity matrix for downstream steps (ensures presence even after MSstats)
  ov_qs_save(data.proc, "int.mat.qs");
  
  if (org != 'NA' & idType != 'NA'){
    feature.vec <- rownames(data.proc);
    #msg("[Annot] Attempting ID mapping for ", length(feature.vec), " features; head IDs: ", paste(utils::head(feature.vec, 5), collapse=", "))
    
    # PROTEOMICS: Convert TO UniProt as central ID (not Entrez)
    anot.id <- .doProteomicsAnnotation(feature.vec, idType, paramSet);
    # Preserve names so downstream lookup (e.g., comp.res export) can align by original IDs
    if (is.null(names(anot.id))) {
      names(anot.id) <- feature.vec
    }

    # Get symbol mapping for UniProt IDs
    if(idType %in% c("s2f", "generic", "ko")){
      symbol.map <- .doGeneIDMapping(anot.id, idType, paramSet, "matrix");
    }else{
      # anot.id now contains UniProt IDs; get their symbols
      symbol.map <- .getUniprotSymbolMap(anot.id, paramSet);
    }
    # Filter to only mapped IDs
    if ("uniprot" %in% colnames(symbol.map)) {
      symbol.map <- symbol.map[which(symbol.map$uniprot %in% anot.id),];
    } else {
      symbol.map <- symbol.map[which(symbol.map$gene_id %in% anot.id),];
    }

    # PROTEOMICS: symbol.map already has uniprot column from .getUniprotSymbolMap()
    # Only add if missing (for special ID types like s2f, ko)
    if (!("uniprot" %in% colnames(symbol.map))) {
      #cat(sprintf("[PerformDataAnnot] DEBUG: Adding uniprot column for special ID type\n"))
      symbol.map <- .addUniprotColumn(symbol.map, idType, feature.vec, anot.id, paramSet);
    }

    #cat(sprintf("[PerformDataAnnot] DEBUG: symbol.map columns: %s\n",
    #            paste(colnames(symbol.map), collapse=", ")))
    #cat(sprintf("[PerformDataAnnot] DEBUG: symbol.map rows: %d\n", nrow(symbol.map)))
    if ("uniprot" %in% colnames(symbol.map)) {
      #cat(sprintf("[PerformDataAnnot] DEBUG: Sample uniprot IDs: %s\n",
      #            paste(head(symbol.map$uniprot, 5), collapse=", ")))
      #cat(sprintf("[PerformDataAnnot] DEBUG: Non-NA uniprot count: %d/%d\n",
      #            sum(!is.na(symbol.map$uniprot)), nrow(symbol.map)))
    }

    saveDataQs(symbol.map, "symbol.map.qs", paramSet$anal.type, dataName);
    #cat(sprintf("[PerformDataAnnot] DEBUG: Saved symbol.map.qs to: %s\n", getwd()))
    
    ov_qs_save(anot.id, "annotation.qs");
    
    hit.inx <- !is.na(anot.id);
    #msg("[Annot] Mapping hits: ", sum(hit.inx), "/", length(anot.id), " (", round(sum(hit.inx)/length(anot.id)*100, 2), "%); first mapped IDs: ", paste(utils::head(anot.id[hit.inx], 5), collapse=", "))
    matched.len <- sum(hit.inx);
    perct <- round(matched.len/length(feature.vec),3)*100;
    thresh <- 0.1 # previous value of 0.25 is causing challenges 
    #for datasets like Ppromelas with low annotation quality
    if (matched.len < length(feature.vec)*thresh){
      current.msg <- paste('Only ', perct, '% ID were matched. You may want to choose another ID type or use default.', sep=""); 
    } else {
      current.msg <- paste("ID annotation: ", "Total [", length(anot.id),
                           "] Matched [", matched.len, "] Unmatched [", sum(!hit.inx),"]", collapse="\n");

      # PROTEOMICS: All data types use UniProt as central ID
      # Phosphosites keep their site-specific IDs, others summarize to protein level
      is_phospho <- (!is.null(paramSet$data.type) && paramSet$data.type == "phospho")

      if (lvlOpt != 'NA' && !is_phospho){
        # Summarize to UniProt protein level (not Entrez gene level)
        matched.uniprot <- anot.id[hit.inx];
        data.anot <- data.proc[hit.inx,];
        rownames(data.anot) <- matched.uniprot;
        current.msg <- paste(current.msg, "Data is now mapped to UniProt ID.");
        paramSet$lvl.opt <- lvlOpt;
        res <- RemoveDuplicates(data.anot, lvlOpt, quiet=F, paramSet, msgSet);
        data.anot <- res[[1]];
        msgSet <- res[[2]];
        dataSet$id.current <- "uniprot";
        dataSet$annotated <- T;
        #msg("[Annot] Summarized to ", nrow(data.anot), " unique UniProt IDs after duplicate handling (lvlOpt=", lvlOpt, ")")
      } else if (is_phospho) {
        # Keep phosphosite IDs as-is, do NOT summarize to protein level
        dataSet$annotated <- FALSE
        #msg("[Annot] PHOSPHO DATA: Keeping phosphosite IDs (not summarizing to protein)")
        current.msg <- paste(current.msg, "Phosphosite IDs preserved (not mapped to protein level).")
      } else {
        #current.msg <- paste(current.msg, "No protein level summarization was performed.");
      }
    }
  } else { # no conversion will be performed
    feature.vec <- rownames(data.proc);
    anot.id <- feature.vec
    perct <- 100;
    hit.inx <- !is.na(anot.id);
    matched.len <- length(feature.vec); # dummies
    minLvl <- 1;
    current.msg <- paste("No annotation was performed."); 
    #msg("[Annot] Skipping annotation (org/idType is NA) - keeping ", length(feature.vec), " original IDs")
  }
  # need to save the ids (mixed gene annotation and original id) 
  # in case, users needs to keep unannotated features
  # this need to be updated to gether with data from now on
  dataSet$data.norm <- data.anot;
  
  ov_qs_save(data.anot, file="orig.data.anot.qs"); # immutable baseline after annotation
  ov_qs_save(data.anot, file="norm.input.anot.qs"); # current normalization input, reset on re-annotation
  col.sum <- colSums(dataSet$data.norm);
  totalCount <-  sum(col.sum);
  avgCount <- totalCount / ncol(dataSet$data.norm);
  minCount <- min(col.sum);
  maxCount <- max(col.sum);

  # Build metadata column names string (not factor levels)
  grp_names = "";
  if(!is.null(dataSet$meta.info) && ncol(dataSet$meta.info) > 0){
    grp_names = paste(names(dataSet$meta.info), collapse = "; ");
  }

  # Build metadata factors info (for reference, if needed)
  lvls = "";
  if(any(dataSet$disc.inx.orig)){
    disc = paste(names(dataSet$disc.inx.orig)[which(dataSet$disc.inx.orig)],collapse = ", ")
    lvls = paste0(lvls,length(which(dataSet$disc.inx.orig))," discrete factors: ",disc,"; ")
  }
  if(any(dataSet$cont.inx.orig)){
    cont = paste(names(dataSet$cont.inx.orig)[which(dataSet$cont.inx.orig)],collapse = ", ")
    lvls = paste0(lvls,length(which(dataSet$cont.inx.orig))," continuous factors: ",cont,".")
  }

  # Use the original missing value count stored BEFORE imputation
  # (stored in paramSet by ReadTabExpressData before replacing NAs with minVal/2)
  if(!is.null(paramSet$original.missing.count)){
    missNum_count <- paramSet$original.missing.count;
  } else {
    # Fallback for older data (shouldn't happen with new code)
    missNum <- which(is.na(dataSet$data.norm)|dataSet$data.norm=="NA"|dataSet$data.norm=="");
    missNum_count <- length(missNum);
  }

  msgSet$current.msg <- current.msg;
  msgSet$summaryVec <- c(matched.len, perct, length(anot.id), sum(!hit.inx), ncol(dataSet$data.norm), ncol(dataSet$meta.info), sprintf("%4.2e", signif(totalCount ,3)), sprintf("%4.2e",signif(avgCount, 3)), sprintf("%4.2e",signif(minCount, 3)), sprintf("%4.2e",signif(maxCount,3)), grp_names, missNum_count)  
  #if(length(missNum)>0){
  #  RemoveMissingPercent(dataSet$name, 0.5)
  ##  ImputeMissingVar(dataSet$name, method="min")
  #}else{
    ov_qs_save(data.anot, file="data.missed.qs");
  #}
  data.anot <- sanitizeSmallNumbers(data.anot);
  fast.write(sanitizeSmallNumbers(data.anot), file="data_annotated.csv");
  .save.annotated.data(data.anot);
  saveSet(paramSet, "paramSet");
  saveSet(msgSet, "msgSet");
  return(dataSet);
}


# Annotating features to internal database
AnnotateGeneData <- function(dataName, org, lvlOpt, idtype){
  paramSet <- readSet(paramSet, "paramSet");
  msgSet <- readSet(msgSet, "msgSet");
  dataSet <- readDataset(dataName);
  
  if(org == "NA"){
    msgSet$current.msg <- "Invalid organism!"
    saveSet(msgSet, "msgSet");
    return(1)
  }
  
  data.raw <- readDataQs("data.raw.qs", paramSet$anal.type, dataName);
  gene.vec <- rownames(data.raw);
  
  #record the info
  paramSet$data.org <- org
  dataSet$q.type.gene <- idtype;
  
  dataSet$gene.org <- org;
  dataSet$gene <- gene.vec;
  if(idtype == "NA"){
    enIDs <- gene.vec;
    symbol.map <- data.frame(gene_id=enIDs, symbol=enIDs, name=enIDs);
  }else{
    if(idtype %in% c("entrez", "symbol", "refseq", "gb", "embl_gene","embl_protein", "embl_transcript", "orf", "tair", "wormbase", "ko", "custom", "s2f")){
      enIDs <- .doGeneIDMapping(gene.vec, idtype, paramSet, "vec");
    }else{
      enIDs <- .doProbeMapping(gene.vec, idtype, paramSet);
      names(enIDs) <- gene.vec;
    }
  
  
  tblNm <- getEntrezTableName(org, "entrez");
  symbol.map <- queryGeneDB(tblNm, org);
  symbol.map <- symbol.map[which(symbol.map$gene_id %in% enIDs),];
  }

  # Add uniprot column
  symbol.map <- .addUniprotColumn(symbol.map, idtype, gene.vec, enIDs, paramSet);

  saveDataQs(symbol.map, "symbol.map.qs", paramSet$anal.type, dataName);

  
  if(idtype == "kos"){
    kos <- enIDs$kos;
    enIDs <- enIDs$entrezs;
    dataSet$kos.name.map <- kos
  }
  
  # Handle case when only KOs are mapped with no corresponding entrez id
  na.inx <- is.na(enIDs);
  
  if(sum(!na.inx) == 0 && idtype == "kos"){
    na.inx <- is.na(kos);
  }
  
  dataSet$gene.name.map <- list(
    hit.values=enIDs,
    match.state = ifelse(is.na(enIDs), 0, 1)
  );
  
  hit.inx <- which(!is.na(enIDs));
  matched.len <- length(hit.inx);
  if(matched.len > 1){
    data.norm <- data.raw[hit.inx,];
    matched.entrez <- enIDs[hit.inx];
    
    
    # now, deal with duplicated entrez id
    # first, average duplicate rows
    dataSet$lvl.opt <- lvlOpt;
    res <- RemoveDuplicates(data.norm, lvlOpt, quiet=F, paramSet, msgSet);
    int.mat <- res[[1]];
    msgSet <- res[[2]];
    # update
    
    data.annotated <-int.mat;
    rownames(int.mat) <- matched.entrez
    if(idtype %in% c("mir_id", "mir_acc", "mirnet")){
      rownames(data.annotated) <- rownames(int.mat);
    }else{
      rownames(data.annotated) <- matched.entrez
    }
    dataSet$enrich_ids = rownames(int.mat);
    names(dataSet$enrich_ids) = doEntrez2SymbolMapping(rownames(int.mat), paramSet$data.org, paramSet$data.idType)
    
    dataSet$id.type <- "entrez";
    
  }else{
    data.annotated <- data.raw;
    dataSet$enrich_ids = rownames(data.annotated)
    dataSet$id.type <- "none";
  }

  if(idtype != "NA"){
    # OPTIMIZED: Calculate mapping ratio once instead of twice
    mapping_ratio <- length(unique(enIDs))/length(gene.vec)
    if(mapping_ratio < 0.3){
      msg <- paste("Less than ", round(mapping_ratio * 100, 2), "% features were mapped");
      msgSet$current.msg <- msg;
      saveSet(msgSet, "msgSet");
      return(0)
    }else{
      msg <- paste("A total of ", length(unique(enIDs)), "unique features were mapped");
    }
  }else{
    msg <- paste("There is a total of ", length(unique(gene.vec)), "unique features.");
  }
  
  saveDataQs(data.annotated, "data.annotated.qs", paramSet$anal.type, dataName);
  msgSet$current.msg <- msg;
  saveSet(msgSet, "msgSet");
  saveSet(paramSet, "paramSet");
  return(RegisterData(dataSet));
}

# PROTEOMICS: Convert a vector of IDs to UniProt IDs (central ID for proteomics)
.doProteomicsAnnotation <- function(feature.vec, idType, paramSet){
  org <- paramSet$data.org;

  # If already UniProt, keep as-is
  if (idType == "uniprot") {
    anot.id <- feature.vec
    names(anot.id) <- feature.vec
    return(anot.id)
  }

  # For other ID types, convert to Entrez first, then Entrez to UniProt
  if(idType %in% c("entrez", "symbol", "refseq", "gb", "embl_gene","embl_protein", "embl_transcript", "orf", "tair", "wormbase", "ko", "custom", "cds", "s2f", "string", "transcript")){
    # Step 1: Convert to Entrez
    entrez.ids <- .doGeneIDMapping(feature.vec, idType, paramSet, "vec");

    # Step 2: Convert Entrez to UniProt
    anot.id <- .convertEntrezToUniprot(entrez.ids, org);

    # Preserve original names
    names(anot.id) <- feature.vec
  }else{
    # Probe IDs: convert to Entrez first
    entrez.ids <- .doProbeMapping(feature.vec, idType, paramSet);

    # Then Entrez to UniProt
    anot.id <- .convertEntrezToUniprot(entrez.ids, org);
    names(anot.id) <- feature.vec;
  }
  return(anot.id);
}

# Helper: Convert Entrez IDs to UniProt IDs
.convertEntrezToUniprot <- function(entrez.vec, org){
  # Query the entrez_uniprot mapping table
  db.map <- queryGeneDB("entrez_uniprot", org);

  if (is.null(db.map) || !is.data.frame(db.map) || nrow(db.map) == 0) {
    warning("[Annot] UniProt mapping table not available for ", org, "; keeping Entrez IDs")
    return(entrez.vec)
  }

  # Match Entrez IDs to get UniProt accessions
  hit.inx <- match(entrez.vec, db.map[, "gene_id"]);
  uniprot.ids <- db.map[hit.inx, "accession"];

  # For unmapped Entrez IDs, keep them as-is
  na.inx <- is.na(uniprot.ids);
  uniprot.ids[na.inx] <- entrez.vec[na.inx];

  return(uniprot.ids)
}

# Helper: Get symbol mapping for UniProt IDs
.getUniprotSymbolMap <- function(uniprot.vec, paramSet){
  org <- paramSet$data.org;

  # Step 1: UniProt → Entrez (via entrez_uniprot table)
  uniprot.map <- queryGeneDB("entrez_uniprot", org);

  if (is.null(uniprot.map) || !is.data.frame(uniprot.map) || nrow(uniprot.map) == 0) {
    # Return minimal map if table not available
    return(data.frame(
      uniprot = uniprot.vec,
      gene_id = uniprot.vec,
      symbol = uniprot.vec,
      name = rep("NA", length(uniprot.vec)),
      stringsAsFactors = FALSE
    ))
  }

  # Match UniProt to get Entrez IDs
  hit.inx <- match(uniprot.vec, uniprot.map[, "accession"]);
  entrez.ids <- uniprot.map[hit.inx, "gene_id"];

  # Step 2: Entrez → Symbol/Name (via entrez table)
  entrez.db <- queryGeneDB("entrez", org);
  entrez.hit <- match(entrez.ids, entrez.db[, "gene_id"]);

  # Build result dataframe
  result <- data.frame(
    uniprot = uniprot.vec,
    gene_id = entrez.ids,
    symbol = entrez.db[entrez.hit, "symbol"],
    name = entrez.db[entrez.hit, "name"],
    stringsAsFactors = FALSE
  )

  # Fill NAs with UniProt ID
  na.inx <- is.na(result$symbol);
  result$symbol[na.inx] <- uniprot.vec[na.inx];
  na.inx <- is.na(result$gene_id);
  result$gene_id[na.inx] <- uniprot.vec[na.inx];
  na.inx <- is.na(result$name);
  result$name[na.inx] <- "NA";

  return(result)
}

#Convert a vector of ids to vector of entrez ids (LEGACY - kept for backward compatibility)
.doAnnotation <- function(feature.vec, idType, paramSet){
  if(idType %in% c("entrez", "symbol", "refseq", "gb", "embl_gene","embl_protein","uniprot", "embl_transcript", "orf", "tair", "wormbase", "ko", "custom", "cds", "s2f", "string", "transcript")){
    anot.id <- .doGeneIDMapping(feature.vec, idType, paramSet, "vec");
  }else{
    anot.id <- .doProbeMapping(feature.vec, idType, paramSet);
    names(anot.id) <- feature.vec;
  }
  return(anot.id);
}

.doGeneIDMapping <- function(feature.vec, idType, paramSet, outputType = "vec", keepNA = F) {
  org <- paramSet$data.org;
  
  if (is.null(feature.vec)) {
    db.map <- queryGeneDB("entrez", org);
    feature.vec <- db.map[, "gene_id"];
    idType <- "entrez";
    # print(".doGeneIDMapping, empty feature.vec, get whole table");
  }
  
  col.nm <- "";
  db.nm <- "";
  
  if (org == "zhangshugang" || org == "cro" || org == "dimmitis"|| org == "hpolygyrus") {
      q.mat <- do.call(rbind, strsplit(feature.vec, "\\."));
      feature.vec <- q.mat[, 1];
  }

  if (idType %in% c("s2f", "ko") || paramSet$data.idType %in% c("s2f", "ko")) {
    col.nm <- "gene_id";
    db.nm <- paste0("entrez_", idType);
  } else if (idType == "symbol") {
    col.nm <- "symbol";
    db.nm <- "entrez";
  } else if (idType == "entrez") {
    col.nm <- "gene_id";
    db.nm <- "entrez";
  } else if (idType == "custom") {
    db.nm <- "custom";
    col.nm <- "gene_id";
  } else if (idType == "cds") {
    col.nm <- "accession";
    db.nm <- "entrez_cds";
  } else if (idType == "uniprot") {
    # For UniProt IDs, keep them as-is (no conversion to Entrez)
    # This ensures compatibility with UniProt-based PPI databases (Rolland, HuRI, IntAct, IRefIndex)
    col.nm <- "accession";
    db.nm <- "uniprot"; # Use UniProt validation table (if available) or keep as-is
    # No strsplit on feature.vec as UniProt primary accessions generally don't have dots.
    # If isoforms with dots are present, they should be handled by the mapping table itself.
  } else {
    if (!(idType == "refseq" && org == "fcd") && !(idType == "string" && org == "cel")) {
      q.mat <- do.call(rbind, strsplit(feature.vec, "\\."));
      feature.vec <- q.mat[, 1];
    }
    col.nm <- "accession";
    if (idType == "tair") {
      db.nm <- "tair";
    } else {
      db.nm <- paste0("entrez_", idType);
    }
  }
  
  # Special handling for UniProt IDs - keep them as-is without conversion
  if (idType == "uniprot") {
    # Normalize UniProt IDs (remove isoforms, phosphosites, etc.)
    normalized_ids <- sub("-\\d+$", "", feature.vec)  # Remove isoform suffixes
    normalized_ids <- sub("_[A-Z]_\\d+$", "", normalized_ids)  # Remove phosphosite annotations
    normalized_ids <- trimws(normalized_ids)

    if (outputType == "vec") {
      # Return normalized UniProt IDs
      return(normalized_ids)
    } else {
      # Return table format with UniProt IDs in both accession and gene_id columns
      return(data.frame(
        accession = normalized_ids,
        gene_id   = normalized_ids,  # Keep as UniProt (not Entrez)
        orig      = feature.vec,
        unmapped  = FALSE,  # All valid UniProt IDs are considered "mapped"
        stringsAsFactors = FALSE
      ))
    }
  }

  db.map <- queryGeneDB(db.nm, org);
  # Gracefully handle missing or malformed mapping tables (avoid match() dimension errors)
  if (is.null(db.map) || !is.data.frame(db.map) || nrow(db.map) == 0) {
    warning("[Annot] Mapping table not available for ", db.nm, " (", org, "); returning NA mappings.")
    if (outputType == "vec") {
      na.vec <- rep(NA_character_, length(feature.vec))
      names(na.vec) <- feature.vec
      return(na.vec)
    } else {
      return(data.frame(
        accession = feature.vec,
        gene_id   = NA_character_,
        orig      = feature.vec,
        unmapped  = TRUE,
        stringsAsFactors = FALSE
      ))
    }
  }
  if (org == "smm" && idType == "symbol") {
    q.mat <- do.call(rbind, strsplit(feature.vec, "\\."));
    feature.vec <- q.mat[, 1];
    q.mat <- do.call(rbind, strsplit(db.map[, col.nm], "\\."));
    db.map[, col.nm] <- q.mat[, 1];
  }
  
  if (!col.nm %in% colnames(db.map)) {
    warning("[Annot] Column ", col.nm, " missing from mapping table ", db.nm, "; returning NA mappings.")
    if (outputType == "vec") {
      na.vec <- rep(NA_character_, length(feature.vec))
      names(na.vec) <- feature.vec
      return(na.vec)
    } else {
      return(data.frame(
        accession = feature.vec,
        gene_id   = NA_character_,
        orig      = feature.vec,
        unmapped  = TRUE,
        stringsAsFactors = FALSE
      ))
    }
  }

  hit.inx <- match(feature.vec, db.map[, col.nm]);
  
  if(outputType == "vec"){
    entrezs <- db.map[hit.inx, "gene_id"];
    
    mode(entrezs) <- "character";
    rm(db.map, feature.vec); gc();
    return(entrezs);
  } else {
    entrezs <- db.map[hit.inx, ];
    
    entrezs <- entrezs[,c(2,1)];
    unmapped_flag <- is.na(entrezs[, 1]) | is.na(entrezs[, 2]);  # Flag unmapped entries
    
    if (!keepNA) {
      na.inx <- is.na(entrezs[, 1]);
      entrezs[, 1][na.inx] <- feature.vec[na.inx];
      na.inx <- is.na(entrezs[, 2]);
      entrezs[, 2][na.inx] <- feature.vec[na.inx];
    } else {
      unmapped.inx <- is.na(entrezs$gene_id);
      if (any(unmapped.inx)) {
        entrezs$accession[unmapped.inx] <- feature.vec[unmapped.inx];
      }
    }
    
    colnames(entrezs) <- c("accession", "gene_id");
    entrezs <- cbind( entrezs,orig = feature.vec, unmapped = unmapped_flag);  # Add orig and unmapped columns
    return(entrezs);
  }
}


# from probe ID to entrez ID 
.doProbeMapping <- function(probe.vec, platform, paramSet){
  org <- paramSet$data.org;

  if(exists("api.lib.path")){
    lib.path <- api.lib.path;
  }else{
    lib.path <- paramSet$lib.path;
  }
  platform.path <- paste(lib.path, org, "/", platform, ".rds", sep="");

  if(!paramSet$on.public.web && !file.exists(platform.path)){
    nmdb <- basename(platform.path);
    download.file(platform.path, destfile = nmdb, method="libcurl", mode = "wb");
    platform.path <- nmdb;
  }
  

  probe.map <- readRDS(platform.path);
  if(is.null(probe.vec)){
    entrez <- probe.map[, "entrez"];
  }else{
    hit.inx <- match(probe.vec, probe.map[, "probe"]);
    entrez <- probe.map[hit.inx, "entrez"];
  }
  return(entrez);
}


# Helper function to add uniprot column to symbol.map
# If original idType is "uniprot", preserves the original uniprot IDs
.addUniprotColumn <- function(symbol.map, idType, feature.vec, mapped.ids, paramSet) {
  # cat(sprintf("[.addUniprotColumn] DEBUG: Called with idType=%s\n", idType))
  # Check if symbol.map has gene_id column
  if (is.null(symbol.map) || !is.data.frame(symbol.map) || !"gene_id" %in% colnames(symbol.map)) {
    # cat(sprintf("[.addUniprotColumn] DEBUG: Invalid symbol.map structure, returning as-is\n"))
    # Return as-is if invalid structure
    return(symbol.map);
  }

  org <- paramSet$data.org;
  # cat(sprintf("[.addUniprotColumn] DEBUG: Organism=%s\n", org))

  # If original ID type is uniprot, preserve the original uniprot IDs
  if (!is.null(idType) && idType == "uniprot") {
    # cat(sprintf("[.addUniprotColumn] DEBUG: Input type is uniprot, preserving original IDs\n"))
    # Check if symbol.map already has accession column (from .doGeneIDMapping with matrix output)
    if ("accession" %in% colnames(symbol.map)) {
      # cat(sprintf("[.addUniprotColumn] DEBUG: Using accession column for uniprot\n"))
      # Use the accession column which contains the original uniprot IDs
      symbol.map$uniprot <- symbol.map$accession;
    } else {
      # cat(sprintf("[.addUniprotColumn] DEBUG: Creating uniprot mapping from named vector\n"))
      # Create mapping from entrez IDs back to original uniprot IDs
      # mapped.ids contains the entrez IDs (named vector with original IDs as names)
      if (!is.null(names(mapped.ids))) {
        # Create a lookup: entrez -> original uniprot
        entrez_to_uniprot <- data.frame(
          gene_id = as.character(mapped.ids),
          uniprot = names(mapped.ids),
          stringsAsFactors = FALSE
        );
        # Remove duplicates, keeping first occurrence
        entrez_to_uniprot <- entrez_to_uniprot[!duplicated(entrez_to_uniprot$gene_id), ];

        # Match symbol.map gene_ids to get uniprot
        hit.inx <- match(as.character(symbol.map$gene_id), entrez_to_uniprot$gene_id);
        symbol.map$uniprot <- entrez_to_uniprot$uniprot[hit.inx];
      } else {
        # Fallback: query database
        symbol.map <- .queryUniprotFromDB(symbol.map, org);
      }
    }
    return(symbol.map);
  }

  # For other ID types, query the uniprot mapping database
  # cat(sprintf("[.addUniprotColumn] DEBUG: Querying database for uniprot mapping\n"))
  # Skip if organism is NA, custom, or uniprot mapping not applicable
  if (is.null(org) || org %in% c("NA", "na", "custom")) {
   # cat(sprintf("[.addUniprotColumn] DEBUG: Organism is NA/custom, setting uniprot to NA\n"))
    symbol.map$uniprot <- NA_character_;
    return(symbol.map);
  }

  symbol.map <- .queryUniprotFromDB(symbol.map, org);
  # cat(sprintf("[.addUniprotColumn] DEBUG: Finished querying database\n"))
  return(symbol.map);
}

# Helper to query uniprot from database
.queryUniprotFromDB <- function(symbol.map, org) {
  # cat(sprintf("[.queryUniprotFromDB] DEBUG: Querying entrez_uniprot table for org=%s\n", org))
  # Query uniprot mapping table
  uniprot.map <- tryCatch({
    queryGeneDB("entrez_uniprot", org);
  }, error = function(e) {
    cat(sprintf("[.queryUniprotFromDB] DEBUG: Error querying database: %s\n", e$message))
    NULL;
  });

  # If uniprot table exists and has required columns, do the mapping
  if (!is.null(uniprot.map) && is.data.frame(uniprot.map) &&
      "gene_id" %in% colnames(uniprot.map) && "accession" %in% colnames(uniprot.map)) {
    # cat(sprintf("[.queryUniprotFromDB] DEBUG: Found uniprot mapping table with %d rows\n", nrow(uniprot.map)))

    # Flatten vectors to avoid dimension errors
    gene_ids <- as.vector(symbol.map$gene_id);
    dim(gene_ids) <- NULL;
    uniprot_gene_ids <- as.vector(uniprot.map$gene_id);
    dim(uniprot_gene_ids) <- NULL;

    # Match gene_id to get uniprot accessions
    hit.inx <- match(gene_ids, uniprot_gene_ids);
    symbol.map$uniprot <- uniprot.map$accession[hit.inx];
    # cat(sprintf("[.queryUniprotFromDB] DEBUG: Mapped %d/%d genes to uniprot\n",
    #            sum(!is.na(symbol.map$uniprot)), nrow(symbol.map)))
  } else {
    # cat(sprintf("[.queryUniprotFromDB] DEBUG: Uniprot mapping table not found or invalid\n"))
    # If table doesn't exist or doesn't have required columns, add NA column
    symbol.map$uniprot <- NA_character_;
  }

  return(symbol.map);
}

queryGeneDB <- function(db.nm, org){
  if(org == "na"){
      print("Not available when organism is not specified");
      return(0);
  }
  paramSet <- readSet(paramSet, "paramSet");    
  if(db.nm == "custom" || org == "custom"){
    db.map <- ov_qs_read("anot_table.qs");
  }else{
    require('RSQLite');
    
    db.path <- paste(paramSet$sqlite.path, org, "_genes.sqlite", sep="")
    if(!PrepareSqliteDB(db.path, paramSet$on.public.web)){
      stop("Sqlite database is missing, please check your internet connection!");
    }
    conv.db <- dbConnect(SQLite(), db.path); 
    tbls <- dbListTables(conv.db)
    if(!db.nm %in% tbls){
        return(0);
    }
    db.map <- dbReadTable(conv.db, db.nm);
    dbDisconnect(conv.db); cleanMem();
  }
  
  return(db.map)
}

getEntrezTableName <- function(data.org, data.idType){
    if(data.org == "generic"){
        tblNm <- paste0("entrez_", data.idType);
    }else{
        tblNm <- "entrez";
    }
    return(tblNm);
}


doEntrez2SymbolMapping <- function(entrez.vec,data.org="NA", data.idType="NA"){

  if(data.idType == "symbol"){
    return(entrez.vec);
  }

  if(data.org == "NA"){
    doEntrez2SymbolMapping(entrez.vec);
  }

  # Map Entrez IDs to symbols via entrez database
  entrez.db <- queryGeneDB("entrez", data.org);
  entrez.db[] <- lapply(entrez.db, as.character)

  hit.inx <- match(entrez.vec, entrez.db[, "gene_id"]);
  symbols <- entrez.db[hit.inx, "symbol"];

  # If no symbol found, use Entrez ID itself
  na.inx <- is.na(symbols);
  symbols[na.inx] <- entrez.vec[na.inx];

  return(symbols);
}

# PROTEOMICS: UniProt-centric ID mapping
doUniprot2SymbolMapping <- function(uniprot.vec, data.org="NA", data.idType="NA"){

  if(data.idType == "symbol"){
    return(uniprot.vec); # nothing to do
  }

  if(data.org == "NA" && data.idType=="NA"){
    paramSet <- readSet(paramSet, "paramSet");
    data.org <- paramSet$data.org;
    data.idType <- paramSet$data.idType;

    # Check if phospho data - strip site info before mapping
    is_phospho <- (!is.null(paramSet$data.type) && paramSet$data.type == "phospho")
    if (is_phospho) {
      # For phospho data, strip site information to get base UniProt ID
      original.ids <- uniprot.vec
      uniprot.vec <- sapply(strsplit(as.character(uniprot.vec), "_"), function(x) x[1])
    }
  }

  if(data.org == "na"){
    return(uniprot.vec); # nothing to do
  }

  # Step 1: UniProt → Entrez (via entrez_uniprot table)
  uniprot.map <- queryGeneDB("entrez_uniprot", data.org);
  uniprot.map[] <- lapply(uniprot.map, as.character)

  hit.inx <- match(uniprot.vec, uniprot.map[, "accession"]);
  entrez.ids <- uniprot.map[hit.inx, "gene_id"];

  # Step 2: Entrez → Symbol (via entrez table)
  entrez.db <- queryGeneDB("entrez", data.org);
  entrez.db[] <- lapply(entrez.db, as.character)

  entrez.hit <- match(entrez.ids, entrez.db[, "gene_id"]);
  symbols <- entrez.db[entrez.hit, "symbol"];

  # If no symbol found, use UniProt ID itself
  na.inx <- is.na(symbols);
  symbols[na.inx] <- uniprot.vec[na.inx];
  return(symbols);
}
# note, entrez.vec could contain NA / NULL – cannot rely on rownames
doEntrezIDAnot <- function(entrez.vec,
                           data.org   = "NA",
                           data.idType = "NA") {
  # NOTE: Despite the name, this now handles UniProt IDs (central ID for proteomics)
  # Kept function name for backward compatibility
  return(doUniprotIDAnot(entrez.vec, data.org, data.idType));
}

# PROTEOMICS: UniProt-centric ID annotation
doUniprotIDAnot <- function(uniprot.vec,
                            data.org   = "NA",
                            data.idType = "NA") {

  # Store original IDs for phospho data
  original.ids <- uniprot.vec

  if (data.org == "NA" && data.idType == "NA") {
    paramSet  <- readSet(paramSet, "paramSet")
    data.org  <- paramSet$data.org
    data.idType <- paramSet$data.idType

    # Check if phospho data - strip site info before annotation
    is_phospho <- (!is.null(paramSet$data.type) && paramSet$data.type == "phospho")
    if (is_phospho) {
      # For phospho data, strip site information to get base UniProt ID
      uniprot.vec <- sapply(strsplit(as.character(uniprot.vec), "_"), function(x) x[1])
    }
  }

  if (tolower(data.org) == "na") {
    # No organism context – fall back to a dummy annotation frame
    return(data.frame(
      gene_id = original.ids,
      symbol  = original.ids,
      name    = rep("NA", length(original.ids)),
      stringsAsFactors = FALSE
    ))
  }

  ## ── 2 · database lookup ---------------------------------------------------
  # Flatten vectors to avoid dimension errors in match()
  uniprot.vec <- as.vector(uniprot.vec)
  dim(uniprot.vec) <- NULL

  # Step 1: UniProt → Entrez (via entrez_uniprot table)
  uniprot.map <- queryGeneDB("entrez_uniprot", data.org)
  uniprot.map[] <- lapply(uniprot.map, as.character)

  accession_col <- as.vector(uniprot.map[, "accession"])
  dim(accession_col) <- NULL

  hit.inx   <- match(uniprot.vec, accession_col)
  entrez.ids <- uniprot.map[hit.inx, "gene_id"]

  # Step 2: Entrez → Symbol/Name (via entrez table)
  entrez.db <- queryGeneDB("entrez", data.org)
  entrez.db[] <- lapply(entrez.db, as.character)

  entrez.hit <- match(entrez.ids, entrez.db[, "gene_id"])

  # Build annotation matrix
  anot.mat <- data.frame(
    gene_id = uniprot.vec,  # Use UniProt as gene_id (central ID)
    symbol = entrez.db[entrez.hit, "symbol"],
    name = entrez.db[entrez.hit, "name"],
    stringsAsFactors = FALSE
  )

  # Fill NAs with UniProt ID
  na.inx <- is.na(anot.mat$symbol)
  anot.mat$symbol[na.inx] <- uniprot.vec[na.inx]
  na.inx <- is.na(anot.mat$name)
  anot.mat$name[na.inx] <- "NA"

  ## ── 3 · fill NAs with original IDs ---------------------------------------
  na.inx <- is.na(hit.inx)
  anot.mat[na.inx, "gene_id"] <- original.ids[na.inx]
  anot.mat[na.inx, "symbol"]  <- original.ids[na.inx]
  anot.mat[na.inx, "name"]    <- "NA"

  # Use original IDs (including site info) as gene_id for phospho data
  # This ensures visual analytics can highlight the correct phosphosite nodes
  anot.mat$gene_id <- original.ids

  rownames(anot.mat) <- NULL
  return(anot.mat)
}


doIdMappingGeneric <- function(orig.vec, gene.map, colNm1, colNm2, type="vec"){

  # Flatten vectors to avoid dimension errors in match()
  orig.vec <- as.vector(orig.vec)
  dim(orig.vec) <- NULL
  col_vec <- as.vector(gene.map[, colNm1])
  dim(col_vec) <- NULL

  hit.inx <- match(orig.vec, col_vec);
  if(colNm2 =="symbol"){
    if(!colNm2 %in% colnames(gene.map)){
        colnames(gene.map)[which(colnames(gene.map) == "accession")] <- "symbol";
    }
  }
  if(type == "vec"){
  result.vec <- gene.map[hit.inx, colNm2];
  
  # if not gene symbol, use id by itself
  na.inx <- is.na(result.vec);
  result.vec[na.inx] <- orig.vec[na.inx];
  return(result.vec);
  }else{
  na.inx <- is.na(hit.inx);
  anot.mat <- gene.map[hit.inx,];
  anot.mat[na.inx, "symbol"] <- orig.vec[na.inx];
  anot.mat[na.inx, "name"] <- orig.vec[na.inx];
  return(anot.mat);
  }
}

##########################################
############# private utility methods #### 
##########################################

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
