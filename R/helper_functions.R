##################################################
## R script for ProteoAnalyst
## Description: Functions related to web interface
## Author: Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
##################################################

GetSigGeneCount <- function(){
  analSet <- readSet(analSet, "analSet");
  return(analSet$sig.gene.count);
}

GetSigGeneCountTotal <- function(){
  analSet <- readSet(analSet, "analSet");
  return(analSet$sig.gene.count.total);
}


CheckRawDataAlreadyNormalized <- function(dataName=""){
  dataSet <- readDataset(dataName);
  #data <- dataSet$data.anot;
  data <- .get.annotated.data();
  if(sum(data > 100, na.rm = TRUE) > 100){ # now we think it is raw counts
    return(0);
  }else{
    return(1);
  }
}

GetMetaCol<- function(dataName=""){
  dataSet <- readDataset(dataName);
  paramSet <- readSet(paramSet, "paramSet");
  anal.type <- paramSet$anal.type;
  if(anal.type == "onedata"){
    colNms <- colnames(dataSet$comp.res);
    if (dataSet$de.method=="limma" || dataSet$de.method=="deqms"){
      inx <- match("AveExpr", colNms)
      return(names(dataSet$comp.res.list));
    } else if(dataSet$de.method=="wtt"){
      inx <- match("t", colNms)
      return(names(dataSet$comp.res.list));
  } else if (dataSet$de.method=="deseq2"){
      inx <- match("baseMean", colNms)
      return(names(dataSet$comp.res.list));
    } else {
      inx <- match("logCPM", colNms)
      return(names(dataSet$comp.res.list));
    }
    resT <- dataSet$comp.res;
    if(inx > 2){
      resT <- resT[,1:inx-1];
      nms <- gsub("logFC.", "logFC_", colnames(resT));
      
      # if there are decimals, we don't want to replace them with "vs"
      # find number of decimals in each column name
      num.dec <- lengths(regmatches(nms, gregexpr("\\.", nms)))
      
      # when only one ".", it's easy
      nms[num.dec == 1] <- gsub("\\.", " vs ", nms[num.dec == 1])
      
      # for three ".", replace only the middle
      nms[num.dec == 3] <- sapply(strsplit(nms[num.dec == 3], "\\."), function(x) {
        g <- seq_along(x)
        g[g < 2] <- 2
        g[g > 2 + 1] <- 2+1
        paste(tapply(x, g, paste, collapse = "."), collapse = " vs ")
      })
      nms <- unlist(nms)
      
      # for two ".", difficult to know which one - just leave as is
      
      return(as.vector(nms));
    }else{
      return(dataSet$par1);
    }
  }else{
    nms <- paste(unique(dataSet$cls), collapse=" vs ");
    return(nms);
  }
}

GetSummaryData <- function(){
  msgSet <- readSet(msgSet, "msgSet");
#print(msgSet$summaryVec);
  return(msgSet$summaryVec);
}

GetMetaColLength <- function(dataName = "") {

  dataSet <- readDataset(dataName)

  # bail out early if no comparison results
  if (is.null(dataSet$comp.res) && is.null(dataSet$comp.res.list)) {
    return(0L)
  }

  method <- tolower(dataSet$de.method)

  if (method == "limma" || method == "deqms") {
    return(length(dataSet$comp.res.list))
  } else if (method == "deseq2") {
    return(length(dataSet$comp.res.list))
  } else {           
    return(length(dataSet$comp.res.list))
  }

  # if the marker column isn't present → length 0
  #if (is.na(inx) || inx <= 1) {
  #  return(0L)
  #}

  #resT <- dataSet$comp.res[, seq_len(inx - 1), drop = FALSE]
  #length(colnames(resT))
}

GetInitLib <- function(){
  paramSet <- readSet(paramSet, "paramSet");
  init.lib <- paramSet$init.lib;
  return(init.lib)
}

GetMetaDatasets<- function(){
  paramSet <- readSet(paramSet, "paramSet");
  mdata.all <- paramSet$mdata.all;
  sel.nms <- names(mdata.all)[mdata.all==1];
  return(sel.nms);
}

SetSelMetaData<- function(selNm){
    paramSet <- readSet(paramSet, "paramSet");
    paramSet$selDataNm <- selNm;
    paramSet$jsonNms$dataName <- selNm;
    saveSet(paramSet, "paramSet");
}

# only for switching single expression data results
SetCurrentData <- function(nm){
  dataSet <- readDataset(nm);
  return(1);
}

GetOmicsDataDims <- function(dataName){
  dataSet <- readDataset(dataName);
  paramSet <- readSet(paramSet, "paramSet");
  if(paramSet$anal.type == "proteinlist"){
  dm <- c(nrow(dataSet$prot.mat), 0);
  naNum <- 0;
  }else{
  data.mat <- dataSet$data.norm;
  dm <- dim(data.mat);
  if (is.null(dm) || length(dm) < 2) {
    if (!is.null(data.mat) && length(data.mat) > 0) {
      dm <- c(length(data.mat), 1);
    } else {
      dm <- c(0, 0);
    }
  }
  naNum <- sum(is.na(data.mat));
  }

  return(c(dm, naNum));
} 

# given dataSet Name, sample name, and class name, do update
# note, for multiple #class, this set which one to use in the subsequent steps
# last one wins

# read in the data and perform
# gene ID mapping using built in libraries
# matchMin is minimal matched probe (%)
# return the total matched gene number

# obtain sample names and their class labels
GetSampleInfo <- function(dataName, clsLbl){
    dataSet <- readDataset(dataName);
    grpInfo <- dataSet$meta.info[[clsLbl]];
    grpLbls <- paste(levels(grpInfo), collapse="\n");
    smplInfo <- paste(Sample = colnames(dataSet$data.orig), "\t", Class=grpInfo, collapse="\n");
    return(c(grpLbls, smplInfo));
}

#for metadata
GetMetaSummaryData<- function(){
    paramSet <- readSet(paramSet, "paramSet");
    inmex.meta <- try(qs::qread("inmex_meta.qs"), silent = TRUE);
    if(inherits(inmex.meta, "try-error")){
      warning("Failed to read inmex_meta.qs");
      return(NULL);
    }
    sel.nms <- unique(inmex.meta$data.lbl)
    sel.nms <- paste(sel.nms, collapse="; ")
    cls.lbls <- unique(inmex.meta$cls.lbl)
    cls.lbls <- paste(cls.lbls, collapse="; ")
    paramSet$summaryVecMeta <- c(length(colnames(inmex.meta$data)),nrow(inmex.meta$data), sel.nms, cls.lbls);
    saveSet(paramSet, "paramSet");
    return(paramSet$summaryVecMeta)
}

GetDatasetNamesString <- function(){
    inmex.meta <- qs::qread("inmex_meta.qs");
    paste(unique(inmex.meta$data.lbl), collapse="||");
}

##Single matrix
GetSampleNumber <-function(){
  data.orig <- qs::qread("data.raw.qs");
  return(ncol(data.orig));
}


GetFilesToBeSaved <-function(naviString){
  paramSet <- readSet(paramSet, "paramSet");
  return(unique(paramSet$partialToBeSaved));
}

GetMetaInfo <- function(dataName=""){
  paramSet <- readSet(paramSet, "paramSet");
  #print(paste0("metainfo==dataname=", dataName));
  if(paramSet$anal.type == "metadata"){
  metaNms<-setdiff(colnames(paramSet$dataSet$meta.info),dataSet$rmMetaCol)
  }else{
  dataSet <- readDataset(dataName);
  metaNms<-setdiff(colnames(dataSet$meta.info),dataSet$rmMetaCol)
  }
  return(metaNms);
}

GetExpressResultfeaturesymbols<-function(){
  analSet <- readSet(analSet, "analSet");
  return(analSet$comp.features.symbols);
}

# Display names for DE result table.
# For phosphosite IDs like "Q9D1F4_T_247", replace the UniProt prefix with gene
# symbol when phospho_symbol_map.qs is available, while preserving site suffix.
GetExpressResultDisplayNames <- function(dataName=""){
  dataSet <- readDataset(dataName);
  ids <- rownames(dataSet$comp.res)
  ids <- as.character(ids)
  dim(ids) <- NULL

  # Default display labels from existing symbol vector (if present)
  labels <- ids
  analSet <- try(readSet(analSet, "analSet"), silent = TRUE)
  if (!inherits(analSet, "try-error") && !is.null(analSet$comp.features.symbols)) {
    syms <- as.character(analSet$comp.features.symbols)
    if (length(syms) == length(ids)) {
      labels <- syms
    }
  }

  # Phospho display enhancement: swap UniProt prefix for gene symbol using
  # precomputed phospho_symbol_map.qs
  # The phospho_symbol_map already contains full display names with isoform and site suffixes
  # (e.g., "GENESYMBOL-2_S_123"), so we use them directly
  phospho.map <- try(readDataQs("phospho_symbol_map.qs",
                                readSet(paramSet, "paramSet")$anal.type,
                                dataName), silent = TRUE)
  if (!inherits(phospho.map, "try-error") && !is.null(phospho.map) && nrow(phospho.map) > 0) {
    phospho.ids <- intersect(ids, rownames(phospho.map))
    if (length(phospho.ids) > 0 && "symbol" %in% colnames(phospho.map)) {
      map.syms <- as.character(phospho.map[phospho.ids, "symbol", drop = TRUE])
      for (i in seq_along(ids)) {
        id <- ids[i]
        if (!(id %in% phospho.ids)) next
        sym <- map.syms[which(phospho.ids == id)[1]]
        if (is.na(sym) || sym == "" || sym == "NA" || sym == id) next
        # Use the symbol directly - it already contains the full display name
        # with isoform and site suffix (e.g., "DOCK10-2_S_12")
        labels[i] <- sym
      }
    }
  }
  return(labels)
}

GetExpressResultGeneIDLinks <- function(dataName=""){
  dataSet <- readDataset(dataName);
  paramSet <- readSet(paramSet, "paramSet");
  ids <- rownames(dataSet$comp.res);

  if(paramSet$data.org == "generic"){
    if(paramSet$data.idType == "ko"){
      annots <- paste0("<a href='https://www.genome.jp/dbget-bin/www_bget?", ids, "' target='_blank'>KEGG</a>");
    } else if(paramSet$data.idType == "s2f"){
      annots <- paste0("<a href='https://www.ecoomicsdb.ca/#/query?ortho=", ids, "' target='_blank'>EODB</a>");
    } else {
      annots <- ids; # Keep as-is
    }
  } else if (paramSet$data.org == "custom"){
    annots <- ids; # Keep as-is
  } else {
    # For standard organisms, need to create appropriate links
    # Handle phosphosite IDs (e.g., Q8N3V7-2_S_123) and isoforms (e.g., P12345-2)
    # Strategy: Display gene symbol (with isoform) but link to base UniProt

    # Get display names (gene symbols) - this handles phospho mapping
    display_names <- GetExpressResultDisplayNames(dataName)

    # Extract base protein IDs from phosphosite IDs if present
    clean_ids <- ids

    # If IDs contain phosphosite suffix (protein_AA_position), extract protein part
    is_phospho <- grepl("_[STY]_\\d+$", ids, perl = TRUE)
    if (any(is_phospho)) {
      clean_ids[is_phospho] <- sub("_[STY]_\\d+$", "", ids[is_phospho])
    }

    # Remove isoform suffixes (e.g., -2, -3) for linking
    clean_ids_no_isoform <- sub("-\\d+$", "", clean_ids)

    # Check if IDs look like UniProt format (alphanumeric, typically 6-10 chars)
    looks_like_uniprot <- grepl("^[A-Z0-9]{6,10}$", clean_ids_no_isoform, perl = TRUE)

    # Initialize link URLs
    annots <- character(length(ids))

    if (any(looks_like_uniprot)) {
      # For UniProt IDs, create UniProt links (more reliable than NCBI for protein isoforms)
      annots[looks_like_uniprot] <- paste0(
        "<a href='https://www.uniprot.org/uniprot/",
        clean_ids_no_isoform[looks_like_uniprot],
        "' target='_blank'>",
        display_names[looks_like_uniprot],
        "</a>"
      )
    }

    # For non-UniProt IDs (likely Entrez), use NCBI links
    if (any(!looks_like_uniprot)) {
      annots[!looks_like_uniprot] <- paste0(
        "<a href='http://www.ncbi.nlm.nih.gov/gene?term=",
        clean_ids_no_isoform[!looks_like_uniprot],
        "' target='_blank'>",
        display_names[!looks_like_uniprot],
        "</a>"
      )
    }
  }
  return(annots);  # Ensure this is a character vector, NOT a list
}


GetExpressResultColNames<-function(){
  resT <- qs::qread("express.de.res.qs");
  colnames(resT);
}

GetExpressResultGeneIDs<-function(dataName=""){
    dataSet <- readDataset(dataName);
    ids <- rownames(dataSet$comp.res)
    ids <- as.character(ids)
    dim(ids) <- NULL
    return(ids)
}

GetExpressGeneIDType<-function(dataName=""){
  dataSet <- readDataset(dataName);
  return(dataSet$id.current);
}

GetExpressResultMatrix <- function(dataName = "", inxt) {
    dataSet  <- readDataset(dataName);
    paramSet <- readSet(paramSet, "paramSet");
    inxt     <- as.numeric(inxt)
    #msg("[GetExpressResultMatrix] method=", dataSet$de.method, " inxt=", inxt, " comp.res.list length=", length(dataSet$comp.res.list))

    # choose base columns and comparison-specific slice
    if (dataSet$de.method == "deseq2") {
        inx <- match("baseMean", colnames(dataSet$comp.res))
        res <- dataSet$comp.res.list[[inxt]];
    } else if (dataSet$de.method=="edger") {
        inx <- match("logCPM", colnames(dataSet$comp.res))
        res <- dataSet$comp.res.list[[inxt]];
    } else if (dataSet$de.method=="limma" || dataSet$de.method=="deqms" || dataSet$de.method=="wtt") {
        inx <- match("AveExpr", colnames(dataSet$comp.res))
        res <- dataSet$comp.res.list[[inxt]];
    } else {
        if (dataSet$de.method == "wtt") {
            inx <- match("t", colnames(dataSet$comp.res))
        }
        res <- dataSet$comp.res
        res <- res[, -(1:(inx - 1)), drop = FALSE]                    # ← fixed slice
        res <- res[rownames(dataSet$comp.res), , drop = FALSE]        # ← align rows
        res <- cbind(dataSet$comp.res[, inxt], res)
        colnames(res)[1] <- colnames(dataSet$comp.res)[inxt]
    }

    # reorder/significant handling only if we have rows
    if (!is.null(dataSet$comp.res) && nrow(dataSet$comp.res) > 0) {
        o <- with(dataSet$comp.res, order(P.Value, -abs(logFC), na.last = TRUE))
        dataSet$comp.res <- dataSet$comp.res[o, , drop = FALSE]
        dataSet$comp.res <- dataSet$comp.res[
            !(rownames(dataSet$comp.res) %in% rownames(dataSet$sig.mat)), ]
        dataSet$comp.res <- rbind(dataSet$sig.mat, dataSet$comp.res)
        dataSet$comp.res <- dataSet$comp.res[complete.cases(dataSet$comp.res), ]

        # Keep per-contrast results aligned to the reordered comp.res
        if (!is.null(dataSet$comp.res.list) && length(dataSet$comp.res.list) > 0) {
            target_order <- rownames(dataSet$comp.res)
            for (i in seq_along(dataSet$comp.res.list)) {
                res_i <- dataSet$comp.res.list[[i]]
                res_i <- res_i[target_order, colnames(res_i), drop = FALSE]
                dataSet$comp.res.list[[i]] <- res_i
            }
        }
    }

    ## --- now extract the column(s) for the return value -------
    if (dataSet$de.method %in% c("limma", "deqms", "edger", "deseq2", "wtt")) {
      # prefer comp.res.list entry; fall back to comp.res if populated
      if (!is.null(dataSet$comp.res.list) && length(dataSet$comp.res.list) >= inxt &&
          nrow(dataSet$comp.res.list[[inxt]]) > 0) {
          #msg("[GetExpressResultMatrix] using comp.res.list[[", inxt, "]] rows=", nrow(dataSet$comp.res.list[[inxt]]))
          res <- dataSet$comp.res.list[[inxt]]
          dataSet$comp.res <- res  # keep dataset in sync for downstream gene IDs
      } else if (!is.null(dataSet$comp.res) && nrow(dataSet$comp.res) > 0) {
          #msg("[GetExpressResultMatrix] comp.res.list empty; using comp.res rows=", nrow(dataSet$comp.res))
          res <- dataSet$comp.res
      }
    } else {
      res <- dataSet$comp.res[ , c(inxt, (inx+1):ncol(dataSet$comp.res)), drop = FALSE]
      res <- res[order(res$P.Value), ]
      colnames(res)[1] <- colnames(dataSet$comp.res)[inxt]
    }
    # Hard-align rows to comp.res when dimensions permit to keep IDs synchronized
    if (!is.null(dataSet$comp.res) && nrow(dataSet$comp.res) > 0) {
      crn <- rownames(dataSet$comp.res)
      dim(crn) <- NULL
      if (nrow(res) == length(crn)) {
        res <- res[crn, , drop = FALSE]
        rownames(res) <- crn
      } else {
        common <- intersect(crn, rownames(res))
        res <- res[common, , drop = FALSE]
        dataSet$comp.res <- dataSet$comp.res[common, , drop = FALSE]
        if (!is.null(dataSet$comp.res.list) && length(dataSet$comp.res.list) >= inxt) {
          dataSet$comp.res.list[[inxt]] <- dataSet$comp.res.list[[inxt]][common, , drop = FALSE]
        }
      }
    }
    # Global ordering: rank by P.Value then |logFC| for display consistency
    if (!is.null(res) && nrow(res) > 0 && "P.Value" %in% colnames(res) && "logFC" %in% colnames(res)) {
      o <- order(res$P.Value, -abs(res$logFC), na.last = TRUE)
      res <- res[o, , drop = FALSE]
      dataSet$comp.res <- res
      if (!is.null(dataSet$comp.res.list) && length(dataSet$comp.res.list) >= inxt) {
        dataSet$comp.res.list[[inxt]] <- res
      }
    }
    if (is.null(res)) {
      #msg("[GetExpressResultMatrix][WARN] result object is NULL")
      return(matrix(numeric(0), nrow = 0, ncol = 0))
    }
    if (nrow(res) == 0) {
      #msg("[GetExpressResultMatrix][WARN] result has 0 rows after processing")
    }

    # Final alignment: ensure the numeric matrix rows follow the master comp.res order
    if (!is.null(dataSet$comp.res) && nrow(dataSet$comp.res) > 0) {
      target_order <- rownames(dataSet$comp.res)
      res <- res[target_order, colnames(res), drop = FALSE]
      # Keep dataset/list in sync for downstream ID retrieval
      dataSet$comp.res <- res
      if (!is.null(dataSet$comp.res.list) && length(dataSet$comp.res.list) >= inxt) {
        dataSet$comp.res.list[[inxt]] <- res
      }
    }
    #msg("[GetExpressResultMatrix] returning matrix dims=", paste(dim(res), collapse="x"))

    RegisterData(dataSet)
    # Shadow save: both qs (backward compat) and Arrow (zero-copy Java access)
    shadow_save(res, "express.de.res.qs")
    return(head(signif(as.matrix(res), 5),1000))
}


###Gene list
GetNumOfLists <- function(){
  paramSet <- readSet(paramSet, "paramSet");
  return(paramSet$numOfLists)
}

GetNumOffeatures <- function(){
  paramSet <- readSet(paramSet, "paramSet");
  return(paramSet$genenum)
}


GetNumOfAnnofeatures <- function(){
  paramSet <- readSet(paramSet, "paramSet");
  return(paramSet$annonum)
}


GetMetaSigGeneCount <- function(){
  analSet <- readSet(analSet, "analSet");
  return(nrow(analSet$meta.mat));
}

GetCurrentJson <- function(type) {
  paramSet <- readSet(paramSet, "paramSet")
  
  # Check if the list paramSet$jsonNms contains the key 'type'
  if (!is.null(paramSet$jsonNms[[type]])) {
    return(paramSet$jsonNms[[type]])
  } else {
    return("NA")  # or return a default value if appropriate
  }
}


SelectDataSet <- function(){
  paramSet <- readSet(paramSet, "paramSet");
  if(!exists('nm.vec')){
    AddErrMsg("No dataset is selected for analysis!");
    return(0);
  }
  mdata.all <- paramSet$mdata.all
  all.nms <- names(mdata.all);
  for(nm in all.nms){
    if(nm %in% nm.vec){
      mdata.all[[nm]] <- 1;
    }else{
      mdata.all[[nm]] <- 0;
    }
  }
  
  if("meta_dat" %in% nm.vec){
    meta.selected <<- TRUE;
  }else{
    meta.selected <<- FALSE;
  }
  
  rm('nm.vec', envir = .GlobalEnv);

  paramSet$mdata.all <- mdata.all
  saveSet(paramSet, "paramSet");
  return(1);
  
}


GetFeatureNum <- function(dataName) {

  dataSet <- readDataset(dataName, quiet = TRUE)

  if (is.null(dataSet) || is.null(dataSet$data.norm)) {
    return(0L)                               # nothing loaded → report zero
  }

  nrow(dataSet$data.norm)
}

# get qualified inx with at least number of replicates
GetDiscreteInx <- function(my.dat, min.rep=2){
  good.inx <- apply(my.dat, 2, function(x){
    x <- x[x!="NA"]
    good1.inx <- length(x) > length(unique(x));
    good2.inx <- min(table(x)) >= min.rep;
    return (good1.inx & good2.inx);
  });
  return(good.inx);
}

GetNumbericalInx <- function(my.dat){
  suppressWarnings({
  good.inx <- apply(my.dat, 2, function(x){
    isNum = as.numeric(as.character(x[x!="NA"]))
    return(all(!is.na(as.numeric(as.character(isNum)))));
  });
  })
  return(good.inx);
}

.set.dataSet <- function(dataSetObj=NA){
  RegisterData(dataSetObj);
  return (1);
}

# remove data object, the current dataSet will be the last one by default 
RemoveData <- function(dataName){
  paramSet <- readSet(paramSet, "paramSet");
  mdata.all <- paramSet$mdata.all;
  if(!is.null(paramSet$mdata.all[[dataName]])){
    paramSet$mdata.all[[dataName]] <- NULL;
  }
  saveSet(paramSet, "paramSet");
}


GetCovSigFileName <-function(dataName){
  dataSet <- readDataset(dataName);
  dataSet$analSet$cov$sig.nm;
}

GetCovSigMat<-function(dataName){
  dataSet <- readDataset(dataName);
  drops <- c("ids","label")
  return(CleanNumber(as.matrix(dataSet$analSet$cov$sig.mat[, !(names(dataSet$analSet$cov$sig.mat) %in% drops)])));
}

GetCovSigIds<-function(dataName){
  dataSet <- readDataset(dataName);
  dataSet$analSet$cov$sig.mat$ids;
}

GetCovSigSymbols<-function(dataName){
  dataSet <- readDataset(dataName);
  dataSet$analSet$cov$sig.mat$label
}

GetCovSigColNames<-function(dataName){
  dataSet <- readDataset(dataName);
  drops <- c("ids","label");
  colnames(dataSet$analSet$cov$sig.mat[,!(names(dataSet$analSet$cov$sig.mat) %in% drops)]);
}

GetCovDENums <- function(dataName){
    deNum <- nrow(dataSet$analSet$cov$sig.mat);
    nonDeNum <- nrow(dataSet$comp.res) - deNum;
    return(c(deNum, nonDeNum));
}


#'Replace infinite numbers
#'@description Replace -Inf, Inf to 99999 and -99999
#'@param bdata Input matrix to clean numbers
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'
CleanNumber <-function(bdata){
  if(sum(bdata==Inf)>0){
    inx <- bdata == Inf;
    bdata[inx] <- NA;
    bdata[inx] <- 999999;
  }
  if(sum(bdata==-Inf)>0){
    inx <- bdata == -Inf;
    bdata[inx] <- NA;
    bdata[inx] <- -999999;
  }
  bdata;
}

GetMetaMethodPVal <-function(){
  paramSet <- readSet(paramSet, "paramSet");
  return(paramSet$BHth);
}

#for enrichment analysis
SetUniverseOpt <- function(universe.opt){
  paramSet <- readSet(paramSet, "paramSet");
  paramSet$universe.opt <- universe.opt;
  if(paramSet$universe.opt == "uploaded"){
    paramSet$universe.opt.readable <- "Uploaded Data";
  }else{
    paramSet$universe.opt.readable <- "Gene Set Library";
  }
    saveSet(paramSet, "paramSet");
}

SetInitLib <-function(library){
  paramSet <- readSet(paramSet, "paramSet");
  paramSet$init.lib <- library;
    saveSet(paramSet, "paramSet");
}
