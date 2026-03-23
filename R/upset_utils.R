##################################################
## R script for ProteoAnalyst
## Description: Compute upset diagram
## Authors: 
## G. Zhou, guangyan.zhou@mail.mcgill.ca
###################################################


#'Prepare data for Upset diagram
#'@param dataSetObj Input the name of the created dataSetObj (see Init.Data)
#'@param fileNm file name of the json file output 
#'@export
PrepareUpsetData <- function(fileNm){
  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");

  mdata.all <- paramSet$mdata.all;
  anal.type <- paramSet$anal.type;
  
  newDat <- list();

  # selected dataset or comparisons for onedata (single gene expression matrix)
  if(anal.type == "metadata"){
  hit.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[hit.inx];
  }else if(anal.type == "onedata"){
  
  }else{
  sel.nms <- names(mdata.all)
  }
  sel.dats <- list();
  

    if(!exists("analSet$inmex.de")){
      analSet$inmex.de <- list();
    }

  # populate gene lists for upset plot based on selected names
  for(nm in sel.nms){
    if(anal.type == "metadata"){
      dataSet <- readDataset(nm);
      sel.dats[[nm]] <- dataSet$sig.mat$ids;
    }else if(anal.type == "onedata"){

    }else{
      dataSet <- readDataset(nm);
      gene.mat <- dataSet$prot.mat;
      
      # convert to entrez
      expr.val <- gene.mat[,1];
      en.ids <- rownames(gene.mat);
      
      names(expr.val) <- en.ids;
      analSet$inmex.de[[nm]] <- en.ids;
      sel.dats[[nm]] <- en.ids;
    }

  }

  if(anal.type == "metadata" & paramSet$meta.selected){
    sel.dats[["meta_dat"]] <- as.character(rownames(analSet$meta.mat));
  }
  
  if(length(sel.dats) == 0){
    AddErrMsg("No signficant features for any dataset!");
    return(0);
  }
  
  sums <- unlist(lapply(sel.dats, length));
  names <- unlist(lapply(sel.dats, paste, collapse = ", "));
  if(anal.type == "metadata"){
    metasum <- length(analSet$meta.stat$idd);
    metaname <- paste(analSet$meta.stat$idd, collapse = ", ");
    allsums <- c(sums, metasum);
    allnames <- c(names, metaname);
  }else{
    allsums <- c(sums);
    allnames <- c(names);
  }


  require(reshape2)
  df <- reshape::melt(sel.dats, value.name="id")
  colnames(df) <- c("name", 'set')
  uniq.nms <- unique(df$name)
  new.df <- dcast(df, name ~ set, value.var='set', fill=0)
  rownames(new.df) <- new.df[,1]
  new.df <- new.df[,-1, drop=F]
  
  gene.map <-  queryGeneDB("entrez", paramSet$data.org);
  gene.map[] <- lapply(gene.map, as.character)
  
  json.list <- list()
  for(i in 1:nrow(new.df)){
    json.list[[i]] <- list()
    json.list[[i]][["sets"]] <- new.df[i,][new.df[i,] != 0]
    entrez.vec <- rownames(new.df)[i];
    hit.inx <- match(entrez.vec, gene.map[, "gene_id"]);
    symbols <- gene.map[hit.inx, "symbol"];
    
    # if not gene symbol, use id by itself
    na.inx <- is.na(symbols);
    symbols[na.inx] <- entrez.vec[na.inx];
    json.list[[i]][["name"]] <- symbols;
    json.list[[i]][["entrez"]] <- entrez.vec;
  }
  
  col.vec <-gg_color_hue(length(sel.dats));
  
  jsonNm <- paste0(fileNm, ".json")
  json.mat <- RJSONIO::toJSON(list(json.list, col.vec));
  sink(jsonNm);
  cat(json.mat);
  sink();
  
  return(1); 
}

#'Prepare data for Upset diagram using DE comparisons (onedata)
#'@param fileNm file name of the json file output
#'@param dataName optional dataset name (defaults to paramSet$dataName)
#'@export
PrepareUpsetDataFromComparisons <- function(fileNm, dataName = ""){
  paramSet <- readSet(paramSet, "paramSet");
  dataSet <- readDataset(dataName);

  if (is.null(dataSet) || is.null(dataSet$comp.res.list) ||
      length(dataSet$comp.res.list) == 0) {
    AddErrMsg("No DE comparison results found.");
    return(0);
  }

  comp_list <- dataSet$comp.res.list
  comp_names <- names(comp_list)
  if (is.null(comp_names) || length(comp_names) == 0) {
    comp_names <- paste0("Comparison_", seq_along(comp_list))
  }
  if (!is.null(paramSet$upset.comp.nms) && length(paramSet$upset.comp.nms) > 0) {
    keep <- comp_names %in% paramSet$upset.comp.nms
    comp_list <- comp_list[keep]
    comp_names <- names(comp_list)
  }
  if (is.null(comp_names) || length(comp_names) == 0) {
    AddErrMsg("No comparisons selected.");
    return(0);
  }

  p.lvl <- paramSet$BHth
  if (is.null(p.lvl) || !is.finite(p.lvl)) {
    p.lvl <- 0.05
  }
  fc.lvl <- dataSet$fc.lvl
  if (is.null(fc.lvl) || !is.finite(fc.lvl)) {
    fc.lvl <- paramSet$fc.thresh
  }
  if (is.null(fc.lvl) || !is.finite(fc.lvl)) {
    fc.lvl <- 0
  }

  use_fdr <- if (is.null(paramSet$use.fdr)) TRUE else isTRUE(paramSet$use.fdr)

  msg("[PrepareUpsetDataFromComparisons] Processing ", length(comp_list), " comparisons")
  msg("[PrepareUpsetDataFromComparisons] Thresholds - p-value: ", p.lvl, ", FC: ", fc.lvl, ", use FDR: ", use_fdr)

  sel.dats <- list()
  for (i in seq_along(comp_list)) {
    rt <- comp_list[[i]]
    if (is.null(rt) || nrow(rt) == 0) {
      sel.dats[[comp_names[[i]]]] <- character(0)
      msg("[PrepareUpsetDataFromComparisons] Comparison ", comp_names[[i]], ": empty result table")
      next
    }

    msg("[PrepareUpsetDataFromComparisons] Comparison ", comp_names[[i]], ": ", nrow(rt), " features")
    msg("[PrepareUpsetDataFromComparisons] Sample rownames: ", paste(head(rownames(rt), 5), collapse = ", "))

    pcol <- NULL
    if (use_fdr && "adj.P.Val" %in% names(rt)) {
      pcol <- "adj.P.Val"
    } else if ("P.Value" %in% names(rt)) {
      pcol <- "P.Value"
    } else if ("padj" %in% names(rt)) {
      pcol <- "padj"
    } else if ("pvalue" %in% names(rt)) {
      pcol <- "pvalue"
    }

    lfc_col <- if ("logFC" %in% names(rt)) {
      "logFC"
    } else if ("log2FoldChange" %in% names(rt)) {
      "log2FoldChange"
    } else {
      NULL
    }

    msg("[PrepareUpsetDataFromComparisons] Using p-value column: ", ifelse(is.null(pcol), "NONE", pcol))
    msg("[PrepareUpsetDataFromComparisons] Using LFC column: ", ifelse(is.null(lfc_col), "NONE", lfc_col))

    deg_pass <- if (!is.null(pcol)) {
      is.finite(rt[[pcol]]) & (rt[[pcol]] <= p.lvl)
    } else {
      rep(TRUE, nrow(rt))
    }
    lfc_pass <- if (!is.null(lfc_col)) {
      is.finite(rt[[lfc_col]]) & (abs(rt[[lfc_col]]) >= fc.lvl)
    } else {
      rep(TRUE, nrow(rt))
    }

    sig_ids <- as.character(rownames(rt)[deg_pass & lfc_pass])
    sel.dats[[comp_names[[i]]]] <- sig_ids
    msg("[PrepareUpsetDataFromComparisons] Comparison ", comp_names[[i]], ": ", length(sig_ids), " significant features")
    if (length(sig_ids) > 0) {
      msg("[PrepareUpsetDataFromComparisons] Sample sig IDs: ", paste(head(sig_ids, 5), collapse = ", "))
    }
  }

  sel.dats <- sel.dats[sapply(sel.dats, length) > 0]
  if (length(sel.dats) == 0) {
    AddErrMsg("No significant features for any comparison.");
    return(0);
  }
  msg("[PrepareUpsetDataFromComparisons] Total comparisons with significant features: ", length(sel.dats))

  require(reshape2)
  df <- reshape::melt(sel.dats, value.name="id")
  colnames(df) <- c("name", 'set')

  msg("[PrepareUpsetDataFromComparisons] Melted dataframe has ", nrow(df), " rows")
  msg("[PrepareUpsetDataFromComparisons] Sample IDs in melted df: ", paste(head(df$name, 10), collapse = ", "))

  new.df <- dcast(df, name ~ set, value.var='set', fill=0)
  rownames(new.df) <- new.df[,1]
  new.df <- new.df[,-1, drop=FALSE]

  msg("[PrepareUpsetDataFromComparisons] Cast dataframe has ", nrow(new.df), " unique features")
  msg("[PrepareUpsetDataFromComparisons] Sample rownames: ", paste(head(rownames(new.df), 10), collapse = ", "))

  gene.map <- queryGeneDB("entrez", paramSet$data.org);
  gene.map[] <- lapply(gene.map, as.character)
  msg("[PrepareUpsetDataFromComparisons] Gene map has ", nrow(gene.map), " entries for organism ", paramSet$data.org)

  json.list <- list()
  unmapped_count <- 0
  for(i in 1:nrow(new.df)){
    json.list[[i]] <- list()
    json.list[[i]][["sets"]] <- new.df[i,][new.df[i,] != 0]
    entrez.vec <- rownames(new.df)[i]
    hit.inx <- match(entrez.vec, gene.map[, "gene_id"])
    symbols <- gene.map[hit.inx, "symbol"]

    na.inx <- is.na(symbols)
    if (na.inx) {
      unmapped_count <- unmapped_count + 1
      if (unmapped_count <= 10) {
        msg("[PrepareUpsetDataFromComparisons] Warning: No symbol found for ID: ", entrez.vec)
      }
    }
    symbols[na.inx] <- entrez.vec[na.inx]
    json.list[[i]][["name"]] <- symbols
    json.list[[i]][["entrez"]] <- entrez.vec
  }

  if (unmapped_count > 0) {
    msg("[PrepareUpsetDataFromComparisons] Total unmapped features: ", unmapped_count, " / ", nrow(new.df))
  }

  col.vec <- gg_color_hue(length(sel.dats))

  jsonNm <- paste0(fileNm, ".json")
  json.mat <- RJSONIO::toJSON(list(json.list, col.vec))
  sink(jsonNm)
  cat(json.mat)
  sink()

  return(1)
}

GetComparisonSigCounts <- function(dataName = "") {
  paramSet <- readSet(paramSet, "paramSet");
  dataSet <- readDataset(dataName);

  if (is.null(dataSet) || is.null(dataSet$comp.res.list) ||
      length(dataSet$comp.res.list) == 0) {
    return(integer(0))
  }

  comp_list <- dataSet$comp.res.list
  comp_names <- names(comp_list)
  if (is.null(comp_names) || length(comp_names) == 0) {
    comp_names <- paste0("Comparison_", seq_along(comp_list))
  }

  p.lvl <- paramSet$BHth
  if (is.null(p.lvl) || !is.finite(p.lvl)) {
    p.lvl <- 0.05
  }
  fc.lvl <- dataSet$fc.lvl
  if (is.null(fc.lvl) || !is.finite(fc.lvl)) {
    fc.lvl <- paramSet$fc.thresh
  }
  if (is.null(fc.lvl) || !is.finite(fc.lvl)) {
    fc.lvl <- 0
  }

  use_fdr <- if (is.null(paramSet$use.fdr)) TRUE else isTRUE(paramSet$use.fdr)

  counts <- integer(length(comp_list))
  for (i in seq_along(comp_list)) {
    rt <- comp_list[[i]]
    if (is.null(rt) || nrow(rt) == 0) {
      counts[[i]] <- 0
      next
    }

    pcol <- NULL
    if (use_fdr && "adj.P.Val" %in% names(rt)) {
      pcol <- "adj.P.Val"
    } else if ("P.Value" %in% names(rt)) {
      pcol <- "P.Value"
    } else if ("padj" %in% names(rt)) {
      pcol <- "padj"
    } else if ("pvalue" %in% names(rt)) {
      pcol <- "pvalue"
    }

    lfc_col <- if ("logFC" %in% names(rt)) {
      "logFC"
    } else if ("log2FoldChange" %in% names(rt)) {
      "log2FoldChange"
    } else {
      NULL
    }

    deg_pass <- if (!is.null(pcol)) {
      is.finite(rt[[pcol]]) & (rt[[pcol]] <= p.lvl)
    } else {
      rep(TRUE, nrow(rt))
    }
    lfc_pass <- if (!is.null(lfc_col)) {
      is.finite(rt[[lfc_col]]) & (abs(rt[[lfc_col]]) >= fc.lvl)
    } else {
      rep(TRUE, nrow(rt))
    }

    counts[[i]] <- sum(deg_pass & lfc_pass, na.rm = TRUE)
  }

  names(counts) <- comp_names
  counts
}

SetUpsetComparisons <- function(){
  paramSet <- readSet(paramSet, "paramSet");
  if (!exists("comp.vec")) {
    paramSet$upset.comp.nms <- NULL
    saveSet(paramSet, "paramSet")
    return(0)
  }
  paramSet$upset.comp.nms <- comp.vec
  rm("comp.vec", envir = .GlobalEnv)
  saveSet(paramSet, "paramSet")
  return(1)
}

#Record upset intersection mode for report
SetUpsetMode <- function(mode){
      paramSet <- readSet(paramSet, "paramSet");
      paramSet$upsetMode <- mode;
  saveSet(paramSet, "paramSet");
}
