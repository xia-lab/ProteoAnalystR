##################################################
## R script for ProteoAnalyst
## Description: Functions for enrichment analysis (GSEA and ORA)
## Authors:
## G. Zhou, guangyan.zhou@mail.mcgill.ca
## Jeff Xia, jeff.xia@mcgill.ca
###################################################  
.performEnrichAnalysis <- function(dataSet, file.nm, fun.type, ora.vec, vis.type){

  dataSet <<- dataSet;

  msgSet <- readSet(msgSet, "msgSet");
  paramSet <- readSet(paramSet, "paramSet");
  require(dplyr)

  # prepare lib
  setres <- .loadEnrichLib(fun.type, paramSet)
  current.featureset <- setres$current.featureset;
  current.setids <<- setres$current.setids

  # prepare query
  ora.nms <- names(ora.vec);

  if(is.null(ora.nms)){
    ora.nms <- ora.vec;
    names(ora.vec) <- ora.vec;
  }

  # PROTEOMICS: Convert UniProt IDs to Entrez IDs for enrichment
  # Store original UniProt IDs for reverse mapping later
  # NOTE: After upload and ID mapping, ProteoAnalyst always uses UniProt IDs internally
  # So we ALWAYS perform this conversion (paramSet$data.idType refers to the initial upload type)
  original.uniprot.vec <- ora.vec;
  original.ora.nms <- ora.nms;

  # Check if ora.vec already contains Entrez IDs (all numeric)
  # Filter out NA values before checking (grepl returns FALSE for NA inputs)
  # Use >= 0.8 threshold to be more permissive (allows up to 20% NA/non-numeric)
  non_na_vec <- ora.vec[!is.na(ora.vec)]
  is.entrez.like <- length(non_na_vec) > 0 && mean(grepl("^[0-9]+$", non_na_vec)) >= 0.8

  if (is.entrez.like) {
    # IDs are already Entrez, no conversion needed
    entrez.vec <- ora.vec
    sym.vec <- doEntrez2SymbolMapping(entrez.vec, paramSet$data.org, "entrez")

    # For reverse mapping, try to get UniProt IDs from Entrez
    uniprot.map <- queryGeneDB("entrez_uniprot", paramSet$data.org)
    hit.inx <- match(entrez.vec, uniprot.map[, "gene_id"])
    uniprot.vec <- uniprot.map[hit.inx, "accession"]

    # Create reverse mapping (Entrez → UniProt)
    uniprot.to.entrez.map <- entrez.vec
    names(uniprot.to.entrez.map) <- ifelse(is.na(uniprot.vec), entrez.vec, uniprot.vec)
    original.uniprot.vec <- ifelse(is.na(uniprot.vec), entrez.vec, uniprot.vec)

    # Use Entrez IDs for matching
    ora.vec <- entrez.vec
    ora.nms <- sym.vec

  } else {
    # IDs are UniProt, convert to Entrez

    # Normalize UniProt IDs (remove phosphosite annotations, isoforms, etc.)
    # e.g., Q9D1F4_T_247 → Q9D1F4, P12345-1 → P12345
    normalized.ora.vec <- ora.vec
    normalized.ora.vec <- sub("_[A-Z]_\\d+$", "", normalized.ora.vec)  # Remove phosphosites
    normalized.ora.vec <- sub("-\\d+$", "", normalized.ora.vec)         # Remove isoforms
    normalized.ora.vec <- trimws(normalized.ora.vec)

    # Query UniProt → Entrez mapping
    uniprot.map <- queryGeneDB("entrez_uniprot", paramSet$data.org)
    hit.inx <- match(normalized.ora.vec, uniprot.map[, "accession"])
    entrez.vec <- uniprot.map[hit.inx, "gene_id"]

    # Get symbols for display names
    sym.vec <- doEntrez2SymbolMapping(entrez.vec, paramSet$data.org, "entrez");

    na.inx <- is.na(entrez.vec)

    # Create UniProt → Entrez mapping dictionary for reverse conversion
    uniprot.to.entrez.map <- entrez.vec
    names(uniprot.to.entrez.map) <- ora.vec

    # Use Entrez IDs for matching, named by symbols
    ora.vec <- entrez.vec
    ora.nms <- sym.vec

    # Remove NAs (unmapped UniProts)
    ora.vec <- ora.vec[!na.inx]
    ora.nms <- ora.nms[!na.inx]
    original.uniprot.vec <- original.uniprot.vec[!na.inx]
    uniprot.to.entrez.map <- uniprot.to.entrez.map[!na.inx]

  }

  # Store mapping in analSet for downstream use (e.g., network analysis)
  analSet <- readSet(analSet, "analSet");
  analSet$uniprot_to_entrez_map <- uniprot.to.entrez.map
  saveSet(analSet, "analSet");
  
  if(paramSet$universe.opt == "library"){
    current.universe <- unique(unlist(current.featureset));
  }else{
    # cut to the universe to uploaded features
    if(paramSet$anal.type == "onedata"){
      data.anot <- .get.annotated.data();
      current.universe <- rownames(data.anot);
    }else if(paramSet$anal.type == "metadata"){
      inmex <- ov_qs_read("inmex_meta.qs");
      current.universe <- rownames(inmex$data);
    }else{
      if(!is.null(paramSet$backgroundUniverse)){
        current.universe <- paramSet$backgroundUniverse;
      }else{
        current.universe <- unique(unlist(current.featureset));
      }
    }
  }
  
  # also make sure pathways only contain features measured in experiment
  #if(!is.null(dataSet$data.anot)){
   if(file.exists("data.anot.qs")){
    current.featureset <- lapply(current.featureset, function(x){x[x %in% current.universe]})
    inds <- lapply(current.featureset, length) > 0
    current.featureset <- current.featureset[inds]
  }

  # prepare for the result table
  set.size<-length(current.featureset);
  res.mat<-matrix(0, nrow=set.size, ncol=5);
  rownames(res.mat)<-names(current.featureset);
  colnames(res.mat)<-c("Total", "Expected", "Hits", "Pval", "FDR");

  q.size<-length(ora.vec);

  # get the matched query for each pathway
  hits.query <- lapply(current.featureset,
                       function(x) {
                         ora.nms[ora.vec%in%unlist(x)];
                       }
  );

  # PROTEOMICS: Convert hits back to UniProt IDs (reverse mapping)
  if(!is.null(uniprot.to.entrez.map)){

    # Create symbol → UniProt mapping
    symbol.to.uniprot.map <- original.uniprot.vec
    names(symbol.to.uniprot.map) <- ora.nms

    # Convert each hit list from symbols to UniProt IDs
    hits.query <- lapply(hits.query, function(symbol.list) {
      # Map symbols to UniProt IDs
      uniprot.list <- symbol.to.uniprot.map[symbol.list]
      # Remove NAs (shouldn't happen but be safe)
      uniprot.list <- uniprot.list[!is.na(uniprot.list)]
      # Remove names for clean JSON output (JavaScript expects array of IDs, not named vector)
      unname(uniprot.list)
    })

  }

  ov_qs_save(hits.query, "hits_query.qs");

  names(hits.query) <- names(current.featureset);
  hit.num<-unlist(lapply(hits.query, function(x){length(unique(x))}), use.names=FALSE);
  
  gene.vec <- current.universe;
  sym.vec <- doEntrez2SymbolMapping(gene.vec, paramSet$data.org, paramSet$data.idType);
  gene.nms <- sym.vec;

  current.featureset.symb <- lapply(current.featureset, 
                       function(x) {
                         gene.nms[gene.vec%in%unlist(x)];
  }
  );

  # total unique gene number
  uniq.count <- length(current.universe);
  
  # unique gene count in each pathway
  set.size <- unlist(lapply(current.featureset, length));
  
  res.mat[,1]<-set.size;
  res.mat[,2]<-q.size*(set.size/uniq.count);
  res.mat[,3]<-hit.num;
  
  # use lower.tail = F for P(X>x)
  raw.pvals <- phyper(hit.num-1, set.size, uniq.count-set.size, q.size, lower.tail=F);
  # Replace NaN values with 1
  raw.pvals[is.nan(raw.pvals)] <- 1
  
  res.mat[,4]<- raw.pvals;
  res.mat[,5] <- p.adjust(raw.pvals, "fdr");
  
  # now, clean up result, synchronize with hit.query
  res.mat <- res.mat[hit.num>0,,drop = F];
  hits.query <- hits.query[hit.num>0];

  if(nrow(res.mat)> 1){
    # order by p value
    ord.inx<-order(res.mat[,4]);
    res.mat <- signif(res.mat[ord.inx,],3);
    hits.query <- hits.query[ord.inx];
    
    res.mat.all <- as.data.frame(res.mat);
    res.mat.all$Pathway <- rownames(res.mat);
    res.mat.all$features <- rep("NA",nrow(res.mat))
    # Iterate through the list and add comma-separated values to the data frame
    for (name in names(hits.query)) {
      if (name %in% res.mat.all$Pathway) {
        res.mat.all[which(res.mat.all$Pathway == name), "features"] <- paste(hits.query[[name]], collapse = ",")
      }
    }
    
    res.mat.all <- res.mat.all[which(res.mat.all$features != "NA"), ];
    res.mat.all$Pathway <- NULL;
    pws <- rownames(res.mat[which(res.mat.all$features != "NA"), ]) 
    fun.ids2 <- as.vector(setres$current.setids[pws]) 
    resTable.all <- data.frame(Pathway = pws, ID = fun.ids2, res.mat.all)

    csv.nm <- paste(file.nm, ".csv", sep="");
    #print(paste(csv.nm, "=====", "enrichAnalysis"));
    write.csv(resTable.all, file=csv.nm, row.names=F);

    imp.inx <- res.mat[,4] <= 0.05;
    imp.inx[is.na(imp.inx)] <- F

    if(sum(imp.inx) < 10){ # too little left, give the top ones
      topn <- ifelse(nrow(res.mat) > 10, 10, nrow(res.mat));
      res.mat <- res.mat[1:topn,];
      hits.query <- hits.query[1:topn];
    }else{
      res.mat <- res.mat[imp.inx,];
      hits.query <- hits.query[imp.inx];
      if(sum(imp.inx) > 120){
        # now, clean up result, synchronize with hit.query
        res.mat <- res.mat[1:120,];
        hits.query <- hits.query[1:120];
      }
    }
  }else{
    # msgSet$current.msg <- "No overlap between queried features and pathway library!"
    saveSet(msgSet, "msgSet");
    return(0);
  }
  
  # Check for and handle duplicate row names in enr.mat
  if(any(duplicated(rownames(res.mat)))) {
    res.mat <- res.mat[!duplicated(rownames(res.mat)), ]
    hits.query <- hits.query[match(rownames(res.mat), names(hits.query))]
    print("Duplicates in enr.mat were removed.")
  } else {
    res.mat <- res.mat
  }

  resTable <- data.frame(Pathway=rownames(res.mat), res.mat);

  # msg("[EnrichAnalysis] Saving enr.mat.qs to: ", getwd())
  # msg("[EnrichAnalysis] res.mat dimensions: ", nrow(res.mat), " x ", ncol(res.mat))
  ov_qs_save(res.mat, "enr.mat.qs");

  # Verify file was written
  if (file.exists("enr.mat.qs")) {
  } else {
  }


  # write json
  fun.anot <- hits.query; 
  total <- resTable$Total; if(length(total) ==1) { total <- matrix(total) };
  fun.pval <- resTable$Pval; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval) };
  fun.padj <- resTable$FDR; if(length(fun.padj) ==1) { fun.padj <- matrix(fun.padj) };
  #print(resTable$Hits);
  hit.num <- paste0(resTable$Hits,"/",resTable$Total); if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
  fun.ids <- as.vector(setres$current.setids[resTable$Pathway]); 
  
  resTable$IDs <- fun.ids;
  if(length(fun.ids) ==1) { fun.ids <- matrix(fun.ids) };
  json.res <- list(
    fun.link = setres$current.setlink[1],
    fun.anot = fun.anot,
    fun.ids = fun.ids,
    fun.pval = fun.pval,
    fun.padj = fun.padj,
    hit.num = hit.num,
    total= total
  );
  json.mat <- rjson::toJSON(json.res);
  json.nm <- paste(file.nm, ".json", sep="");
  
  sink(json.nm)
  cat(json.mat);
  sink();
  
  # write csv
  fun.hits <<- hits.query;
  fun.pval <<- fun.pval;
  hit.num <<- resTable$Hits;
  #csv.nm <- paste(file.nm, ".csv", sep="");  
  #print(csv.nm)  
  #fast.write(resTable, file=csv.nm, row.names=F);
  paramSet$partialToBeSaved <- c(paramSet$partialToBeSaved, c(json.nm))
  
  imgSet <- readSet(imgSet, "imgSet");
  rownames(resTable) <- NULL;
  imgSet$enrTables[[vis.type]] <- list()
  imgSet$enrTables[[vis.type]]$table <- resTable;
  imgSet$enrTables[[vis.type]]$library <- fun.type
  imgSet$enrTables[[vis.type]]$algo <- "Overrepresentation Analysis"

    imgSet$enrTables[[vis.type]]$current.featureset <- current.featureset;
    imgSet$enrTables[[vis.type]]$hits.query <- hits.query;
    imgSet$enrTables[[vis.type]]$current.setids <- current.setids;
    imgSet$enrTables[[vis.type]]$res.mat<- res.mat;
    imgSet$enrTables[[vis.type]]$current.featureset.symb <- current.featureset.symb;
  
    saveSet(imgSet, "imgSet");
  saveSet(paramSet, "paramSet");
  
  saveSet(msgSet, "msgSet");
  return(1);
}

.loadEnrichLib <- function(fun.type, paramSet){
  #if custom return here.
  if(fun.type == "custom"){
    return(.loadCustomEnrichLib(fun.type, paramSet));
  }

  if(paramSet$data.org == "generic"){
    folderNm <- paramSet$data.idType;
  }else{
    folderNm <- paramSet$data.org;
  }

  # For KEGG, prefer the jointpa genetic library from MetaboAnalyst when available
  if(fun.type == "kegg" && !is.null(paramSet$jointpa.lib.path)){
    jointpa.file <- paste0(paramSet$jointpa.lib.path, folderNm, ".qs");
    if(file.exists(jointpa.file)){
      return(.loadJointpaGeneticLib(jointpa.file));
    }
  }

  if(exists("api.lib.path")){
    lib.path <- api.lib.path;
  }else{
    lib.path <- paramSet$lib.path;
  }

  my.path <- paste(lib.path, folderNm, "/", fun.type, ".rds", sep="");
  if(!paramSet$on.public.web && !file.exists(platform.path)){
    nmdb <- basename(my.path);
    download.file(my.path, destfile = nmdb, method="libcurl", mode = "wb");
    my.path <- nmdb;
  }

  if (!file.exists(my.path)) {
    AddErrMsg(paste0("[EnrichLib] Library file not found: ", my.path));
    return(0);
  }
  my.lib <- tryCatch(
    readRDS(my.path),
    error = function(e) {
      stop("[EnrichLib] Failed to read library file: ", my.path, " | ", conditionMessage(e))
    }
  );

  if(substr(fun.type, 0, 2)=="go"){
    if(is.null(names(my.lib))){ # some go lib does not give names
      names(my.lib) <- c("link", "term", "sets");
    }
  }

  current.featureset <- my.lib$sets;
  #remove empty pathways
  keep.inx <- lapply(current.featureset,length)>0
  current.featureset <- current.featureset[keep.inx]
  my.lib$term <- my.lib$term[keep.inx]
  set.ids<- names(current.featureset);
  names(set.ids) <- names(current.featureset) <- my.lib$term;
  
  if(substr(fun.type, 0, 2)=="go"){
    names(current.featureset) = firstup(names(current.featureset))
    names(current.featureset) = gsub("-", "_", names(current.featureset))
    names(set.ids) = firstup(names(set.ids));
    names(set.ids) = gsub("-", "_", names(set.ids));
  }
  ov_qs_save(current.featureset, "current_featureset.qs");
  res <- list();
  res$current.setlink <- my.lib$link;
  res$current.setids <- set.ids;
  res$current.featureset <- current.featureset;
  return(res);
}

# Load MetaboAnalyst jointpa genetic .qs file and return the same structure as .loadEnrichLib()
# jointpa mset.list: keyed by pathway IDs, values = space-separated "hsa:1234" gene ID strings
# jointpa path.ids: named by pathway names, values = pathway IDs
.loadJointpaGeneticLib <- function(file.path) {
  jlib <- tryCatch(
    ov_qs_read(file.path),
    error = function(e) {
      stop("[JointpaLib] Failed to read: ", file.path, " | ", conditionMessage(e))
    }
  );

  # Reverse path.ids: pathway_ID -> pathway_name
  rev.ids <- setNames(names(jlib$path.ids), jlib$path.ids);

  # Parse mset.list: each entry is a character vector of space-separated KEGG ID groups
  # Strip "org:" prefix, filter out "cpd:" entries and NAs, deduplicate
  sets <- lapply(jlib$mset.list, function(x) {
    ids <- unlist(strsplit(x, " ", fixed = TRUE));
    ids <- ids[!is.na(ids)];
    ids <- ids[!grepl("^cpd:", ids, perl = TRUE)];
    ids <- unique(gsub("^[a-z]+:", "", ids, perl = TRUE));
    ids[nchar(ids) > 0]
  });
  # sets is keyed by pathway IDs; term = pathway names (parallel)
  term <- rev.ids[names(sets)];

  keep.inx <- sapply(sets, length) > 0;
  sets <- sets[keep.inx];
  term <- term[keep.inx];

  set.ids <- names(sets);  # pathway IDs
  names(set.ids) <- names(sets) <- term;  # re-key by pathway names

  ov_qs_save(sets, "current_featureset.qs");
  res <- list();
  res$current.setlink <- "https://www.genome.jp/dbget-bin/www_bget?pathway+";
  res$current.setids <- set.ids;
  res$current.featureset <- sets;
  return(res);
}

GetRidgePlot <- function(dataName, imgNm = "abc", dpi=96, format="png", fun.type = "kegg", ridgeType = "ora", ridgeColor = "teal", gseaRankOpt="", sigLevel = 0.05, pwNum=20, inx = 1){
    dataSet <- readDataset(dataName);
    return(compute.ridgeline(dataSet, imgNm, dpi, format, fun.type, ridgeType, ridgeColor,gseaRankOpt, sigLevel, pwNum, inx));
}

PerformUpsetORA <- function(dataName="", file.nm, fun.type, IDs){
  paramSet <- readSet(paramSet, "paramSet");
  dataSet <- readDataset(dataName);
  gene.vec <- unlist(strsplit(IDs, "; "));


  # Check if phospho data and use appropriate functions
  is_phospho <- (!is.null(paramSet$data.type) && paramSet$data.type == "phospho")

  if (is_phospho) {
    # For phospho data, gene.vec contains phosphosite IDs (uniprot+site)
    res <- .performEnrichAnalysisPhospho(dataSet, file.nm, fun.type, gene.vec, "upset");
  } else {
    # Regular data - use standard enrichment
    sym.vec <- doEntrez2SymbolMapping(gene.vec, paramSet$data.org, paramSet$data.idType);
    if (sum(!is.na(sym.vec)) > 0) {
    }
    names(gene.vec) <- sym.vec;
    res <- .performEnrichAnalysis(dataSet, file.nm, fun.type, gene.vec, "upset");
  }

  return(res);
}

PerformGSEA<- function(dataName, file.nm, fun.type, netNm, mType, selectedFactorInx=1, mode = "multi",rankOpt=""){
    return(my.perform.gsea(dataName, file.nm, fun.type, netNm, mType, selectedFactorInx, mode,rankOpt));
}

ComputeRankedVec <- function(data, opt, inx = 1){
   return(my.compute.ranked.vec(data, opt, inx));
}

PlotGSView <-function(cmpdNm, format="png", dpi=96, width=NA){
   return(plot.gs.view(cmpdNm, format, dpi, width));

}

PlotGSViewNew <-function(cmpdNm, format="png", dpi=96, imgName){
   return(plot.gs.view(cmpdNm, format, dpi, NA, imgName));

}


.loadCustomEnrichLib <- function(fun.type, paramSet){
  
  # Determine folder name based on paramSet information
  if(paramSet$data.org == "generic"){
    folderNm <- paramSet$data.idType;
  }else{
    folderNm <- paramSet$data.org;
  }

  
  # Load the custom gene set library
  my.lib <- ov_qs_read("custom_lib.qs")

  # Extract the specific gene set based on the function type (e.g., cell line)
  current.featureset <- my.lib
  
  # Remove any empty pathways (or cell lines with no features)
  keep.inx <- lapply(current.featureset, length) > 0
  current.featureset <- current.featureset[keep.inx]

  # Get the names of the gene sets (e.g., cell line names)
  set.ids <- names(current.featureset)
  names(set.ids) <- names(current.featureset)

  # If the function type pertains to GO terms (or similar), normalize the names
  if(substr(fun.type, 0, 2) == "go") {
    names(current.featureset) <- firstup(names(current.featureset)) # Capitalize first letter
    names(current.featureset) <- gsub("-", "_", names(current.featureset)) # Replace hyphen with underscore
    names(set.ids) <- firstup(names(set.ids))
    names(set.ids) <- gsub("-", "_", names(set.ids))
  }

  # Save the processed gene set to a new file
  ov_qs_save(current.featureset, "current_featureset.qs")
  
  # Create the result object to return
  res <- list()
  res$current.setlink <- "" # Empty placeholder for set link
  res$current.setids <- set.ids # Names (IDs) of the gene sets (cell lines)
  res$current.featureset <- current.featureset # The actual gene set data

  return(res)
}

PerformDefaultEnrichment <- function(dataName="", file.nm, fun.type){
  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");
  #save.image("defaultenr.RData");
  anal.type <- paramSet$anal.type;
    if(anal.type=="onedata"){
      dataSet <- readDataset(dataName); #instead of function parameter
      if(nrow(dataSet$sig.mat) == 0){
        print("No DE features were identified, can not perform ORA analysis!");
        return(0);
      }

      gene.vec <- rownames(dataSet$sig.mat);
    }else if(anal.type=="metadata"){
      gene.vec <- rownames(analSet$meta.mat);
    }else{
      gene.vec <- rownames(paramSet$all.ent.mat);
    }

  # Check if phospho data and use appropriate functions
  is_phospho <- (!is.null(paramSet$data.type) && paramSet$data.type == "phospho")

  if (is_phospho) {
    # For phospho data, gene.vec contains phosphosite IDs (uniprot+site)
    res <- .performEnrichAnalysisPhospho(dataSet, file.nm, fun.type, gene.vec, "default");
  } else {
    # Regular data - use standard enrichment
    sym.vec <- doEntrez2SymbolMapping(gene.vec, paramSet$data.org, paramSet$data.idType);
    names(gene.vec) <- sym.vec;
    res <- .performEnrichAnalysis(dataSet, file.nm, fun.type, gene.vec, "default");
  }


  return(1);
}

GetSigSetCount <- function(enrType, type, pval=0.05){
  pval <- as.numeric(pval);
  imgSet <- readSet(imgSet, "imgSet");
    
  tbl <- imgSet$enrTables[[enrType]]$table
  count <- 0;
  if(type == "raw"){
   count<-sum(tbl$Pval<pval);
  }else{
    count<-sum(tbl$FDR<pval);
  }
  return(count);
}


GetSetIDLinks <- function(type=""){
  imgSet <- readSet(imgSet, "imgSet");
  fun.type <- imgSet$enrTables[[type]]$library;

  paramSet <- readSet(paramSet, "paramSet");
  ids <- imgSet$enrTables[[type]]$table$IDs
  pathways <- imgSet$enrTables[[type]]$table$Pathway
    if(fun.type %in% c("go_bp", "go_mf", "go_cc")){
        annots <- paste("<a href='https://www.ebi.ac.uk/QuickGO/term/", ids, "' target='_blank'>", pathways, "</a>", sep="");
    }else if(fun.type %in% c("go_panthbp", "go_panthmf", "go_panthcc")){
        annots <- paste("<a href='https://www.pantherdb.org/panther/categoryList.do?searchType=basic&fieldName=all&organism=all&fieldValue=", ids, "&listType=5' target='_blank'>", pathways, "</a>", sep="");
    }else if(fun.type == "kegg"){
        annots <- paste("<a href='https://www.genome.jp/dbget-bin/www_bget?pathway+", ids, "' target='_blank'>", pathways, "</a>", sep="");
    }else if(fun.type == "reactome"){
        annots <- paste("<a href='https://reactome.org/content/query?q=", ids, "' target='_blank'>", pathways, "</a>", sep="");
    }else{
        annots <- ids;
    }
  
  return(annots);
}

GetHTMLPathSet <- function(type, setNm){
  imgSet <- readSet(imgSet, "imgSet");
  current.featureset <- imgSet$enrTables[[type]]$current.featureset.symb;
  hits.query <- imgSet$enrTables[[type]]$hits.query;
  set <- current.featureset[[setNm]]; 
  
  #set <- cur.setids[[setNm]];
  
  hits <- hits.query
  
  # highlighting with different colors
  red.inx <- which(set %in% hits[[setNm]]);
  
  # use actual cmpd names
  #nms <- names(set);
  nms <- doEntrez2SymbolMapping(set);
  nms[red.inx] <- paste("<font color=\"red\">", "<b>", nms[red.inx], "</b>", "</font>",sep="");

  return(cbind(setNm, paste(unique(nms), collapse="; ")));
}


GetEnrResultMatrix <-function(type){
  imgSet <- readSet(imgSet, "imgSet");
  #print(names(imgSet$enrTables[[type]]));
  res <- imgSet$enrTables[[type]]$res.mat
  return(signif(as.matrix(res), 5));
}

GetEnrResultColNames<-function(type){
  imgSet <- readSet(imgSet, "imgSet");
  res <- imgSet$enrTables[[type]]$res.mat
  colnames(res);
}

GetEnrResSetIDs<-function(type){
  imgSet <- readSet(imgSet, "imgSet");
  res <- imgSet$enrTables[[type]]$table;
  return(res$IDs);
}

GetEnrResSetNames<-function(type){
  imgSet <- readSet(imgSet, "imgSet");
  res <- imgSet$enrTables[[type]]$table;
  return(res$Pathway);
}


GetGseaSetCount <- function(type, pval=0.05){
  pval <- as.numeric(pval);
  imgSet <- readSet(imgSet, "imgSet");
    
  tbl <- imgSet$enrTables[["gsea"]]$table
  count <- 0;
  if(type == "raw"){
   count<-sum(tbl$Pval<pval);
  }else{
    count<-sum(tbl$Padj<pval);
  }
  return(count);
}


GetGseaIDLinks <- function(dataName=""){
  dataSet <- readDataset(dataName);
  imgSet <- readSet(imgSet, "imgSet");
  fun.type <- imgSet$enrTables[["gsea"]]$library;

  paramSet <- readSet(paramSet, "paramSet");
  ids <- imgSet$enrTables[["gsea"]]$table$IDs
  pathways <- imgSet$enrTables[["gsea"]]$table$Name

    if(fun.type %in% c("go_bp", "go_mf", "go_cc")){
        annots <- paste("<a href='https://www.ebi.ac.uk/QuickGO/term/", ids, "' target='_blank'>", pathways, "</a>", sep="");
    }else if(fun.type %in% c("go_panthbp", "go_panthmf", "go_panthcc")){
        annots <- paste("<a href='https://www.pantherdb.org/panther/categoryList.do?searchType=basic&fieldName=all&organism=all&fieldValue=", ids, "&listType=5' target='_blank'>", pathways, "</a>", sep="");
    }else if(fun.type == "kegg"){
        annots <- paste("<a href='https://www.genome.jp/dbget-bin/www_bget?pathway+", ids, "' target='_blank'>", pathways, "</a>", sep="");
    }else if(fun.type == "reactome"){
        annots <- paste("<a href='https://reactome.org/content/query?q=", ids, "' target='_blank'>", pathways, "</a>", sep="");
    }else{
        annots <- ids;
    }
  
  return(annots);
}

GetGseaHTMLPathSet <- function(setNm){
  imgSet <- readSet(imgSet, "imgSet");
  current.featureset <- imgSet$enrTables[["gsea"]]$current.featureset.symb;
  hits.query <- imgSet$enrTables[["gsea"]]$hits.query;
  set <- current.featureset[[setNm]]; 
    
  hits <- hits.query
  
  # highlighting with different colors
  red.inx <- which(set %in% hits[[setNm]]);
  
  # use actual cmpd names
  #nms <- names(set);
  nms <- set;
  nms[red.inx] <- paste("<font color=\"red\">", "<b>", nms[red.inx], "</b>", "</font>",sep="");

  return(cbind(setNm, paste(unique(nms), collapse="; ")));
}


GetGseaResultMatrix <-function(){
  imgSet <- readSet(imgSet, "imgSet");
  res <- imgSet$enrTables[["gsea"]]$res.mat

  return(signif(as.matrix(res), 5));
}

GetGseaResultColNames <-function(){
  imgSet <- readSet(imgSet, "imgSet");
  res <- imgSet$enrTables[["gsea"]]$res.mat
  colnames(res);
}

GetGseaResSetIDs <-function(){
  imgSet <- readSet(imgSet, "imgSet");
  res <- imgSet$enrTables[["gsea"]]$table;

  return(res$IDs);
}

GetGseaResSetNames<-function(){
  imgSet <- readSet(imgSet, "imgSet");
  res <- imgSet$enrTables[["gsea"]]$table;
  return(res$Name);
}

GetEnrResTypes<-function(){
  imgSet <- readSet(imgSet, "imgSet");
  nms <- names(imgSet$enrTables);
  nms <- setdiff(nms, "default");
  nms <- setdiff(nms, "gsea")
  return(nms);
}

GetEnrResLibrary<-function(type){
  imgSet <- readSet(imgSet, "imgSet");
  summary <- imgSet$enrTables[[type]];
  res <- summary$library
  return(res);
}


PerformNetEnrichment <- function(dataName="", file.nm, fun.type, IDs){
  #dataName <<- dataName;
  #file.nm <<- file.nm;
  #fun.type <<- fun.type;
  #IDs <<- IDs;
  #save.image("PerformNetEnrichment.RData");
  dataSet <- readDataset(dataName);
  paramSet <- readSet(paramSet, "paramSet");
  data.org <- paramSet$data.org;
  # prepare query
  ora.vec <- NULL;
  idtype <- "entrez";
  
    ora.vec <- unlist(strsplit(IDs, "; "));
    names(ora.vec) <- as.character(ora.vec);


  if(fun.type %in% c("trrust", "encode", "jaspar", "mirnet", "met", "drugbank", "disease")){
    res <- PerformRegEnrichAnalysis(dataSet, file.nm, fun.type, ora.vec, "inverse");
  }else{
    # Check if phospho data and use appropriate functions
    is_phospho <- (!is.null(paramSet$data.type) && paramSet$data.type == "phospho")

    if (is_phospho) {
      # For phospho data, ora.vec contains phosphosite IDs (uniprot+site)
      res <- .performEnrichAnalysisPhospho(dataSet, file.nm, fun.type, ora.vec, "coexp");
    } else {
      # Regular data - use standard enrichment
      res <- .performEnrichAnalysis(dataSet, file.nm, fun.type, ora.vec, "coexp");
    }
  }

  return(res);
}

PerformRegEnrichAnalysis <- function(dataSet, file.nm, fun.type, ora.vec, netInv){
    return(my.reg.enrich(dataSet, file.nm, fun.type, ora.vec, netInv));
}


FindCommunities <- function(method = "walktrap",
                            use.weight = FALSE,
                            component = "largest",#c("largest", "all", "seed"),
                            min_features = 5,
                            require_query_hit = TRUE) {
  # save.image("find.Rdata");
  paramSet    <- readSet(paramSet, "paramSet")
  analSet    <- readSet(analSet, "analSet")

  seed.expr   <- paramSet$seed.expr
  ppi.comps <- analSet$ppi.comps
  current.net <- ppi.comps[[current.net.nm]]
  
  if (igraph::vcount(current.net) < 2L) return("NA||Graph too small")
  
  # ---- choose component(s) ----------------------------------------------------
  pick_largest <- function(g) {
    comp <- igraph::components(g)
    igraph::induced_subgraph(g, which(comp$membership == which.max(comp$csize)))
  }
  
  graph_list <- switch(
    component,
    "largest" = list(pick_largest(current.net)),
    "all" = {
      comp <- igraph::components(current.net)
      split(seq_len(igraph::vcount(current.net)), comp$membership) |>
        lapply(function(vs) igraph::induced_subgraph(current.net, vs))
    },
    "seed" = {
      seeds <- intersect(seed.proteins, igraph::V(current.net)$name)
      if (length(seeds) == 0) return("NA||No seeds in network!")
      comp <- igraph::components(current.net)
      # keep components that contain at least one seed
      keep_comp <- unique(comp$membership[match(seeds, igraph::V(current.net)$name)])
      vs <- which(comp$membership %in% keep_comp)
      list(igraph::induced_subgraph(current.net, vs))
    }
  )
  
  # ---- helper: run detection on one graph ------------------------------------
  run_one <- function(g) {
    if (!igraph::is_connected(g)) g <- pick_largest(g)
    if (igraph::vcount(g) < 2L) return(list(vec = character(0), tbl = NULL))
    
    # -- FIX 1: ensure vertex names exist; fallback to other attrs/indices -----
    vnames <- igraph::V(g)$name
    if (is.null(vnames)) vnames <- rep(NA_character_, igraph::vcount(g))
    if (all(is.na(vnames))) {
      alt <- NULL
      if (!is.null(igraph::V(g)$Id))    alt <- as.character(igraph::V(g)$Id)
      if (!is.null(igraph::V(g)$Label)) alt <- as.character(igraph::V(g)$Label)
      vnames <- if (!is.null(alt)) alt else as.character(seq_len(igraph::vcount(g)))
      igraph::V(g)$name <- vnames
    }
    
    # symbol mapping with fallback to name
    hit.x <- match(vnames, ppi.net[["node.data"]][, 1])
    sybls <- ppi.net[["node.data"]][hit.x, 2]
    sybls[is.na(sybls)] <- vnames[is.na(sybls)]
    names(sybls) <- vnames
    
    # -- FIX 2: robust weight assignment (no invalid indexing) -----------------
    weights <- igraph::E(g)$weight
    # community detection (weights used when supported)
    fc <- switch(
      method,
      "walktrap"  = igraph::cluster_walktrap(g, weights = weights),
      "infomap"   = igraph::cluster_infomap(g, e.weights = weights),
      "labelprop" = igraph::cluster_label_prop(g, weights = weights),
      { return(list(err = "NA||Unknown method!")) }
    )
    
    comm_list <- igraph::communities(fc)
    if (length(comm_list) == 0 || igraph::modularity(fc, weights = weights) == 0) {
      return(list(vec = character(0), tbl = NULL))
    }
    
    # iterate communities
    # OPTIMIZED: Use list to avoid rbind in loop
    community.vec  <- character(0)
    gene_community_list <- vector("list", length(comm_list))
    rowcount <- 0L

    for (i in seq_along(comm_list)) {
      vids       <- comm_list[[i]]     # vertex indices
      comm_size  <- length(vids)
      if (comm_size < min_features) next

      comm_names  <- vnames[vids]
      # ensure no NA names in printable path
      comm_names[is.na(comm_names)] <- as.character(vids[is.na(comm_names)])

      # label strings (symbols if available)
      comm_labels <- sybls[comm_names]
      comm_labels[is.na(comm_labels)] <- comm_names

      # query hits
      qnums <- comm_size

      if (require_query_hit && qnums == 0) next

      # in/out degree test
      subg    <- igraph::induced_subgraph(g, vids)
      in.deg  <- igraph::degree(subg)
      out.deg <- igraph::degree(g, vids) - in.deg
      ppval   <- suppressWarnings(
        wilcox.test(in.deg, out.deg, exact = FALSE, alternative = "two.sided")$p.value
      )
      ppval   <- signif(ppval, 3)

      # record
      rowcount <- rowcount + 1L
      pids     <- paste(comm_labels, collapse = "->")

      com.mat <- cbind(Id     = comm_names,
                       Label  = comm_labels,
                       Module = as.character(i))  # keep as before
      gene_community_list[[rowcount]] <- com.mat

      community.vec[rowcount] <- paste(c(comm_size, qnums, ppval, pids), collapse = ";")
    }

    # OPTIMIZED: Combine once at the end
    gene.community <- if(rowcount > 0) do.call(rbind, gene_community_list[1:rowcount]) else NULL
    
    if (length(community.vec) == 0) return(list(vec = character(0), tbl = NULL))
    
    # order: size desc, p-value asc (same as before)
    community_data <- do.call(rbind, lapply(community.vec, function(x) {
      parts <- strsplit(x, ";")[[1]]
      data.frame(size = as.numeric(parts[1]), p_value = as.numeric(parts[3]))
    }))
    ord <- with(community_data, order(-size, p_value))
    
    list(vec = community.vec[ord], tbl = gene.community)
  }
  # run and combine
  out_list <- lapply(graph_list, run_one)
  if (any(vapply(out_list, function(x) !is.null(x$err), logical(1)))) { 
    return(out_list[[which(vapply(out_list, function(x) !is.null(x$err), logical(1)))] ]$err)
  }
  
  vecs <- unlist(lapply(out_list, `[[`, "vec"), use.names = FALSE)
  tbls <- do.call(rbind, lapply(out_list, `[[`, "tbl"))
  
  if (length(vecs) == 0) return("NA||No communities were detected!")
  
  all.communities <- paste(vecs, collapse = "||")
  if (!is.null(tbls) && nrow(tbls) > 0) {
    colnames(tbls) <- c("Id", "Label", "Module")
    fast.write(tbls, file = "module_table.csv", row.names = FALSE)
    return(all.communities)
  } else {
    return("NA")
  }
}


community.significance.test <- function(graph, vs, ...) {
  subgraph <- induced_subgraph(graph, vs)
  in.degrees <- degree(subgraph)
  out.degrees <- degree(graph, vs) - in.degrees
  wilcox.test(in.degrees, out.degrees, ...)
}


# from to should be valid nodeIDs
GetShortestPaths <- function(from, to){
  current.net <- ppi.comps[[current.net.nm]];
  paths <- all_shortest_paths(current.net, from, to)$res;
  if(length(paths) == 0){
    return (paste("No connection between the two nodes!"));
  }
  
  path.vec <- vector(mode="character", length=length(paths));
  for(i in 1:length(paths)){
    path.inx <- paths[[i]]; 
    path.ids <- V(current.net)$name[path.inx];
    path.sybls <- path.ids;
    pids <- paste(path.ids, collapse="->");
    psbls <- paste(path.sybls, collapse="->");
    path.vec[i] <- paste(c(pids, psbls), collapse=";")
  }
  
  if(length(path.vec) > 50){
    path.vec <- path.vec[1:50];
  }
  
  all.paths <- paste(path.vec, collapse="||");
  return(all.paths);
}

DecomposeGraph <- function(gObj,analSet, minNodeNum = 3, jsonBool = F){
  # now decompose to individual connected subnetworks
    if(jsonBool == "netjson"){
        comps <-list(gObj)
    }else{
        comps <-decompose(gObj, min.vertices=minNodeNum);
    }
  if(length(comps) == 0){
    current.msg <<- paste("No subnetwork was identified with at least", minNodeNum, "nodes!");
    return(NULL);
  }
  
  # first compute subnet stats
  net.stats <- ComputeSubnetStats(comps);
  ord.inx <- order(net.stats[,1], decreasing=TRUE);
  net.stats <- net.stats[ord.inx,];
  comps <- comps[ord.inx];
  names(comps) <- rownames(net.stats) <- paste("subnetwork", 1:length(comps), sep="");
  
  # note, we report stats for all nets (at least 3 nodes);
  hit.inx <- net.stats$Node >= minNodeNum;
  comps <- comps[hit.inx];
  #overall <- list();
  #overall[["overall"]] <- g
  #ppi.comps <- append(overall, ppi.comps, after=1);
  
  # now record
  net.stats <<- net.stats;
  sub.stats <- unlist(lapply(comps, vcount)); 
  analSet$ppi.comps <- comps;
  analSet$net.stats <- net.stats;
  analSet$substats <- sub.stats;
  return(analSet);
}

regEnrichment <- function(dataName, file.nm, fun.type, IDs, netInv){

  regBool <<- "false";
  res <- PerformNetEnrichment(dataName, file.nm, fun.type, IDs);
}

# ── Server-side PNG: Enrichment Network (igraph) ──
PlotEnrichNetworkPNG <- function(dataName, imgName, format="png", dpi=150, width=NA) {
  tryCatch({
    require(igraph); require(reshape)
    cat("[PlotEnrichNetworkPNG] wd=", getwd(), "\n")
    cat("[PlotEnrichNetworkPNG] enr.mat.qs exists:", file.exists("enr.mat.qs"), "\n")
    qs_files <- list.files(pattern = "\\.qs$")
    cat("[PlotEnrichNetworkPNG] qs files in wd:", paste(qs_files, collapse=", "), "\n")
    enr.mat <- qs::qread("enr.mat.qs")
    hits.query <- qs::qread("hits_query.qs")
    if (is.null(enr.mat) || nrow(enr.mat) == 0) return(0)
    if ("FDR" %in% colnames(enr.mat)) {
      ord.inx <- order(enr.mat[, "FDR"])
      enr.mat <- enr.mat[ord.inx[1:min(20, nrow(enr.mat))], , drop = FALSE]
    }
    hits <- as.numeric(enr.mat[, "Hits"]); id <- rownames(enr.mat)
    names(hits) <- id; hits.query <- hits.query[id]
    n <- nrow(enr.mat); w <- matrix(0, n, n); colnames(w) <- rownames(w) <- id
    for (i in 1:n) for (j in i:n) w[i,j] <- overlap_ratio(hits.query[id[i]], hits.query[id[j]], "mixed")
    w[lower.tri(w)] <- t(w)[lower.tri(w)]
    wd <- reshape::melt(w); wd <- wd[wd[,1]!=wd[,2] & !is.na(wd[,3]),]
    g <- graph_from_data_frame(wd[,-3], directed=FALSE)
    g <- delete_edges(g, E(g)[wd[,3] < 0.3])
    if (vcount(g) == 0) return(0)
    cnt2 <- hits[V(g)$name]
    V(g)$size <- if (all(cnt2==cnt2[1])) rep(12, length(cnt2)) else scales::rescale(log(cnt2+1,10), to=c(8,24))
    V(g)$color <- "orange"; V(g)$frame.color <- "white"
    V(g)$label <- V(g)$name; V(g)$label.cex <- 0.6; V(g)$label.color <- "black"
    E(g)$arrow.mode <- 0; E(g)$color <- "#cccccc"
    l <- layout_with_graphopt(g)
    imgPath <- paste0(imgName, ".", format)
    w.val <- if (is.na(width)) 8 else width/dpi
    png(imgPath, width=w.val, height=w.val*0.75, units="in", res=dpi)
    par(mar=c(1,1,2,1)); plot(g, layout=l, main="Enrichment Network (KEGG)"); dev.off()
    return(1)
  }, error = function(e) { message("PlotEnrichNetworkPNG error: ", e$message); return(0) })
}

# ── Server-side PNG: Gene-Pathway Enrichment Heatmap (ProteoAnalyst version) ──
PlotEnrichHeatmapPNG <- function(dataName, imgName, format="png", dpi=150, width=NA) {
  tryCatch({
    enr.mat <- qs::qread("enr.mat.qs")
    current.geneset <- if (file.exists("current_featureset.qs")) qs::qread("current_featureset.qs") else NULL
    if (is.null(enr.mat) || nrow(enr.mat) < 2 || is.null(current.geneset)) return(0)

    # PA stores UniProt IDs in prot.mat; use analSet$uniprot_to_entrez_map for Entrez IDs
    user.entrez <- NULL
    analSet <- tryCatch(readSet(analSet, "analSet"), error = function(e) NULL)
    if (!is.null(analSet) && !is.null(analSet$uniprot_to_entrez_map)) {
      umap <- analSet$uniprot_to_entrez_map
      user.entrez <- unique(unname(umap[!is.na(umap)]))
    }
    if (is.null(user.entrez) || length(user.entrez) == 0) {
      dataSet <- readDataset(dataName)
      if (!is.null(dataSet) && !is.null(dataSet$prot.mat)) user.entrez <- rownames(dataSet$prot.mat)
    }
    if (is.null(user.entrez) || length(user.entrez) == 0) return(0)

    if ("FDR" %in% colnames(enr.mat)) {
      enr.mat <- enr.mat[order(as.numeric(enr.mat[,"FDR"]))[1:min(15,nrow(enr.mat))], , drop=FALSE]
    }
    matched.pw <- intersect(rownames(enr.mat), names(current.geneset))
    if (length(matched.pw) < 2) return(0)
    hit.genes.per.pw <- lapply(matched.pw, function(pw) intersect(user.entrez, current.geneset[[pw]]))
    names(hit.genes.per.pw) <- matched.pw
    all.hit.genes <- unique(unlist(hit.genes.per.pw))
    if (length(all.hit.genes) < 2) return(0)
    gp.mat <- matrix(0, length(all.hit.genes), length(matched.pw))
    rownames(gp.mat) <- all.hit.genes; colnames(gp.mat) <- matched.pw
    for (pw in matched.pw) { g <- hit.genes.per.pw[[pw]]; if (length(g)>0) gp.mat[g,pw] <- 1 }
    tryCatch({
      ps <- readSet(paramSet, "paramSet")
      syms <- doEntrez2SymbolMapping(rownames(gp.mat), ps$data.org, ps$data.idType)
      if (!is.null(syms) && length(syms)==nrow(gp.mat)) rownames(gp.mat) <- make.unique(syms)
    }, error=function(e) {})
    if (nrow(gp.mat) > 30) gp.mat <- gp.mat[order(rowSums(gp.mat), decreasing=TRUE)[1:30], , drop=FALSE]
    gp.mat <- gp.mat[rowSums(gp.mat)>0, colSums(gp.mat)>0, drop=FALSE]
    if (nrow(gp.mat)<2 || ncol(gp.mat)<2) return(0)
    colnames(gp.mat) <- make.unique(substr(colnames(gp.mat),1,40))
    imgPath <- paste0(imgName, ".", format)
    w.val <- if (is.na(width)) max(7, ncol(gp.mat)*0.6+3) else width/dpi
    h.val <- max(5, nrow(gp.mat)*0.25+2)
    png(imgPath, width=w.val, height=h.val, units="in", res=dpi, bg="white")
    par(mar=c(1, 8, max(4, max(nchar(colnames(gp.mat)))*0.3), 1))
    nr <- nrow(gp.mat); nc <- ncol(gp.mat)
    plot(NA, xlim=c(0,nc), ylim=c(0,nr), xaxt="n", yaxt="n", xlab="", ylab="", bty="n", asp=NA)
    for (i in 1:nr) for (j in 1:nc) if (gp.mat[i,j]==1) rect(j-1,nr-i,j-0.05,nr-i+0.95, col="#F4837D", border="white", lwd=0.5)
    for (i in 0:nr) abline(h=i, col="#e0e0e0", lwd=0.3); for (j in 0:nc) abline(v=j, col="#e0e0e0", lwd=0.3)
    mtext(rownames(gp.mat), side=2, at=nr:1-0.5, las=1, line=0.5, cex=0.7, adj=1)
    text(x=1:nc-0.5, y=nr+0.3, labels=colnames(gp.mat), srt=45, adj=0, xpd=TRUE, cex=0.65)
    mtext("Input Genes", side=2, line=6.5, cex=0.9, font=2)
    mtext("Enriched Terms", side=3, line=max(2, max(nchar(colnames(gp.mat)))*0.15), cex=0.9, font=2)
    dev.off(); return(1)
  }, error = function(e) { message("PlotEnrichHeatmapPNG error: ", e$message); tryCatch(dev.off(), error=function(x){}); return(0) })
}
