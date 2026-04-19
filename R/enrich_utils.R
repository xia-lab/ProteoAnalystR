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
    # msg("[EnrichAnalysis] ora.vec has no names, using values as names")
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

    # msg("[EnrichAnalysis] Found ", sum(!is.na(uniprot.vec)), " UniProt mappings for Entrez IDs")
    # msg("[EnrichAnalysis] Final query size: ", length(ora.vec))
  } else {
    # IDs are UniProt, convert to Entrez
    msg("[EnrichAnalysis] Converting UniProt IDs to Entrez IDs for enrichment matching...")
    msg("[EnrichAnalysis] Input query size: ", length(ora.vec))

    # Normalize UniProt IDs (remove phosphosite annotations, isoforms, etc.)
    # e.g., Q9D1F4_T_247 → Q9D1F4, P12345-1 → P12345
    normalized.ora.vec <- ora.vec
    normalized.ora.vec <- sub("_[A-Z]_\\d+$", "", normalized.ora.vec)  # Remove phosphosites
    normalized.ora.vec <- sub("-\\d+$", "", normalized.ora.vec)         # Remove isoforms
    normalized.ora.vec <- trimws(normalized.ora.vec)

    n_modified <- sum(normalized.ora.vec != ora.vec, na.rm = TRUE)
    if (!is.na(n_modified) && n_modified > 0) {
      msg("[EnrichAnalysis] Normalized ", n_modified, " IDs with phosphosite/isoform annotations")
      modified_idx <- which(normalized.ora.vec != ora.vec & !is.na(normalized.ora.vec) & !is.na(ora.vec))
      if (length(modified_idx) > 0) {
        isoform_examples <- ora.vec[modified_idx]
        normalized_examples <- normalized.ora.vec[modified_idx]
        msg("[EnrichAnalysis] Sample original IDs: ", paste(head(isoform_examples, 5), collapse = ", "))
        msg("[EnrichAnalysis] Sample normalized IDs: ", paste(head(normalized_examples, 5), collapse = ", "))
      }
    } else {
      msg("[EnrichAnalysis] No isoform/phosphosite annotations detected")
    }

    # Query UniProt → Entrez mapping
    msg("[EnrichAnalysis] Querying UniProt→Entrez database for organism: ", paramSet$data.org)
    uniprot.map <- queryGeneDB("entrez_uniprot", paramSet$data.org)
    msg("[EnrichAnalysis] Database contains ", nrow(uniprot.map), " UniProt→Entrez mappings")
    hit.inx <- match(normalized.ora.vec, uniprot.map[, "accession"])
    entrez.vec <- uniprot.map[hit.inx, "gene_id"]

    # Get symbols for display names
    sym.vec <- doEntrez2SymbolMapping(entrez.vec, paramSet$data.org, "entrez");

    # Count conversions
    na.inx <- is.na(entrez.vec)
    msg("[EnrichAnalysis] Mapped ", sum(!na.inx), " UniProt IDs to Entrez IDs (",
        round(100*sum(!na.inx)/length(ora.vec), 1), "%)")
    msg("[EnrichAnalysis] Failed to map ", sum(na.inx), " UniProt IDs")
    if (sum(na.inx) > 0) {
      failed_ids <- ora.vec[na.inx]
      msg("[EnrichAnalysis] Sample failed IDs: ", paste(head(failed_ids, 10), collapse = ", "))
      normalized_failed <- normalized.ora.vec[na.inx]
      msg("[EnrichAnalysis] Normalized failed IDs: ", paste(head(normalized_failed, 10), collapse = ", "))
    }
    if (sum(!na.inx) > 0) {
      msg("[EnrichAnalysis] Sample Entrez IDs: ", paste(head(entrez.vec[!na.inx], 10), collapse = ", "))
      msg("[EnrichAnalysis] Sample gene symbols: ", paste(head(sym.vec[!na.inx], 10), collapse = ", "))
    }

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

    # msg("[EnrichAnalysis] Final query size after removing unmapped: ", length(ora.vec))
  }

  # Store mapping in analSet for downstream use (e.g., network analysis)
  analSet <- readSet(analSet, "analSet");
  analSet$uniprot_to_entrez_map <- uniprot.to.entrez.map
  saveSet(analSet, "analSet");
  # msg("[EnrichAnalysis] Saved UniProt→Entrez mapping to analSet")
  
  # msg("[EnrichAnalysis] Setting up universe (universe.opt: ", paramSet$universe.opt, ")")
  if(paramSet$universe.opt == "library"){
    current.universe <- unique(unlist(current.featureset));
    # msg("[EnrichAnalysis] Using library universe")
  }else{
    # cut to the universe to uploaded features
    # msg("[EnrichAnalysis] Using uploaded data universe (anal.type: ", paramSet$anal.type, ")")
    if(paramSet$anal.type == "onedata"){
      data.anot <- .get.annotated.data();
      current.universe <- rownames(data.anot);
      # msg("[EnrichAnalysis] Got universe from onedata")
    }else if(paramSet$anal.type == "metadata"){
      inmex <- ov_qs_read("inmex_meta.qs");
      current.universe <- rownames(inmex$data);
      # msg("[EnrichAnalysis] Got universe from metadata")
    }else{
      if(!is.null(paramSet$backgroundUniverse)){
        current.universe <- paramSet$backgroundUniverse;
        # msg("[EnrichAnalysis] Got universe from paramSet$backgroundUniverse")
      }else{
        current.universe <- unique(unlist(current.featureset));
        # msg("[EnrichAnalysis] Got universe from feature set (fallback)")
      }
    }
  }
  # msg("[EnrichAnalysis] Universe size: ", length(current.universe))
  
  # also make sure pathways only contain features measured in experiment
  #if(!is.null(dataSet$data.anot)){
   if(file.exists("data.anot.qs")){
    # msg("[EnrichAnalysis] Filtering pathways to universe (data.anot.qs exists)")
    current.featureset <- lapply(current.featureset, function(x){x[x %in% current.universe]})
    inds <- lapply(current.featureset, length) > 0
    current.featureset <- current.featureset[inds]
    # msg("[EnrichAnalysis] After filtering: ", length(current.featureset), " pathways remain")
  }

  # prepare for the result table
  set.size<-length(current.featureset);
  res.mat<-matrix(0, nrow=set.size, ncol=5);
  rownames(res.mat)<-names(current.featureset);
  colnames(res.mat)<-c("Total", "Expected", "Hits", "Pval", "FDR");

  q.size<-length(ora.vec);
  # msg("[EnrichAnalysis] Query size (q.size): ", q.size)

  # get the matched query for each pathway
  # msg("[EnrichAnalysis] Calculating hits for each pathway...")
  hits.query <- lapply(current.featureset,
                       function(x) {
                         ora.nms[ora.vec%in%unlist(x)];
                       }
  );

  # PROTEOMICS: Convert hits back to UniProt IDs (reverse mapping)
  if(!is.null(uniprot.to.entrez.map)){
    # msg("[EnrichAnalysis] Converting hit lists back to UniProt IDs...")

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

    # msg("[EnrichAnalysis] Converted hit lists to UniProt IDs")
  }

  ov_qs_save(hits.query, "hits_query.qs");

  names(hits.query) <- names(current.featureset);
  hit.num<-unlist(lapply(hits.query, function(x){length(unique(x))}), use.names=FALSE);
  # msg("[EnrichAnalysis] Total pathways with hits: ", sum(hit.num > 0), " out of ", length(hit.num))
  # msg("[EnrichAnalysis] Max hits in a pathway: ", max(hit.num))
  # msg("[EnrichAnalysis] Min hits in a pathway: ", min(hit.num))
  
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
  # msg("[EnrichAnalysis] Filtering results - keeping only pathways with hits...")
  res.mat <- res.mat[hit.num>0,,drop = F];
  hits.query <- hits.query[hit.num>0];
  # msg("[EnrichAnalysis] After filtering: ", nrow(res.mat), " pathways with hits remain")

  if(nrow(res.mat)> 1){
    # msg("[EnrichAnalysis] Sorting results by p-value...")
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
    # msg("[EnrichAnalysis] Saved full results to: ", csv.nm)

    imp.inx <- res.mat[,4] <= 0.05;
    imp.inx[is.na(imp.inx)] <- F
    # msg("[EnrichAnalysis] Pathways with p <= 0.05: ", sum(imp.inx))

    if(sum(imp.inx) < 10){ # too little left, give the top ones
      topn <- ifelse(nrow(res.mat) > 10, 10, nrow(res.mat));
      # msg("[EnrichAnalysis] Too few significant results (", sum(imp.inx), "), taking top ", topn, " pathways")
      res.mat <- res.mat[1:topn,];
      hits.query <- hits.query[1:topn];
    }else{
      # msg("[EnrichAnalysis] Filtering to significant pathways (p <= 0.05)")
      res.mat <- res.mat[imp.inx,];
      hits.query <- hits.query[imp.inx];
      if(sum(imp.inx) > 120){
        # msg("[EnrichAnalysis] Too many results (", sum(imp.inx), "), limiting to top 120")
        # now, clean up result, synchronize with hit.query
        res.mat <- res.mat[1:120,];
        hits.query <- hits.query[1:120];
      }
    }
  }else{
    # msg("[EnrichAnalysis] ERROR: Only ", nrow(res.mat), " pathway(s) with hits - cannot proceed!")
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
    # msg("[EnrichAnalysis] Successfully saved enr.mat.qs (", file.size("enr.mat.qs"), " bytes)")
  } else {
    warning("[EnrichAnalysis] WARNING: enr.mat.qs was not created!")
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
    stop("[EnrichLib] Library file not found: ", my.path)
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

GetRidgePlot <- function(dataName, imgNm = "abc", dpi=96, format="png", fun.type = "kegg", ridgeType = "ora", ridgeColor = "teal", gseaRankOpt="", sigLevel = 0.05, pwNum=20, inx = 1){
    dataSet <- readDataset(dataName);
    if(!exists("compute.ridgeline")){ # public web on same user dir
        rc.path <- paste0(resource.dir, "rscripts/ProteoAnalystR/R/utils_ridgeline.Rc")
        r.path <- paste0(resource.dir, "rscripts/ProteoAnalystR/R/utils_ridgeline.R")
        if (file.exists(rc.path)) {
          ok <- tryCatch({
            compiler::loadcmp(rc.path)
            TRUE
          }, error = function(e) {
            FALSE
          })
          if (!ok) {
            source(r.path, local = .GlobalEnv)
          }
        } else {
          source(r.path, local = .GlobalEnv)
        }
    }
    return(compute.ridgeline(dataSet, imgNm, dpi, format, fun.type, ridgeType, ridgeColor,gseaRankOpt, sigLevel, pwNum, inx));
}

PerformUpsetORA <- function(dataName="", file.nm, fun.type, IDs){
  paramSet <- readSet(paramSet, "paramSet");
  dataSet <- readDataset(dataName);
  gene.vec <- unlist(strsplit(IDs, "; "));

  msg("[UpsetORA] Starting upset enrichment analysis")
  msg("[UpsetORA] Input IDs string length: ", nchar(IDs))
  msg("[UpsetORA] Parsed gene.vec length: ", length(gene.vec))
  msg("[UpsetORA] Sample input IDs: ", paste(head(gene.vec, 10), collapse = ", "))
  msg("[UpsetORA] Organism: ", paramSet$data.org)
  msg("[UpsetORA] ID type: ", paramSet$data.idType)
  msg("[UpsetORA] Function type: ", fun.type)

  # Check if phospho data and use appropriate functions
  is_phospho <- (!is.null(paramSet$data.type) && paramSet$data.type == "phospho")
  msg("[UpsetORA] Is phospho data: ", is_phospho)

  if (is_phospho) {
    # For phospho data, gene.vec contains phosphosite IDs (uniprot+site)
    msg("[UpsetORA] Using phospho enrichment analysis pathway")
    res <- .performEnrichAnalysisPhospho(dataSet, file.nm, fun.type, gene.vec, "upset");
  } else {
    # Regular data - use standard enrichment
    msg("[UpsetORA] Using standard enrichment analysis pathway")
    msg("[UpsetORA] Converting IDs to symbols for naming...")
    sym.vec <- doEntrez2SymbolMapping(gene.vec, paramSet$data.org, paramSet$data.idType);
    msg("[UpsetORA] Symbol conversion: ", sum(!is.na(sym.vec)), " successful, ", sum(is.na(sym.vec)), " failed")
    if (sum(!is.na(sym.vec)) > 0) {
      msg("[UpsetORA] Sample symbols: ", paste(head(sym.vec[!is.na(sym.vec)], 10), collapse = ", "))
    }
    names(gene.vec) <- sym.vec;
    msg("[UpsetORA] Calling .performEnrichAnalysis with ", length(gene.vec), " genes")
    res <- .performEnrichAnalysis(dataSet, file.nm, fun.type, gene.vec, "upset");
  }

  msg("[UpsetORA] Enrichment analysis completed, result: ", res)
  return(res);
}

PerformGSEA<- function(dataName, file.nm, fun.type, netNm, mType, selectedFactorInx=1, mode = "multi",rankOpt=""){
    if(!exists("my.perform.gsea")){ 
        compiler::loadcmp(paste0(resource.dir, "rscripts/ProteoAnalystR/R/utils_gsea.Rc"));    
    }
    return(my.perform.gsea(dataName, file.nm, fun.type, netNm, mType, selectedFactorInx, mode,rankOpt));
}

ComputeRankedVec <- function(data, opt, inx = 1){
   if(!exists("my.compute.ranked.vec")){ 

        compiler::loadcmp(paste0(resource.dir, "rscripts/ProteoAnalystR/R/utils_gsea.Rc"));    
        
    }
   return(my.compute.ranked.vec(data, opt, inx));
}

PlotGSView <-function(cmpdNm, format="png", dpi=96, width=NA){
   if(!exists("plot.gs.view")){ 
        compiler::loadcmp(paste0(resource.dir, "rscripts/ProteoAnalystR/R/utils_gsea.Rc"));    
   }
   return(plot.gs.view(cmpdNm, format, dpi, width));

}

PlotGSViewNew <-function(cmpdNm, format="png", dpi=96, imgName){
   if(!exists("plot.gs.view")){ 
        compiler::loadcmp(paste0(resource.dir, "rscripts/ProteoAnalystR/R/utils_gsea.Rc"));    
   }
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
    if(!exists("my.reg.enrich")){ # public web on same user dir
        compiler::loadcmp(paste0(resource.dir, "rscripts/ProteoAnalystR/R/_utils_regenrich.Rc"));    
    }
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
