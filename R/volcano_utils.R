  ##################################################
## R scripts for ProteoAnalyst 
## Functions related to volcano plot
## Author: Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

#'Prepare data for volcano plot visualization
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: MIT
#'@export
#'
Volcano.Anal <- function(dataName="", fileNm="name", paired=FALSE, fcthresh=0, threshp=0.05, analType="NA", inx=1, dpi=96, format="png"){
  #save.image('volc.RData');

  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");

  anal.type <- paramSet$anal.type;

  inx <- 1
  #print("Prepare volcano anal");
  limit_fc <- T; #whether limit to -10 to 10 fc
  if(anal.type == "metadata"){
    if(paramSet$selDataNm == "meta_default"){
       if(is.null(paramSet$fc.thresh)){
         paramSet$fc.thresh <- 0; 
       }

      data <- ov_qs_read("allMeta.mat.qs")    
      p.value <- data[, 2]
      data <- cbind(unname(analSet$meta.avgFC[rownames(data)]), data);
      fcthresh <- paramSet$fc.thresh;
      threshp <- paramSet$BHth;
      inx <- 1;
    }else{
      dataSet <- readDataset(paramSet$selDataNm);
      
      data <- as.matrix(analSet$inmex.ind[paramSet$selDataNm][[1]])
      p.value <- data[, "Pval"]
    }
    limit_fc <-F;
  }else{
    dataSet <- readDataset(dataName);
    data <- as.matrix(dataSet$comp.res);
    if(is.null(paramSet$use.fdr) || paramSet$use.fdr){
        p.value <- data[, "adj.P.Val"];
    }else{
        p.value <- data[, "P.Value"];
    }
  }

  paramSet$fcthreshu <- fcthresh
  inx.p <- p.value <= threshp;
  
  zero.inx <- unname(p.value) == 0;
  p.value[zero.inx] <- min(p.value[!zero.inx])/10;
 
  p.log <- -log10(p.value);
  if(paramSet$data.org == "omk"){
    anot.id <- rownames(data);
    gene.anot <- data.frame(gene_id=anot.id, symbol=anot.id, name=anot.id, stringsAsFactors=FALSE)
  }else if (anal.type == "metadata" || dataSet$annotated ){ # annotated to entrez
    anot.id <- rownames(data);
    gene.anot <- doEntrezIDAnot(anot.id, paramSet$data.org, paramSet$data.idType);
  }else{
    anot.id <- rownames(data);
    gene.anot <- data.frame(gene_id=anot.id, symbol=anot.id, name=anot.id, stringsAsFactors=FALSE)
    paramSet$init.lib <- "NA"
  }

  # For phosphoproteomics, use symbol-based phosphosite labels for display
  # (e.g., GeneSymbol_T_247) while preserving raw phosphosite IDs in anot.id/gene_id.
  is_phospho <- (!is.null(paramSet$data.type) && paramSet$data.type == "phospho")
  if (is_phospho && !is.null(gene.anot) && nrow(gene.anot) > 0) {
    phospho.map <- try(readDataQs("phospho_symbol_map.qs", paramSet$anal.type, dataName), silent = TRUE)
    if (!inherits(phospho.map, "try-error") && !is.null(phospho.map) && nrow(phospho.map) > 0 &&
        "symbol" %in% colnames(phospho.map)) {
      row_ids <- as.character(rownames(data))
      disp_labels <- row_ids
      hit.ids <- intersect(row_ids, rownames(phospho.map))
      if (length(hit.ids) > 0) {
        map.syms <- as.character(phospho.map[hit.ids, "symbol", drop = TRUE])
        for (ii in seq_along(row_ids)) {
          pid <- row_ids[ii]
          if (!(pid %in% hit.ids)) next
          sym <- map.syms[which(hit.ids == pid)[1]]
          if (is.na(sym) || sym == "" || sym == "NA" || sym == pid) next
          # Use symbol directly - it already contains the full display name
          # with isoform and site suffix (e.g., "DOCK10-2_S_12")
          disp_labels[ii] <- sym
        }
        # Keep same semantics as other modules: symbol/name are display labels.
        if ("symbol" %in% colnames(gene.anot) && length(gene.anot$symbol) == length(disp_labels)) {
          gene.anot$symbol <- disp_labels
        }
        if ("name" %in% colnames(gene.anot) && length(gene.anot$name) == length(disp_labels)) {
          gene.anot$name <- disp_labels
        }
      }
    }
  }
  
  #gene symbol to be used for boxplot   
  
  # create a named matrix of sig vars for display
  fc.log <- data[, inx];
  if(limit_fc){
    hit.maxPos <- (which(fc.log> 10) )
    hit.maxNeg <- (which(fc.log< -10) )
    fc.log[hit.maxPos] <- 10;
    fc.log[hit.maxNeg] <- 10;
  }
  #fc.all <- res$fc.all;
  
  if(fcthresh != 0){
    inx.up <- fc.log > fcthresh & p.value < threshp;
    inx.down <- fc.log < -fcthresh & p.value < threshp;
  }else{
    inx.up <- fc.log > 0 & p.value < threshp;
    inx.down <- fc.log < 0 & p.value < threshp;
  }
  
  # create named sig table for display
  inx.imp <- (inx.up | inx.down) & inx.p;
  sig.var <- cbind(fc.log[inx.imp,drop=F], p.value[inx.imp,drop=F], p.log[inx.imp,drop=F]);
  colnames(sig.var) <- c("log2(FC)", "p.value", "-log10(p)");
  # first order by log(p), then by log(FC)
  ord.inx <- order(sig.var[,3], abs(sig.var[,1]), decreasing=T);
  sig.var <- sig.var[ord.inx,, drop=F];
  
  sig.var <- signif (sig.var, 5);
  sig.var1 <- sig.var;
  sig.var1 <- cbind(rownames(sig.var), sig.var);
  colnames(sig.var1) <- c("name", "log2(FC)", "p.value", "-log10(p)");
  
  ###########################
  ## for Volcano data
  ##########################
  
  if(paramSet$init.lib != "NA"){
    saveSet(paramSet, "paramSet");
    PerformVolcanoEnrichment(dataName, "enrichment_result", paramSet$init.lib, "null", "all", inx)
    paramSet <- readSet(paramSet, "paramSet");
    msgSet <- readSet(msgSet, "msgSet");
  }
  
  fileName <- "volcano.csv";
  jsonNm <- "volcano.json";
  json.obj <- rjson::toJSON(sig.var1);
  sink(jsonNm);
  cat(json.obj);
  sink();
  fast.write(signif (sig.var,5),file=fileName);
  colnames(gene.anot)[1] <- "anot.id"
  
  volcano <- list (
    raw.threshx = fcthresh,
    raw.threshy = threshp,
    paired = paired,
    thresh.y = -log10(threshp),
    fc.symb =rownames(data),
    fc.log = fc.log,
    fc.log.uniq = jitter(fc.log),
    inx.up = inx.up,
    inx.down = inx.down,
    p.log = p.log,
    p.raw = p.value,
    inx.p = inx.p,
    sig.mat = sig.var,
    conv = gene.anot,
    analType = anal.type,
    org=paramSet$data.org,
    dat.opt = paramSet$selDataNm,
    naviString = "Volcano Plot"
  );
  
  analSet$volcano <- volcano;
  saveSet(analSet, "analSet");
  sigDownIds <- GetVolcanoUpLftIDs();
  sigUpIds <- GetVolcanoUpRgtIDs();
  nonSigIds <- GetVolcanoDnIDs();
  
  volcano[["sigDownIds"]] <- sigDownIds;
  volcano[["sigUpIds"]] <- sigUpIds;
  volcano[["nonSigIds"]] <- nonSigIds;
  
  jsonNm <- paste0(fileNm, ".json");
  json.obj <- rjson::toJSON(volcano);
  sink(jsonNm);
  cat(json.obj);
  sink();
  
  if(paramSet$init.lib == "NA"){
    enr.mat <- "NA"
  }else{
    enr.mat <- ov_qs_read("enr.mat.qs");
    #fast.write(enr.mat, file="enrichment_result.csv", row.names=T);
  }
  sink("enrichment_result.json");
  cat(json.obj);
  sink();
  #paramSet$partialToBeSaved <- c(paramSet$partialToBeSaved, c(jsonNm, "enrichment_result.csv"))
  paramSet$jsonNms["volcano"] <- fileNm;

    # Generate volcano_data
    volcano_data <- data.frame(
        gene = gene.anot$symbol,                  # Gene names
        log2FoldChange = fc.log,                  # Log fold change values
        pValue = p.value,                         # Raw p-values
        negLog10PValue = p.log,                   # -log10 p-value values
        significant = ifelse(inx.up & inx.p, "upregulated", 
                             ifelse(inx.down & inx.p, "downregulated", "nonsignificant")) # Determine significance
    )


  # Append hover_text after defining volcano_data
  volcano_data$hover_text <- with(volcano_data, paste("Gene: ", gene, 
                                                      "<br>Log2 FC: ", log2FoldChange, 
                                                      "<br>P-value: ", pValue))
    library(ggplot2)
    library(plotly)

    # Create a ggplot
    gg_volcano <- ggplot(volcano_data, aes(x = log2FoldChange, y = negLog10PValue, 
                                           color = significant, 
                                           text = paste("Gene: ", gene, "<br>Log2 FC: ", log2FoldChange, 
                                                        "<br>P-value: ", format(pValue, scientific = TRUE)))) +
      geom_point(alpha = 0.6) +  # Adjust point transparency
      scale_color_manual(values = c("upregulated" = "red", "downregulated" = "blue", "nonsignificant" = "grey")) +
      labs(x = "Log2 Fold Change", 
           y = "-Log10 P-value") +
      theme_minimal()

    # Convert to ggplotly for interactive plot, including tooltips
    pwidget <- ggplotly(gg_volcano, tooltip = "text")

    # Customize the layout to optimize hover interaction
    pwidget <- pwidget %>% layout(hovermode = 'closest')

    # Print the plot
    pwidget

  imgSet <- readSet(imgSet, "imgSet");
  widgetNm <- paste0(fileNm, ".rda");
  imgSet$volcanoPlotly <- widgetNm;

  save(pwidget, file = widgetNm);

  imgSet$volcanoPlot <- paste0(fileNm, ".png");

  Cairo::Cairo(file = imgSet$volcanoPlot, unit="px", dpi=dpi, width=1000, height=800, type=format, bg="white");
  print(gg_volcano)
  dev.off()

  saveSet(imgSet, "imgSet");
  saveSet(paramSet, "paramSet");
  saveSet(analSet, "analSet");
  print("Volcano OK");
  
  return(1);
}


GetVolcanoDnMat <- function(){
  analSet <- readSet(analSet, "analSet");
  vcn <- analSet$volcano;
  imp.inx <- (vcn$inx.up | vcn$inx.down) & vcn$inx.p;
  blue.inx <- which(!imp.inx);
  
  if(sum(blue.inx)>0){
    xs <- vcn$fc.log[blue.inx]
    ys <- vcn$p.log[blue.inx];
    return(as.matrix(cbind(xs, ys)));
  }else{
    return(as.matrix(cbind(-1, -1)));
  }
}


GetVolcanoUpLftMat <- function(){
  analSet <- readSet(analSet, "analSet");
  vcn <- analSet$volcano;
  imp.inx <- vcn$inx.down & vcn$inx.p;
  red.inx <- which(imp.inx);
  if(sum(red.inx)>0){
    xs <- vcn$fc.log[red.inx]
    ys <- vcn$p.log[red.inx];
    return(as.matrix(cbind(xs, ys)));
  }else{
    return(as.matrix(cbind(-1, -1)));
  }
}

GetVolcanoUpRgtMat <- function(){
  analSet <- readSet(analSet, "analSet");
  vcn <- analSet$volcano;
  imp.inx <- vcn$inx.up & vcn$inx.p;
  red.inx <- which(imp.inx);
  if(sum(red.inx)>0){
    xs <- vcn$fc.log[red.inx]
    ys <- vcn$p.log[red.inx];
    return(as.matrix(cbind(xs, ys)));
  }else{
    return(as.matrix(cbind(-1, -1)));
  }
}

GetVolcanoUpLftIDs <- function(){
  analSet <- readSet(analSet, "analSet");
  vcn <- analSet$volcano;
  imp.inx <- vcn$inx.down & vcn$inx.p;
  red.inx <- which(imp.inx);
  if(sum(red.inx)>0){
    return(names(vcn$fc.log)[red.inx]);
  }else{
    return("NA");
  }
}

GetVolcanoUpRgtIDs <- function(){
  analSet <- readSet(analSet, "analSet");
  vcn <- analSet$volcano;
  imp.inx <- vcn$inx.up & vcn$inx.p;
  red.inx <- which(imp.inx);
  if(sum(red.inx)>0){
    return(names(vcn$fc.log)[red.inx]);
  }else{
    return("NA");
  }
}

GetVolcanoDnIDs <- function(){
  analSet <- readSet(analSet, "analSet");
  vcn <- analSet$volcano;
  imp.inx <- (vcn$inx.up | vcn$inx.down) & vcn$inx.p;
  blue.inx <- which(!imp.inx);
  if(sum(blue.inx)>0){
    return(names(vcn$fc.log)[blue.inx]);
  }else{
    return("NA");
  }
}


PerformVolcanoEnrichment<-function(dataName="", file.nm, fun.type, IDs, type, inx){
  dataSet <- readDataset(dataName);
  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");

  fcthreshu <- suppressWarnings(as.numeric(paramSet$fcthreshu));
  if (length(fcthreshu) == 0 || is.na(fcthreshu)) {
    fcthreshu <- 0
    cat("[Volcano Enrichment] fcthreshu missing/invalid; defaulting to 0\n")
  }
  print(paste0(fcthreshu, "===fcthresh"));
  anal.type <- paramSet$anal.type;
  inx <- as.numeric(inx)
  if(anal.type == "onedata"){
    if(dataSet$type == "array"){
      sigmat <- dataSet$sig.mat
    } else {
      sigmat <- dataSet$sig.mat
    }
  }else{
    if(paramSet$selDataNm == "meta_default"){
      sigmat <- analSet$meta.mat
      sigmat <- cbind(unname(analSet$meta.avgFC[rownames(sigmat)]), sigmat);
      inx <- 1;
    }else{
      sigmat <- analSet$inmex.ind[paramSet$selDataNm][[1]][which(analSet$inmex.ind[paramSet$selDataNm][[1]][,'Pval'] < as.numeric(paramSet$pvalu)),];
    }
  }
  
  if(type == "focus"){
    gene.vec <- unlist(strsplit(IDs, "; "));
  }else if(type == "all"){
    gene.vecup <- rownames(sigmat[which(sigmat[,inx] > fcthreshu),]);
    gene.vecdown <- rownames(sigmat[which(sigmat[,inx] < -fcthreshu),]);
    gene.vec <- c(gene.vecup, gene.vecdown);
  }else if(type == "up"){
    gene.vec <- rownames(sigmat[which(sigmat[,inx] > fcthreshu),]);
  }else{
    gene.vec <- rownames(sigmat[which(sigmat[,inx] < -fcthreshu),]);
  }

  # Check if phospho data and use appropriate functions
  is_phospho <- (!is.null(paramSet$data.type) && paramSet$data.type == "phospho")

  if (is_phospho) {
    # For phospho data, gene.vec contains phosphosite IDs (uniprot+site)
    # Use phospho-aware enrichment
    cat(sprintf("[Volcano Enrichment] Phospho data detected - using phosphosite-aware enrichment\n"))
    res <- .performEnrichAnalysisPhospho(dataSet, file.nm, fun.type, gene.vec, "volcano")
  } else {
    # PROTEOMICS: gene.vec contains UniProt IDs (central ID)
    # Need to convert to Entrez IDs for enrichment library matching
    # Step 1: UniProt → Entrez (for matching enrichment libraries)
    uniprot.map <- queryGeneDB("entrez_uniprot", paramSet$data.org)
    hit.inx <- match(gene.vec, uniprot.map[, "accession"])
    entrez.vec <- uniprot.map[hit.inx, "gene_id"]

    # Step 2: Get symbols for display
    sym.vec <- doEntrez2SymbolMapping(gene.vec, paramSet$data.org, paramSet$data.idType)

    # For enrichment: use Entrez IDs (values) named by symbols (names)
    entrez.vec.named <- entrez.vec
    names(entrez.vec.named) <- sym.vec

    # Remove NAs (unmapped UniProts)
    entrez.vec.named <- entrez.vec.named[!is.na(entrez.vec.named)]

    res <- .performEnrichAnalysis(dataSet, file.nm, fun.type, entrez.vec.named, "volcano")
  }

  return(res);
}


# note: hit.query, resTable must synchronize
# ora.vec should contains entrez ids, named by their gene symbols
PerformVolcanoBatchEnrichment <- function(dataName="", file.nm, fun.type, IDs, inx){
  dataSet <- readDataset(dataName);
  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");

  msgSet <- readSet(msgSet, "msgSet");
  anal.type <- paramSet$anal.type;
  # prepare lib
  inx <- as.numeric(inx);
  if(anal.type == "onedata"){
    if(dataSet$type == "array"){
      sigmat <- dataSet$sig.mat
    } else {
      sigmat <- dataSet$sig.mat
    }
  }else{
    sigmat <- analSet$inmex.ind[paramSet$selDataNm][[1]][which(analSet$inmex.ind[paramSet$selDataNm][[1]][,'Pval'] < as.numeric(paramSet$pvalu)),];
  }
  
  one.path.vec <- unlist(strsplit(IDs, "; "));
  
  fcthreshu <- suppressWarnings(as.numeric(paramSet$fcthreshu))
  if (length(fcthreshu) == 0 || is.na(fcthreshu)) {
    fcthreshu <- 0
  }
  gene.vecup <- rownames(sigmat[which(sigmat[,inx] > fcthreshu),]);
  gene.vecdown <- rownames(sigmat[which(sigmat[,inx] < -fcthreshu),]);
  ora.vec <- c(gene.vecup, gene.vecdown);

  # PROTEOMICS: ora.vec contains UniProt IDs, need to convert to Entrez for enrichment
  uniprot.map <- queryGeneDB("entrez_uniprot", paramSet$data.org)
  hit.inx <- match(ora.vec, uniprot.map[, "accession"])
  entrez.vec <- uniprot.map[hit.inx, "gene_id"]

  # Get symbols for display
  sym.vec <- doEntrez2SymbolMapping(ora.vec, paramSet$data.org, paramSet$data.idType);

  # Use Entrez IDs for matching, named by symbols
  ora.vec <- entrez.vec
  names(ora.vec) <- sym.vec

  # Remove NAs (unmapped)
  na.inx <- is.na(ora.vec)
  ora.vec <- ora.vec[!na.inx]
  sym.vec <- sym.vec[!na.inx]
  
  current.featureset <- list()
  current.featureset[["Set"]] <- one.path.vec
  current.featureset[["Set2"]] <- one.path.vec
  
  # prepare query
  ora.nms <- names(ora.vec);
  
  # prepare for the result table
  set.size<-length(current.featureset);
  res.mat<-matrix(0, nrow=set.size, ncol=5);
  rownames(res.mat)<-names(current.featureset);
  colnames(res.mat)<-c("Total", "Expected", "Hits", "Pval", "FDR");
  
  # need to cut to the universe covered by the pathways, not all features
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
  
  hits.inx <- ora.vec %in% current.universe;
  ora.vec <- ora.vec[hits.inx];
  ora.nms <- ora.nms[hits.inx];
  
  q.size<-length(ora.vec);
  
  # get the matched query for each pathway
  hits.query <- lapply(current.featureset, 
                       function(x) {
                         ora.nms[ora.vec%in%unlist(x)];
                       }
  );
  
  ov_qs_save(hits.query, "hits_query.qs");
  
  names(hits.query) <- names(current.featureset);
  hit.num<-unlist(lapply(hits.query, function(x){length(unique(x))}), use.names=FALSE);
  
  # total unique gene number
  #uniq.count <- length(current.universe);
  uniq.count <- nrow(dataSet$data.norm);
  
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
  res.mat[,5] <- raw.pvals;
  
  # now, clean up result, synchronize with hit.query
  res.mat <- res.mat[hit.num>0,,drop = F];
  hits.query <- hits.query[hit.num>0];
  if(nrow(res.mat)> 0){
    # order by p value
    ord.inx<-order(res.mat[,4]);
    #res.mat <- signif(res.mat[ord.inx,],3);
    hits.query <- hits.query[ord.inx];
    
    imp.inx <- res.mat[,4] <= 0.05;
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
    return(0);
  }
  
  #get gene symbols
  resTable <- data.frame(Pathway=rownames(res.mat), res.mat);
  res.mat[,"Hits"] <- res.mat[,"Hits"]

  # Check for and handle duplicate row names in enr.mat
  if(any(duplicated(rownames(res.mat)))) {
    res.mat <- res.mat[!duplicated(rownames(res.mat)), ]
    hits.query <- hits.query[match(rownames(res.mat), names(hits.query))]

    print("Duplicates in enr.mat were removed.")
  } else {
    res.mat <- res.mat
  }

  ov_qs_save(res.mat, "enr.mat.qs");
  msgSet$current.msg <- "Functional enrichment analysis was completed";
  
  # write json
  fun.anot <- hits.query; 
  total <- resTable$Total; if(length(total) ==1) { total <- matrix(total) };
  fun.pval <- resTable$Pval; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval) };
  fun.padj <- resTable$FDR; if(length(fun.padj) ==1) { fun.padj <- matrix(fun.padj) };
  hit.num <- paste0(resTable$Hits,"/",resTable$Total); if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
  fun.ids <- as.vector(names(fun.anot)); 
  if(length(fun.ids) ==1) { fun.ids <- matrix(fun.ids) };
  json.res <- list(
    fun.link = "",
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
  fun.pval <<- resTable$Pval;
  hit.num <<- resTable$Hits;
  csv.nm <- paste(file.nm, ".csv", sep="");    
  fast.write(resTable, file=csv.nm, row.names=F);
  
  saveSet(msgSet, "msgSet");
  return(1);
}
