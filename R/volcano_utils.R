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

  comp.result <- .getVolcanoCompartments(gene.anot$anot.id, rownames(data), paramSet, gene.anot$symbol)
  comp.vec <- comp.result$compartment
  comp.raw.vec <- comp.result$raw

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
    compartment = as.list(comp.vec),
    compartmentRaw = as.list(comp.raw.vec),
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
  
  if(paramSet$init.lib == "NA" || !file.exists("enr.mat.qs")){
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


.getVolcanoCompartments <- function(anot.ids, raw.ids, paramSet, symbols = NULL) {
  comp.map <- rep("Unknown", length(raw.ids))
  raw.map  <- rep("", length(raw.ids))
  names(comp.map) <- raw.ids
  names(raw.map)  <- raw.ids
  diag.log <- if (exists(".paProteinDiagLog")) {
    .paProteinDiagLog
  } else {
    function(...) {
      txt <- paste0(...)
      try(message(txt), silent = TRUE)
      try(cat(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", txt, "\n"),
              file = "protein_localization_diagnostics.log", append = TRUE), silent = TRUE)
    }
  }
  max.rest.fallback <- if (exists("pa.uniprot.rest.max.ids", envir = .GlobalEnv)) {
    get("pa.uniprot.rest.max.ids", envir = .GlobalEnv)
  } else {
    50
  }

  if (is.null(symbols) || length(symbols) != length(raw.ids)) {
    diag.log("[VolcanoCompartment] symbols unavailable or length mismatch; raw.ids=",
             length(raw.ids), "; symbols=", ifelse(is.null(symbols), "NULL", length(symbols)))
    symbols <- rep(NA_character_, length(raw.ids))
  }

  org <- paramSet$data.org
  diag.log("[VolcanoCompartment] lookup start; org=", org, "; features=", length(raw.ids),
           "; sample raw=", paste(head(raw.ids, 5), collapse = ","),
           "; sample anot=", paste(head(anot.ids, 5), collapse = ","),
           "; sample symbols=", paste(head(symbols, 5), collapse = ","))

  if (is.null(org) || org == "omk" || org == "NA") {
    diag.log("[VolcanoCompartment] skipped because organism is unavailable: ", org)
    return(list(compartment=comp.map, raw=raw.map))
  }

  lib.path <- if (exists("api.lib.path")) api.lib.path else paramSet$lib.path

  loc.path <- paste0(lib.path, org, "/", org, "_localization.qs")
  if (!file.exists(loc.path)) loc.path <- paste0(lib.path, org, "/localization.qs")
  if (!file.exists(loc.path)) {
    diag.log("[VolcanoCompartment] localization file missing for org=", org, "; lib.path=", lib.path)
    return(list(compartment=comp.map, raw=raw.map))
  }

  loc.data <- try(ov_qs_read(loc.path), silent = TRUE)
  if (inherits(loc.data, "try-error") || is.null(loc.data)) {
    diag.log("[VolcanoCompartment] failed to read localization file: ", loc.path,
             "; error=", ifelse(inherits(loc.data, "try-error"), as.character(loc.data), "NULL"))
    return(list(compartment=comp.map, raw=raw.map))
  }
  if (!all(c("EntrezID", "Broad.category") %in% colnames(loc.data))) {
    diag.log("[VolcanoCompartment] localization file missing required columns; columns=",
             paste(colnames(loc.data), collapse = ","))
    return(list(compartment=comp.map, raw=raw.map))
  }
  diag.log("[VolcanoCompartment] localization loaded; rows=", nrow(loc.data),
           "; path=", loc.path)

  # anot.ids are UniProt IDs (doEntrezIDAnot/doUniprotIDAnot returns UniProt as gene_id).
  # Map UniProt -> Entrez before hitting the localization table. Some accessions
  # can be missing from stale entrez_uniprot tables, so recover through symbols
  # and a small UniProt primary-gene fallback where possible.
  clean.ids <- if (exists(".normalizeUniProtAccessionsForLookup")) {
    .normalizeUniProtAccessionsForLookup(anot.ids)
  } else if (exists(".paNormalizeUniprotIds")) {
    .paNormalizeUniprotIds(anot.ids)
  } else {
    sub("_.*$", "", as.character(anot.ids))
  }
  uniprot.db <- try(queryGeneDB("entrez_uniprot", org), silent = TRUE)
  if (!inherits(uniprot.db, "try-error") && !is.null(uniprot.db) &&
      "accession" %in% colnames(uniprot.db) && "gene_id" %in% colnames(uniprot.db)) {
    up2entrez <- setNames(as.character(uniprot.db$gene_id), as.character(uniprot.db$accession))
    entrez.ids <- up2entrez[clean.ids]
    diag.log("[VolcanoCompartment] entrez_uniprot lookup complete; table rows=", nrow(uniprot.db),
             "; mapped=", sum(!is.na(entrez.ids) & nzchar(as.character(entrez.ids))), "/", length(entrez.ids),
             "; missing sample=", paste(head(clean.ids[is.na(entrez.ids) | !nzchar(as.character(entrez.ids))], 5), collapse = ","))
  } else {
    entrez.ids <- as.character(anot.ids)
    diag.log("[VolcanoCompartment] entrez_uniprot unavailable; using annotation IDs directly; error=",
             ifelse(inherits(uniprot.db, "try-error"), as.character(uniprot.db), "not a valid table"))
  }

  loc.symbols <- if ("Gene.name" %in% colnames(loc.data)) {
    valid <- !is.na(loc.data$Gene.name) & nzchar(as.character(loc.data$Gene.name))
    setNames(as.character(loc.data$EntrezID[valid]), toupper(as.character(loc.data$Gene.name[valid])))
  } else {
    character(0)
  }

  missing <- which(is.na(entrez.ids) | !nzchar(as.character(entrez.ids)) |
                     !(as.character(entrez.ids) %in% as.character(loc.data$EntrezID)))
  if (length(missing) > 0 && length(loc.symbols) > 0) {
    clean.symbols <- toupper(trimws(as.character(symbols)))
    by.symbol <- loc.symbols[clean.symbols[missing]]
    fill <- !is.na(by.symbol) & nzchar(by.symbol)
    if (any(fill)) {
      entrez.ids[missing[fill]] <- by.symbol[fill]
    }
    diag.log("[VolcanoCompartment] symbol fallback using existing annotation symbols: ",
             sum(fill), "/", length(missing),
             "; sample symbols=", paste(head(clean.symbols[missing], 5), collapse = ","))
  } else if (length(missing) > 0) {
    diag.log("[VolcanoCompartment] symbol fallback skipped; missing=", length(missing),
             "; loc.symbols=", length(loc.symbols))
  }

  missing <- which(is.na(entrez.ids) | !nzchar(as.character(entrez.ids)) |
                     !(as.character(entrez.ids) %in% as.character(loc.data$EntrezID)))
  if (length(missing) > max.rest.fallback) {
    diag.log("[VolcanoCompartment] UniProt primary-symbol fallback skipped because unresolved batch is too large: ",
             length(missing), " > ", max.rest.fallback,
             ". Refresh entrez_uniprot SQLite instead.")
  } else if (length(missing) > 0 && length(loc.symbols) > 0 && exists(".fetchUniprotPrimarySymbols")) {
    fetched.symbols <- try(.fetchUniprotPrimarySymbols(clean.ids[missing]), silent = TRUE)
    if (!inherits(fetched.symbols, "try-error") && length(fetched.symbols) > 0) {
      fetched <- toupper(trimws(as.character(fetched.symbols[clean.ids[missing]])))
      by.symbol <- loc.symbols[fetched]
      fill <- !is.na(by.symbol) & nzchar(by.symbol)
      if (any(fill)) {
        entrez.ids[missing[fill]] <- by.symbol[fill]
      }
      diag.log("[VolcanoCompartment] UniProt primary-symbol fallback: ",
               sum(fill), "/", length(missing),
               "; sample fetched symbols=", paste(head(fetched, 5), collapse = ","))
    } else {
      diag.log("[VolcanoCompartment] UniProt primary-symbol fallback failed; missing=", length(missing),
               "; error=", ifelse(inherits(fetched.symbols, "try-error"), as.character(fetched.symbols), "empty"))
    }
  } else if (length(missing) > 0) {
    diag.log("[VolcanoCompartment] UniProt primary-symbol fallback skipped; missing=", length(missing),
             "; helper.exists=", exists(".fetchUniprotPrimarySymbols"),
             "; loc.symbols=", length(loc.symbols))
  }

  loc.entrez <- as.character(loc.data$EntrezID)
  hits <- which(!is.na(entrez.ids) & entrez.ids %in% loc.entrez)
  diag.log("[VolcanoCompartment] final localization hits=", length(hits), "/", length(raw.ids),
           "; still missing=", length(raw.ids) - length(hits),
           "; hit sample=", paste(head(raw.ids[hits], 5), collapse = ","))
  if (length(hits) > 0) {
    for (i in hits) {
      loc.rows <- loc.data[as.character(loc.data$EntrezID) == as.character(entrez.ids[i]), , drop = FALSE]
      broad <- paste(unique(as.character(loc.rows$Broad.category)), collapse = "; ")
      main <- if ("Main.location" %in% colnames(loc.rows)) paste(unique(as.character(loc.rows$Main.location)), collapse = "; ") else NA_character_
      resolved <- .paPrimaryCompartment(broad, main)
      comp.map[i] <- resolved$primary
      raw.map[i] <- resolved$all_categories
    }
  }
  diag.log("[VolcanoCompartment] assigned compartment counts: ",
           paste(names(table(comp.map)), as.integer(table(comp.map)), sep = "=", collapse = "; "))
  list(compartment=comp.map, raw=raw.map)
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
      sigmat <- analSet$inmex.ind[paramSet$selDataNm][[1]][which(analSet$inmex.ind[paramSet$selDataNm][[1]][,'Pval'] < as.numeric(paramSet$pvalu)), , drop=FALSE];
    }
  }

  if(!is.matrix(sigmat)) sigmat <- as.matrix(sigmat);

  if(type == "focus"){
    gene.vec <- unlist(strsplit(IDs, "; "));
  }else if(type == "all"){
    gene.vecup <- rownames(sigmat[which(sigmat[,inx] > fcthreshu), , drop=FALSE]);
    gene.vecdown <- rownames(sigmat[which(sigmat[,inx] < -fcthreshu), , drop=FALSE]);
    gene.vec <- c(gene.vecup, gene.vecdown);
  }else if(type == "up"){
    gene.vec <- rownames(sigmat[which(sigmat[,inx] > fcthreshu), , drop=FALSE]);
  }else{
    gene.vec <- rownames(sigmat[which(sigmat[,inx] < -fcthreshu), , drop=FALSE]);
  }

  # Check if phospho data and use appropriate functions
  is_phospho <- (!is.null(paramSet$data.type) && paramSet$data.type == "phospho")

  if (is_phospho) {
    res <- .performEnrichAnalysisPhospho(dataSet, file.nm, fun.type, gene.vec, "volcano")
  } else {
    sym.vec <- doEntrez2SymbolMapping(gene.vec, paramSet$data.org, paramSet$data.idType);
    names(gene.vec) <- sym.vec;
    res <- .performEnrichAnalysis(dataSet, file.nm, fun.type, gene.vec, "volcano");
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
    sigmat <- analSet$inmex.ind[paramSet$selDataNm][[1]][which(analSet$inmex.ind[paramSet$selDataNm][[1]][,'Pval'] < as.numeric(paramSet$pvalu)), , drop=FALSE];
  }

  if(!is.matrix(sigmat)) sigmat <- as.matrix(sigmat);

  one.path.vec <- unlist(strsplit(IDs, "; "));

  fcthreshu <- suppressWarnings(as.numeric(paramSet$fcthreshu))
  if (length(fcthreshu) == 0 || is.na(fcthreshu)) {
    fcthreshu <- 0
  }
  gene.vecup <- rownames(sigmat[which(sigmat[,inx] > fcthreshu), , drop=FALSE]);
  gene.vecdown <- rownames(sigmat[which(sigmat[,inx] < -fcthreshu), , drop=FALSE]);
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
    if(sum(imp.inx) < 20){ # too little left, give the top ones
      topn <- ifelse(nrow(res.mat) > 20, 20, nrow(res.mat));
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

GetVolcanoPathwayHeatmapJSON <- function(dataName="", IDs=""){
  dataSet <- readDataset(dataName);
  analSet <- readSet(analSet, "analSet");

  ids <- trimws(unlist(strsplit(IDs, ";", fixed = TRUE)));
  ids <- ids[ids != "" & !is.na(ids)];
  if (length(ids) == 0 || is.null(dataSet$data.norm)) {
    return(rjson::toJSON(list(status = "empty", rows = list(), cols = list(), values = list())));
  }

  dat <- dataSet$data.norm;
  if (!is.null(dataSet$rmidx) && length(dataSet$rmidx) > 0 && ncol(dat) >= max(dataSet$rmidx)) {
    dat <- dat[, -dataSet$rmidx, drop = FALSE];
  }
  dat <- as.matrix(dat);

  hit.ids <- ids[ids %in% rownames(dat)];
  if (length(hit.ids) == 0 && !is.null(analSet$volcano$conv)) {
    conv <- analSet$volcano$conv;
    if ("symbol" %in% colnames(conv) && "anot.id" %in% colnames(conv)) {
      hit.inx <- match(ids, as.character(conv$symbol));
      hit.ids <- as.character(conv$anot.id[hit.inx]);
      hit.ids <- hit.ids[!is.na(hit.ids) & hit.ids %in% rownames(dat)];
    }
  }
  hit.ids <- unique(hit.ids);
  if (length(hit.ids) == 0) {
    return(rjson::toJSON(list(status = "empty", rows = list(), cols = colnames(dat), values = list())));
  }

  mat <- dat[hit.ids, , drop = FALSE];
  row.labels <- hit.ids;
  if (!is.null(analSet$volcano$conv)) {
    conv <- analSet$volcano$conv;
    if ("symbol" %in% colnames(conv) && "anot.id" %in% colnames(conv)) {
      hit.inx <- match(hit.ids, as.character(conv$anot.id));
      labels <- as.character(conv$symbol[hit.inx]);
      ok <- !is.na(labels) & labels != "" & labels != "NA";
      row.labels[ok] <- labels[ok];
    }
  }

  values <- lapply(seq_len(nrow(mat)), function(i) {
    vals <- suppressWarnings(as.numeric(mat[i, ]));
    vals[is.infinite(vals)] <- NA;
    vals
  });

  annotations <- list();
  if (!is.null(dataSet$meta.info) && ncol(dataSet$meta.info) > 0) {
    meta <- data.frame(dataSet$meta.info, check.names = FALSE);
    smp.inx <- match(colnames(mat), rownames(meta));
    keep.cols <- seq_len(min(4, ncol(meta)));
    annotations <- lapply(keep.cols, function(j) {
      vals <- as.character(meta[smp.inx, j]);
      vals[is.na(vals)] <- "NA";
      meta.nm <- colnames(meta)[j];
      meta.type <- "discrete";
      if (!is.null(dataSet$disc.inx) && meta.nm %in% names(dataSet$disc.inx) && !isTRUE(dataSet$disc.inx[meta.nm])) {
        meta.type <- "continuous";
      }
      list(
        name = meta.nm,
        type = meta.type,
        values = vals
      )
    });
  }

  row.rank <- .getVolcanoHeatmapClusterRanks(mat, "row");
  col.rank <- .getVolcanoHeatmapClusterRanks(mat, "col");

  rjson::toJSON(list(
    status = "ok",
    rows = as.vector(row.labels),
    ids = as.vector(hit.ids),
    cols = as.vector(colnames(mat)),
    values = values,
    annotations = annotations,
    feature.cluster = row.rank,
    sample.cluster = col.rank
  ));
}

.getVolcanoHeatmapClusterRanks <- function(mat, margin = "row") {
  n <- if (margin == "row") nrow(mat) else ncol(mat);
  base <- seq_len(n);
  res <- list(pval = base, ward = base, average = base, single = base, complete = base);
  if (n <= 1) {
    return(res);
  }

  dat <- if (margin == "row") mat else t(mat);
  dat <- t(apply(dat, 1, function(x) {
    x <- suppressWarnings(as.numeric(x));
    if (all(is.na(x))) {
      return(rep(0, length(x)));
    }
    x[is.na(x)] <- median(x, na.rm = TRUE);
    sx <- stats::sd(x);
    if (is.na(sx) || sx == 0) {
      return(rep(0, length(x)));
    }
    as.numeric(scale(x));
  }));

  dist.obj <- try(stats::dist(dat), silent = TRUE);
  if (inherits(dist.obj, "try-error")) {
    return(res);
  }

  methods <- c(ward = "ward.D", average = "average", single = "single", complete = "complete");
  for (nm in names(methods)) {
    hc <- try(stats::hclust(dist.obj, method = methods[[nm]]), silent = TRUE);
    if (!inherits(hc, "try-error")) {
      ord <- hc$order;
      rk <- match(seq_len(n), ord);
      res[[nm]] <- rk;
    }
  }
  res;
}
