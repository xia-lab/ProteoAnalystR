##################################################
## R script for ProteoAnalyst
## Description: Compute Ridgeline plot
## Authors: 
## Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
## Jessica Ewald, jessica.ewald@mail.mcgill.ca
##################################################

compute.ridgeline <- function(dataSet, imgNm = "abc", dpi=96, format="png", fun.type = "kegg", ridgeType = "ora", ridgeColor = "teal",rankOpt="fc", sigLevel = 0.05, pwNum=20, inx = 1){
  
  #save.image("ridge.RData");
  paramSet <- readSet(paramSet, "paramSet");
  msgSet <- readSet(msgSet, "msgSet");
  analSet <- readSet(analSet, "analSet");
  imageName <- paste0(imgNm, "dpi" , dpi, ".", format);
  anal.type <- paramSet$anal.type;
  require("dplyr");
  require("fgsea");
  require("reshape");
  require("ggplot2");
  require("ggridges");
  safe_qread <- function(path) {
    if (!file.exists(path)) {
      stop("[Ridgeline] Required file missing: ", path)
    }
    tryCatch(
      ov_qs_read(path),
      error = function(e) {
        stop("[Ridgeline] Failed to read ", path, " | ", conditionMessage(e))
      }
    )
  }
  # process colors
  if(ridgeColor == "teal"){
    high.col = "#3C7C60";
    low.col = "#C9E1D6";
  } else if (ridgeColor == "orange"){
    high.col = "#ffa34d";
    low.col = "#994a00";    
  } else {
    high.col = "#00da00";
    low.col = "#005000";   
  }
  
  #get pw library
  setres <- .loadEnrichLib(fun.type, paramSet);
  current.featureset <- setres$current.featureset;
  # get DEGs
  if(anal.type == "proteinlist"){
    if(paramSet$numOfLists > 1){
      dataSet <- readDataset(paramSet$selDataNm);
    }
    sigmat <- as.data.frame(dataSet$prot.mat)
    msg("[Ridgeline] For proteinlist - sigmat has ", nrow(sigmat), " rows")
    msg("[Ridgeline] sigmat rownames (first 5): ", paste(head(rownames(sigmat), 5), collapse=", "))

    sigmat$entrez <- rownames(sigmat);
    universe <- unique(unlist(current.featureset));
    expr.vec <- sigmat[,1];

    if(sum(expr.vec) == 0){
      msgSet$current.msg <- "Uploaded gene list needs to contain fold-change value to perform Ridgeline analysis!";
      saveSet(msgSet, "msgSet");
      return(-2);
    }
  }else if(anal.type == "onedata"){
    if(ridgeType == "ora"){
      sigmat <- dataSet$sig.mat;
      # Check if there are any significant features
      if(is.null(sigmat) || nrow(sigmat) == 0){
        msgSet <- readSet(msgSet, "msgSet");
        msgSet$current.msg <- "No significant features found for enrichment analysis. Please relax the differential expression threshold (e.g., increase p-value cutoff or decrease fold-change threshold) to obtain more features for pathway enrichment.";
        saveSet(msgSet, "msgSet");
        return(0);
      }
    }else{
      sigmat <- dataSet$comp.res;
      allmat <- dataSet$comp.res;
    }
    sigmat$entrez <- rownames(sigmat);
    universe <- rownames(dataSet$data.norm);
  }else{
    meta.avgFC <- analSet$meta.avgFC;
    inx <- 1;
    if(paramSet$selDataNm == "meta_default"){
        if(ridgeType == "ora"){
          sigmat <- analSet$meta.mat;
        # Check if there are any significant features
        if(is.null(sigmat) || nrow(sigmat) == 0){
          msgSet <- readSet(msgSet, "msgSet");
          msgSet$current.msg <- "No significant features found for enrichment analysis. Please relax the differential expression threshold (e.g., increase p-value cutoff or decrease fold-change threshold) to obtain more features for pathway enrichment.";
          saveSet(msgSet, "msgSet");
          return(0);
        }
        allmat <- safe_qread("meta.resTable.qs");
        sigmat <- cbind(unname(meta.avgFC[rownames(sigmat)]), sigmat);

      }else{
        allmat <- safe_qread("meta.resTable.qs");
        sigmat <- allmat;
        sigmat <- cbind(unname(meta.avgFC[rownames(sigmat)]), sigmat);

      }
      allmat$logFC <- unname(meta.avgFC[rownames(allmat)]);
      universe <- rownames(allmat);
    }else{
      dataSet <- readDataset(paramSet$selDataNm);
      if(ridgeType == "ora"){
        sigmat <- dataSet$sig.mat;
        # Check if there are any significant features
        if(is.null(sigmat) || nrow(sigmat) == 0){
          msgSet <- readSet(msgSet, "msgSet");
          msgSet$current.msg <- "No significant features found for enrichment analysis. Please relax the differential expression threshold (e.g., increase p-value cutoff or decrease fold-change threshold) to obtain more features for pathway enrichment.";
          saveSet(msgSet, "msgSet");
          return(0);
        }
      }else{
        sigmat <- dataSet$comp.res;
        allmat <- dataSet$comp.res;
      }
      sigmat <- as.data.frame(sigmat);
      sigmat$entrez <- rownames(sigmat);
      universe <- rownames(dataSet$data.norm);
    }
  }
  if(ridgeType == "ora"){
    gene.vec <- rownames(sigmat);

    # Check if phospho data and use appropriate functions
    is_phospho <- (!is.null(paramSet$data.type) && paramSet$data.type == "phospho")

    msg("[Ridgeline] About to perform enrichment analysis")
    msg("[Ridgeline] Current working directory: ", getwd())
    msg("[Ridgeline] Number of genes: ", length(gene.vec))
    msg("[Ridgeline] Is phospho data: ", is_phospho)
    msg("[Ridgeline] Fun type: ", fun.type)

    # Check if there are any significant genes
    if(length(gene.vec) == 0){
      msgSet <- readSet(msgSet, "msgSet");
      msgSet$current.msg <- "No significant features found for enrichment analysis. Please relax the differential expression threshold (e.g., increase p-value cutoff or decrease fold-change threshold) to obtain more features for pathway enrichment.";
      saveSet(msgSet, "msgSet");
      return(0);
    }

    if (is_phospho) {
      # For phospho data, gene.vec contains phosphosite IDs (uniprot+site)
      msg("[Ridgeline] Calling .performEnrichAnalysisPhospho")
      .performEnrichAnalysisPhospho(dataSet, imgNm, fun.type, gene.vec, "ridgeline")
    } else {
      # Regular data - use standard enrichment
      msg("[Ridgeline] Calling .performEnrichAnalysis")
      sym.vec <- doEntrez2SymbolMapping(gene.vec, paramSet$data.org, paramSet$data.idType);
      names(gene.vec) <- sym.vec;
      .performEnrichAnalysis(dataSet, imgNm, fun.type, gene.vec, "ridgeline")
    }

    msg("[Ridgeline] Enrichment analysis completed")
    msg("[Ridgeline] Files in working directory:")
    msg(paste(list.files(pattern = "*.qs"), collapse = ", "))

    # Check if file exists before reading
    if (!file.exists("enr.mat.qs")) {
      stop("[Ridgeline] ERROR: enr.mat.qs file not found after enrichment analysis. Working directory: ", getwd())
    }

    msg("[Ridgeline] Reading enr.mat.qs file")
    res <- safe_qread("enr.mat.qs");
    colnames(res) <- c("size", "expected", "overlap", "pval", "padj");

    res <- res[,c(4,5,3,1,2)]
    res <- as.data.frame(cbind(pathway=rownames(res), res));
    res$padj <- as.numeric(res$padj)
    res$pval <- as.numeric(res$pval)

    # Load hits.query saved by enrichment analysis (contains UniProt IDs after conversion)
    msg("[Ridgeline] Loading hits_query.qs file")
    if (!file.exists("hits_query.qs")) {
      stop("[Ridgeline] ERROR: hits_query.qs file not found. Working directory: ", getwd())
    }
    saved.hits.query <- safe_qread("hits_query.qs");
    msg("[Ridgeline] Loaded hits.query with ", length(saved.hits.query), " pathways")
  } else {
    # GSEA branch
    msg("[Ridgeline] GSEA analysis - preparing ranked vector")

    rankedVec<- ComputeRankedVec(dataSet, rankOpt, paramSet$selectedFactorInx);
    msg("[Ridgeline] RankedVec length: ", length(rankedVec))
    msg("[Ridgeline] Sample rankedVec names: ", paste(head(names(rankedVec), 5), collapse = ", "))

    # For proteomics data, names(rankedVec) are UniProt IDs, need to convert to Entrez
    # Universe for proteinlist is Entrez IDs from library, for onedata is UniProt IDs from data
    if(anal.type == "proteinlist"){
      # Universe already contains Entrez IDs from current.featureset
      gene.vec <- universe;
      msg("[Ridgeline] Proteinlist: Universe already contains Entrez IDs (", length(gene.vec), " genes)")
    } else {
      # Universe contains UniProt IDs, need to convert to Entrez
      msg("[Ridgeline] Converting universe UniProt IDs to Entrez IDs...")

      # Normalize UniProt IDs (remove phosphosite annotations, isoforms)
      normalized.universe <- sub("_[A-Z]_\\d+$", "", universe)  # Remove phosphosites
      normalized.universe <- sub("-\\d+$", "", normalized.universe)  # Remove isoforms

      uniprot.map <- queryGeneDB("entrez_uniprot", paramSet$data.org)
      hit.inx <- match(normalized.universe, uniprot.map[, "accession"])
      entrez.vec <- uniprot.map[hit.inx, "gene_id"]
      gene.vec <- entrez.vec[!is.na(entrez.vec)]
      msg("[Ridgeline] Converted ", length(gene.vec), " UniProt IDs to Entrez IDs")
    }

    # Convert Entrez IDs to symbols for fgsea matching
    sym.vec <- doEntrez2SymbolMapping(gene.vec, paramSet$data.org, "entrez");
    gene.nms <- sym.vec;
    msg("[Ridgeline] Converted ", length(gene.nms), " Entrez IDs to symbols")

    # Create symbol-based gene sets from Entrez-based current.featureset
    current.featureset.symb <- lapply(current.featureset,
                         function(x) {
                           gene.nms[gene.vec%in%unlist(x)];
    }
    );

    # Convert rankedVec names from UniProt to symbols (via Entrez)
    msg("[Ridgeline] Converting rankedVec names from UniProt to symbols...")

    # Normalize UniProt IDs (remove phosphosite annotations, isoforms)
    normalized.names <- sub("_[A-Z]_\\d+$", "", names(rankedVec))  # Remove phosphosites
    normalized.names <- sub("-\\d+$", "", normalized.names)        # Remove isoforms

    uniprot.map <- queryGeneDB("entrez_uniprot", paramSet$data.org)
    hit.inx <- match(normalized.names, uniprot.map[, "accession"])
    entrez.vec <- uniprot.map[hit.inx, "gene_id"]
    sym.vec <- doEntrez2SymbolMapping(entrez.vec, paramSet$data.org, "entrez")
    names(rankedVec) <- sym.vec
    # Remove NAs
    rankedVec <- rankedVec[!is.na(names(rankedVec))]
    msg("[Ridgeline] RankedVec after conversion: ", length(rankedVec), " genes with symbols")


    if(fun.type %in% c("go_bp", "go_mf", "go_cc")){
      res <- fgsea::fgsea(pathways = current.featureset.symb, 
                          stats    = rankedVec,
                          minSize  = 5,
                          maxSize = 500,
                          scoreType = "std",
                          nperm=10000)    
    }else{
      res <- fgsea::fgsea(pathways = current.featureset.symb, 
                          stats    = rankedVec,
                          minSize  = 5,
                          maxSize = 500,
                          scoreType = "std")   
      
    }


  }
  res <- .signif_df(res, 4);
  res <- res[order(res$pval),];
  resTable <- res;
  # process results;
  res <- res[,c(1,2,3)];
  colnames(res) <- c("name", "pval", "adj.pval");
  res.sig <- res;
  totalSigPws <- dim(res[res$pval < sigLevel, ])[1];
  
  if(pwNum != -1){
    if(dim(res.sig)[1] > pwNum){ # limit size if too many sig results
      res.sig <- res.sig[1:pwNum, ]
    }
  }
  # Get hits per pathway early (needed for gs.plot creation)
  # For ORA, use the saved hits.query from enrichment analysis (contains UniProt IDs)
  # For GSEA, compute from current.featureset and convert to UniProt
  if(ridgeType == "ora"){
    # Use saved hits.query (already has UniProt IDs after conversion)
    hits.query <- saved.hits.query[resTable$pathway];
    msg("[Ridgeline] Using saved hits.query for ORA (", length(hits.query), " pathways)")
  } else {
    # GSEA: compute hits.query from current.featureset (Entrez IDs)
    # Use gene.vec (Entrez) as the universe for matching
    hits.query <- lapply(current.featureset,
                         function(x) {
                           x[x %in% gene.vec];
                         }
    );
    hits.query <- hits.query[resTable$pathway];

    # Convert Entrez IDs back to UniProt IDs for consistency
    msg("[Ridgeline] Converting GSEA hits.query from Entrez to UniProt...")
    uniprot.map <- queryGeneDB("entrez_uniprot", paramSet$data.org)

    hits.query <- lapply(hits.query, function(entrez.list) {
      # Find UniProt IDs for these Entrez IDs
      hit.inx <- match(entrez.list, uniprot.map[, "gene_id"])
      uniprot.list <- uniprot.map[hit.inx, "accession"]
      uniprot.list[!is.na(uniprot.list)]
    })

    msg("[Ridgeline] Converted GSEA hits.query to UniProt IDs (", length(hits.query), " pathways)")
  }

  # prepare data for plotting
  # For phospho data, rownames(sigmat) are phosphosite IDs that need to be converted to entrez
  is_phospho <- (!is.null(paramSet$data.type) && paramSet$data.type == "phospho")

  if (is_phospho && !is.null(paramSet$phospho.mapping)) {
    # Convert phosphosite IDs to entrez IDs using the mapping
    phospho.to.entrez <- paramSet$phospho.mapping$phospho.to.entrez
    phosphosite.ids <- rownames(sigmat)

    # Map phosphosites to entrez (keep only mapped ones)
    entrez.ids <- phospho.to.entrez[phosphosite.ids]
    entrez.ids <- entrez.ids[!is.na(entrez.ids)]

    # Get fold changes for mapped phosphosites
    phosphosite.with.entrez <- names(entrez.ids)
    fc.values <- sigmat[phosphosite.with.entrez, inx]

    # Create data frame with entrez IDs and fold changes
    degs.plot <- data.frame(entrez = as.character(entrez.ids),
                           log2FC = fc.values,
                           stringsAsFactors = FALSE)
    degs.plot <- reshape::melt(degs.plot);
    colnames(degs.plot)[1] <- "entrez";
  } else {
    # Regular data - use rownames as-is
    degs.plot <- data.frame(entrez = rownames(sigmat), log2FC = sigmat[,inx]);
    degs.plot <- reshape::melt(degs.plot);
    colnames(degs.plot)[1] <- "entrez";
  }

  msg("[Ridgeline] Created degs.plot with ", nrow(degs.plot), " rows")
  msg("[Ridgeline] degs.plot$entrez (first 5): ", paste(head(degs.plot$entrez, 5), collapse=", "))

  # For protein lists, degs.plot$entrez contains Entrez IDs (from dataSet$prot.mat rownames)
  # But hits.query contains UniProt IDs, so we need to convert back
  if(anal.type == "proteinlist" && paramSet$data.idType == "uniprot") {
    symbol.map <- readSet(symbol.map, "symbol.map")
    if("uniprot" %in% colnames(symbol.map)) {
      # Create Entrez -> UniProt mapping
      entrez.to.uniprot <- setNames(symbol.map$uniprot, as.character(symbol.map$gene_id))
      # Convert degs.plot$entrez from Entrez to UniProt
      degs.plot$entrez <- entrez.to.uniprot[as.character(degs.plot$entrez)]
      # Remove NAs (genes that couldn't be mapped back)
      degs.plot <- degs.plot[!is.na(degs.plot$entrez), ]
      msg("[Ridgeline] Converted degs.plot IDs from Entrez to UniProt for matching")
      msg("[Ridgeline] After conversion - degs.plot has ", nrow(degs.plot), " rows")
      msg("[Ridgeline] After conversion - degs.plot$entrez (first 5): ", paste(head(degs.plot$entrez, 5), collapse=", "))
    } else {
      msg("[Ridgeline] WARNING: symbol.map does not have uniprot column, cannot convert IDs")
    }
  }

  # Use hits.query (which has UniProt IDs for both ORA and GSEA) instead of current.featureset (which has Entrez IDs)
  # This ensures gs.plot matches with degs.plot
  gs.plot <- reshape::melt(hits.query);
  colnames(gs.plot) <- c("entrez", "name");
  msg("[Ridgeline] Created gs.plot from hits.query (", ridgeType, " mode)")
  msg("[Ridgeline] gs.plot has ", nrow(gs.plot), " rows")
  msg("[Ridgeline] gs.plot$entrez (first 5): ", paste(head(gs.plot$entrez, 5), collapse=", "))
  msg("[Ridgeline] gs.plot$name (first 5): ", paste(head(gs.plot$name, 5), collapse=", "))

  # For phospho data in ORA mode, gs.plot$entrez contains phosphosite IDs
  # Need to convert to Entrez IDs to match degs.plot
  if (is_phospho && ridgeType == "ora") {
    msg("[Ridgeline] Converting gs.plot phosphosite IDs to Entrez IDs...")

    # Strip phosphosite IDs to UniProt IDs and remove isoforms
    uniprot_ids <- sapply(strsplit(as.character(gs.plot$entrez), "_"), function(x) x[1])
    uniprot_ids <- sub("-\\d+$", "", uniprot_ids)

    # Convert UniProt to Entrez
    uniprot.map <- queryGeneDB("entrez_uniprot", paramSet$data.org)
    hit.inx <- match(uniprot_ids, uniprot.map[, "accession"])
    entrez.ids <- uniprot.map[hit.inx, "gene_id"]

    # Update gs.plot with Entrez IDs
    gs.plot$entrez <- as.character(entrez.ids)

    # Remove rows with NA Entrez IDs
    gs.plot <- gs.plot[!is.na(gs.plot$entrez), ]

    msg("[Ridgeline] After conversion - gs.plot has ", nrow(gs.plot), " rows")
    msg("[Ridgeline] After conversion - gs.plot$entrez (first 5): ", paste(head(gs.plot$entrez, 5), collapse=", "))
  }

  df <- merge(res.sig, gs.plot, by = "name", all.x = TRUE, all.y = FALSE);
  df <- merge(df, degs.plot, by = "entrez", all.x = TRUE, all.y = FALSE);
  df <- na.omit(df)

  msg("[Ridgeline] After merges - df has ", nrow(df), " rows and ", ncol(df), " columns")
  msg("[Ridgeline] df column names: ", paste(colnames(df), collapse=", "))
  if(nrow(df) > 0) {
    msg("[Ridgeline] df$name has ", length(levels(df$name)), " unique levels")
  }

  # Check if we have any data to plot
  if (nrow(df) == 0) {
    msgSet <- readSet(msgSet, "msgSet");
    msgSet$current.msg <- "No data available for ridgeline plot. Enrichment results may not have overlapping genes with expression data.";
    saveSet(msgSet, "msgSet");
    return(0);
  }

  # calculate the mean fold change to order the pathways in the plot
  means <- aggregate(df$value, by = list(df$name), mean);
  means <- means[order(means$x, decreasing = FALSE), ];
  df$name <- factor(df$name, levels = means$Group.1);
  
  # make the plot
  rp <- ggplot(df, aes(x = value, y = name, fill = adj.pval)) +
    geom_density_ridges(
      jittered_points = TRUE, point_shape = "|", point_size = 5, point_color = "#898A89",
      color = "white",
      scale = 1.5, rel_min_height = .02, size = 0.25,
      position = position_points_jitter(height = 0)) +
    geom_vline(xintercept = 0, color = "red") +
    scale_y_discrete(expand = c(0, 0), name = "Gene Set") + 
    scale_x_continuous(expand = c(0, 0), name = "log2FC") +
    scale_fill_gradient("adj. pval",
                        low = high.col, high = low.col) + 
    coord_cartesian(clip = "off") +
    theme_ridges(center = TRUE) +
    theme(legend.position = "right",
          text = element_text(size=12, color = "black"),
          axis.title = element_text(size=12, face = "bold"),
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(size=12,color = "black"))
  
  Cairo::Cairo(file=imageName, width=10, height=8, type=format, bg="white", dpi=dpi, unit="in");
  print(rp);
  dev.off();
  
  ##interative ridge json data
  ridge_bw <- rp$layers[[1]]$computed_stat_params$bandwidth;
  jsonNm <- paste0(imgNm, ".json");

  # Map IDs to symbols
  # For both ORA and GSEA, df$entrez now contains UniProt IDs (after conversion)
  symb <- doEntrez2SymbolMapping(df$entrez, paramSet$data.org, "uniprot");
  msg("[Ridgeline] Mapped ", length(symb), " UniProt IDs to symbols (", ridgeType, " mode)")
  df$symbol <- symb;
  
  data.list <- list();
  gene.list <- list();
  pval.list <- list();
  col.list <- list();
  
  for(i in 1:length(levels(df$name))){
    nm <- as.character(levels(df$name)[i]);
    data.list[[ nm ]] <- as.vector(unlist(df[which(df$name == nm), "value"]));
    gene.list[[ nm ]] <- as.vector(unlist(df[which(df$name == nm), "symbol"]));
    pval.list[[ nm ]] <- unname(unlist(res[which(res$name == nm), "pval"]));
  }

  msg("[Ridgeline] After loop - data.list length: ", length(data.list))
  msg("[Ridgeline] After loop - gene.list length: ", length(gene.list))
  msg("[Ridgeline] After loop - pval.list length: ", length(pval.list))

  minFc <- min(df$value);
  maxFc <- max(df$value);
  minPval <- min(df$pval);
  maxPval <- max(df$pval);

  # enr result for table display (hits.query already assigned earlier)
  fun.anot <- hits.query
  if(ridgeType == "ora"){
    total <- resTable[,5]; if(length(total) ==1) { total <- matrix(total) };
    fun.pval <- resTable[,"pval"]; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval) };
    fun.padj <- resTable[,"padj"]; if(length(fun.padj) ==1) { fun.padj <- matrix(fun.padj) };
    if(ridgeType == "ora"){
      hit.num <- resTable[,4]; if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
    }else{
      hit.num <- resTable[,"size"]; if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
    }  
  }else{
    total <- as.list(resTable[,5])[[1]]; if(length(total) ==1) { total <- matrix(total) };
    fun.pval <- as.list(resTable[,"pval"])[[1]]; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval) };
    fun.padj <- as.list(resTable[,"padj"])[[1]]; if(length(fun.padj) ==1) { fun.padj <- matrix(fun.padj) };
    if(ridgeType == "ora"){
      hit.num <- as.list(resTable[,4])[[1]]; if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
    }else{
      hit.num <- as.list(resTable[,"size"])[[1]]; if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
    }  
    
  }

  # Get pathway IDs - handle NULL cases
  if(!is.null(setres$current.setids) && length(names(fun.anot)) > 0) {
    fun.ids <- as.vector(setres$current.setids[names(fun.anot)]);
    # Replace NAs with empty strings
    fun.ids[is.na(fun.ids)] <- ""
  } else {
    # Fallback: use pathway names as IDs
    fun.ids <- names(fun.anot)
  }
  if(length(fun.ids) ==1) { fun.ids <- matrix(fun.ids) };

  msg("[Ridgeline] Creating enr.res with ", length(fun.anot), " pathways")
  msg("[Ridgeline] fun.anot is: ", class(fun.anot), " with length ", length(fun.anot))
  msg("[Ridgeline] fun.ids length: ", length(fun.ids))

  enr.res <- list(
    fun.anot = fun.anot,
    fun.ids = fun.ids,
    fun.pval = fun.pval,
    fun.padj = fun.padj,
    hit.num = hit.num,
    total= total
  );
  
  if(ridgeType == "gsea"){
    enr.res[["ES"]] <- unname(unlist(resTable[,"ES"]));
  }

  msg("[Ridgeline] About to create res.list")
  msg("[Ridgeline] - data.list is ", class(data.list), " with length ", length(data.list))
  msg("[Ridgeline] - gene.list is ", class(gene.list), " with length ", length(gene.list))
  msg("[Ridgeline] - df is ", class(df), " with ", nrow(df), " rows")
  msg("[Ridgeline] - enr.res is ", class(enr.res), " with ", length(enr.res), " elements")
  msg("[Ridgeline] - bandwidth: ", ridge_bw)

  res.list <- list(data=data.list,
                   proteinlist=gene.list,
                   df=df, minPval=minPval,
                   maxPval=maxPval,
                   min=minFc,
                   max=maxFc,
                   minPval = min(res.sig$pval),
                   maxPval = max(res.sig$pval),
                   bandwidth=ridge_bw,
                   pathwayPvals = pval.list,
                   pathwayCols = col.list,
                   enrRes = enr.res,
                   dat.opt = paramSet$selDataNm,
                   naviString="ridge");

  msg("[Ridgeline] Created res.list with ", length(res.list), " elements")
  msg("[Ridgeline] res.list keys: ", paste(names(res.list), collapse=", "))

  if(ridgeType == "gsea"){
  csv.nm <- paste0(imgNm, ".csv");
  fast.write(resTable, file=csv.nm);
  }

  analSet$ridgeline <- res.list;
  saveSet(analSet, "analSet");

  msg("[Ridgeline] Before JSON serialization - res.list structure:")
  msg("[Ridgeline] - names: ", paste(names(res.list), collapse=", "))
  msg("[Ridgeline] Writing JSON to: ", jsonNm)

  json.obj <- rjson::toJSON(res.list);
  sink(jsonNm);
  cat(json.obj);
  sink();

  msg("[Ridgeline] JSON file written successfully")
  
  #for link sharing
  paramSet$jsonNms$ridge <- jsonNm
  paramSet$partialToBeSaved <- c( paramSet$partialToBeSaved, c(jsonNm));
  saveSet(paramSet, "paramSet");
  
  imgSet <- readSet(imgSet, "imgSet");
  rownames(resTable) <- NULL;
  
  imgSet$compute.ridgeline <- imageName;
    saveSet(imgSet, "imgSet");
  
  return(totalSigPws)
}
