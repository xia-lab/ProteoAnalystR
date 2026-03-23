

my.enrich.net<-function(dataSet, netNm="abc", type="list", overlapType="mixed", analSet){
  enr.mat <- qs:::qread("enr.mat.qs");

  hits <-  enr.mat[,"Hits"];
  pvals <- enr.mat[,"Pval"];
  
  pvalue <- pvals;
  id <- names(pvalue);
  
  paramSet <- readSet(paramSet, "paramSet");
  anal.type <- paramSet$anal.type;
  
  if(is.null(enr.mat)){
    return(0);
  }
  
  require(igraph);
  require(reshape);

  current.featureset <- qs::qread("current_featureset.qs");
  hits.query <- qs::qread("hits_query.qs")
  hits.query <- hits.query[rownames(enr.mat)];
  featuresets <- hits.query;
  n <- nrow(enr.mat);
  w <- matrix(NA, nrow=n, ncol=n);
  colnames(w) <- rownames(w) <- id;
  for (i in 1:n) {
    for (j in i:n) {
      w[i,j] <- overlap_ratio(featuresets[id[i]], featuresets[id[j]], overlapType)
    }
  }
  wd <- reshape::melt(w);
  wd <- wd[wd[,1] != wd[,2],];
  wd <- wd[!is.na(wd[,3]),];
  
  g <- graph_from_data_frame(wd[,-3], directed=F);
  if(type == "list"){
    g <- delete_edges(g, E(g)[wd[,3] < 0.3]);
  }else{
    g <- delete_edges(g, E(g)[wd[,3] < 0.3]);
  }
  idx <- unlist(sapply(V(g)$name, function(x) which(x == id)));
  
  # define local function
  my.normalize <- function(x){
    return((x- min(x)) /(max(x)-min(x)))
  }
  my.rescale <- function(x, from, to){
    (x - min(x)) / max(x - min(x)) * (to - from) + from
  }
  

  # Compute rank-based color gradient instead of p-value based
  # Rank pathways by p-value (lower p-value = better rank = higher color intensity)
  pvalue_ranks <- rank(pvalue, ties.method = "first")
  # Normalize ranks to 0-1 range
  normalized_ranks <- (pvalue_ranks - 1) / (length(pvalue_ranks) - 1)

  # Ensure that you compute colors only for existing vertices in the graph
  existing_vertices <- V(g)$name
  vertex_colors <- ComputeColorGradient(normalized_ranks[names(normalized_ranks) %in% existing_vertices], "black", F, F)
  vertex_colorsw <- ComputeColorGradient(normalized_ranks[names(normalized_ranks) %in% existing_vertices], "white", F, F)

  # Assign colors only to existing vertices
  V(g)$color <- vertex_colors
  V(g)$colorw <- vertex_colorsw
  
  cnt <- hits;
  names(cnt) <- id;
  cnt2 <- cnt[V(g)$name];
  
  if (all(cnt2 == cnt2[1])){
    V(g)$size <- rep(16, length(cnt2))
  }else{
    V(g)$size <- my.rescale(log(cnt2+1, base=10), 8, 32);
  }
  
  # layout
  pos.xy <- layout_nicely(g);
  
  # now create the json object
  nodes <- vector(mode="list");
  node.nms <- V(g)$name;
  node.sizes <- V(g)$size;
  node.cols <- V(g)$color;
  node.colsw <- V(g)$colorw;
  
  for(i in 1:length(node.sizes)){
    nodes[[i]] <- list(
      id = node.nms[i],
      label=node.nms[i],
      size = node.sizes[i],
      true_size=node.sizes[i], 
      molType="set",
      colorb=node.cols[i],
      colorw=node.colsw[i],
      posx = pos.xy[i,1],
      posy = pos.xy[i,2]
    );
  }
  
  edge.mat <- as_edgelist(g);
  edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2]);
  
  # covert to json
  bedges <- stack(hits.query);
  b.mat <- matrix(NA, nrow=nrow(bedges), ncol=2);
  b.mat[,1] <- bedges[,"values"];
  b.mat[,2] <- as.character(bedges[,"ind"]);
  b.mat <- b.mat[complete.cases(b.mat),]
  b.mat <- unique(b.mat)

  # For phospho data, values are phosphosite IDs; map labels to corresponding node IDs in bg
  # Here source should be feature nodes, target is the enriched pathway/set
  colnames(b.mat) <- c("source", "target");
  bg <- graph_from_data_frame(b.mat, directed=F);
  idx <- unlist(sapply(V(bg)$name, function(x) which(x == id)));
  cols <- color_scale("red", "#E5C494");
  
  V(bg)$color[V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(-log(pvalue), "black", F, F);
  V(bg)$colorw[V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(-log(pvalue), "white", F, F);
  node.nms <- V(bg)$name;

  # Helper function: Convert UniProt IDs to gene symbols for proteomics data
  convert.uniprot.to.symbols <- function(uniprot.ids, org) {
    # Check if IDs are Entrez-like or UniProt
    is.entrez.like <- mean(grepl("^[0-9]+$", uniprot.ids)) > 0.9

    if(is.entrez.like) {
      # Already Entrez IDs
      return(doEntrez2SymbolMapping(uniprot.ids, org, "entrez"))
    } else {
      # UniProt IDs - convert to Entrez first, then to symbols
      normalized.ids <- sub("_[A-Z]_\\d+$", "", uniprot.ids)  # Remove phosphosites
      normalized.ids <- sub("-\\d+$", "", normalized.ids)      # Remove isoforms
      normalized.ids <- trimws(normalized.ids)

      uniprot.map <- queryGeneDB("entrez_uniprot", org)
      hit.inx <- match(normalized.ids, uniprot.map[, "accession"])
      entrez.ids <- uniprot.map[hit.inx, "gene_id"]

      return(doEntrez2SymbolMapping(entrez.ids, org, "entrez"))
    }
  }

  # Initialize expvals as empty named vector (will be populated in conditionals below)
  expvals <- setNames(rep(0, length(node.nms)), node.nms)

  if(anal.type == "onedata"){
    tbl <- dataSet$comp.res
    # Gene node names in bg are UniProt IDs (or Entrez IDs from enrichment conversion)
    # Match directly using rownames(tbl) against V(bg)$name
    gene.node.nms <- V(bg)$name[!V(bg)$name %in% rownames(enr.mat)]

    # Try direct match first (rownames are UniProt IDs matching node names)
    keep.inx <- which(rownames(tbl) %in% gene.node.nms)

    if(length(keep.inx) == 0) {
      # Fallback: normalize both sides (remove phosphosite/isoform suffixes)
      norm.tbl.ids <- sub("_[A-Z]_\\d+$", "", rownames(tbl))
      norm.tbl.ids <- sub("-\\d+$", "", norm.tbl.ids)
      norm.node.nms <- sub("_[A-Z]_\\d+$", "", gene.node.nms)
      norm.node.nms <- sub("-\\d+$", "", norm.node.nms)
      keep.inx <- which(norm.tbl.ids %in% norm.node.nms)
    }

    if(length(keep.inx) > 0) {
      tbl.sub <- tbl[keep.inx, , drop=FALSE]
      expr.val <- as.numeric(tbl.sub[,paramSet$selectedFactorInx])
      # Key by the matching node names (UniProt IDs)
      matched.node.ids <- gene.node.nms[match(rownames(tbl.sub), gene.node.nms)]
      if(all(is.na(matched.node.ids))) {
        # Use normalized matching
        norm.tbl.ids <- sub("_[A-Z]_\\d+$", "", rownames(tbl.sub))
        norm.tbl.ids <- sub("-\\d+$", "", norm.tbl.ids)
        norm.node.nms <- sub("_[A-Z]_\\d+$", "", gene.node.nms)
        norm.node.nms <- sub("-\\d+$", "", norm.node.nms)
        matched.node.ids <- gene.node.nms[match(norm.tbl.ids, norm.node.nms)]
      }
      names(expr.val) <- matched.node.ids
      expr.val <- expr.val[!is.na(names(expr.val))]
      # Merge into expvals (keyed by node names)
      expvals[names(expr.val)] <- expr.val

      if(length(expr.val) > 0) {
        # Color gene nodes by expression values (fold change)
        gene.expvals <- expvals[gene.node.nms]
        gene.expvals <- gene.expvals[!is.na(gene.expvals)]
        V(bg)$color[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(unname(gene.expvals), "black", T, T);
        V(bg)$colorw[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(unname(gene.expvals), "black", T, T);
      } else {
        V(bg)$color[!V(bg)$name %in% rownames(enr.mat)] <- "#00FFFF";
        V(bg)$colorw[!V(bg)$name %in% rownames(enr.mat)] <- "#668B8B"
      }
    } else {
      V(bg)$color[!V(bg)$name %in% rownames(enr.mat)] <- "#00FFFF";
      V(bg)$colorw[!V(bg)$name %in% rownames(enr.mat)] <- "#668B8B"
    }
  }else if(anal.type == "proteinlist" && sum(as.numeric(paramSet$all.prot.mat[,1])) != 0){
    tbl <- paramSet$all.prot.mat
    # tbl[,1] = fold change values, tbl[,2] = gene symbols
    # Gene node names in bg are UniProt IDs - need to convert to match
    gene.node.nms <- V(bg)$name[which(!V(bg)$name %in% rownames(enr.mat))]

    # Convert gene node UniProt IDs to symbols for matching against tbl[,2]
    gene.node.symbols <- convert.uniprot.to.symbols(gene.node.nms, paramSet$data.org)
    gene.node.symbols[is.na(gene.node.symbols)] <- ""

    # Match tbl symbols against gene node symbols
    keep.inx <- which(tbl[,2] %in% gene.node.symbols)

    if(length(keep.inx) > 0) {
      tbl.sub <- tbl[keep.inx, , drop=FALSE]
      expr.val <- as.numeric(tbl.sub[,1])
      # Map back from symbol to UniProt node name
      matched.ids <- gene.node.nms[match(tbl.sub[,2], gene.node.symbols)]
      names(expr.val) <- matched.ids
      expr.val <- expr.val[!is.na(names(expr.val))]
      expvals[names(expr.val)] <- expr.val

      if(length(expr.val) > 0) {
        gene.expvals <- expvals[gene.node.nms]
        gene.expvals <- gene.expvals[!is.na(gene.expvals)]
        V(bg)$color[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(unname(gene.expvals), "black", T, T);
        V(bg)$colorw[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(unname(gene.expvals), "black", T, T);
      } else {
        V(bg)$color[!V(bg)$name %in% rownames(enr.mat)] <- "#00FFFF";
        V(bg)$colorw[!V(bg)$name %in% rownames(enr.mat)] <- "#668B8B"
      }
    } else {
      V(bg)$color[!V(bg)$name %in% rownames(enr.mat)] <- "#00FFFF";
      V(bg)$colorw[!V(bg)$name %in% rownames(enr.mat)] <- "#668B8B"
    }

  }else if(anal.type =="metadata"){
    gene.node.nms <- V(bg)$name[!V(bg)$name %in% rownames(enr.mat)]
    if(paramSet$selDataNm == "meta_default"){
      tbl <- analSet$meta.mat.all
      # Match rownames(tbl) directly against gene node UniProt IDs
      keep.inx <- which(rownames(tbl) %in% gene.node.nms)
      if(length(keep.inx) == 0) {
        norm.tbl.ids <- sub("_[A-Z]_\\d+$", "", rownames(tbl))
        norm.tbl.ids <- sub("-\\d+$", "", norm.tbl.ids)
        norm.node.nms <- sub("_[A-Z]_\\d+$", "", gene.node.nms)
        norm.node.nms <- sub("-\\d+$", "", norm.node.nms)
        keep.inx <- which(norm.tbl.ids %in% norm.node.nms)
      }

      if(length(keep.inx) > 0) {
        tbl.sub <- tbl[keep.inx, , drop=FALSE]
        expr.val <- analSet$meta.avgFC[rownames(tbl)[keep.inx]]
        matched.ids <- gene.node.nms[match(rownames(tbl.sub), gene.node.nms)]
        if(all(is.na(matched.ids))) {
          norm.ids <- sub("-\\d+$", "", sub("_[A-Z]_\\d+$", "", rownames(tbl.sub)))
          norm.nms <- sub("-\\d+$", "", sub("_[A-Z]_\\d+$", "", gene.node.nms))
          matched.ids <- gene.node.nms[match(norm.ids, norm.nms)]
        }
        names(expr.val) <- matched.ids
        expr.val <- expr.val[!is.na(names(expr.val))]
        expvals[names(expr.val)] <- expr.val

        if(length(expr.val) > 0) {
          gene.expvals <- expvals[gene.node.nms]
          gene.expvals <- gene.expvals[!is.na(gene.expvals)]
          V(bg)$color[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(unname(gene.expvals), "black", T, T);
          V(bg)$colorw[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(unname(gene.expvals), "black", T, T);
        } else {
          V(bg)$color[!V(bg)$name %in% rownames(enr.mat)] <- "#00FFFF";
          V(bg)$colorw[!V(bg)$name %in% rownames(enr.mat)] <- "#668B8B"
        }
      } else {
        V(bg)$color[!V(bg)$name %in% rownames(enr.mat)] <- "#00FFFF";
        V(bg)$colorw[!V(bg)$name %in% rownames(enr.mat)] <- "#668B8B"
      }
    }else{
      dataSet <- readDataset(paramSet$selDataNm);
      tbl <- dataSet$comp.res;
      # Match rownames(tbl) directly against gene node UniProt IDs
      keep.inx <- which(rownames(tbl) %in% gene.node.nms)
      if(length(keep.inx) == 0) {
        norm.tbl.ids <- sub("_[A-Z]_\\d+$", "", rownames(tbl))
        norm.tbl.ids <- sub("-\\d+$", "", norm.tbl.ids)
        norm.node.nms <- sub("_[A-Z]_\\d+$", "", gene.node.nms)
        norm.node.nms <- sub("-\\d+$", "", norm.node.nms)
        keep.inx <- which(norm.tbl.ids %in% norm.node.nms)
      }

      if(length(keep.inx) > 0) {
        tbl.sub <- tbl[keep.inx, , drop=FALSE]
        expr.val <- as.numeric(tbl.sub[,"logFC"])
        matched.ids <- gene.node.nms[match(rownames(tbl.sub), gene.node.nms)]
        if(all(is.na(matched.ids))) {
          norm.ids <- sub("-\\d+$", "", sub("_[A-Z]_\\d+$", "", rownames(tbl.sub)))
          norm.nms <- sub("-\\d+$", "", sub("_[A-Z]_\\d+$", "", gene.node.nms))
          matched.ids <- gene.node.nms[match(norm.ids, norm.nms)]
        }
        names(expr.val) <- matched.ids
        expr.val <- expr.val[!is.na(names(expr.val))]
        expvals[names(expr.val)] <- expr.val

        if(length(expr.val) > 0) {
          gene.expvals <- expvals[gene.node.nms]
          gene.expvals <- gene.expvals[!is.na(gene.expvals)]
          inx <- !V(bg)$name %in% rownames(enr.mat);
          V(bg)$color[inx] <- ComputeColorGradient(unname(gene.expvals), "black", T, T);
          V(bg)$colorw[inx] <- ComputeColorGradient(unname(gene.expvals), "black", T, T);
        } else {
          V(bg)$color[!V(bg)$name %in% rownames(enr.mat)] <- "#00FFFF";
          V(bg)$colorw[!V(bg)$name %in% rownames(enr.mat)] <- "#668B8B"
        }
      } else {
        V(bg)$color[!V(bg)$name %in% rownames(enr.mat)] <- "#00FFFF";
        V(bg)$colorw[!V(bg)$name %in% rownames(enr.mat)] <- "#668B8B"
      }
    }
  }else{
    # expvals already initialized above - just assign default colors
    V(bg)$color[!V(bg)$name %in% rownames(enr.mat)] <- "#00FFFF";
    V(bg)$colorw[!V(bg)$name %in% rownames(enr.mat)] <- "#668B8B"
  }

  node.dgr2 <- as.numeric(degree(bg));
  V(bg)$size <- my.rescale(log(node.dgr2, base=10), 8, 24); 
  
  # layout
  pos.xy <- layout_nicely(bg);
  
  # now create the json object
  bnodes <- vector(mode="list");
  node.sizes <- V(bg)$size;
  node.cols <- V(bg)$color;
  node.colsw <- V(bg)$colorw;
  
  shapes <- rep("circle", length(node.nms));
  hit.inx <- node.nms %in% b.mat[,"source"];
  shapes[hit.inx] <- "gene";

  # OPTIMIZED: Vectorized ID conversion to avoid per-node database queries
  # Separate pathway nodes from gene nodes
  pathway.inx <- node.nms %in% rownames(enr.mat)
  gene.ids <- node.nms[!pathway.inx]

  # For phospho enrichment networks, preserve phosphosite suffix while replacing
  # UniProt prefix with gene symbol using precomputed phospho_symbol_map.
  phospho.display.labels <- NULL
  if(length(gene.ids) > 0 && file.exists("phospho_symbol_map.qs")) {
    phospho.map <- try(readDataQs("phospho_symbol_map.qs", paramSet$anal.type, paramSet$dataName), silent = TRUE)
    if(!inherits(phospho.map, "try-error") && !is.null(phospho.map) && nrow(phospho.map) > 0 &&
       "symbol" %in% colnames(phospho.map)) {
      phospho.display.labels <- gene.ids
      hit.ids <- intersect(gene.ids, rownames(phospho.map))
      if(length(hit.ids) > 0) {
        map.syms <- as.character(phospho.map[hit.ids, "symbol", drop = TRUE])
        for(j in seq_along(gene.ids)) {
          gid <- gene.ids[j]
          if(!(gid %in% hit.ids)) next
          sym <- map.syms[which(hit.ids == gid)[1]]
          if(is.na(sym) || sym == "" || sym == "NA" || sym == gid) next
          # Use symbol directly - it already contains the full display name
          # with isoform and site suffix (e.g., "DOCK10-2_S_12")
          phospho.display.labels[j] <- sym
        }
      }
    }
  }

  # Convert all gene IDs at once (single database query)
  if(!is.null(phospho.display.labels)) {
    gene.symbols <- phospho.display.labels
  } else if(length(gene.ids) > 0) {
    gene.symbols <- convert.uniprot.to.symbols(gene.ids, paramSet$data.org)
    # Handle NAs - keep original ID if conversion failed
    gene.symbols[is.na(gene.symbols)] <- gene.ids[is.na(gene.symbols)]
  } else {
    gene.symbols <- character(0)
  }

  # Build node labels vector
  node.lbls <- character(length(node.nms))
  node.lbls[pathway.inx] <- node.nms[pathway.inx]  # Pathways keep their names
  node.lbls[!pathway.inx] <- gene.symbols           # Genes get symbol labels
  for(i in 1:length(node.sizes)){
    bnodes[[i]] <- list(
      id = node.nms[i],
      label=node.lbls[i], 
      size=node.sizes[i], 
      colorb=node.cols[i],
      colorw=node.colsw[i],
      true_size=node.sizes[i], 
      molType=shapes[i],
      exp= unname(expvals[node.nms[i]]),
      posx = pos.xy[i,1],
      posy = pos.xy[i,2]
    );
  }
  
  ppi.comps <- vector(mode="list");
  paramSet$current.net.nm <- netNm
  ppi.comps[[netNm]] <- bg;
  analSet$ppi.comps <- ppi.comps
  
  bedge.mat <- as_edgelist(bg);
  bedge.mat <- cbind(id=paste0("b", 1:nrow(bedge.mat)), source=bedge.mat[,1], target=bedge.mat[,2]);
  # analSet$list.features contains gene IDs (Entrez or UniProt), convert to symbols
  initsbls <- convert.uniprot.to.symbols(analSet$list.features, paramSet$data.org)
  names(initsbls) <- analSet$list.features
  
  #for rjson generation
  edge.mat <- apply(edge.mat, 1, as.list)
  bedge.mat <- apply(bedge.mat, 1, as.list)
  enr.mat <- apply(enr.mat, 1, as.list)
  
  #paramSet$current.sigmajsR.nm <- paste0(netNm, ".rda");
  #save(nodes, edge.mat, file = paramSet$current.sigmajsR.nm);

  netData <- list(nodes=nodes, 
                  edges=edge.mat, 
                  bnodes=bnodes, 
                  bedges=bedge.mat, 
                  enr=unname(enr.mat), 
                  id=names(enr.mat), 
                  sizes=analSet$listSizes, 
                  hits=hits.query, 
                  proteinlist=initsbls, 
                  analType=anal.type, 
                  org=paramSet$data.org, 
                  backgroundColor=list("#514F6A", "#222222"),
                  dat.opt = paramSet$selDataNm,
                  naviString = "Enrichment Network");
  
  netName <- paste0(netNm, ".json");
  paramSet$partialToBeSaved <- c( paramSet$partialToBeSaved, c(netName));
  paramSet$jsonNms$network <- netName;
  saveSet(paramSet, "paramSet");
  
  analSet$enrichNet <- netData;
  saveSet(analSet, "analSet");
  sink(netName);
  cat(rjson::toJSON(netData));
  sink();
  return(analSet);
}
