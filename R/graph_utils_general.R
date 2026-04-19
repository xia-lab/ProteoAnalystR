##################################################
## R script for ProteoAnalyst
## Description: General graph manipulation functions 
## Authors: 
## G. Zhou, guangyan.zhou@mail.mcgill.ca
## J. Xia, jeff.xia@mcgill.ca
###################################################

GetColorSchema <- function(my.grps){
  # Use consistent color palette for all group sizes
  my.grps <- as.factor(my.grps);
  grp.num <- length(levels(my.grps));

  if(grp.num > 9){
    pal12 <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
               "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
               "#FFFF99", "#B15928");
    dist.cols <- colorRampPalette(pal12)(grp.num);
    lvs <- levels(my.grps);
    colors <- vector(mode="character", length=length(my.grps));
    for(i in 1:length(lvs)){
      colors[my.grps == lvs[i]] <- dist.cols[i];
    }
  }else{
    # Use custom palette with better contrast (replaced cyan #00FFFF with magenta #FF00FF for 4th color)
    # This ensures 4th group (magenta) contrasts well with 1st group (red)
    pal9 <- c("#E31A1C", "#33A02C", "#1F78B4", "#FF00FF", "#FF7F00",
              "#6A3D9A", "#B15928", "#FFFF99", "#A6CEE3");
    lvs <- levels(my.grps);
    colors <- vector(mode="character", length=length(my.grps));
    for(i in 1:length(lvs)){
      colors[my.grps == lvs[i]] <- pal9[i];
    }
  }
  return (colors);
}

# new range [a, b]
rescale2NewRange <- function(qvec, a, b){
  q.min <- min(qvec);
  q.max <- max(qvec);
  if(length(qvec) < 50){
    a <- a*2;
  }
  if(q.max == q.min){
    new.vec <- rep(8, length(qvec));
  }else{
    coef.a <- (b-a)/(q.max-q.min);
    const.b <- b - coef.a*q.max;
    new.vec <- coef.a*qvec + const.b;
  }
  return(new.vec);
}


ComputeColorGradient <- function(nd.vec, background="black", centered, colorblind){
  require("RColorBrewer");
  
  minval <- min(nd.vec, na.rm=TRUE);
  maxval <- max(nd.vec, na.rm=TRUE);
  res <- maxval-minval;
  
  if(res == 0){
    return(rep("#FF0000", length(nd.vec)));
  }
  color <- GetColorGradient(background, centered, colorblind);
  breaks <- generate_breaks(nd.vec, length(color), center = centered);
  return(scale_vec_colours(nd.vec, col = color, breaks = breaks));
}

generate_breaks <- function(x, n, center = F){
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(seq(0, 1, length.out = n + 1))
  }
  if(center){
    m <- max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
    if (!is.finite(m) || m == 0) {
      res <- seq(0, 1, length.out = n + 1)
    } else {
      res <- seq(-m, m, length.out = n + 1)
    }
  }
  else{
    mn <- min(x, na.rm = T)
    mx <- max(x, na.rm = T)
    if (!is.finite(mn) || !is.finite(mx) || mn == mx) {
      res <- seq(0, 1, length.out = n + 1)
    } else {
      res <- seq(mn, mx, length.out = n + 1)
    }
  }
  return(res)
}

scale_vec_colours <- function(x, col = rainbow(10), breaks = NA){
  breaks <- sort(unique(breaks));
  return(col[as.numeric(cut(x, breaks = breaks, include.lowest = T))])
}

GetColorGradient <- function(background, center, colorblind=F) {
  if (background == "black") {
    if (center) {
      if (colorblind) {
        return(c(colorRampPalette(c("#3182bd", "#6baed6", "#bdd7e7"))(50), colorRampPalette(c("#FF7783", "#E32636", "#BD0313"))(50)))
      } else {
        return(c(colorRampPalette(c("#31A231", "#5BC85B", "#90EE90"))(50), colorRampPalette(c("#FF7783", "#E32636", "#BD0313"))(50)))
      }
    } else {
      if (colorblind) {
        return(c(colorRampPalette(c("#3182bd", "#bbfdff"))(50), colorRampPalette(c("#FF9DA6", "#FF7783", "#E32636", "#BD0313"))(50)))
      } else {
        return(colorRampPalette(rev(heat.colors(9)))(100))
      }
    }
  } else {
    if (center) {
      if (colorblind) {
        return(c(colorRampPalette(c("#3182bd", "#bbfdff"))(50), colorRampPalette(c("#FF9DA6", "#FF7783", "#E32636", "#BD0313"))(50)))
      } else {
        return(c(colorRampPalette(c("#137B13", "#31A231", "#5BC85B", "#90EE90"))(50), colorRampPalette(c("#FF7783", "#E32636", "#BD0313", "#96000D"))(50)))
      }
    } else {
      return(colorRampPalette(hsv(h = seq(0.72, 1, 0.035), s = 0.72, v = 1))(100))
    }
  }
}

UpdateNetworkLayout <- function(algo, filenm, focus=""){
  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");

  current.net.nm <- paramSet$current.net.nm;
  current.net <- analSet$ppi.comps[[current.net.nm]];
  pos.xy <- PerformLayOut(current.net.nm, algo, focus);
  nms <- V(current.net)$name;

  # IMPORTANT: Keep node IDs as Entrez to match the network structure
  # The network uses Entrez IDs as primary identifiers (node.id field)
  # UniProt IDs are stored separately in node.uniprot field
  # This ensures layout updates can match nodes correctly
  node.ids <- as.character(nms);

  nodes <- vector(mode="list");
  if(algo %in% c("fr", "kk")){
    for(i in 1:length(node.ids)){
      nodes[[i]] <- list(
        id=node.ids[i],
        x=pos.xy[i,1],
        y=pos.xy[i,2],
        z=pos.xy[i,3]
      );
    }
  }else{
    for(i in 1:length(node.ids)){
      nodes[[i]] <- list(
        id=node.ids[i], 
        x=pos.xy[i,1], 
        y=pos.xy[i,2]
      );
    }
  }
  # now only save the node pos to json
  netData <- list(nodes=nodes);
  sink(filenm);
  cat(RJSONIO::toJSON(netData));
  sink();
  return(filenm);
}


ExtractModule<- function(nodeids, type="enr"){
  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");

  set.seed(8574);
  nodes <- strsplit(nodeids, ";")[[1]];

  msg("[ExtractModule] Input nodeids: ", nodeids)
  msg("[ExtractModule] Parsed nodes (", length(nodes), "): ", paste(head(nodes, 10), collapse=", "))
  msg("[ExtractModule] Current network name: ", paramSet$current.net.nm)

  g <- analSet$ppi.comps[[paramSet$current.net.nm]];

  if(is.null(g)) {
    msg("[ExtractModule] ERROR: Network not found: ", paramSet$current.net.nm)
    return("NA")
  }

  msg("[ExtractModule] Network has ", vcount(g), " vertices and ", ecount(g), " edges")
  msg("[ExtractModule] Network vertex names (first 10): ", paste(head(V(g)$name, 10), collapse=", "))

  # Check if input nodes are UniProt IDs (contain letters) vs Entrez IDs (all numeric)
  # Graph vertices are stored as Entrez IDs
  first_node <- nodes[1]
  is_uniprot <- grepl("[A-Za-z]", first_node)

  if (is_uniprot) {
    msg("[ExtractModule] Input IDs appear to be UniProt format, converting to Entrez...")
    org <- if (!is.null(paramSet$data.org)) paramSet$data.org else "hsa"

    # Convert UniProt to Entrez using database
    tryCatch({
      uniprot.map <- queryGeneDB("entrez_uniprot", org)

      if (!is.null(uniprot.map) && is.data.frame(uniprot.map) &&
          "gene_id" %in% colnames(uniprot.map) && "accession" %in% colnames(uniprot.map)) {

        # Match UniProt accessions to get Entrez IDs
        hit.inx <- match(nodes, as.character(uniprot.map$accession))
        entrez.ids <- uniprot.map$gene_id[hit.inx]

        # Remove NAs
        valid.inx <- !is.na(entrez.ids)
        nodes <- as.character(entrez.ids[valid.inx])

        msg("[ExtractModule] Converted ", sum(valid.inx), " UniProt IDs to Entrez IDs")
        msg("[ExtractModule] Converted IDs (first 10): ", paste(head(nodes, 10), collapse=", "))
      } else {
        msg("[ExtractModule] ERROR: Could not load UniProt mapping database")
        return("NA")
      }
    }, error = function(e) {
      msg("[ExtractModule] ERROR: Failed to convert UniProt IDs: ", e$message)
      return("NA")
    })
  }

  # try to see if the nodes themselves are already connected
  hit.inx <- V(g)$name %in% nodes;
  msg("[ExtractModule] Found ", sum(hit.inx), " of ", length(nodes), " nodes in network")

  if(sum(hit.inx) == 0) {
    msg("[ExtractModule] ERROR: None of the requested nodes found in network")
    return("NA")
  }

  gObj <- induced_subgraph(g, V(g)$name[hit.inx]);
  msg("[ExtractModule] Induced subgraph has ", vcount(gObj), " vertices")

  # now find connected components
  comps <-decompose(gObj, min.vertices=1);
  msg("[ExtractModule] Found ", length(comps), " connected component(s)")

  if(length(comps) == 1){ # nodes are all connected
    g <- comps[[1]];
    msg("[ExtractModule] All nodes connected in single component")
  }else{
    # extract modules
    msg("[ExtractModule] Multiple components, computing shortest paths...")

    # Filter nodes to only those present in the graph
    nodes_in_graph <- nodes[nodes %in% V(g)$name]
    msg("[ExtractModule] Using ", length(nodes_in_graph), " nodes that exist in graph")

    if(length(nodes_in_graph) < 2) {
      msg("[ExtractModule] ERROR: Need at least 2 nodes in graph for path computation")
      return("NA")
    }

    paths.list <-list();
    sd.len <- length(nodes_in_graph);
    for(pos in 1:sd.len){
      from_node <- nodes_in_graph[pos]
      to_nodes <- nodes_in_graph[-(1:pos)]

      if(length(to_nodes) > 0) {
        msg("[ExtractModule] Computing paths from '", from_node, "' to ", length(to_nodes), " other nodes")
        paths.list[[pos]] <- shortest_paths(g, from_node, to_nodes)$vpath;
      }
    }
    nds.inxs <- unique(unlist(paths.list));
    msg("[ExtractModule] Found ", length(nds.inxs), " unique vertices in shortest paths")

    if(length(nds.inxs) == 0) {
      msg("[ExtractModule] ERROR: No connecting paths found between nodes")
      return("NA")
    }

    nodes2rm <- V(g)$name[-nds.inxs];
    g <- simplify(delete_vertices(g, nodes2rm));
    msg("[ExtractModule] After path extraction: ", vcount(g), " vertices, ", ecount(g), " edges")
  }
  nodeList <- igraph::as_data_frame(g, "vertices");
  msg("[ExtractModule] Extracted ", nrow(nodeList), " nodes for module")

  if(nrow(nodeList) < 3){
    msg("[ExtractModule] ERROR: Module too small (", nrow(nodeList), " nodes, need >= 3)")
    return ("NA");
  }

  if(ncol(nodeList) == 1){
    nodeList <- data.frame(Id=nodeList[,1], Label=nodeList[,1]);
  }

  paramSet$module.count <- paramSet$module.count + 1;
  module.nm <- paste("module", paramSet$module.count, sep="");
  msg("[ExtractModule] Creating module: ", module.nm)

  colnames(nodeList) <- c("Id", "Label");
  ndFileNm <- paste(module.nm, "_node_list.csv", sep="");
  fast.write(nodeList, file=ndFileNm, row.names=F);
  msg("[ExtractModule] Wrote node list to: ", ndFileNm)

  edgeList <- igraph::as_data_frame(g, "edges");
  edgeList <- cbind(rownames(edgeList), edgeList);
  colnames(edgeList) <- c("Id", "Source", "Target");
  edgFileNm <- paste(module.nm, "_edge_list.csv", sep="");
  fast.write(edgeList, file=edgFileNm, row.names=F);
  msg("[ExtractModule] Wrote edge list to: ", edgFileNm)

  # Add layout coordinates and styling attributes using same logic as compartment network
  msg("[ExtractModule] Computing layout and node attributes...")
  vc <- vcount(g)
  nms <- V(g)$name

  # Calculate network metrics
  node.dgr <- igraph::degree(g)
  node.btw <- igraph::betweenness(g, directed = FALSE, normalized = FALSE)

  # Get expression values if available
  node.expr <- rep(0, vc)
  if (!is.null(V(g)$expr)) {
    node.expr <- V(g)$expr
  }

  # Calculate node sizes based on degree (same logic as compartment network)
  if (vc > 500) {
    min.size <- 2
  } else if (vc > 200) {
    min.size <- 3
  } else {
    min.size <- 4
  }
  node.sizes <- rescale2NewRange((log(node.dgr + 1))^2, min.size, 12)

  # Compute expression-based colors (centered gradient) - store as expcolb/expcolw/expcolc
  centered <- TRUE
  exp.colsb <- ComputeColorGradient(node.expr, "black", centered, FALSE)
  exp.colsw <- ComputeColorGradient(node.expr, "white", centered, FALSE)
  exp.colsc <- ComputeColorGradient(node.expr, "colorblind", centered, TRUE)

  # Topology colors based on betweenness
  topo.val <- log(node.btw + 1)
  topo.colsb <- ComputeColorGradient(topo.val, "black", FALSE, FALSE)
  topo.colsw <- ComputeColorGradient(topo.val, "white", FALSE, FALSE)
  topo.colsc <- ComputeColorGradient(topo.val, "colorblind", FALSE, TRUE)

  # Compute layout
  if (vc > 500) {
    layout.coords <- igraph::layout_with_fr(g, niter = 100)
  } else {
    layout.coords <- igraph::layout_with_fr(g, niter = 500)
  }

  # Normalize coordinates to [0, 100] (like compartment network)
  pos.x <- (layout.coords[, 1] - min(layout.coords[, 1])) /
           (max(layout.coords[, 1]) - min(layout.coords[, 1])) * 100
  pos.y <- (layout.coords[, 2] - min(layout.coords[, 2])) /
           (max(layout.coords[, 2]) - min(layout.coords[, 2])) * 100

  # Check if compartment information was preserved from original graph
  if (!is.null(V(g)$broad_category) && !is.null(V(g)$main_location)) {
    n_with_compartment <- sum(V(g)$broad_category != "Unknown")
    msg("[ExtractModule] Preserved compartment info for ", n_with_compartment, "/", vc, " nodes")

    # Log compartment distribution
    comp_table <- table(V(g)$broad_category)
    msg("[ExtractModule] Compartment distribution: ", paste(names(comp_table), "=", comp_table, collapse=", "))
  } else {
    msg("[ExtractModule] WARNING: No compartment information found on extracted nodes")
    msg("[ExtractModule] Setting default compartment attributes...")

    # Set default compartment attributes if not present
    V(g)$broad_category <- "Unknown"
    V(g)$main_location <- "Unknown"
  }

  # Get compartment colors (matching PPI network export logic)
  category.colors <- list(
    "Nucleus" = "#e41a1c",
    "Cell surface & adhesion" = "#377eb8",
    "Cytoskeleton" = "#4daf4a",
    "Endomembrane" = "#984ea3",
    "Mitochondria & metabolic organelles" = "#ff7f00",
    "Mitochondrial & metabolic organelles" = "#17becf",
    "Cytosol" = "#a65628",
    "Extracellular" = "#f781bf",
    "Unknown" = "#999999"
  )

  # Assign compartment colors as default (like PPI network)
  comp.colors <- sapply(V(g)$broad_category, function(cat) {
    col <- category.colors[[cat]]
    if (is.null(col)) col <- "#999999"
    return(col)
  })

  # Add attributes to graph
  # Default colors are compartment colors (colorb/colorw for light/dark backgrounds)
  # Expression colors stored separately as expcolb/expcolw/expcolc for alternative coloring
  V(g)$posx <- pos.x
  V(g)$posy <- pos.y
  V(g)$size <- node.sizes
  V(g)$color <- comp.colors      # Default: compartment color
  V(g)$colorw <- comp.colors     # Default: compartment color
  V(g)$expcolb <- exp.colsb      # Alternative: expression color (black bg)
  V(g)$expcolw <- exp.colsw      # Alternative: expression color (white bg)
  V(g)$expcolc <- exp.colsc      # Alternative: expression color (colorblind)
  V(g)$topocolb <- topo.colsb    # Alternative: topology color (black bg)
  V(g)$topocolw <- topo.colsw    # Alternative: topology color (white bg)
  V(g)$topocolc <- topo.colsc    # Alternative: topology color (colorblind)
  V(g)$expr <- node.expr

  msg("[ExtractModule] Layout computed: ", vc, " nodes positioned with degree-based sizing")
  msg("[ExtractModule] Default colors set to compartment colors")

  # record the module
  analSet$ppi.comps[[module.nm]] <- g;
  msg("[ExtractModule] Stored module in analSet$ppi.comps: ", module.nm)
  msg("[ExtractModule] Module graph has ", vcount(g), " vertices and ", ecount(g), " edges")

  # Save analSet BEFORE calling UpdateSubnetStats so it can read the updated ppi.comps
  saveSet(analSet, "analSet");

  UpdateSubnetStats();

  saveSet(paramSet, "paramSet");

  msg("[ExtractModule] Calling CorrIgraph2SigmaJS to generate JSON...")
  ##if(type == "enr"){
    #convertIgraph2JSON(module.nm, filenm);
  #}else{
    # Use g directly instead of re-reading from analSet
    # Pass module.nm (without .json) - CorrIgraph2SigmaJS will add .json
    CorrIgraph2SigmaJS(g, module.nm, paramSet, analSet)
  #}

  # CorrIgraph2SigmaJS creates the file with .json extension
  filenm <- paste(module.nm, ".json", sep="");
  msg("[ExtractModule] SUCCESS: Generated JSON file: ", filenm)
  return (filenm);
}

# also save to GraphML
ExportNetwork <- function(fileName){
  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");
  current.net <- analSet$ppi.comps[[paramSet$current.net.nm]];
  write.graph(current.net, file=fileName, format="graphml");
}

##################################################
## R script for ProteoAnalyst
## Description: Perform network layout
## Author: 
## Jeff Xia, jeff.xia@mcgill.ca
## G. Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

PerformLayOut <- function(net.nm, algo, focus=""){
  analSet <- readSet(analSet, "analSet");

  g <- analSet$ppi.comps[[net.nm]];
  vc <- vcount(g);
  if(algo == "Default"){
    if(vc > 5000) {
      pos.xy <- layout_with_lgl(g);
    }else if(vc < 100){
      pos.xy <- layout_with_kk(g);
    }else{
      pos.xy <- layout_with_fr(g);
    }
  }else if(algo == "FrR"){
    # igraph API changed across versions; prefer modern signature and
    # fall back to legacy "area" argument when needed.
    pos.xy <- tryCatch(
      layout_with_fr(g, niter = 1000, grid = "nogrid"),
      error = function(e) {
        layout_with_fr(g, area = 34 * vc^2)
      }
    );
  }else if(algo == "random"){
    pos.xy <- layout_randomly (g);
  }else if(algo == "lgl"){
    pos.xy <- layout_with_lgl(g);
  }else if(algo == "gopt"){
    # this is a slow one
    if(vc > 3000) {
      maxiter <- 50;
    }else if(vc > 2000) {
      maxiter <- 100;
    }else if(vc > 1000) {
      maxiter <- 200;
    }else{
      maxiter <- 500;
    }
    pos.xy <- layout_with_graphopt(g, niter=maxiter);
  }else if(algo == "fr"){
    pos.xy <- layout_with_fr(g, dim=3, niter=500);
  }else if(algo == "kk"){
    pos.xy <- layout_with_kk(g, dim=3, maxiter=500);
  }else if(algo == "tripartite"){
    l <- layout_with_sugiyama(g, layers = V(g)$layers*(vc/4));
    pos.xy <- -l$layout[,2:1];
  }else if(algo == "concentric"){
    require(graphlayouts);
    # the fist element in the list for concentric is the central node.
    if(focus==""){
      inx <- 1;
    }else{
      inx <- which(V(g)$name == focus);
    }
    coords <- layout_with_focus(g,inx);
    pos.xy <- coords$xy;
  }else if(algo == "backbone"){
    require(graphlayouts);
    if(length(V(g)$name)<2000){
      coords <- layout_with_stress(g);
      pos.xy <- coords;
    }else{
      coords <- layout_with_sparse_stress(g,pivots=100);
      pos.xy <- coords;
    }
    
  }
  pos.xy;
}

##################################################
## Shared network helpers
##################################################

GetNetsName <- function(){
  if (exists("net.stats") && !is.null(net.stats)) {
    return(rownames(net.stats))
  }
  character(0)
}

GetNetsNodeNum <- function(){
  if (exists("net.stats") && !is.null(net.stats)) {
    return(as.numeric(net.stats$Node))
  }
  integer(0)
}

GetNetsEdgeNum <- function(){
  if (exists("net.stats") && !is.null(net.stats)) {
    return(as.numeric(net.stats$Edge))
  }
  integer(0)
}

GetNetsQueryNum <- function(){
  if (exists("net.stats") && !is.null(net.stats)) {
    return(as.numeric(net.stats$Query))
  }
  integer(0)
}

##################################################
## Cellular Localization Network Visualization
##################################################

#' PrepareLocalizationNetwork
#'
#' Generates a cytoscape.js JSON file for cellular localization network viewer
#' Uses Broad.category from localization.qs to create compound nodes representing
#' cellular compartments (Nucleus, Mitochondria, Cell surface, etc.)
#'
#' @param fileName Character, base name for output file (without .json extension)
#' @param runCompartmentLayout Logical, if TRUE compute the expensive compartment-aware FR layout;
#' if FALSE use fast non-compartment random coordinates.
#' @return Integer, 1 if successful, 0 otherwise
#'
PrepareLocalizationNetwork <- function(fileName = "localization_network",
                                       runCompartmentLayout = FALSE) {

  paramSet <- readSet(paramSet, "paramSet")
  analSet <- readSet(analSet, "analSet")

  # Get the PPI network graph - CreateGraph stores it in analSet$ppi.comps
  overall.graph <- NULL

  # Try to get from analSet$ppi.comps (preferred - where CreateGraph stores it)
  if (!is.null(analSet$ppi.comps) && length(analSet$ppi.comps) > 0) {
    # Use the current network name if available
    if (!is.null(paramSet$current.net.nm) && paramSet$current.net.nm %in% names(analSet$ppi.comps)) {
      overall.graph <- analSet$ppi.comps[[paramSet$current.net.nm]]
      #cat(sprintf("[Localization] Using network: %s\n", paramSet$current.net.nm))
    } else if (!is.null(analSet$current.compartment.net) && analSet$current.compartment.net %in% names(analSet$ppi.comps)) {
      overall.graph <- analSet$ppi.comps[[analSet$current.compartment.net]]
      #cat(sprintf("[Localization] Using compartment network: %s\n", analSet$current.compartment.net))
    } else {
      # Use the first network
      overall.graph <- analSet$ppi.comps[[1]]
      #cat(sprintf("[Localization] Using network: %s\n", names(analSet$ppi.comps)[1]))
    }
  } else if (exists("overall.graph", envir = .GlobalEnv) && !is.null(get("overall.graph", envir = .GlobalEnv))) {
    # Fallback to global variable if it exists
    overall.graph <- get("overall.graph", envir = .GlobalEnv)
    msg("[Localization] Using global overall.graph\n")
  }

  # Check if we have a network graph
  if (is.null(overall.graph)) {
    AddErrMsg("No network graph found. Please build a PPI network first.")
    return(0)
  }

  #cat(sprintf("[Localization] Network has %d nodes, %d edges\n", igraph::vcount(overall.graph), igraph::ecount(overall.graph)))
  flush.console()

  # Get organism
  org <- paramSet$data.org
  if (is.null(org) || org == "") {
    org <- "hsa"  # default to human
  }
  # print(paste("org====", org));

  # Load localization data
  loc.path <- paste0(paramSet$lib.path, org, "/", org, "_localization.qs")
  if (!file.exists(loc.path)) {
    AddErrMsg(paste0("Localization data not found for organism: ", org))
    return(0)
  }

  loc.data <- ov_qs_read(loc.path)
  id.type <- if (!is.null(paramSet$data.idType)) paramSet$data.idType else "entrez"

  # Get network nodes and edges
  require(igraph)
  original.nodes <- V(overall.graph)$name  # Keep original for UniProt mapping later
  nodes <- V(overall.graph)$name
  edges.df <- igraph::as_data_frame(overall.graph, what = "edges")

  # CRITICAL: Ensure nodes are Entrez IDs for consistent JSON export
  # The graph might have UniProt IDs, but we need to work with Entrez IDs
  # Get Entrez IDs for localization lookup and as primary node identifiers
  entrez.ids <- nodes  # Default: assume nodes are Entrez IDs

  # Check if graph has entrez attribute (for UniProt-based graphs)
  if (!is.null(V(overall.graph)$entrez)) {
    entrez.ids <- V(overall.graph)$entrez
    nodes <- entrez.ids  # CRITICAL: Update nodes to be Entrez IDs
    cat(sprintf("[Localization] Using graph entrez attribute (converted to Entrez IDs)\n"))
  } else {
    # Detect if nodes are UniProt IDs (not purely numeric)
    # UniProt IDs typically start with letters: P12345, Q9Y2Q0, etc.
    if (length(nodes) > 0 && !grepl("^[0-9]+$", nodes[1])) {
      cat("[Localization] Detected UniProt IDs, converting to Entrez...\n")

      # Convert UniProt to Entrez using database (pattern from enrich_utils.R)
      tryCatch({
        # Query UniProt → Entrez mapping
        uniprot.map <- queryGeneDB("entrez_uniprot", org)

        if (!is.null(uniprot.map) && is.data.frame(uniprot.map) &&
            "gene_id" %in% colnames(uniprot.map) && "accession" %in% colnames(uniprot.map)) {

          # Match UniProt accessions to get Entrez IDs
          hit.inx <- match(nodes, uniprot.map$accession)
          entrez.ids <- uniprot.map$gene_id[hit.inx]

          # Keep original UniProt ID where no Entrez mapping found
          na.inx <- is.na(entrez.ids)
          entrez.ids[na.inx] <- nodes[na.inx]

          # CRITICAL: Update nodes to be Entrez IDs for all downstream processing
          nodes <- entrez.ids

          mapped_count <- sum(!na.inx)
          cat(sprintf("[Localization] Converted %d/%d UniProt IDs to Entrez (nodes now use Entrez)\n",
                      mapped_count, length(nodes)))
        } else {
          cat("[Localization] Warning: UniProt to Entrez mapping table not available\n")
        }
      }, error = function(e) {
        cat(sprintf("[Localization] Error converting UniProt to Entrez: %s\n", e$message))
        cat("[Localization] Using original node IDs\n")
      })
    } else {
      cat(sprintf("[Localization] Using node names for lookup (sample: %s)\n", nodes[1]))
    }
  }

  # Map nodes to localization data (EntrezID)
  loc.map <- loc.data[match(as.character(entrez.ids), as.character(loc.data$EntrezID)), ]

  # Keep edge endpoints aligned with the node IDs used downstream.
  # When the original graph uses UniProt vertex names but localization switches
  # nodes to Entrez IDs, igraph edge data still carries the original names.
  # Remap edges through the same vertex-level ID mapping before any layout or
  # collapsed graph construction.
  vertex.id.map <- setNames(as.character(nodes), as.character(original.nodes))
  edges.df$from <- unname(vertex.id.map[as.character(edges.df$from)])
  edges.df$to <- unname(vertex.id.map[as.character(edges.df$to)])

  invalid.edge.inx <- is.na(edges.df$from) | is.na(edges.df$to) |
    !nzchar(edges.df$from) | !nzchar(edges.df$to)
  if (any(invalid.edge.inx)) {
    msg(sprintf(
      "[Localization] Dropping %d edges with endpoints not present in remapped vertices\n",
      sum(invalid.edge.inx)
    ))
    edges.df <- edges.df[!invalid.edge.inx, , drop = FALSE]
  }

  # Debug: Check how many nodes matched
  matched_count <- sum(!is.na(loc.map$EntrezID))
  cat(sprintf("[Localization] Matched %d/%d nodes to localization data\n", matched_count, length(nodes)))
  if (matched_count == 0) {
    cat("[Localization] WARNING: NO nodes matched! Sample nodes:\n")
    cat(paste("  ", head(nodes, 5), collapse="\n"))
    cat("\n")
    cat("[Localization] Entrez IDs used for lookup:\n")
    cat(paste("  ", head(as.character(entrez.ids), 5), collapse="\n"))
    cat("\n")
    cat("[Localization] Sample EntrezIDs from loc.data:\n")
    cat(paste("  ", head(as.character(loc.data$EntrezID), 5), collapse="\n"))
    cat("\n")
  }

  # Count proteins per broad category
  loc.map$Broad.category[is.na(loc.map$Broad.category)] <- "Unknown"
  category.counts <- table(loc.map$Broad.category)
  cat(sprintf("[Localization] Category distribution: %s\n",
              paste(names(category.counts), "=", category.counts, collapse=", ")))

  if (isTRUE(runCompartmentLayout)) {
    # Strategy: Add temporary hub-spoke edges within each compartment to force clustering
    # For each compartment, connect all nodes to the highest-degree node (hub)
    # This creates a star topology that clusters nodes together efficiently
    msg("[Localization] Adding temporary hub-spoke edges for layout...\n")
    flush.console()

    # Create edge list for layout (original + temporary)
    # OPTIMIZED: Collect temporary edges in list, then combine once (much faster)
    temp.edges.list <- vector("list", length(nodes))  # Overestimate size
    temp.edge.count <- 0

    # For each compartment, find hub node and connect all others to it
    for (cat in names(category.counts)) {
      comp.nodes <- nodes[loc.map$Broad.category == cat & !is.na(loc.map$Broad.category)]
      if (length(comp.nodes) <= 1) next

      # Find the node with highest degree in this compartment
      comp.degrees <- numeric(length(comp.nodes))
      for (i in seq_along(comp.nodes)) {
        node <- comp.nodes[i]
        # Count edges connected to this node
        comp.degrees[i] <- sum(edges.df$from == node | edges.df$to == node)
      }
      hub.idx <- which.max(comp.degrees)
      hub.node <- comp.nodes[hub.idx]

      # Add temporary edges from all other nodes to the hub
      for (i in seq_along(comp.nodes)) {
        if (i == hub.idx) next  # Skip hub itself

        target.node <- comp.nodes[i]

        # Check if edge already exists in original edges
        exists <- any((edges.df$from == hub.node & edges.df$to == target.node) |
                        (edges.df$from == target.node & edges.df$to == hub.node))

        if (!exists) {
          temp.edge.count <- temp.edge.count + 1
          temp.edges.list[[temp.edge.count]] <- data.frame(
            from = hub.node,
            to = target.node,
            stringsAsFactors = FALSE
          )
        }
      }
    }

    # Combine original edges with temporary edges (single rbind operation)
    if (temp.edge.count > 0) {
      temp.edges <- do.call(rbind, temp.edges.list[1:temp.edge.count])
      edges.for.layout <- edges.df[, c("from", "to"), drop = FALSE]
      layout.edges <- rbind(edges.for.layout, temp.edges)
    } else {
      layout.edges <- edges.df[, c("from", "to"), drop = FALSE]
    }

    # Create layout graph with all edges (original + temporary)
    layout.graph <- graph_from_data_frame(layout.edges, directed = FALSE, vertices = nodes)

    # Set edge weights: keep strong intra-compartment pull but avoid extreme spacing
    edge.weights <- numeric(ecount(layout.graph))
    w.temp.intra <- 80      # temporary hub edges within a compartment
    w.orig.intra <- 15      # original PPI edges within a compartment
    w.inter <- 1            # edges between compartments

    for (i in 1:ecount(layout.graph)) {
      # get.edge() is defunct in igraph >= 2.1.0; ends() returns the edge endpoints
      edge <- as.vector(igraph::ends(layout.graph, es = i, names = TRUE))
      from.idx <- match(edge[1], nodes)
      to.idx <- match(edge[2], nodes)
      from.comp <- loc.map$Broad.category[from.idx]
      to.comp <- loc.map$Broad.category[to.idx]

      # Check if this is a temporary edge (same compartment, added for clustering)
      is.original <- any((edges.df$from == edge[1] & edges.df$to == edge[2]) |
                           (edges.df$from == edge[2] & edges.df$to == edge[1]))

      if (!is.na(from.comp) && !is.na(to.comp) && from.comp == to.comp && !is.original) {
        edge.weights[i] <- w.temp.intra
      } else if (!is.na(from.comp) && !is.na(to.comp) && from.comp == to.comp) {
        edge.weights[i] <- w.orig.intra
      } else {
        edge.weights[i] <- w.inter
      }
    }
    E(layout.graph)$weight <- edge.weights

    # Compute layout using Fruchterman-Reingold with weighted edges
    msg("[Localization] Computing FR layout with temporary edges...\n")
    flush.console()

    set.seed(42)  # Reproducible layout
    layout.coords <- layout_with_fr(
      layout.graph,
      weights = E(layout.graph)$weight,
      niter = 8000,
      start.temp = sqrt(vcount(layout.graph)) * 2,
      grid = "nogrid"
    )

    msg("[Localization] Layout complete, temporary edges will not be exported\n")
    flush.console()

    # Normalize coordinates with wider aspect ratio [0, 1600] x [0, 1000] for better display
    x.coords <- layout.coords[, 1]
    y.coords <- layout.coords[, 2]

    x.min <- min(x.coords)
    x.max <- max(x.coords)
    y.min <- min(y.coords)
    y.max <- max(y.coords)

    if (is.finite(x.max - x.min) && (x.max - x.min) > 0) {
      x.coords <- (x.coords - x.min) / (x.max - x.min) * 1600
    } else {
      x.coords <- rep(800, length(nodes))
    }
    if (is.finite(y.max - y.min) && (y.max - y.min) > 0) {
      y.coords <- (y.coords - y.min) / (y.max - y.min) * 1000
    } else {
      y.coords <- rep(500, length(nodes))
    }
  } else {
    # Fast initial render: no expensive compartment-aware force layout.
    # Use igraph's layout_nicely as requested.
    set.seed(42)
    vc <- vcount(overall.graph);
    if(vc > 5000) {
      pos.xy <- layout_with_lgl(overall.graph);
    }else if(vc < 100){
      pos.xy <- layout_with_kk(overall.graph);
    }else{
      pos.xy <- layout_with_fr(overall.graph);
    }

    x.coords <- pos.xy[, 1]
    y.coords <- pos.xy[, 2]

    msg("[Localization] Skipping compartment-aware layout for initial load (using layout_nicely).\n")
    flush.console()
  }

  # Define colors for each broad category
  category.colors <- list(
    "Nucleus" = "#e41a1c",
    "Cell surface & adhesion" = "#377eb8",
    "Cytoskeleton" = "#4daf4a",
    "Endomembrane" = "#984ea3",
    "Mitochondria & metabolic organelles" = "#ff7f00",
    "Mitochondrial & metabolic organelles" = "#17becf",  # Teal/cyan color
    "Cytosol" = "#a65628",
    "Extracellular" = "#f781bf",
    "Unknown" = "#999999"
  )

  # Get expression values if available
  expr.vals <- rep(0, length(nodes))
  if (!is.null(V(overall.graph)$expr)) {
    expr.vals <- V(overall.graph)$expr
  } else if (exists("expr.vec") && !is.null(expr.vec)) {
    expr.vals <- expr.vec[match(nodes, names(expr.vec))]
    expr.vals[is.na(expr.vals)] <- 0
  }

  # Build CorrIgraph2SigmaJS-style payload -------------------------------
  degs <- igraph::degree(overall.graph)
  rescale_size <- function(x, from = 8, to = 20) {
    rng <- range(x)
    if (rng[2] - rng[1] == 0)
      return(rep((from + to) / 2, length(x)))
    (x - rng[1]) / (rng[2] - rng[1]) * (to - from) + from
  }
  sizes <- rescale_size(log10(degs + 1))

  # Map Entrez IDs to gene symbols
  symVec <- doEntrez2SymbolMapping(nodes, org, "entrez")
  initsbls <- symVec
  names(initsbls) <- nodes
  initsbls[is.na(initsbls) | !nzchar(initsbls)] <- nodes[is.na(initsbls) | !nzchar(initsbls)]

  uniprot.vec <- rep(NA_character_, length(nodes))
  has_stored_mapping <- FALSE

  # PRIORITY 1: If original graph nodes were UniProt IDs, use them directly
  if (exists("original.nodes") && length(original.nodes) > 0 &&
      !grepl("^[0-9]+$", original.nodes[1])) {
    # original.nodes are UniProt IDs, use them directly for best quality
    uniprot.vec <- as.character(original.nodes)
    has_stored_mapping <- TRUE
    msg("[Localization] Using original graph node IDs as UniProt (direct mapping)\n")
  } else if (!is.null(analSet$uniprot_to_entrez_map)) {
    msg("[Localization] Using stored UniProt→Entrez mapping from PPI network building\n")
    # Reverse the mapping: Entrez → UniProt
    entrez_to_uniprot <- names(analSet$uniprot_to_entrez_map)
    names(entrez_to_uniprot) <- as.character(analSet$uniprot_to_entrez_map)

    # Map nodes (Entrez IDs) back to UniProt
    for (i in seq_along(nodes)) {
      entrez.id <- as.character(nodes[i])
      if (entrez.id %in% names(entrez_to_uniprot)) {
        uniprot.vec[i] <- entrez_to_uniprot[entrez.id]
      }
    }

    mapped_count <- sum(!is.na(uniprot.vec))
    # cat(sprintf("[Localization] Mapped %d/%d nodes using stored mapping\n",
    #            mapped_count, length(nodes)))
    has_stored_mapping <- TRUE
  }

  # Fallback: If no stored mapping or some nodes unmapped, use database lookup
  # This can create duplicates because Entrez→UniProt is many-to-one!
  if (!has_stored_mapping || any(is.na(uniprot.vec))) {
    msg("[Localization] Warning: Using database lookup for Entrez→UniProt (may create duplicates)\n")

    tryCatch({
      # Query entrez_uniprot database table
      uniprot.map <- queryGeneDB("entrez_uniprot", org)

      if (!is.null(uniprot.map) && is.data.frame(uniprot.map) &&
          "gene_id" %in% colnames(uniprot.map) && "accession" %in% colnames(uniprot.map)) {
        # Only map unmapped nodes
        for (i in which(is.na(uniprot.vec))) {
          hit.inx <- match(as.character(nodes[i]), as.character(uniprot.map$gene_id))
          if (!is.na(hit.inx)) {
            uniprot.vec[i] <- uniprot.map$accession[hit.inx]
          }
        }

        # cat(sprintf("[Localization] Total mapped: %d/%d nodes\n",
        #             sum(!is.na(uniprot.vec)), length(nodes)))
      } else {
        cat(sprintf("[Localization] Warning: entrez_uniprot table not available for organism %s\n", org))
      }
    }, error = function(e) {
      cat(sprintf("[Localization] Warning: Could not load UniProt mapping: %s\n", e$message))
    })
  }

  base.ids <- as.character(nodes)  # Always use Entrez IDs as primary node ID

  entrez.to.finalid <- base.ids
  names(entrez.to.finalid) <- as.character(nodes)
  node.groups <- split(seq_along(nodes), base.ids)
  collapsed.nodes <- names(node.groups)
  collapsed.idx <- vapply(node.groups, function(idx) idx[1], integer(1))

  gene.nodes <- lapply(collapsed.nodes, function(base.id) {
    idx <- node.groups[[base.id]]
    i <- idx[1]
    entrez.id <- as.character(nodes[i])  # primary Entrez ID for display
    category <- loc.map$Broad.category[i]
    if (is.na(category)) category <- "Unknown"
    comp.id <- gsub("[^A-Za-z0-9_]", "_", category)
    gene.name <- loc.map$Gene.name[i]
    if (is.na(gene.name) || gene.name == "") gene.name <- entrez.id
    if (!is.na(symVec[i]) && nzchar(symVec[i])) gene.name <- symVec[i]
    main.loc <- loc.map$Main.location[i]
    if (is.na(main.loc)) main.loc <- "Unknown"
    col <- category.colors[[category]]
    if (is.null(col)) col <- "#999999"

    uniprot.values <- unique(na.omit(uniprot.vec[idx]))
    uniprot.values <- uniprot.values[nzchar(uniprot.values)]
    uniprot.primary <- if (length(uniprot.values) > 0) {
      as.character(uniprot.values[1])
    } else {
      ""  # Empty string if no UniProt mapping (base.id is now always Entrez)
    }
    uniprot.ids <- if (length(uniprot.values) > 1) {
      as.character(uniprot.values)
    } else {
      character(0)
    }

    list(
      id        = base.id,
      label     = gene.name,
      uniprot   = uniprot.primary,
      uniprot_ids = uniprot.ids,
      entrez    = entrez.id,
      size      = sizes[i],
      true_size = sizes[i],
      molType   = "gene",
      colorb    = col,
      colorw    = col,
      exp       = round(expr.vals[i], 3),
      posx      = round(x.coords[i], 2),
      posy      = round(y.coords[i], 2),
      compartment = comp.id,
      type      = "gene",
      location  = main.loc
    )
  })

  # Compartment meta-nodes
  comp.nodes <- lapply(names(category.counts), function(cat) {
    col <- category.colors[[cat]]
    if (is.null(col)) col <- "#999999"
    ids.in.cat <- which(loc.map$Broad.category == cat)
    posx <- mean(x.coords[ids.in.cat])
    posy <- mean(y.coords[ids.in.cat])
    comp.id <- gsub("[^A-Za-z0-9_]", "_", cat)
    list(
      id        = comp.id,
      label     = cat,
      size      = 0,
      true_size = 0,
      molType   = "compartment",
      colorb    = col,
      colorw    = col,
      exp       = 0,
      posx      = ifelse(is.nan(posx), 0, posx),
      posy      = ifelse(is.nan(posy), 0, posy),
      compartment = comp.id,
      type      = "compartment"
    )
  })

  # Edges with width scaling
  # Use the entrez.to.finalid mapping built during node creation
  # This ensures edges reference the exact same deduplicated node IDs
  if (!"weight" %in% colnames(edges.df)) {
    edges.df$weight <- 1
  }
  edges.df$from_final <- entrez.to.finalid[as.character(edges.df$from)]
  edges.df$to_final <- entrez.to.finalid[as.character(edges.df$to)]
  invalid.final.inx <- is.na(edges.df$from_final) | is.na(edges.df$to_final)
  if (any(invalid.final.inx)) {
    msg(sprintf(
      "[Localization] Dropping %d collapsed edges with unmapped node IDs\n",
      sum(invalid.final.inx)
    ))
    edges.df <- edges.df[!invalid.final.inx, , drop = FALSE]
  }
  edges.df <- edges.df[edges.df$from_final != edges.df$to_final, , drop = FALSE]
  if (nrow(edges.df) == 0) {
    AddErrMsg("No interactions remain after collapsing duplicate UniProt nodes.")
    return(0)
  }
  wMin <- min(edges.df$weight)
  wMax <- max(edges.df$weight)
  rescale_w <- function(w, from = 0.5, to = 2.5) {
    if (wMax == wMin) return((from + to) / 2)
    (w - wMin) / (wMax - wMin) * (to - from) + from
  }
  edges.list <- lapply(seq_len(nrow(edges.df)), function(i) {
    w <- as.numeric(edges.df$weight[i])
    # Map Entrez IDs to final deduplicated node IDs
    # Use the SAME mapping we created during node generation
    from.finalid <- as.character(edges.df$from_final[i])
    to.finalid <- as.character(edges.df$to_final[i])

    list(
      id     = paste0("e", i),
      source = from.finalid,
      target = to.finalid,
      weight = w,
      size   = rescale_w(w)
    )
  })

  # Compartments summary
  compartment.info <- lapply(names(category.counts), function(cat) {
    col <- category.colors[[cat]]
    if (is.null(col)) col <- "#999999"
    comp.id <- gsub("[^A-Za-z0-9_]", "_", cat)
    list(
      id = comp.id,
      label = cat,
      color = col,
      geneCount = as.integer(category.counts[[cat]])
    )
  })
  names(compartment.info) <- names(category.counts)

  # Node table for grid (degree/betweenness/expr/location)
  collapsed.graph <- igraph::graph_from_data_frame(
    edges.df[, c("from_final", "to_final")],
    directed = FALSE,
    vertices = collapsed.nodes
  )
  degs.collapsed <- igraph::degree(collapsed.graph)
  if (igraph::vcount(collapsed.graph) < 1000) {
    btws <- igraph::betweenness(collapsed.graph)
  } else {
    btws <- igraph::betweenness(collapsed.graph, cutoff = 3)
  }
  expr.collapsed <- expr.vals[collapsed.idx]
  loc.collapsed <- loc.map$Broad.category[collapsed.idx]
  node.table <- data.frame(
    id          = as.character(collapsed.nodes),
    label       = sapply(gene.nodes, function(n) n$label),
    uniprot     = sapply(gene.nodes, function(n) n$uniprot),
    degree      = as.numeric(degs.collapsed),
    betweenness = as.numeric(btws),
    expr        = as.numeric(expr.collapsed),
    location    = as.character(loc.collapsed),
    stringsAsFactors = FALSE
  )
  node.table <- node.table[order(node.table$degree, decreasing = TRUE), ]

  metadata <- list(
    compartmentCount = length(category.counts),
    geneCount = length(collapsed.nodes),
    edgeCount = nrow(edges.df),
    organism = org,
    compartments = compartment.info
  )

  network.json <- list(
    nodes           = c(gene.nodes),
    edges           = edges.list,
    backgroundColor = list("#f5f5f5", "#0066CC"),
    naviString      = "Compartment Network",
    org             = org,
    hasPeptideData  = HasPeptideLevelData(),
    nodeTable       = node.table,
    compartments    = compartment.info,
    metadata        = metadata
  )

  output.file <- paste0(fileName, ".json")
  flush.console()

  write_result <- tryCatch({
    jsonlite::write_json(network.json, output.file, auto_unbox = TRUE, pretty = TRUE)
    TRUE
  }, error = function(e) {
    AddErrMsg(paste0("Failed to write JSON file: ", e$message))
    FALSE
  })

  if (!write_result) {
    return(0)
  }

  # Expose node mapping for downstream routines (e.g., communities)
  ppi.net <- list()
  ppi.net[["node.data"]] <- data.frame(
    Id = as.character(collapsed.nodes),
    Label = sapply(gene.nodes, function(n) n$label),
    stringsAsFactors = FALSE
  )
  ppi.net <<- ppi.net


  return(1)
}
