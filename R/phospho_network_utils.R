##################################################
## R script for ProteoAnalyst
## Description: Kinase-substrate network visualization (ClueGO-style)
## Creates hub-and-spoke networks where kinases are hubs and substrates are spokes
## Authors:
## Jeff Xia, jeff.xia@mcgill.ca
## Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

#' Initialize kinase-substrate network for visualization
#' Creates a ClueGO-style network where kinases are hubs and phosphosites are spokes
#'
#' @param dataName Dataset name
#' @param database Kinase database to use
#' @param fdr_cutoff FDR cutoff for significant kinases (default: 0.05)
#' @return 1 on success, 0 on failure
InitKinaseSubstrateNetwork <- function(dataName, database = "phosphositeplus", fdr_cutoff = 0.05) {
  dataSet <- readDataset(dataName)
  paramSet <- readSet(paramSet, "paramSet")
  analSet <- readSet(analSet, "analSet")
  msgSet <- readSet(msgSet, "msgSet")

  msg("[InitKinaseSubstrateNetwork] Starting kinase-substrate network initialization...")
  msg("[InitKinaseSubstrateNetwork] Database: ", database)
  msg("[InitKinaseSubstrateNetwork] FDR cutoff: ", fdr_cutoff)

  # Check if kinase enrichment results exist (prefer in-memory analSet, fallback to file)
  kinase_results <- NULL
  if (!is.null(analSet$kinase.enrich) && nrow(analSet$kinase.enrich) > 0) {
    kinase_results <- analSet$kinase.enrich
    msg("[InitKinaseSubstrateNetwork] Loaded ", nrow(kinase_results), " kinase enrichment results from analSet")
  } else {
    kinase_results_file <- "kinase_enrichment_results.csv"
    if (file.exists(kinase_results_file)) {
      kinase_results <- read.csv(kinase_results_file, stringsAsFactors = FALSE)
      msg("[InitKinaseSubstrateNetwork] Loaded ", nrow(kinase_results), " kinase enrichment results from file")
    } else {
      msgSet$current.msg <- "No kinase enrichment results found. Please run kinase enrichment analysis first."
      saveSet(msgSet, "msgSet")
      return(0)
    }
  }

  # Filter significant kinases
  sig_kinases <- kinase_results[kinase_results$FDR <= fdr_cutoff, ]
  msg("[InitKinaseSubstrateNetwork] Found ", nrow(sig_kinases), " significant kinases (FDR <= ", fdr_cutoff, ")")

  # If too few significant kinases (< 3), use top N approach to ensure meaningful network
  min_kinases <- 3
  if (nrow(sig_kinases) < min_kinases) {
    # If fewer than minimum significant kinases, take top 10 by FDR/p-value
    msg("[InitKinaseSubstrateNetwork] Only ", nrow(sig_kinases), " significant kinases at FDR <= ", fdr_cutoff,
            ". Using top 10 kinases by FDR instead to create meaningful network.")

    # Sort by FDR (or p-value if FDR not available)
    if ("FDR" %in% colnames(kinase_results)) {
      kinase_results <- kinase_results[order(kinase_results$FDR), ]
    } else if ("P_value" %in% colnames(kinase_results)) {
      kinase_results <- kinase_results[order(kinase_results$P_value), ]
    }

    # Take top 10 (or fewer if less than 10 total)
    n_kinases <- min(10, nrow(kinase_results))
    sig_kinases <- kinase_results[1:n_kinases, ]

    msgSet$current.msg <- paste0("Only ", nrow(sig_kinases), " significant kinases at FDR <= ", fdr_cutoff,
                                  ". Showing top ", n_kinases, " kinases by FDR.")
    saveSet(msgSet, "msgSet")

    msg("[InitKinaseSubstrateNetwork] Selected top ", n_kinases, " kinases (FDR range: ",
            round(min(sig_kinases$FDR), 4), " - ", round(max(sig_kinases$FDR), 4), ")")
  }

  # Load kinase-substrate database
  ks_db <- .loadKinaseSubstrateDB(database, paramSet)
  db_empty <- (is.data.frame(ks_db) && nrow(ks_db) == 0) || (is.list(ks_db) && length(ks_db) == 0)
  if (is.null(ks_db) || db_empty) {
    msgSet$current.msg <- "Failed to load kinase-substrate database"
    saveSet(msgSet, "msgSet")
    return(0)
  }
  db_size <- if (is.data.frame(ks_db)) nrow(ks_db) else length(ks_db)
  msg("[InitKinaseSubstrateNetwork] Loaded kinase-substrate database: ", db_size, " relationships")

  # Get user's phosphosite data
  if (is.null(dataSet$sig.mat) || nrow(dataSet$sig.mat) == 0) {
    msgSet$current.msg <- "No significant phosphosites found. Please run differential expression first."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  sig_sites <- rownames(dataSet$sig.mat)
  all_sites <- rownames(dataSet$data.norm)
  msg("[InitKinaseSubstrateNetwork] User has ", length(sig_sites), " significant sites, ",
          length(all_sites), " total sites")

  # Build kinase-substrate network
  network_data <- .buildKinaseSubstrateNetwork(
    sig_kinases = sig_kinases,
    ks_db = ks_db,
    sig_sites = sig_sites,
    all_sites = all_sites,
    dataSet = dataSet,
    paramSet = paramSet
  )

  if (is.null(network_data) || network_data$node_count == 0) {
    msgSet$current.msg <- "Failed to build kinase-substrate network. No matching substrates found."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Store network data in analSet
  analSet$kinase.network <- network_data

  # Also store kinase enrichment results for reference
  analSet$kinase.network$enrichment_results <- sig_kinases

  saveSet(analSet, "analSet")

  msg("[InitKinaseSubstrateNetwork] Network created: ",
          network_data$node_count, " nodes, ", network_data$edge_count, " edges")
  msgSet$current.msg <- paste0("Kinase-substrate network initialized: ",
                                network_data$kinase_count, " kinases, ",
                                network_data$substrate_count, " substrates")
  saveSet(msgSet, "msgSet")

  return(1)
}

#' Build kinase-substrate network structure
#' @keywords internal
.buildKinaseSubstrateNetwork <- function(sig_kinases, ks_db, sig_sites, all_sites, dataSet, paramSet) {

  # Convert site IDs to substrate format for matching
  # User data: "Q80XQ2_S_578" -> need to match with database format
  site_info <- .parsePhosphositeIDs(sig_sites)
  all_site_info <- .parsePhosphositeIDs(all_sites)

  # Load phospho_symbol_map if available for gene symbols
  phospho_map <- NULL
  if (file.exists("phospho_symbol_map.qs")) {
    phospho_map <- readDataQs("phospho_symbol_map.qs", paramSet$anal.type, paramSet$dataName)
  }

  # Initialize network lists
  kinase_nodes <- list()
  substrate_nodes <- list()
  edges <- list()
  edge_id <- 1

  # Store kinase-to-substrates mapping (like hits.query in enrichment network)
  kinase_substrates_map <- list()

  # For each significant kinase
  for (i in 1:nrow(sig_kinases)) {
    kinase_name <- sig_kinases$Kinase[i]
    kinase_fdr <- sig_kinases$FDR[i]
    kinase_pval <- sig_kinases$P_value[i]
    kinase_fold <- sig_kinases$Fold_Enrichment[i]

    # Find substrates for this kinase in the database
    # ks_db can be a list (kinase -> substrates) or data frame
    kinase_substrate_ids <- NULL
    if (is.list(ks_db) && !is.data.frame(ks_db)) {
      # List format: kinase names are keys, substrate IDs are values
      if (kinase_name %in% names(ks_db)) {
        kinase_substrate_ids <- ks_db[[kinase_name]]
      }
    } else if (is.data.frame(ks_db)) {
      # Data frame format with Kinase and Substrate columns
      if ("Kinase" %in% colnames(ks_db) && "Substrate" %in% colnames(ks_db)) {
        kinase_substrate_ids <- ks_db$Substrate[ks_db$Kinase == kinase_name]
      }
    }

    if (is.null(kinase_substrate_ids) || length(kinase_substrate_ids) == 0) next

    # Match database substrates to user's phosphosite data
    matched_substrates <- .matchSubstrateIDsToSites(
      substrate_ids = kinase_substrate_ids,
      site_info = site_info,
      all_site_info = all_site_info,
      sig_sites = sig_sites,
      all_sites = all_sites
    )

    if (length(matched_substrates) == 0) next

    # Store kinase-to-substrates mapping
    kinase_substrates_map[[kinase_name]] <- matched_substrates

    # Add kinase node (only if it has matched substrates)
    kinase_node <- list(
      id = paste0("kinase_", kinase_name),
      label = kinase_name,
      type = "kinase",
      fdr = kinase_fdr,
      pval = kinase_pval,
      fold_enrichment = kinase_fold,
      substrate_count = length(matched_substrates),
      # Color by fold enrichment (for ClueGO-style visualization)
      score = log2(kinase_fold),  # Log2 fold enrichment for color scale
      significant = kinase_fdr <= 0.05
    )
    kinase_nodes[[length(kinase_nodes) + 1]] <- kinase_node

    # Add substrate nodes and edges
    for (substrate_id in matched_substrates) {
      # Check if substrate already added
      substrate_exists <- any(sapply(substrate_nodes, function(n) n$id == substrate_id))

      if (!substrate_exists) {
        # Get substrate data
        substrate_fc <- NA
        substrate_pval <- NA
        substrate_symbol <- substrate_id

        if (substrate_id %in% rownames(dataSet$sig.mat)) {
          if ("logFC" %in% colnames(dataSet$sig.mat)) {
            substrate_fc <- dataSet$sig.mat[substrate_id, "logFC"]
          }
          if ("adj.P.Val" %in% colnames(dataSet$sig.mat)) {
            substrate_pval <- dataSet$sig.mat[substrate_id, "adj.P.Val"]
          }
        }

        # Get gene symbol from phospho_map if available
        if (!is.null(phospho_map) && substrate_id %in% rownames(phospho_map)) {
          substrate_symbol <- phospho_map[substrate_id, "symbol"]
        }

        substrate_node <- list(
          id = substrate_id,
          label = substrate_symbol,
          type = "substrate",
          phosphosite = substrate_id,
          logFC = substrate_fc,
          pval = substrate_pval,
          score = substrate_fc,  # Use log2FC for color scale
          significant = !is.na(substrate_pval) && substrate_pval <= 0.05
        )
        substrate_nodes[[length(substrate_nodes) + 1]] <- substrate_node
      }

      # Add edge from kinase to substrate
      edge <- list(
        id = paste0("edge_", edge_id),
        source = paste0("kinase_", kinase_name),
        target = substrate_id,
        type = "kinase_substrate"
      )
      edges[[length(edges) + 1]] <- edge
      edge_id <- edge_id + 1
    }
  }

  # Combine all nodes
  all_nodes <- c(kinase_nodes, substrate_nodes)

  msg("[.buildKinaseSubstrateNetwork] Network summary:")
  msg("  Kinases: ", length(kinase_nodes))
  msg("  Substrates: ", length(substrate_nodes))
  msg("  Edges: ", length(edges))

  return(list(
    nodes = all_nodes,
    edges = edges,
    node_count = length(all_nodes),
    edge_count = length(edges),
    kinase_count = length(kinase_nodes),
    substrate_count = length(substrate_nodes),
    kinase_nodes = kinase_nodes,
    substrate_nodes = substrate_nodes,
    kinase_substrates_map = kinase_substrates_map,  # For hits.query equivalent
    sig_sites = sig_sites  # All significant phosphosites
  ))
}

#' Parse phosphosite IDs into components
#' @keywords internal
.parsePhosphositeIDs <- function(site_ids) {
  # Parse IDs like "Q80XQ2_S_578" or "P12345_S123"
  parsed <- data.frame(
    site_id = site_ids,
    uniprot = sapply(strsplit(site_ids, "_"), function(x) x[1]),
    stringsAsFactors = FALSE
  )

  # Try to extract residue and position
  parsed$residue <- sapply(strsplit(site_ids, "_"), function(x) {
    if (length(x) >= 2) {
      # Could be "S" or "S_578" format
      res <- x[2]
      if (nchar(res) == 1) return(res)  # Single letter residue
      return(substring(res, 1, 1))  # Extract first letter
    }
    return(NA)
  })

  parsed$position <- sapply(strsplit(site_ids, "_"), function(x) {
    if (length(x) == 3) return(x[3])  # "Q80XQ2_S_578" format
    if (length(x) == 2) {
      # "P12345_S123" format - extract numbers
      num <- gsub("^[A-Z]+", "", x[2])
      if (nchar(num) > 0) return(num)
    }
    return(NA)
  })

  return(parsed)
}

#' Match kinase substrate IDs from database to user's phosphosite data
#' Substrate IDs from database are typically in format: "UNIPROT;SITE" (e.g., "Q80XQ2;S578")
#' @keywords internal
.matchSubstrateIDsToSites <- function(substrate_ids, site_info, all_site_info, sig_sites, all_sites) {
  if (length(substrate_ids) == 0) return(character(0))

  # Vectorized parsing of substrate IDs
  # Parse all at once instead of looping
  parsed_list <- strsplit(substrate_ids, ";")

  # Extract uniprot and site for each substrate
  substrate_uniprots <- character(length(substrate_ids))
  substrate_sites <- character(length(substrate_ids))

  for (i in seq_along(parsed_list)) {
    parts <- parsed_list[[i]]
    if (length(parts) == 2) {
      substrate_uniprots[i] <- parts[1]
      substrate_sites[i] <- parts[2]
    } else {
      # Try underscore separator
      parts <- strsplit(substrate_ids[i], "_")[[1]]
      if (length(parts) >= 2) {
        substrate_uniprots[i] <- parts[1]
        substrate_sites[i] <- paste(parts[-1], collapse = "_")
      }
    }
  }

  # Remove hyphens and extract residue/position vectorized
  sites_clean <- gsub("-", "", substrate_sites)
  residues <- substring(sites_clean, 1, 1)
  positions <- gsub("^[A-Z]+", "", sites_clean)

  # Create lookup key from user's data for fast matching
  # Format: "UNIPROT_RESIDUE_POSITION"
  user_keys <- paste0(site_info$uniprot, "_", site_info$residue, "_", site_info$position)
  db_keys <- paste0(substrate_uniprots, "_", residues, "_", positions)

  # Vectorized matching
  matched_indices <- match(db_keys, user_keys)
  matched_indices <- matched_indices[!is.na(matched_indices)]

  matched <- character(0)
  if (length(matched_indices) > 0) {
    matched <- site_info$site_id[matched_indices]
  }

  # For unmatched, try pattern matching (slower, only on remaining)
  unmatched_mask <- is.na(match(db_keys, user_keys))
  if (any(unmatched_mask)) {
    unmatched_uniprots <- substrate_uniprots[unmatched_mask]
    unmatched_residues <- residues[unmatched_mask]
    unmatched_positions <- positions[unmatched_mask]

    # Build all patterns at once
    patterns <- paste0("^", unmatched_uniprots, "_", unmatched_residues, "(_|)", unmatched_positions, "$")
    pattern_str <- paste(patterns, collapse = "|")

    if (nchar(pattern_str) > 0 && nchar(pattern_str) < 10000) {  # Avoid regex that's too long
      loose_matches <- grep(pattern_str, sig_sites, value = TRUE)
      matched <- c(matched, loose_matches)
    }
  }

  return(unique(matched))
}

#' Generate kinase-substrate network visualization in enrichment network JSON format
#'
#' @param dataName Dataset name
#' @param netNm Network name for output files
#' @param fdr_cutoff FDR cutoff for kinase filtering
#' @return 1 on success, 0 on failure
PerformKinaseSubstrateNetworkView <- function(dataName = "", netNm = "kinase_substrate_net", fdr_cutoff = 0.05) {
  dataSet <- readDataset(dataName)
  paramSet <- readSet(paramSet, "paramSet")
  analSet <- readSet(analSet, "analSet")
  msgSet <- readSet(msgSet, "msgSet")

  msg("[PerformKinaseSubstrateNetworkView] Generating network view...")

  # Check if network data exists
  if (is.null(analSet$kinase.network)) {
    msgSet$current.msg <- "Kinase-substrate network not initialized. Please run InitKinaseSubstrateNetwork first."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  network_data <- analSet$kinase.network

  # Generate JSON in enrichment network format
  json_data <- .exportKinaseNetworkToJSON(
    network_data = network_data,
    netNm = netNm,
    paramSet = paramSet
  )

  if (is.null(json_data)) {
    msgSet$current.msg <- "Failed to export network to JSON"
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Write JSON file
  json_file <- paste0(netNm, ".json")
  con <- file(json_file, open = "w", encoding = "UTF-8")
  writeLines(json_data, con = con)
  close(con)
  msg("[PerformKinaseSubstrateNetworkView] Wrote network JSON: ", json_file)

  # Generate summary table
  .generateKinaseNetworkTable(network_data, netNm)

  msgSet$current.msg <- paste0("Kinase-substrate network generated: ",
                                network_data$kinase_count, " kinases, ",
                                network_data$substrate_count, " substrates")
  saveSet(msgSet, "msgSet")

  return(1)
}

#' Export kinase-substrate network to JSON format (same as enrichment network)
#' @keywords internal
.exportKinaseNetworkToJSON <- function(network_data, netNm, paramSet) {
  require(rjson)
  require(igraph)

  # Get enrichment results
  enrichment_results <- network_data$enrichment_results
  kinase_substrates_map <- network_data$kinase_substrates_map
  sig_sites <- network_data$sig_sites

  # 1. Build kinase-kinase graph based on shared substrates (like pathway overlap)
  kinase_nodes_list <- network_data$kinase_nodes
  kinase_names <- sapply(kinase_nodes_list, function(n) gsub("^kinase_", "", n$id))

  # Calculate kinase-kinase similarity based on shared substrates
  n_kinases <- length(kinase_names)
  w <- matrix(NA, nrow=n_kinases, ncol=n_kinases)
  colnames(w) <- rownames(w) <- kinase_names

  for (i in 1:n_kinases) {
    for (j in i:n_kinases) {
      if (i == j) {
        w[i,j] <- 1
      } else {
        # Calculate Jaccard similarity of substrates
        ki <- kinase_names[i]
        kj <- kinase_names[j]
        subs_i <- kinase_substrates_map[[ki]]
        subs_j <- kinase_substrates_map[[kj]]
        if (length(subs_i) > 0 && length(subs_j) > 0) {
          intersection <- length(intersect(subs_i, subs_j))
          union <- length(union(subs_i, subs_j))
          w[i,j] <- intersection / union
          w[j,i] <- w[i,j]
        }
      }
    }
  }

  # Create igraph from similarity matrix (filter edges with similarity >= 0.3)
  wd <- reshape::melt(w)
  wd <- wd[wd[,1] != wd[,2],]
  wd <- wd[!is.na(wd[,3]),]

  # Only keep edges with similarity >= 0.3
  wd <- wd[wd[,3] >= 0.3, ]

  # Create graph - this will only include nodes that appear in edges
  # Need to explicitly add all kinases as vertices (including isolated ones)
  msg("[DEBUG] Total kinases before graph: ", length(kinase_names))
  msg("[DEBUG] Edges after filtering (>=0.3): ", nrow(wd))
  if (nrow(wd) > 0) {
    msg("[DEBUG] Sample edges: ", paste(head(wd[,1]), "->", head(wd[,2]), collapse="; "))
    g <- graph_from_data_frame(wd[,-3], directed=F, vertices=kinase_names)
  } else {
    # No edges pass threshold - create graph with only vertices
    msg("[DEBUG] No edges passed threshold - creating empty graph with vertices")
    g <- make_empty_graph(n=0, directed=F)
    g <- add_vertices(g, length(kinase_names), name=kinase_names)
  }

  msg("[DEBUG] Graph vertices count: ", vcount(g))
  msg("[DEBUG] Graph edges count: ", ecount(g))

  # Layout kinase network
  pos.xy <- layout_nicely(g)

  # 2. Create nodes list (kinase nodes with enrichment data)
  nodes <- vector(mode="list")
  node.nms <- V(g)$name
  msg("[DEBUG] node.nms (kinases in graph): ", paste(head(node.nms, 10), collapse=", "))

  # Rescale function for node sizes
  my.rescale <- function(x, from, to){
    if (length(x) == 0 || all(is.na(x))) return(rep(from, length(x)))
    x_range <- max(x, na.rm=TRUE) - min(x, na.rm=TRUE)
    if (x_range == 0) return(rep(from, length(x)))
    (x - min(x, na.rm=TRUE)) / x_range * (to - from) + from
  }

  # Get node sizes based on substrate counts
  node.sizes <- sapply(node.nms, function(kn) {
    kinase_node <- kinase_nodes_list[[which(kinase_names == kn)]]
    kinase_node$substrate_count
  })
  V(g)$size <- my.rescale(log(node.sizes+1, base=10), 8, 32)

  # Get node colors based on fold enrichment (log2 scores already stored in kinase nodes)
  scores <- sapply(node.nms, function(kn) {
    kinase_node <- kinase_nodes_list[[which(kinase_names == kn)]]
    if (!is.null(kinase_node$score) && is.finite(kinase_node$score)) {
      return(kinase_node$score)
    }
    log2(kinase_node$fold_enrichment)
  })
  scores[!is.finite(scores)] <- 0
  V(g)$color <- ComputeColorGradient(scores, "black", F, F)
  V(g)$colorw <- ComputeColorGradient(scores, "white", F, F)

  # Create nodes list
  for(i in 1:length(V(g)$name)){
    nodes[[i]] <- list(
      id = node.nms[i],
      label = node.nms[i],
      size = V(g)$size[i],
      true_size = V(g)$size[i],
      molType = "kinase",
      colorb = V(g)$color[i],
      colorw = V(g)$colorw[i],
      posx = pos.xy[i,1],
      posy = pos.xy[i,2]
    )
  }

  # 3. Create edges between kinases (pathway-pathway connections)
  edge.mat <- as_edgelist(g)
  msg("[DEBUG] Edges from graph (edge.mat): ", nrow(edge.mat), " edges")
  if (nrow(edge.mat) > 0) {
    msg("[DEBUG] Sample edge.mat: ", paste(head(edge.mat[,1]), "->", head(edge.mat[,2]), collapse="; "))
  }
  edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2])

  # 4. Create background graph (kinase-substrate connections)
  b.mat <- do.call(rbind, lapply(names(kinase_substrates_map), function(kname) {
    substrates <- kinase_substrates_map[[kname]]
    if (length(substrates) > 0) {
      cbind(source = substrates, target = kname)
    }
  }))

  bg <- graph_from_data_frame(b.mat, directed=F)

  # Color background nodes (substrates)
  V(bg)$color <- rep("#00FFFF", length(V(bg)))  # Default cyan for substrates
  V(bg)$colorw <- rep("#668B8B", length(V(bg)))

  # Color kinase nodes in background graph
  # Color kinase nodes in background graph by fold enrichment scores
  V(bg)$color[V(bg)$name %in% kinase_names] <- ComputeColorGradient(scores, "black", F, F)
  V(bg)$colorw[V(bg)$name %in% kinase_names] <- ComputeColorGradient(scores, "white", F, F)

  # Try to color substrate nodes by logFC if available
  analSet <- readSet(analSet, "analSet")
  if (!is.null(analSet$kinase.network$substrate_nodes)) {
    substrate_nodes <- analSet$kinase.network$substrate_nodes
    substrate_ids <- V(bg)$name[!V(bg)$name %in% kinase_names]

    expvals <- rep(0, length(substrate_ids))
    names(expvals) <- substrate_ids

    for (sub_node in substrate_nodes) {
      if (sub_node$id %in% substrate_ids && !is.na(sub_node$logFC)) {
        expvals[sub_node$id] <- sub_node$logFC
      }
    }

    V(bg)$color[!V(bg)$name %in% kinase_names] <- ComputeColorGradient(unname(expvals), "black", T, T)
    V(bg)$colorw[!V(bg)$name %in% kinase_names] <- ComputeColorGradient(unname(expvals), "black", T, T)
  }

  # Size background nodes
  node.dgr2 <- as.numeric(degree(bg))
  V(bg)$size <- my.rescale(log(node.dgr2, base=10), 8, 24)

  # Layout background graph
  pos.xy.bg <- layout_nicely(bg)

  # 5. Create bnodes (background nodes - all nodes including kinases and substrates)
  bnodes <- vector(mode="list")
  bg.node.nms <- V(bg)$name
  node.sizes.bg <- V(bg)$size
  node.cols.bg <- V(bg)$color
  node.colsw.bg <- V(bg)$colorw

  # Determine node shapes
  shapes <- rep("circle", length(bg.node.nms))
  shapes[bg.node.nms %in% kinase_names] <- "kinase"
  shapes[bg.node.nms %in% sig_sites] <- "gene"  # Phosphosites shown as genes

  # Get gene symbols for labels
  phospho_map <- NULL
  if (file.exists("phospho_symbol_map.qs")) {
    phospho_map <- readDataQs("phospho_symbol_map.qs", paramSet$anal.type, paramSet$dataName)
  }

  # Create labels with gene symbol and site info (e.g., "Smarca4;pS695")
  node.lbls <- sapply(bg.node.nms, function(nm) {
    if (nm %in% kinase_names) {
      return(nm)  # Kinases use their name
    } else {
      # For phosphosites, format as "Symbol;pRESIDUEPOSITION"
      # Extract residue and position from phosphosite ID (e.g., "Q3TKT4_S_695")
      parts <- strsplit(nm, "_")[[1]]

      if (length(parts) >= 3) {
        residue <- parts[2]
        position <- parts[3]
        site_str <- paste0("p", residue, position)  # e.g., "pS695"

        if (!is.null(phospho_map) && nm %in% rownames(phospho_map)) {
          symbol <- phospho_map[nm, "symbol"]
          return(paste0(symbol, ";", site_str))  # e.g., "Smarca4;pS695"
        } else {
          return(site_str)  # Fallback to just site if no symbol
        }
      } else {
        return(nm)  # Fallback to full ID if parsing fails
      }
    }
  }, USE.NAMES = FALSE)  # IMPORTANT: Don't create named vector

  # Get expression values for substrates
  expvals.bg <- rep(0, length(bg.node.nms))
  names(expvals.bg) <- bg.node.nms
  if (!is.null(analSet$kinase.network$substrate_nodes)) {
    for (sub_node in substrate_nodes) {
      if (sub_node$id %in% bg.node.nms && !is.na(sub_node$logFC)) {
        expvals.bg[sub_node$id] <- sub_node$logFC
      }
    }
  }

  for(i in 1:length(node.sizes.bg)){
    bnodes[[i]] <- list(
      id = bg.node.nms[i],
      label = node.lbls[i],
      size = node.sizes.bg[i],
      colorb = node.cols.bg[i],
      colorw = node.colsw.bg[i],
      true_size = node.sizes.bg[i],
      molType = shapes[i],
      exp = unname(expvals.bg[bg.node.nms[i]]),
      posx = pos.xy.bg[i,1],
      posy = pos.xy.bg[i,2]
    )
  }

  # 6. Create bedges (background edges - kinase to substrate)
  bedge.mat <- as_edgelist(bg)
  bedge.mat <- cbind(id=paste0("b", 1:nrow(bedge.mat)), source=bedge.mat[,1], target=bedge.mat[,2])

  # 7. Create proteinlist (all significant phosphosites with symbols)
  proteinlist <- rep(NA, length(sig_sites))
  names(proteinlist) <- sig_sites
  if (!is.null(phospho_map)) {
    for (site_id in sig_sites) {
      if (site_id %in% rownames(phospho_map)) {
        proteinlist[site_id] <- phospho_map[site_id, "symbol"]
      } else {
        proteinlist[site_id] <- site_id
      }
    }
  } else {
    proteinlist <- sig_sites
    names(proteinlist) <- sig_sites
  }

  # 8. Create enr matrix (enrichment results for each kinase)
  # IMPORTANT: Only include kinases that are in the final graph (node.nms)
  # Some kinases may be filtered out if they have no connections after Jaccard filtering
  # Note: enrichment_results has columns: Kinase, Substrates_Total, Substrates_Sig, P_value, Fold_Enrichment, FDR
  # Need to match enrichment network format: Total, Expected, Hits, Pval, FDR
  # Plus add Fold_Enrichment for kinase-specific display
  kinases_in_graph <- node.nms  # Only kinases that remain in the graph

  enr.mat <- enrichment_results[enrichment_results$Kinase %in% kinases_in_graph,
                                c("Substrates_Total", "Substrates_Sig", "P_value", "FDR", "Fold_Enrichment"),
                                drop = FALSE]
  rownames(enr.mat) <- enrichment_results$Kinase[enrichment_results$Kinase %in% kinases_in_graph]

  # Reorder to match kinases_in_graph order
  enr.mat <- enr.mat[kinases_in_graph, , drop = FALSE]

  # Add Expected column (for kinase enrichment, this is background rate * total sig sites)
  # Expected = (Total substrates in kinase / Total background sites) * Total significant sites
  # Get total background sites from dataSet
  dataSet <- readDataset(paramSet$dataName)
  total_all_sites <- nrow(dataSet$data.norm)  # Total phosphosites in background
  total_sig_sites <- length(sig_sites)  # Total significant phosphosites
  enr.mat$Expected <- (enr.mat$Substrates_Total / total_all_sites) * total_sig_sites

  # Reorder and rename columns to match enrichment network format exactly
  # Include Fold_Enrichment for kinase network display
  enr.mat <- enr.mat[, c("Substrates_Total", "Expected", "Substrates_Sig", "P_value", "FDR", "Fold_Enrichment")]
  colnames(enr.mat) <- c("Total", "Expected", "Hits", "Pval", "FDR", "Fold_Enrichment")

  # 9. Get sizes (number of substrates per kinase)
  # Only for kinases in the graph
  sizes <- sapply(kinases_in_graph, function(kn) {
    length(kinase_substrates_map[[kn]])
  })
  names(sizes) <- kinases_in_graph

  # Convert to lists for rjson
  edge.mat <- apply(edge.mat, 1, as.list)
  bedge.mat <- apply(bedge.mat, 1, as.list)

  # IMPORTANT: Store names before converting to list, then use unname()
  # This ensures 'id' field gets the kinase names and 'enr' is unnamed array
  enr.mat.names <- rownames(enr.mat)
  enr.mat <- apply(enr.mat, 1, as.list)
  names(enr.mat) <- enr.mat.names  # Restore names for id field

  # Filter hits mapping to only include kinases in the graph
  hits_in_graph <- kinase_substrates_map[kinases_in_graph]

  msg("[DEBUG] Final data structure:")
  msg("[DEBUG]   nodes count: ", length(nodes))
  msg("[DEBUG]   edges count: ", length(edge.mat))
  msg("[DEBUG]   bnodes count: ", length(bnodes))
  msg("[DEBUG]   bedges count: ", length(bedge.mat))
  msg("[DEBUG]   enr count: ", length(enr.mat))
  msg("[DEBUG]   id count: ", length(names(enr.mat)))
  msg("[DEBUG]   hits count: ", length(hits_in_graph))
  msg("[DEBUG]   id values: ", paste(head(names(enr.mat), 10), collapse=", "))

  # 10. Create final netData structure (EXACTLY like enrichment network)
  netData <- list(
    nodes = nodes,
    edges = edge.mat,
    bnodes = bnodes,
    bedges = bedge.mat,
    enr = unname(enr.mat),
    id = names(enr.mat),
    sizes = sizes,
    hits = hits_in_graph,
    proteinlist = proteinlist,
    analType = paramSet$anal.type,
    org = paramSet$data.org,
    backgroundColor = list("#514F6A", "#222222"),
    dat.opt = paramSet$selDataNm,
    naviString = "Kinase-Substrate Network"
  )

  # Convert to JSON using rjson (same as enrichment network)
  json_str <- rjson::toJSON(netData)

  msg("[DEBUG] JSON generated, length: ", nchar(json_str))

  # Save the kinase-kinase graph to analSet for layout updates
  # This allows UpdateNetworkLayout() to access the graph
  analSet <- readSet(analSet, "analSet")
  if (is.null(analSet$ppi.comps)) {
    analSet$ppi.comps <- list()
  }
  analSet$ppi.comps[[netNm]] <- g
  saveSet(analSet, "analSet")

  # Also update paramSet with current network name
  paramSet <- readSet(paramSet, "paramSet")
  paramSet$current.net.nm <- netNm
  saveSet(paramSet, "paramSet")

  return(json_str)
}

#' Get node color based on score and type
#' @keywords internal
.getNodeColor <- function(node) {
  if (node$type == "kinase") {
    # Color kinases by fold enrichment (score = log2 fold)
    score <- ifelse(is.null(node$score), 0, node$score)
    if (score > 1) {
      return("#d73027")  # Red for highly enriched
    } else if (score > 0.5) {
      return("#fc8d59")  # Orange for moderately enriched
    } else {
      return("#fee090")  # Yellow for weakly enriched
    }
  } else {
    # Color substrates by logFC
    score <- ifelse(is.null(node$score) || is.na(node$score), 0, node$score)
    if (score > 1) {
      return("#4575b4")  # Blue for upregulated
    } else if (score < -1) {
      return("#91bfdb")  # Light blue for downregulated
    } else {
      return("#e0f3f8")  # Very light blue for not changed
    }
  }
}

#' Generate summary table for kinase-substrate network
#' @keywords internal
.generateKinaseNetworkTable <- function(network_data, netNm) {
  # Create summary table with kinase information
  kinase_nodes <- Filter(function(n) n$type == "kinase", network_data$nodes)

  kinase_table <- do.call(rbind, lapply(kinase_nodes, function(n) {
    data.frame(
      Kinase = n$label,
      FDR = n$fdr,
      P_value = n$pval,
      Fold_Enrichment = n$fold_enrichment,
      Substrates = n$substrate_count,
      stringsAsFactors = FALSE
    )
  }))

  # Sort by FDR
  kinase_table <- kinase_table[order(kinase_table$FDR), ]

  # Write to CSV
  csv_file <- paste0(netNm, "_kinases.csv")
  write.csv(kinase_table, csv_file, row.names = FALSE)
  msg("[.generateKinaseNetworkTable] Wrote kinase table: ", csv_file)

  return(kinase_table)
}

#' Regenerate kinase-substrate network with a different database
#'
#' Runs complete kinase enrichment pipeline with a new database selection
#'
#' @param dataName Dataset name
#' @param database Database name (e.g., "phosphositeplus", "networkin", "all")
#' @param netNm Network name for output file
#' @param fdr_cutoff FDR cutoff (default: 0.05)
#' @return 1 on success, 0 on failure
SwitchKinaseDatabase <- function(dataName = "", database = "phosphositeplus", netNm = "kinase_substrate_net", fdr_cutoff = 0.05) {

  paramSet <- readSet(paramSet, "paramSet")
  msgSet <- readSet(msgSet, "msgSet")

  msg("[SwitchKinaseDatabase] Regenerating kinase network with database: ", database)

  # Step 1: Perform kinase enrichment analysis
  msg("[SwitchKinaseDatabase] Step 1/3: Running kinase enrichment...")
  enrich_result <- PerformKinaseEnrichment(dataName, database)
  if (enrich_result <= 0) {
    msgSet$current.msg <- "Kinase enrichment failed. No significant kinases found."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Step 2: Initialize kinase-substrate network
  msg("[SwitchKinaseDatabase] Step 2/3: Initializing network...")
  init_result <- InitKinaseSubstrateNetwork(dataName, database, fdr_cutoff)
  if (init_result == 0) {
    msgSet$current.msg <- "Kinase network initialization failed."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Step 3: Generate network visualization
  msg("[SwitchKinaseDatabase] Step 3/3: Generating visualization...")
  view_result <- PerformKinaseSubstrateNetworkView(dataName, netNm, fdr_cutoff)
  if (view_result == 0) {
    msgSet$current.msg <- "Kinase network visualization failed."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  msgSet$current.msg <- paste0("Kinase network regenerated with database '", database,
                                "'. Found ", enrich_result, " significant kinases.")
  saveSet(msgSet, "msgSet")

  msg("[SwitchKinaseDatabase] Complete! Network file: ", netNm, ".json")
  return(1)
}
