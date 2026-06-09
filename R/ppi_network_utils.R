##################################################
## R helpers to build ProteoAnalyst PPI networks
## Based on the NetworkAnalyst workflow (SearchNetDB/CreateGraph)
## The implementation queries the local SQLite PPI database
## and materializes the result in analSet$ppi.comps so Java clients
## can interrogate it via existing utility bridges.
##################################################

.ppi_query_cache <- list(edges = NULL, seeds = NULL, table = NULL, order = 0, minScore = 0)

#' ComputeSubnetStats
#'
#' Helper function to compute basic statistics for network components
#' This is needed by DecomposeGraph which is called during PPI network building
#'
#' @param comps List of igraph components
#' @return Data frame with Node, Edge, Query counts for each component
#'
ComputeSubnetStats <- function(comps){
  library(igraph);
  net.stats <- as.data.frame(matrix(0, ncol = 3, nrow = length(comps)));
  colnames(net.stats) <- c("Node", "Edge", "Query");
  for(i in 1:length(comps)){
    g <- comps[[i]];
    num_nodes <- vcount(g);
    num_edges <- ecount(g);
    # Count query nodes (from original upload)
    num_query <- if (!is.null(V(g)$is_query)) {
      sum(V(g)$is_query, na.rm = TRUE)
    } else {
      cat(sprintf("[PPI] Module %d: is_query attribute is NULL!\n", i))
      0  # If is_query attribute doesn't exist, default to 0
    }
    net.stats[i,] <- c(num_nodes, num_edges, num_query);
  }
  return(net.stats);
}

.collectSeedProteins <- function() {
  paramSet <- readSet(paramSet, "paramSet")
  seeds <- character(0)

  if (!is.null(paramSet$mdata.all)) {
    for (nm in names(paramSet$mdata.all)) {
      ds <- readDataset(nm, quiet = TRUE)
      if (!is.null(ds)) {
        if (!is.null(ds$seeds.proteins)) {
          seeds <- c(seeds, ds$seeds.proteins)
        } else if (!is.null(ds$sig.mat)) {
          seeds <- c(seeds, rownames(ds$sig.mat))
        }
      }
    }
  }

  if (length(seeds) == 0 && exists("dataSets")) {
    for (nm in names(dataSets)) {
      ds <- dataSets[[nm]]
      if (!is.null(ds$seeds.proteins)) {
        seeds <- c(seeds, ds$seeds.proteins)
      } else if (!is.null(ds$sig.mat)) {
        seeds <- c(seeds, rownames(ds$sig.mat))
      }
    }
  }

  if (length(seeds) == 0 ){
      ds <- readDataset(paramSet$dataName)
      if (!is.null(ds$seeds.proteins)) {
        seeds <- c(seeds, ds$seeds.proteins)
      } else if (!is.null(ds$sig.mat)) {
        seeds <- c(seeds, rownames(ds$sig.mat))
      }
  }

  seeds <- unique(as.character(na.omit(seeds)))
  #msg(sprintf("[PPI] collected %d unique seed(s)\n", length(seeds)))
  if (length(seeds) > 0) {
    # msg(sprintf("[PPI] Sample collected seeds: %s\n", paste(head(seeds, 5), collapse=", ")))
  }
  # flush.console()
  seeds
}

.ensurePpiList <- function() {
  analSet <- readSet(analSet, "analSet")
  if (is.null(analSet$ppi.comps)) {
    analSet$ppi.comps <- list()
  }
  analSet
}

.roundScore <- function(confThresh) {
  if (length(confThresh) == 0) {
    return(0)
  }
  score <- as.numeric(confThresh)[1]
  if (is.na(score)) {
    return(0)
  }
  if (score <= 1) {
    score <- score * 1000
  }
  round(score, 4)
}

.normalizeUniprotIds <- function(ids) {
  if (length(ids) == 0) {
    return(ids)
  }
  cleaned <- trimws(as.character(ids))
  cleaned <- sub("^.*\\|([^|]+)\\|.*$", "\\1", cleaned)
  cleaned <- sub("^([^;]+);.*$", "\\1", cleaned)
  cleaned <- sub("^([^,]+),.*$", "\\1", cleaned)
  cleaned <- sub("-\\d+$", "", cleaned)
  # Remove phosphosite annotations (e.g., Q9D1F4_T_247 → Q9D1F4)
  # Pattern: _[A-Z]_\d+ (underscore, amino acid, underscore, position)
  cleaned <- sub("_[A-Z]_\\d+$", "", cleaned)
  cleaned
}

.isEntrezLike <- function(ids) {
  if (length(ids) == 0) {
    return(FALSE)
  }
  mean(grepl("^[0-9]+$", ids)) > 0.9
}

.collectUniprotMappingFromDatasets <- function() {
  paramSet <- readSet(paramSet, "paramSet")
  mappings <- character(0)

  append_mapping <- function(ds) {
    if (is.null(ds) || is.null(ds$GeneAnotDB)) {
      return(character(0))
    }
    gene.db <- ds$GeneAnotDB
    if (!is.data.frame(gene.db) || !"gene_id" %in% colnames(gene.db)) {
      return(character(0))
    }
    key.col <- NULL
    if ("orig" %in% colnames(gene.db)) {
      key.col <- "orig"
    } else if ("accession" %in% colnames(gene.db)) {
      key.col <- "accession"
    }
    if (is.null(key.col)) {
      return(character(0))
    }
    keys <- as.character(gene.db[[key.col]])
    vals <- as.character(gene.db$gene_id)
    valid <- !(is.na(keys) | is.na(vals) | !nzchar(keys) | !nzchar(vals))
    if (!any(valid)) {
      return(character(0))
    }
    setNames(vals[valid], keys[valid])
  }

  if (!is.null(paramSet$mdata.all)) {
    for (nm in names(paramSet$mdata.all)) {
      ds <- readDataset(nm, quiet = TRUE)
      mappings <- c(mappings, append_mapping(ds))
    }
  }

  if (length(mappings) == 0 && exists("dataSets")) {
    for (nm in names(dataSets)) {
      mappings <- c(mappings, append_mapping(dataSets[[nm]]))
    }
  }

  if (length(mappings) == 0 && !is.null(paramSet$dataName)) {
    ds <- readDataset(paramSet$dataName)
    mappings <- c(mappings, append_mapping(ds))
  }

  if (length(mappings) == 0) {
    return(NULL)
  }
  mappings <- mappings[!duplicated(names(mappings))]
  mappings
}

.setCurrentMessage <- function(text) {
  msgSet <- readSet(msgSet, "msgSet")
  msgSet$current.msg <- text
  saveSet(msgSet, "msgSet")
  current.msg <<- text
}

.normalizePpiEdges <- function(edges) {
  if (is.null(edges) || !is.data.frame(edges) || nrow(edges) == 0) {
    return(edges)
  }

  if (!all(c("id1", "id2") %in% colnames(edges))) {
    AddErrMsg(paste(
      "PPI query did not return the required edge columns id1/id2.",
      "Available columns:",
      paste(colnames(edges), collapse = ", ")
    ));
    return(0);
  }

  if (!"combined_score" %in% colnames(edges)) {
    score.candidates <- c("score", "confidence", "conf", "weight", "experimental")
    score.col <- score.candidates[score.candidates %in% colnames(edges)][1]

    if (!is.na(score.col) && nzchar(score.col)) {
      edges$combined_score <- suppressWarnings(as.numeric(edges[[score.col]]))
      edges$combined_score[is.na(edges$combined_score)] <- 1000
    } else {
      edges$combined_score <- 1000
    }
  }

  edges[, c("id1", "id2", "combined_score"), drop = FALSE]
}

SearchNetDB <- function(dummy = NA, dbType = "ppi", dbName = "NA", requireExp = FALSE, confThresh = 0.5, order = 0) {
  paramSet <- readSet(paramSet, "paramSet")
  analSet <- readSet(analSet, "analSet")
  seeds <- .collectSeedProteins()
  if (length(seeds) == 0) {
    .setCurrentMessage("No seed features or proteins are available. Please upload data or select a list before running the PPI builder.")
    return(c(0, 0, 0))
  }

  if (dbType != "ppi") {
    .setCurrentMessage(paste("Unsupported database type:", dbType))
    return(c(0, 0, 0))
  }

  sqlite.path <- paramSet$sqlite.path
  if (!nzchar(sqlite.path)) {
    .setCurrentMessage("Cannot find the PPI SQLite database path on this server.")
    return(c(0, 0, 0))
  }

  min.score <- .roundScore(confThresh)

  # Determine if this database uses UniProt IDs (ppi_uniprot.sqlite) or Entrez IDs (ppi.sqlite).
  # Routing is per-database, not based on file existence.
  UNIPROT_DBS <- c("intact", "huri", "rolland", "irefinx")
  use_uniprot_sqlite <- dbName %in% UNIPROT_DBS &&
                        file.exists(paste0(sqlite.path, "ppi_uniprot.sqlite"))
  msg(sprintf("[PPI] Database: %s | ID mode: %s | ppi_uniprot.sqlite used: %s",
              dbName, ifelse(use_uniprot_sqlite, "UniProt", "Entrez"), use_uniprot_sqlite))
  if (length(seeds) > 0) {
  }

  uniprot_seeds <- seeds
  entrez_seeds <- seeds
  data.idType <- paramSet$data.idType

  if (.isEntrezLike(seeds)) {
    entrez_seeds <- as.character(seeds)
    dataset_uniprot_map <- .collectUniprotMappingFromDatasets()
    if (!is.null(dataset_uniprot_map) && length(dataset_uniprot_map) > 0) {
      analSet$uniprot_to_entrez_map <- dataset_uniprot_map
      saveSet(analSet, "analSet")
    }
  } else if (use_uniprot_sqlite) {

    normalized_seeds <- .normalizeUniprotIds(seeds)
    entrez_seeds <- as.character(normalized_seeds)

    dataset_uniprot_map <- .collectUniprotMappingFromDatasets()
    if (!is.null(dataset_uniprot_map) && length(dataset_uniprot_map) > 0) {
      analSet$uniprot_to_entrez_map <- dataset_uniprot_map
      saveSet(analSet, "analSet")
    }
  } else {
    dataset_uniprot_map <- .collectUniprotMappingFromDatasets()
    org <- if (!is.null(paramSet$data.org)) paramSet$data.org else "hsa"
    uniprot.map <- queryGeneDB("entrez_uniprot", org)

    if (!is.null(uniprot.map) && is.data.frame(uniprot.map)) {
      uniprot_to_entrez <- setNames(as.character(uniprot.map$gene_id), as.character(uniprot.map$accession))
      normalized_seeds <- .normalizeUniprotIds(uniprot_seeds)
      lookup_seeds <- as.character(normalized_seeds)
      entrez_seeds <- uniprot_to_entrez[lookup_seeds]
      matched_count <- sum(!is.na(entrez_seeds))

      if (matched_count == 0) {
        if (!is.null(dataset_uniprot_map) && length(dataset_uniprot_map) > 0) {
          alt_entrez <- dataset_uniprot_map[as.character(uniprot_seeds)]
          alt_matched <- sum(!is.na(alt_entrez))
          if (alt_matched > 0) {
            entrez_seeds <- alt_entrez
            matched_count <- alt_matched
          }
        }
      }

      if (matched_count == 0) {
        names(uniprot_to_entrez) <- toupper(names(uniprot_to_entrez))
        lookup_seeds <- toupper(lookup_seeds)
        entrez_seeds <- uniprot_to_entrez[lookup_seeds]
        matched_count <- sum(!is.na(entrez_seeds))
      }

      entrez_seeds <- entrez_seeds[!is.na(entrez_seeds)]

      if (length(entrez_seeds) > 0) {
        uniprot_to_entrez_map <- uniprot_to_entrez[lookup_seeds]
        names(uniprot_to_entrez_map) <- as.character(uniprot_seeds)
        uniprot_to_entrez_map <- uniprot_to_entrez_map[!is.na(uniprot_to_entrez_map)]
        analSet$uniprot_to_entrez_map <- uniprot_to_entrez_map
        saveSet(analSet, "analSet")
      } else if (.isEntrezLike(uniprot_seeds)) {
        entrez_seeds <- as.character(uniprot_seeds)
        if (!is.null(dataset_uniprot_map) && length(dataset_uniprot_map) > 0) {
          analSet$uniprot_to_entrez_map <- dataset_uniprot_map
          saveSet(analSet, "analSet")
        }
      } else {
        .setCurrentMessage("Failed to convert UniProt IDs to Entrez for PPI database query.")
        return(c(0, 0, 0))
      }
    } else {
    }
  }

  seeds <- as.character(entrez_seeds)
  original_seeds <- seeds  # Store original seeds before neighbor expansion

  if (length(seeds) == 0) {
    .setCurrentMessage("No valid IDs available for PPI database query after conversion.")
    return(c(0, 0, 0))
  }

  org <- paramSet$data.org
  if (is.null(org) || org == "") {
    org <- "hsa"
  }

  loc.path <- paste0(paramSet$lib.path, org, "/localization.qs")
  loc.data <- if (file.exists(loc.path)) ov_qs_read(loc.path) else NULL

  msg(sprintf("[PPI] Querying %s in %s mode with %d seeds (min.score=%s, order=%s)",
              dbName, ifelse(use_uniprot_sqlite, "UniProt", "Entrez"), length(seeds), min.score, order))

  edges <- tryCatch({
    QueryPpiSQLite(sqlite.path, dbName, seeds, requireExp, min.score, use.uniprot = use_uniprot_sqlite)
  }, error = function(e) {
    msg(sprintf("[PPI] ERROR during database query: %s", conditionMessage(e)))
    .setCurrentMessage(paste("Database query failed:", conditionMessage(e)))
    NULL
  })

  if (is.null(edges)) {
    if (!exists("current.msg", inherits = TRUE) || is.null(current.msg) || current.msg == "") {
      .setCurrentMessage("Database query returned NULL.")
    }
    return(c(0, 0, 0))
  } else if (nrow(edges) == 0) {
    msg(sprintf("[PPI] No matches found in %s", dbName))
    if (!exists("current.msg", inherits = TRUE) || is.null(current.msg) || current.msg == "") {
      .setCurrentMessage("No matches found for the given seed features/proteins.")
    }
    return(c(0, 0, 0))
  } else {
    all.db.nodes <- unique(c(edges$id1, edges$id2))
    matched.seeds <- sum(seeds %in% all.db.nodes)
    msg(sprintf("[PPI] Query successful: %d edges found; %d/%d seeds matched (%.1f%%)",
                nrow(edges), matched.seeds, length(seeds), 100 * matched.seeds / max(1, length(seeds))))
  }

  edges <- unique(.normalizePpiEdges(edges))
  invalid.edge.count <- sum(is.na(edges$id1) | is.na(edges$id2) | !nzchar(as.character(edges$id1)) | !nzchar(as.character(edges$id2)))
  if (invalid.edge.count > 0) {
    msg(sprintf("[PPI] Removing %d invalid edges with missing node IDs before network filtering", invalid.edge.count))
    edges <- edges[!(is.na(edges$id1) | is.na(edges$id2) | !nzchar(as.character(edges$id1)) | !nzchar(as.character(edges$id2))), , drop = FALSE]
  }
  edges <- edges[edges$id1 != edges$id2, , drop = FALSE]
  if (nrow(edges) == 0) {
    .setCurrentMessage("Only self-loop edges were returned; nothing to build.")
    return(c(0, 0, 0))
  }

  # Zero-order network should only include interactions between the uploaded proteins
  if (order == 0) {
    edges <- edges[edges$id1 %in% seeds & edges$id2 %in% seeds, , drop = FALSE]
    if (nrow(edges) == 0) {
      .setCurrentMessage("No direct interactions found among the uploaded proteins at zero order.")
      return(c(0, 0, 0))
    }
  }

  # Apply compartment-aware filtering if localization data is available and order >= 1
  if (!is.null(loc.data) && order >= 1) {
    # msg("[PPI] Applying compartment-aware filtering...\n")

    # Create compartment map (EntrezID -> Broad.category)
    comp.map <- setNames(as.character(loc.data$Broad.category), as.character(loc.data$EntrezID))
    comp.map[is.na(comp.map)] <- "Unknown"

    # Map seeds to compartments
    seed.comps <- comp.map[seeds]
    seed.comps[is.na(seed.comps)] <- "Unknown"

    # Identify all nodes from initial query
    all.nodes <- unique(c(edges$id1, edges$id2))
    is.seed <- all.nodes %in% seeds

    # PASS 1: Intra-compartment expansion (add neighbors in same compartment)
    # OPTIMIZED: Use list to collect edges, then combine once (10-100x faster)
    intra.edges.list <- vector("list", nrow(edges))
    intra.count <- 0
    network.nodes <- seeds  # Start with seeds

    for (i in 1:nrow(edges)) {
      id1 <- as.character(edges$id1[i])
      id2 <- as.character(edges$id2[i])

      # Get compartments
      comp1 <- comp.map[id1]
      comp2 <- comp.map[id2]
      if (is.na(comp1)) comp1 <- "Unknown"
      if (is.na(comp2)) comp2 <- "Unknown"

      # Check if at least one node is a seed
      id1.is.seed <- id1 %in% seeds
      id2.is.seed <- id2 %in% seeds

      if (!id1.is.seed && !id2.is.seed) {
        # Neither is a seed, skip for now
        next
      }

      # Check compartment compatibility
      if (comp1 == comp2) {
        # Same compartment: always add (first-order expansion)
        intra.count <- intra.count + 1
        intra.edges.list[[intra.count]] <- edges[i, , drop = FALSE]
        network.nodes <- unique(c(network.nodes, id1, id2))
      } else {
        # Different compartments: only add if both are already in network (seeds)
        if (id1.is.seed && id2.is.seed) {
          intra.count <- intra.count + 1
          intra.edges.list[[intra.count]] <- edges[i, , drop = FALSE]
          network.nodes <- unique(c(network.nodes, id1, id2))
        }
      }
    }

    # Combine all edges at once (much faster than growing with rbind)
    if (intra.count > 0) {
      intra.edges <- do.call(rbind, intra.edges.list[1:intra.count])
    } else {
      intra.edges <- data.frame()
    }

    edges <- if (intra.count > 0) intra.edges else data.frame()

    # For dense databases without a confidence score filter (non-STRING), prune non-seed
    # neighbor nodes that connect to only one seed — these are peripheral and inflate the
    # network. Only retain non-seeds with >= 2 seed connections (linker nodes).
    if (!grepl("string$", dbName) && nrow(edges) > 0) {
      all.edge.nodes <- unique(c(as.character(edges$id1), as.character(edges$id2)))
      non.seeds <- setdiff(all.edge.nodes, seeds)
      if (length(non.seeds) > 0) {
        seed.conn <- vapply(non.seeds, function(nd) {
          partners <- c(as.character(edges$id2[as.character(edges$id1) == nd]),
                        as.character(edges$id1[as.character(edges$id2) == nd]))
          sum(partners %in% seeds)
        }, integer(1))
        keep.nodes <- c(seeds, non.seeds[seed.conn >= 2])
        edges <- edges[as.character(edges$id1) %in% keep.nodes &
                       as.character(edges$id2) %in% keep.nodes, , drop = FALSE]
        network.nodes <- unique(c(seeds, as.character(edges$id1), as.character(edges$id2)))
        msg(sprintf("[PPI] After linker filter: %d nodes, %d edges", length(network.nodes), nrow(edges)))
      }
    }
  }

  nodes <- unique(c(edges$id1, edges$id2))
  .setCurrentMessage(paste("Found", length(nodes), "proteins and", nrow(edges), "interactions."))
  .ppi_query_cache <<- list(edges = edges, seeds = original_seeds, table = dbName, order = order, minScore = min.score)
  msgSet <- readSet(msgSet, "msgSet")
  msgSet$current.msg <- current.msg
  saveSet(msgSet, "msgSet")
  return(c(length(nodes), nrow(edges), length(original_seeds)))
}

CreateGraph <- function(dummy = NA) {
  cache <- .ppi_query_cache
  if (is.null(cache$edges) || nrow(cache$edges) == 0) {
  .setCurrentMessage("No PPI result is cached. Please run SearchNetDB first.")
  return(c(0, 0, 0))
}

  # msg(sprintf("[PPI] creating graph from %d edges\n", nrow(cache$edges)))
  # flush.console()
  invalid.edge.count <- sum(is.na(cache$edges$id1) | is.na(cache$edges$id2) | !nzchar(as.character(cache$edges$id1)) | !nzchar(as.character(cache$edges$id2)))
  if (invalid.edge.count > 0) {
    msg(sprintf("[PPI] Removing %d invalid cached edges before graph creation", invalid.edge.count))
    cache$edges <- cache$edges[!(is.na(cache$edges$id1) | is.na(cache$edges$id2) | !nzchar(as.character(cache$edges$id1)) | !nzchar(as.character(cache$edges$id2))), , drop = FALSE]
  }
  g <- igraph::graph_from_data_frame(cache$edges[, c("id1", "id2")], directed = FALSE)
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

  # Seed proteins whose edges were fully pruned (e.g. by the linker filter) won't
  # appear in the edge list and are therefore missing from the graph. Add them back
  # as isolated vertices so their expression values are preserved in the output.
  if (!is.null(cache$seeds) && length(cache$seeds) > 0) {
    missing.seeds <- setdiff(as.character(cache$seeds), V(g)$name)
    if (length(missing.seeds) > 0) {
      g <- igraph::add_vertices(g, length(missing.seeds), name = missing.seeds)
      msg(sprintf("[CreateGraph] Added %d isolated seed(s) missing from edge list", length(missing.seeds)))
    }
  }

  if (igraph::vcount(g) == 0) {
    current.msg <<- "Failed to construct graph: no nodes remain after simplification."
    return(c(0, 0, 0))
  }

  paramSet <- readSet(paramSet, "paramSet")
  # msg(sprintf("[PPI] graph nodes=%d edges=%d\n", igraph::vcount(g), igraph::ecount(g)))
  # flush.console()

  # Attach expression values (if provided) to nodes for downstream exports.
  # Fill from multiple sources in order and keep unresolved nodes as NA until the end.
  expr.vec <- rep(NA_real_, igraph::vcount(g))
  names(expr.vec) <- V(g)$name
  analSet <- readSet(analSet, "analSet")

  fillExprFromNamedValues <- function(named.vals) {
    if (is.null(named.vals) || length(named.vals) == 0) {
      return(invisible(NULL))
    }
    src.names <- names(named.vals)
    if (is.null(src.names)) {
      return(invisible(NULL))
    }
    src.names <- as.character(src.names)
    src.vals <- suppressWarnings(as.numeric(named.vals))
    keep <- !(is.na(src.names) | !nzchar(src.names) | is.na(src.vals))
    if (!any(keep)) {
      return(invisible(NULL))
    }
    src.names <- src.names[keep]
    src.vals <- src.vals[keep]

    # Direct ID matching.
    direct <- src.vals[match(names(expr.vec), src.names)]
    miss <- is.na(expr.vec) & !is.na(direct)
    expr.vec[miss] <<- direct[miss]

    # Normalized matching for UniProt-like IDs.
    if (any(is.na(expr.vec))) {
      tgt.norm <- .normalizeUniprotIds(names(expr.vec))
      src.norm <- .normalizeUniprotIds(src.names)
      norm.map <- src.vals[match(tgt.norm, src.norm)]
      miss.norm <- is.na(expr.vec) & !is.na(norm.map)
      expr.vec[miss.norm] <<- norm.map[miss.norm]
    }
  }

  # Build FC mapping directly from DE results (preferred source for Expr.).
  collectCompResLogFC <- function() {
    out <- numeric(0)
    append_from_dataset <- function(ds) {
      if (is.null(ds) || is.null(ds$comp.res) || !is.data.frame(ds$comp.res) || nrow(ds$comp.res) == 0) {
        return()
      }
      if (!("logFC" %in% colnames(ds$comp.res))) {
        return()
      }
      vals <- suppressWarnings(as.numeric(ds$comp.res$logFC))
      ids <- rownames(ds$comp.res)
      keep <- !(is.na(ids) | !nzchar(ids) | is.na(vals))
      if (!any(keep)) {
        return()
      }
      ids <- as.character(ids[keep])
      vals <- vals[keep]

      # 1) direct IDs from comp.res rownames
      direct <- stats::setNames(vals, ids)
      out <<- c(out, direct)

      # 2) mapped IDs (UniProt -> Entrez), for PPI graph node names
      up2ent <- analSet$uniprot_to_entrez_map
      if (!is.null(up2ent) && length(up2ent) > 0) {
        up.keys <- as.character(names(up2ent))
        up.vals <- as.character(up2ent)
        valid.map <- !(is.na(up.keys) | !nzchar(up.keys) | is.na(up.vals) | !nzchar(up.vals))
        if (any(valid.map)) {
          up.keys <- up.keys[valid.map]
          up.vals <- up.vals[valid.map]

          mapped.direct <- up.vals[match(ids, up.keys)]
          hit.direct <- !is.na(mapped.direct) & nzchar(mapped.direct)
          if (any(hit.direct)) {
            out <<- c(out, stats::setNames(vals[hit.direct], mapped.direct[hit.direct]))
          }

          # normalized UniProt matching: sp|P12345|NAME, isoforms, phosphosite suffixes, etc.
          id.norm <- .normalizeUniprotIds(ids)
          up.norm <- .normalizeUniprotIds(up.keys)
          mapped.norm <- up.vals[match(id.norm, up.norm)]
          hit.norm <- !is.na(mapped.norm) & nzchar(mapped.norm)
          if (any(hit.norm)) {
            out <<- c(out, stats::setNames(vals[hit.norm], mapped.norm[hit.norm]))
          }
        }
      }
    }

    if (!is.null(paramSet$mdata.all)) {
      for (nm in names(paramSet$mdata.all)) {
        ds <- readDataset(nm, quiet = TRUE)
        append_from_dataset(ds)
      }
    }
    if (length(out) == 0 && exists("dataSets")) {
      for (nm in names(dataSets)) {
        append_from_dataset(dataSets[[nm]])
      }
    }
    if (length(out) == 0 && !is.null(paramSet$dataName)) {
      ds <- readDataset(paramSet$dataName, quiet = TRUE)
      append_from_dataset(ds)
    }

    if (length(out) == 0) {
      return(out)
    }
    # Preserve first occurrence for deterministic behavior
    out <- out[!duplicated(names(out))]
    out
  }

  msg(sprintf("[CreateGraph-EXPR] graph nodes=%d; sample node names: %s",
              igraph::vcount(g), paste(head(V(g)$name, 5), collapse=",")))
  msg(sprintf("[CreateGraph-EXPR] uniprot_to_entrez_map size=%d; sample: %s",
              length(analSet$uniprot_to_entrez_map),
              if (length(analSet$uniprot_to_entrez_map) > 0)
                paste0(head(names(analSet$uniprot_to_entrez_map), 3), "->",
                       head(analSet$uniprot_to_entrez_map, 3), collapse=",") else "EMPTY"))

  de.fc.vals <- collectCompResLogFC()
  msg(sprintf("[CreateGraph-EXPR] collectCompResLogFC returned %d values; sample: %s",
              length(de.fc.vals),
              if (length(de.fc.vals) > 0) paste(head(names(de.fc.vals), 3), collapse=",") else "NONE"))
  if (length(de.fc.vals) > 0) {
    fillExprFromNamedValues(de.fc.vals)
  }

  msg(sprintf("[CreateGraph-EXPR] after DE fill: %d/%d nodes have non-zero expr",
              sum(!is.na(expr.vec) & expr.vec != 0), length(expr.vec)))

  if (!is.null(paramSet$all.ent.mat)) {
    exp.mat <- paramSet$all.ent.mat
    if (!is.null(dim(exp.mat)) && nrow(exp.mat) > 0) {
      cat(sprintf("[CreateGraph-EXPR] all.ent.mat: %d rows, sample rownames: %s, sample values: %s\n",
                  nrow(exp.mat), paste(head(rownames(exp.mat), 3), collapse=","),
                  paste(head(suppressWarnings(as.numeric(exp.mat[,1])), 3), collapse=",")))
      exp.vals <- suppressWarnings(as.numeric(exp.mat[, 1]))
      names(exp.vals) <- rownames(exp.mat)
      fillExprFromNamedValues(exp.vals)
      cat(sprintf("[CreateGraph-EXPR] after all.ent.mat direct fill: %d/%d nodes non-zero\n",
                  sum(!is.na(expr.vec) & expr.vec != 0), length(expr.vec)))
      # If rownames are UniProt-like, re-key via uniprot_to_entrez_map and try again
      up2ent <- analSet$uniprot_to_entrez_map
      if (!is.null(up2ent) && length(up2ent) > 0) {
        norm.up <- .normalizeUniprotIds(rownames(exp.mat))
        norm.keys <- .normalizeUniprotIds(names(up2ent))
        entrez.from.up <- up2ent[match(norm.up, norm.keys)]
        valid.conv <- !is.na(entrez.from.up) & nzchar(entrez.from.up)
        if (any(valid.conv)) {
          conv.vals <- suppressWarnings(as.numeric(exp.mat[valid.conv, 1]))
          names(conv.vals) <- as.character(entrez.from.up[valid.conv])
          fillExprFromNamedValues(conv.vals)
          cat(sprintf("[CreateGraph-EXPR] after UniProt->Entrez conversion fill: %d/%d nodes non-zero (converted %d IDs)\n",
                      sum(!is.na(expr.vec) & expr.vec != 0), length(expr.vec), sum(valid.conv)))
        } else {
          cat("[CreateGraph-EXPR] UniProt->Entrez map present but no overlap with all.ent.mat rownames\n")
        }
      } else {
        cat("[CreateGraph-EXPR] no uniprot_to_entrez_map available\n")
      }
    }
  } else {
    cat("[CreateGraph-EXPR] paramSet$all.ent.mat is NULL - expression values not available\n")
  }

  if (!is.null(paramSet$all.prot.mat)) {
    exp.tbl <- paramSet$all.prot.mat
    if (!is.null(dim(exp.tbl)) && nrow(exp.tbl) > 0 && ncol(exp.tbl) >= 2) {
      msg(sprintf("[CreateGraph-EXPR] all.prot.mat: %d rows, sample col2: %s",
                  nrow(exp.tbl), paste(head(as.character(exp.tbl[, 2]), 3), collapse=",")))
      exp.vals <- suppressWarnings(as.numeric(exp.tbl[, 1]))
      names(exp.vals) <- as.character(exp.tbl[, 2])
      fillExprFromNamedValues(exp.vals)
    }
  }

  cat(sprintf("[CreateGraph-EXPR] FINAL: %d/%d nodes have non-zero expr\n",
              sum(!is.na(expr.vec) & expr.vec != 0), length(expr.vec)))

  expr.vec[is.na(expr.vec)] <- 0
  V(g)$expr <- expr.vec
  V(g)$abundance <- expr.vec  # Also set abundance for compatibility
  expr.vec <<- expr.vec[!is.na(expr.vec)]

  # Mark seed nodes and attach both Entrez and UniProt identifiers. PPI graphs
  # can be stored with UniProt vertex names, while localization and symbols are
  # keyed by Entrez.
  node.names <- V(g)$name
  org <- if (!is.null(paramSet$data.org) && nzchar(paramSet$data.org)) {
    paramSet$data.org
  } else if (!is.null(paramSet$org) && nzchar(paramSet$org)) {
    paramSet$org
  } else {
    "hsa"
  }
  graph.ids <- .paResolveGraphIds(g, org)
  V(g)$entrez <- graph.ids$entrez
  V(g)$uniprot <- graph.ids$uniprot

  seed.ids <- unique(trimws(as.character(cache$seeds)))
  seed.ids <- seed.ids[!is.na(seed.ids) & nzchar(seed.ids)]
  seed.uniprot <- .paNormalizeUniprotIds(seed.ids)
  is.seed <- node.names %in% seed.ids |
    V(g)$entrez %in% seed.ids |
    V(g)$uniprot %in% seed.uniprot
  V(g)$is_query <- is.seed;
  V(g)$size <- .paScaleGraphNodeSizes(g)

  gene.symbols <- tryCatch({
    doEntrez2SymbolMapping(V(g)$entrez, org, "entrez")
  }, error = function(e) {
    rep(NA_character_, length(node.names))
  })
  gene.symbols[is.na(gene.symbols) | !nzchar(gene.symbols)] <- node.names[is.na(gene.symbols) | !nzchar(gene.symbols)]
  V(g)$gene_symbol <- gene.symbols

  # Add compartment/localization information as vertex attributes
  # This preserves compartment data for extracted modules
  loc.path <- paste0(paramSet$lib.path, org, "/", org, "_localization.qs")

  if (file.exists(loc.path)) {
    msg(sprintf("[CreateGraph] Loading localization data from: %s\n", loc.path))
    loc.data <- ov_qs_read(loc.path)
    loc.map <- loc.data[match(as.character(V(g)$entrez), as.character(loc.data$EntrezID)), ]

    comp.res <- .paResolveCompartmentAnnotations(
      ifelse(is.na(loc.map$Broad.category) | loc.map$Broad.category == "",
             "Unknown",
             as.character(loc.map$Broad.category)),
      ifelse(is.na(loc.map$Main.location),
             "Unknown",
             as.character(loc.map$Main.location))
    )

    V(g)$broad_category <- comp.res$primary
    V(g)$compartment_all <- comp.res$all_categories
    V(g)$main_location <- comp.res$all_locations

    msg(sprintf("[CreateGraph] Added compartment info to %d/%d nodes\n",
                sum(V(g)$broad_category != "Unknown"), vcount(g)))
  } else {
    # No localization data available - set defaults
    V(g)$broad_category <- "Unknown"
    V(g)$compartment_all <- "Unknown"
    V(g)$main_location <- "Unknown"
    msg(sprintf("[CreateGraph] Localization data not found, all nodes marked as Unknown\n"))
  }

  # Store the overall graph globally (for compatibility with NetworkAnalyst)
  overall.graph <<- g
  # msg(sprintf("[PPI] Stored overall.graph globally\n"))

  # Decompose graph into connected subnetworks
  # msg(sprintf("[PPI] Decomposing graph into connected components...\n"))
  # flush.console()

  analSet <- .ensurePpiList()

  # Use minNodeNum = 2 for very small networks (allows 2-node networks to be visualized)
  minNodes <- ifelse(vcount(g) < 3, 2, 3)

  analSet <- DecomposeGraph(g, analSet, minNodeNum = minNodes)
  substats <- analSet$substats

  if (is.null(substats) || length(substats) == 0) {
    msg("[PPI] WARNING: No subnetworks found after decomposition\n")
    return(c(0, 0, 0))
  }

  # Store overall graph in analSet
  analSet$overall.graph <- overall.graph

  # Update global variables for compatibility
  net.stats <<- analSet$net.stats
  ppi.comps <<- analSet$ppi.comps

  # Set current network to the first (largest) module
  if (length(ppi.comps) > 0) {
    current.net.nm <<- names(ppi.comps)[1]
    paramSet$current.net.nm <- current.net.nm
    msg(sprintf("[PPI] Set current network to %s with %d nodes\n",
                current.net.nm, vcount(ppi.comps[[current.net.nm]])))
  }

  saveSet(paramSet, "paramSet")
  saveSet(analSet, "analSet")

  # Build summary message
  num_seeds <- length(cache$seeds)
  num_nodes <- igraph::vcount(g)
  num_edges <- igraph::ecount(g)
  num_modules <- length(ppi.comps)

  msg <- sprintf("Built network with %d nodes, %d edges from %d query genes. Decomposed into %d connected subnetwork(s).",
                 num_nodes, num_edges, num_seeds, num_modules)
  .setCurrentMessage(msg)

  # msg(sprintf("[PPI] %s\n", msg))
  # msg(sprintf("[PPI] Subnetwork sizes: %s\n", paste(substats, collapse=", ")))
  # msg(sprintf("[PPI] Network statistics table:\n"))
  print(net.stats)
  flush.console()

  # NOTE: DO NOT clear cache - it's needed for subsequent operations
  # (GetMinConnectedGraphs, BuildSeedProteinNet, ComputePCSFNet)
  # Cache is now preserved to enable network refinement operations
  # .ppi_query_cache <<- list(edges = NULL, seeds = NULL, table = NULL, order = 0, minScore = 0)

  msgSet <- readSet(msgSet, "msgSet")
  msgSet$current.msg <- current.msg
  saveSet(msgSet, "msgSet")

  # Return: num_seeds, num_nodes, num_edges, num_modules, substats
  output <- c(num_seeds, num_nodes, num_edges, num_modules, substats)
  return(output)
}

##################################################
## Network Visualization Export
##################################################

#' PrepareNetwork
#'
#' Exports a PPI network to JSON format for visualization
#' This is the main entry point called from Java for network viewer
#'
#' @param net.nm Name of the network (e.g., "module1", "ppi_1")
#' @param json.nm Output JSON file name
#' @return 1 if successful, 0 otherwise
#'
PrepareNetwork <- function(net.nm, json.nm) {

  # Get the network from analSet or from the first PPI network
  analSet <- readSet(analSet, "analSet")

  # Try to find the network by name
  if (!is.null(analSet$ppi.comps) && net.nm %in% names(analSet$ppi.comps)) {
    g <- analSet$ppi.comps[[net.nm]]
  } else if (!is.null(analSet$ppi.comps) && length(analSet$ppi.comps) > 0) {
    # Use the first network if named network not found
    g <- analSet$ppi.comps[[1]]
    net.nm <- names(analSet$ppi.comps)[1]
  } else if (exists("overall.graph") && !is.null(overall.graph)) {
    # Fallback to overall.graph
    g <- overall.graph
  } else {
    AddErrMsg("No PPI network found. Please build a network first.")
    return(0)
  }

  if (igraph::vcount(g) == 0) {
    AddErrMsg("Network is empty (no nodes).")
    return(0)
  }

  result <- tryCatch({
    msg(sprintf("[PrepareNetwork] Calling .exportNetworkToJSON...\n"))
    .exportNetworkToJSON(g, json.nm, net.nm)
    msg(sprintf("[PrepareNetwork] SUCCESS: Network exported to %s\n", json.nm))
    1
  }, error = function(e) {
    # msg(sprintf("[PrepareNetwork] ERROR in .exportNetworkToJSON:\n"))
    # msg(sprintf("[PrepareNetwork]   Error message: %s\n", e$message))
    # msg(sprintf("[PrepareNetwork]   Error class: %s\n", paste(class(e), collapse=", ")))
    if (!is.null(e$call)) {
      msg(sprintf("[PrepareNetwork]   Error call: %s\n", deparse(e$call)[1]))
    }
    traceback_lines <- capture.output(traceback())
    if (length(traceback_lines) > 0) {
      msg(sprintf("[PrepareNetwork]   Traceback:\n"))
      for (line in traceback_lines) {
        msg(sprintf("[PrepareNetwork]     %s\n", line))
      }
    }
    AddErrMsg(paste("Failed to export network:", e$message))
    0
  })

  msg(sprintf("[PrepareNetwork] Export result: %d (1=success, 0=failure)\n", result))

  if (result == 1) {
    current.net.nm <<- net.nm
  }

  return(result)
}

#' .exportNetworkToJSON
#'
#' Internal function to convert igraph to JSON for cytoscape/sigma.js visualization
#'
#' @param g igraph object
#' @param json.nm Output file name
#' @param net.nm Network name
#'
.exportNetworkToJSON <- function(g, json.nm, net.nm) {

  paramSet <- readSet(paramSet, "paramSet")

  # Validate graph
  if (is.null(g) || !igraph::is.igraph(g)) {
    AddErrMsg("Invalid graph object");
    return(0);
  }

  if (igraph::vcount(g) == 0) {
    AddErrMsg("Graph has no vertices");
    return(0);
  }

  # Get node names
  nms <- V(g)$name
  if (is.null(nms) || length(nms) == 0) {
    nms <- as.character(1:vcount(g))
    V(g)$name <- nms
  }
  n <- length(nms)

  # msg(sprintf("[exportNetworkToJSON] Exporting network '%s' with %d nodes, %d edges\n",
  #            net.nm, n, ecount(g)))

  # Get organism from paramSet
  org <- paramSet$data.org
  if(is.null(org) || !nzchar(org)) {
    warning(".exportNetworkToJSON: organism not found in paramSet, using 'hsa' as default")
    org <- "hsa"
  }
  graph.ids <- .paResolveGraphIds(g, org)
  V(g)$entrez <- graph.ids$entrez
  V(g)$uniprot <- graph.ids$uniprot

  # Calculate network metrics
  node.dgr <- igraph::degree(g)
  node.btw <- igraph::betweenness(g, directed = FALSE, normalized = FALSE)

  # Get expression values if available
  node.expr <- rep(0, n)
  if (!is.null(V(g)$expr)) {
    node.expr <- V(g)$expr
  } else if (exists("expr.vec") && !is.null(expr.vec)) {
    node.expr <- expr.vec[match(nms, names(expr.vec))]
    node.expr[is.na(node.expr)] <- 0
  }

  # Use existing graph sizes when present so switching between extracted and
  # original networks does not silently rescale the same nodes.
  stored.sizes <- if (!is.null(V(g)$size) && length(V(g)$size) == n) {
    suppressWarnings(as.numeric(V(g)$size))
  } else {
    rep(NA_real_, n)
  }
  if (all(is.finite(stored.sizes))) {
    node.sizes <- stored.sizes
  } else {
    node.sizes <- .paScaleGraphNodeSizes(g)
    V(g)$size <- node.sizes
  }

  # Compute colors based on expression (centered gradient)
  centered <- TRUE
  node.colsb <- ComputeColorGradient(node.expr, "black", centered, FALSE)
  node.colsw <- ComputeColorGradient(node.expr, "white", centered, FALSE)
  node.colsc <- ComputeColorGradient(node.expr, "colorblind", centered, TRUE)

  # Topology colors based on betweenness
  topo.val <- log(node.btw + 1)
  topo.colsb <- ComputeColorGradient(topo.val, "black", FALSE, FALSE)
  topo.colsw <- ComputeColorGradient(topo.val, "white", FALSE, FALSE)
  topo.colsc <- ComputeColorGradient(topo.val, "colorblind", FALSE, TRUE)

  # Get layout coordinates
  if (n > 500) {
    # Use faster layout for large networks
    layout.coords <- igraph::layout_with_fr(g, niter = 100)
  } else {
    layout.coords <- igraph::layout_with_fr(g, niter = 500)
  }

  # Normalize coordinates to [0, 100] (like localization network)
  pos.x <- (layout.coords[, 1] - min(layout.coords[, 1])) /
           (max(layout.coords[, 1]) - min(layout.coords[, 1])) * 100
  pos.y <- (layout.coords[, 2] - min(layout.coords[, 2])) /
           (max(layout.coords[, 2]) - min(layout.coords[, 2])) * 100

  # Load localization data for compartment assignment. Localization tables use
  # Entrez IDs, even when the graph itself uses UniProt vertex names.
  loc.path <- paste0(paramSet$lib.path, org, "/", org, "_localization.qs")

  node.labels <- tryCatch({
    doEntrez2SymbolMapping(graph.ids$entrez, org, "entrez")
  }, error = function(e) {
    rep(NA_character_, length(nms))
  })
  if (!is.null(V(g)$gene_symbol) && length(V(g)$gene_symbol) == length(nms)) {
    use.attr.label <- is.na(node.labels) | !nzchar(node.labels)
    node.labels[use.attr.label] <- as.character(V(g)$gene_symbol)[use.attr.label]
  }
  node.labels[is.na(node.labels) | !nzchar(node.labels)] <- nms[is.na(node.labels) | !nzchar(node.labels)]

  loc.map <- NULL
  node.categories <- rep("Unknown", length(nms))
  node.category.all <- rep("Unknown", length(nms))
  node.locations <- rep("Unknown", length(nms))
  if (file.exists(loc.path)) {
    msg(sprintf("[exportNetworkToJSON] Loading localization data from: %s
", loc.path))
    loc.data <- ov_qs_read(loc.path)
    loc.map <- loc.data[match(as.character(graph.ids$entrez), as.character(loc.data$EntrezID)), ]
    raw.categories <- ifelse(is.na(loc.map$Broad.category) | loc.map$Broad.category == "",
                             "Unknown",
                             as.character(loc.map$Broad.category))
    raw.locations <- ifelse(is.na(loc.map$Main.location) | loc.map$Main.location == "",
                            "Unknown",
                            as.character(loc.map$Main.location))
    comp.res <- .paResolveCompartmentAnnotations(raw.categories, raw.locations)
    node.categories <- comp.res$primary
    node.category.all <- comp.res$all_categories
    node.locations <- comp.res$all_locations
    msg(sprintf("[exportNetworkToJSON] Mapped %d/%d nodes to localization data
",
                sum(!is.na(loc.map$Broad.category)), length(nms)))
  } else {
    msg(sprintf("[exportNetworkToJSON] Localization data not found, using single compartment
"))
  }

  if (!is.null(V(g)$broad_category) && length(V(g)$broad_category) == length(nms)) {
    attr.locations <- if (!is.null(V(g)$main_location) && length(V(g)$main_location) == length(nms)) {
      as.character(V(g)$main_location)
    } else {
      rep(NA_character_, length(nms))
    }
    attr.res <- .paResolveCompartmentAnnotations(as.character(V(g)$broad_category), attr.locations)
    keep.attr <- !is.na(attr.res$primary) & nzchar(attr.res$primary) & attr.res$primary != "Unknown"
    node.categories[keep.attr] <- attr.res$primary[keep.attr]
    node.category.all[keep.attr] <- attr.res$all_categories[keep.attr]
    node.locations[keep.attr] <- attr.res$all_locations[keep.attr]
  }
  if (!is.null(V(g)$main_location) && length(V(g)$main_location) == length(nms)) {
    attr.locations <- as.character(V(g)$main_location)
    keep.attr <- !is.na(attr.locations) & nzchar(attr.locations) & attr.locations != "Unknown"
    node.locations[keep.attr] <- attr.locations[keep.attr]
  }

  V(g)$broad_category <- node.categories
  V(g)$compartment_all <- node.category.all
  V(g)$main_location <- node.locations

  uniprot.vec <- graph.ids$uniprot

  phosphosite.map <- list()
  if (!is.null(paramSet$phospho.mapping) && !is.null(paramSet$phospho.mapping$entrez.to.phospho)) {
    phosphosite.map <- paramSet$phospho.mapping$entrez.to.phospho
    msg(sprintf("[exportNetworkToJSON] Found phosphosite mapping for %d proteins
", length(phosphosite.map)))
  }

  category.colors <- .paCompartmentColors()

  gene.nodes <- lapply(seq_along(nms), function(i) {
    node.id <- as.character(nms[i])
    entrez.id <- if (!is.na(graph.ids$entrez[i]) && nzchar(graph.ids$entrez[i])) {
      as.character(graph.ids$entrez[i])
    } else {
      ""
    }
    category <- node.categories[i]
    category.all <- node.category.all[i]
    main.loc <- node.locations[i]

    comp.id <- gsub("[^A-Za-z0-9_]", "_", category)
    comp.color <- category.colors[[category]]
    if (is.null(comp.color)) comp.color <- "#999999"

    uniprot.id <- if (!is.na(uniprot.vec[i]) && !is.null(uniprot.vec[i]) && nzchar(uniprot.vec[i])) {
      as.character(uniprot.vec[i])
    } else {
      ""
    }

    is.query <- if (!is.null(V(g)$is_query)) V(g)$is_query[i] else FALSE

    phosphosites <- NULL
    if (length(phosphosite.map) > 0 && nzchar(entrez.id) && entrez.id %in% names(phosphosite.map)) {
      phosphosites <- phosphosite.map[[entrez.id]]
      if (length(phosphosites) > 0) {
        phosphosites <- paste(phosphosites, collapse = ";")
      } else {
        phosphosites <- NULL
      }
    }

    node.data <- list(
      id = node.id,
      label = node.labels[i],
      uniprot = uniprot.id,
      entrez = entrez.id,
      size = node.sizes[i],
      true_size = node.sizes[i],
      molType = "gene",
      type = "gene",
      colorb = comp.color,
      colorw = comp.color,
      topocolb = topo.colsb[i],
      topocolw = topo.colsw[i],
      topocolc = topo.colsc[i],
      expcolb = node.colsb[i],
      expcolw = node.colsw[i],
      expcolc = node.colsc[i],
      exp = round(node.expr[i], 3),
      posx = round(pos.x[i], 2),
      posy = round(pos.y[i], 2),
      compartment = comp.id,
      location = main.loc,
      main_location = main.loc,
      broad_category = category,
      compartment_all = category.all,
      all_compartments = category.all,
      degree = node.dgr[i],
      between = round(node.btw[i], 2),
      query = if (is.query) "Y" else "N"
    )

    if (!is.null(phosphosites)) {
      node.data$phosphosites <- phosphosites
    }

    node.data
  })

  compartment.map <- list()
  for (i in seq_along(nms)) {
    category <- node.categories[i]
    comp.id <- gsub("[^A-Za-z0-9_]", "_", category)

    if (is.null(compartment.map[[comp.id]])) {
      comp.color <- category.colors[[category]]
      if (is.null(comp.color)) comp.color <- "#999999"
      compartment.map[[comp.id]] <- list(
        id = comp.id,
        label = category,
        color = comp.color,
        node_ids = c()
      )
    }
    compartment.map[[comp.id]]$node_ids <- c(compartment.map[[comp.id]]$node_ids, as.character(nms[i]))
  }

  comp.nodes <- lapply(names(compartment.map), function(comp.id) {
    comp.info <- compartment.map[[comp.id]]
    node.indices <- which(sapply(gene.nodes, function(n) n$compartment == comp.id))
    if (length(node.indices) > 0) {
      comp.x <- mean(pos.x[node.indices])
      comp.y <- mean(pos.y[node.indices])
    } else {
      comp.x <- 50
      comp.y <- 50
    }

    list(
      id = comp.id,
      label = comp.info$label,
      size = 0,
      true_size = 0,
      molType = "compartment",
      type = "compartment",
      colorb = comp.info$color,
      colorw = comp.info$color,
      exp = 0,
      posx = round(comp.x, 2),
      posy = round(comp.y, 2),
      compartment = comp.id
    )
  })

  all.nodes <- c(gene.nodes, comp.nodes)

  edge.mat <- igraph::as_edgelist(g, names = TRUE)
  edges.list <- lapply(seq_len(nrow(edge.mat)), function(i) {
    list(
      id = paste0("e", i),
      source = as.character(edge.mat[i, 1]),
      target = as.character(edge.mat[i, 2]),
      weight = 1,
      size = 1
    )
  })

  # Build compartment metadata (like localization network)
  compartment.info <- lapply(names(compartment.map), function(comp.id) {
    comp.data <- compartment.map[[comp.id]]
    list(
      id = comp.id,
      label = comp.data$label,
      color = comp.data$color,
      nodes = as.list(comp.data$node_ids)
    )
  })

  # Find largest compartment
  largest.comp.id <- names(compartment.map)[which.max(sapply(compartment.map, function(x) length(x$node_ids)))]

  # Build nodeTable for Node Explorer
  # Extract node info from gene.nodes (excluding compartment meta-nodes)
  node.table <- lapply(gene.nodes, function(node) {
    list(
      id = node$id,
      label = node$label,
      uniprot = if (!is.null(node$uniprot) && node$uniprot != "") node$uniprot else "",
      entrez = if (!is.null(node$entrez) && node$entrez != "") node$entrez else "",
      degree = if (!is.null(node$degree)) node$degree else 0,
      betweenness = if (!is.null(node$between)) node$between else 0,
      expr = if (!is.null(node$exp)) node$exp else 0,
      location = if (!is.null(node$main_location)) node$main_location else "Unknown",
      primary_compartment = if (!is.null(node$broad_category)) node$broad_category else "Unknown",
      all_compartments = if (!is.null(node$all_compartments)) node$all_compartments else "Unknown"
    )
  })

  # Build JSON structure matching localization network format
  network.json <- list(
    nodes = all.nodes,
    edges = edges.list,
    compartments = compartment.info,
    nodeTable = node.table,
    idType = graph.ids$id_type,
    org = org,
    hasPeptideData = HasPeptideLevelData(),
    metadata = list(
      largestComponent = largest.comp.id,
      networkName = net.nm
    )
  )

  # Save to file
  output.path <- file.path(paste0(json.nm, ".json"))

  # msg(sprintf("[exportNetworkToJSON] Writing JSON...\n"))
  jsonlite::write_json(network.json, output.path, auto_unbox = TRUE, pretty = FALSE)
  return(invisible(NULL))
}



# zero-order network - create ppi nets from only input (seeds)

#' Create network from only input (seeds)
#'
#' @export
#'
BuildSeedProteinNet <- function(){
    cache <- .ppi_query_cache
    if (is.null(cache$seeds) || length(cache$seeds) == 0) {
        .setCurrentMessage("No seed proteins cached. Please run SearchNetDB first.")
        return(c(0, 0, 0))
    }

    analSet <- .ensurePpiList()

    # Get nodes from overall graph
    nodes <- V(overall.graph)$name

    # Filter to only seed proteins (from cache)
    hit.inx <- nodes %in% cache$seeds
    nodes2rm <- nodes[!hit.inx]

    # msg(sprintf("[PPI] Building zero-order network from %d seed proteins\n", sum(hit.inx)))

    # Delete non-seed nodes
    g <- igraph::simplify(igraph::delete_vertices(overall.graph, nodes2rm))

    if (igraph::vcount(g) == 0) {
        .setCurrentMessage("No seed-only network could be built.")
        return(c(0, 0, 0))
    }

    # Mark all nodes as query (since this is seeds-only)
    V(g)$is_query <- TRUE

    # msg(sprintf("[PPI] Seed-only network: %d nodes, %d edges\n",
    #            igraph::vcount(g), igraph::ecount(g)))

    # Store as overall graph
    overall.graph <<- g
    analSet$overall.graph <- g

    # Decompose into subnetworks
    analSet <- DecomposeGraph(g, analSet, minNodeNum = 2)

    if (is.null(analSet$substats) || length(analSet$substats) == 0) {
        .setCurrentMessage("No subnetworks found in seed-only network.")
        return(c(0, 0, 0))
    }

    # Update global variables
    net.stats <<- analSet$net.stats
    ppi.comps <<- analSet$ppi.comps

    # Set current network to first module
    if (length(ppi.comps) > 0) {
        current.net.nm <<- names(ppi.comps)[1]
    }

    num_seeds <- length(cache$seeds)
    num_nodes <- igraph::vcount(g)
    num_edges <- igraph::ecount(g)
    num_modules <- length(ppi.comps)

    # msg <- sprintf("Built seed-only network with %d nodes, %d edges. Decomposed into %d subnetwork(s).",
    #               num_nodes, num_edges, num_modules)
    .setCurrentMessage(msg)

    # msg(sprintf("[PPI] %s\n", msg))
    # msg(sprintf("[PPI] Subnetwork sizes: %s\n", paste(analSet$substats, collapse=", ")))

    saveSet(analSet, "analSet")

    return(c(num_nodes, num_edges, num_seeds, num_modules))
}

ComputePCSFNet <- function(){
    cache <- .ppi_query_cache
    if (is.null(cache$seeds) || length(cache$seeds) == 0) {
        .setCurrentMessage("No seed proteins cached. Please run SearchNetDB first.")
        return(c(0, 0, 0))
    }

    if (!exists("overall.graph", envir = .GlobalEnv) || is.null(overall.graph)) {
        .setCurrentMessage("No overall graph exists. Please run CreateGraph first.")
        return(c(0, 0, 0))
    }

    paramSet <- readSet(paramSet, "paramSet")
    analSet <- .ensurePpiList()

    # msg("[PPI] Computing Prize-Collecting Steiner Forest (PCSF)...\n")

    # Get edge list from overall graph
    edg <- as.data.frame(igraph::as_edgelist(overall.graph))
    edg$V3 <- rep(1, nrow(edg))
    colnames(edg) <- c("from", "to", "cost")

    # Create igraph from edge list
    node_names <- unique(c(as.character(edg[,1]), as.character(edg[,2])))
    ppi <- igraph::graph_from_data_frame(edg[,1:2], vertices=node_names, directed=FALSE)
    E(ppi)$weight <- as.numeric(edg[,3])
    ppi <- igraph::simplify(ppi)

    # msg(sprintf("[PPI] Input graph for PCSF: %d nodes, %d edges\n",
    #            igraph::vcount(ppi), igraph::ecount(ppi)))

    # Create prize vector: high prizes for seed nodes, low/zero for others
    # This is the key to PCSF - it finds the minimum cost subnetwork connecting high-prize nodes
    all.nodes <- V(ppi)$name
    prizes <- rep(0, length(all.nodes))  # Default prize = 0 for non-seed nodes
    names(prizes) <- all.nodes

    # Assign prizes to seed proteins only
    seed.indices <- all.nodes %in% cache$seeds

    # If expression values exist, use them as prizes for seeds
    if (exists("expr.vec", envir = .GlobalEnv) && !is.null(expr.vec) && sum(abs(expr.vec), na.rm=TRUE) > 0) {
        # msg("[PPI] Using expression values as prizes for seed proteins\n")
        # Map expression values to seed nodes
        for (seed in cache$seeds) {
            if (seed %in% names(expr.vec) && seed %in% all.nodes) {
                prizes[seed] <- abs(expr.vec[seed])  # Use absolute value as prize
            }
        }
        # Normalize to avoid extremely large values
        max.prize <- max(prizes, na.rm = TRUE)
        if (max.prize > 0) {
            prizes <- prizes / max.prize * 10  # Scale to 0-10 range
        }
    } else {
        # msg("[PPI] No expression values found, using uniform prizes for seeds\n")
        prizes[seed.indices] <- 1.0  # Uniform prize for all seeds
    }

    # Ensure seed proteins have non-zero prizes
    prizes[seed.indices][prizes[seed.indices] == 0] <- 1.0

    # msg(sprintf("[PPI] Assigned prizes to %d seed proteins (out of %d total nodes)\n",
    #            sum(seed.indices), length(all.nodes)))
    # msg(sprintf("[PPI] Prize range: %.3f to %.3f\n", min(prizes), max(prizes)))

    # Run Steiner Forest algorithm with seed-specific prizes
    g <- Compute.SteinerForest(ppi, prizes, w = 5, b = 100, mu = 0.0005)

    if (igraph::vcount(g) == 0) {
        .setCurrentMessage("PCSF algorithm returned empty network.")
        return(c(0, 0, 0))
    }

    # msg(sprintf("[PPI] PCSF result: %d nodes, %d edges\n",
    #             igraph::vcount(g), igraph::ecount(g)))

    # Mark query nodes (nodes that were in original seeds)
    node.names <- V(g)$name
    is.seed <- node.names %in% cache$seeds
    V(g)$is_query <- is.seed

    # msg(sprintf("[PPI] Marked %d/%d nodes as query nodes in PCSF network\n",
    #            sum(is.seed), length(node.names)))

    # Store as overall graph
    overall.graph <<- g
    analSet$overall.graph <- g

    # Decompose into subnetworks
    analSet <- DecomposeGraph(g, analSet, minNodeNum = 2)

    if (is.null(analSet$substats) || length(analSet$substats) == 0) {
        .setCurrentMessage("No subnetworks found in PCSF network.")
        return(c(0, 0, 0))
    }

    # Update global variables
    net.stats <<- analSet$net.stats
    ppi.comps <<- analSet$ppi.comps

    # Set current network to first module
    if (length(ppi.comps) > 0) {
        current.net.nm <<- names(ppi.comps)[1]
    }

    num_seeds <- length(cache$seeds)
    num_query_in_net <- sum(is.seed)
    num_nodes <- igraph::vcount(g)
    num_edges <- igraph::ecount(g)
    num_modules <- length(ppi.comps)

    msg <- sprintf("PCSF network built with %d nodes, %d edges from %d seeds. Decomposed into %d subnetwork(s).",
                   num_nodes, num_edges, num_query_in_net, num_modules)
    .setCurrentMessage(msg)

    # msg(sprintf("[PPI] %s\n", msg))
    # msg(sprintf("[PPI] Subnetwork sizes: %s\n", paste(analSet$substats, collapse=", ")))

    saveSet(analSet, "analSet")

    return(c(num_nodes, num_edges, num_query_in_net, num_modules))
}

GetMinConnectedGraphs <- function(max.len = 200){
    tryCatch({
        # cat("[PPI] GetMinConnectedGraphs called\n")
        # flush.console()

        # Handle NA parameter from Java call
        if (is.na(max.len) || is.null(max.len)) {
            max.len <- 200
        }
        # cat(sprintf("[PPI] max.len = %d\n", max.len))
        # flush.console()

        cache <- .ppi_query_cache
        if (is.null(cache$seeds) || length(cache$seeds) == 0) {
            # cat("[PPI] ERROR: No seeds in cache\n")
            # flush.console()
            .setCurrentMessage("No seed proteins cached. Please run SearchNetDB first.")
            return(c(0, 0, 0))
        }
        # cat(sprintf("[PPI] Found %d seeds in cache\n", length(cache$seeds)))
        # flush.console()

        if (!exists("overall.graph", envir = .GlobalEnv)) {
            # cat("[PPI] ERROR: overall.graph does not exist\n")
            # flush.console()
            .setCurrentMessage("No overall graph exists. Please run CreateGraph first.")
            return(c(0, 0, 0))
        }

        if (is.null(overall.graph) || igraph::vcount(overall.graph) == 0) {
            # cat("[PPI] ERROR: overall.graph is null or empty\n")
            # flush.console()
            .setCurrentMessage("Overall graph is empty. Please build a network first.")
            return(c(0, 0, 0))
        }
        # msg(sprintf("[PPI] overall.graph has %d nodes\n", igraph::vcount(overall.graph)))
        # flush.console()

        set.seed(8574)
        analSet <- .ensurePpiList()

        # msg("[PPI] Computing minimum connected network using shortest paths...\n")
        # cat("[PPI] Step 1: Getting seeds from cache\n")
        # flush.console()

        # Get seed proteins from cache, but only keep those present in the graph
        all.graph.nodes <- V(overall.graph)$name
        my.seeds <- cache$seeds[cache$seeds %in% all.graph.nodes]

        if (length(my.seeds) == 0) {
            # cat("[PPI] ERROR: No seed proteins found in current graph\n")
            # flush.console()
            .setCurrentMessage("No seed proteins are present in the current network.")
            return(c(0, 0, 0))
        }

        sd.len <- length(my.seeds)
        # cat(sprintf("[PPI] Step 2: Starting with %d seed proteins (from %d original seeds)\n",
        #            sd.len, length(cache$seeds)))
        # flush.console()

        # First trim overall.graph to remove non-seed nodes of degree 1
        # cat("[PPI] Step 3: Computing node degrees\n")
        # flush.console()
        dgrs <- igraph::degree(overall.graph)
        keep.inx <- dgrs > 1 | (names(dgrs) %in% my.seeds)
        nodes2rm <- V(overall.graph)$name[!keep.inx]

        cat(sprintf("[PPI] Step 4: Trimming graph (removing %d nodes)\n", length(nodes2rm)))
        flush.console()
        if (length(nodes2rm) > 0) {
            g.trimmed <- igraph::simplify(igraph::delete_vertices(overall.graph, nodes2rm))
            cat(sprintf("[PPI] Trimmed %d degree-1 non-seed nodes\n", length(nodes2rm)))
        } else {
            g.trimmed <- overall.graph
        }
        # flush.console()

    # Restrict operation if too many seeds (shortest path computation is expensive)
    if(sd.len > max.len){
        hit.inx <- names(dgrs) %in% my.seeds
        sd.dgrs <- dgrs[hit.inx]
        sd.dgrs <- rev(sort(sd.dgrs))

        # Limit to top max.len seeds by degree
        actual.len <- min(max.len, length(sd.dgrs))
        my.seeds <- names(sd.dgrs)[1:actual.len]
        sd.len <- actual.len

        msg <- sprintf("Minimum connected network computed using top %d seed proteins (by degree).", sd.len)
        .setCurrentMessage(msg)
        # msg(sprintf("[PPI] %s\n", msg))
    }else{
        msg <- sprintf("Minimum connected network computed using all %d seed proteins.", sd.len)
        .setCurrentMessage(msg)
        # msg(sprintf("[PPI] %s\n", msg))
    }

        # Calculate shortest paths between all pairs of seeds
        # cat("[PPI] Step 5: Computing shortest paths\n")
        # flush.console()
        paths.list <- list()
        # msg(sprintf("[PPI] Computing shortest paths between %d seeds...\n", sd.len))

        for(pos in 1:sd.len){
            if (pos < sd.len) {
                # Get paths from current seed to all remaining seeds
                remaining.seeds <- my.seeds[(pos+1):sd.len]
                paths.list[[pos]] <- igraph::shortest_paths(
                    g.trimmed,
                    from = my.seeds[pos],
                    to = remaining.seeds,
                    mode = "all"
                )$vpath
            }
        }

        cat("[PPI] Step 6: Extracting unique nodes from paths\n")
        flush.console()
        # Extract unique node indices from all paths
        nds.inxs <- unique(unlist(paths.list))

        # cat(sprintf("[PPI] Found %d unique node indices in paths\n", length(nds.inxs)))
        flush.console()

        if (length(nds.inxs) == 0) {
            # cat("[PPI] ERROR: No connecting paths found\n")
            flush.console()
            .setCurrentMessage("No connecting paths found between seed proteins.")
            return(c(0, 0, 0))
        }

    # Keep only nodes that are on shortest paths
    nodes2rm <- V(g.trimmed)$name[-nds.inxs]
    g <- igraph::simplify(igraph::delete_vertices(g.trimmed, nodes2rm))

    if (igraph::vcount(g) == 0) {
        .setCurrentMessage("Minimum connected network is empty.")
        return(c(0, 0, 0))
    }

    # msg(sprintf("[PPI] Minimum connected network: %d nodes, %d edges\n",
    #            igraph::vcount(g), igraph::ecount(g)))

    # Mark query nodes (nodes that were in original seeds)
    node.names <- V(g)$name
    is.seed <- node.names %in% cache$seeds
    V(g)$is_query <- is.seed

    # msg(sprintf("[PPI] Marked %d/%d nodes as query nodes\n",
    #            sum(is.seed), length(node.names)))

    # Store as overall graph
    overall.graph <<- g
    analSet$overall.graph <- g

    # Decompose into subnetworks
    analSet <- DecomposeGraph(g, analSet, minNodeNum = 2)

    if (is.null(analSet$substats) || length(analSet$substats) == 0) {
        .setCurrentMessage("No subnetworks found in minimum connected network.")
        return(c(0, 0, 0))
    }

    # Update global variables
    net.stats <<- analSet$net.stats
    ppi.comps <<- analSet$ppi.comps

    # Set current network to first module
    if (length(ppi.comps) > 0) {
        current.net.nm <<- names(ppi.comps)[1]
    }

    num_seeds <- length(cache$seeds)
    num_query_in_net <- sum(is.seed)
    num_nodes <- igraph::vcount(g)
    num_edges <- igraph::ecount(g)
    num_modules <- length(ppi.comps)

    msg <- sprintf("Minimum connected network: %d nodes, %d edges from %d seeds. Decomposed into %d subnetwork(s).",
                   num_nodes, num_edges, num_query_in_net, num_modules)
    .setCurrentMessage(msg)

    # msg(sprintf("[PPI] %s\n", msg))
    # msg(sprintf("[PPI] Subnetwork sizes: %s\n", paste(analSet$substats, collapse=", ")))

        saveSet(analSet, "analSet")

        return(c(num_nodes, num_edges, num_query_in_net, num_modules))

    }, error = function(e) {
        # msg(sprintf("[PPI] ERROR in GetMinConnectedGraphs: %s\n", conditionMessage(e)))
        .setCurrentMessage(paste("Error computing minimum connected network:", conditionMessage(e)))
        return(c(0, 0, 0))
    })
}


# Adapted from PCSF
# https://github.com/IOR-Bioinformatics/PCSF
Compute.SteinerForest <- function(ppi, terminals, w = 2, b = 1, mu = 0.0005, dummies){

  # Gather the terminal genes to be analyzed, and their scores
  terminal_names <- names(terminals)
  terminal_values <- as.numeric(terminals)

  # Incorporate the node prizes
  node_names <- V(ppi)$name
  node_prz <- vector(mode = "numeric", length = length(node_names))
  index <- match(terminal_names, node_names)
  percent <- signif((length(index) - sum(is.na(index)))/length(index)*100, 4)
  if (percent < 5){
    print("Less than 1% of your terminal nodes are matched in the interactome!");
    return(NULL);
  }
  paste0("  ", percent, "% of your terminal nodes are included in the interactome\n");
  terminal_names <- terminal_names[!is.na(index)]
  terminal_values <- terminal_values[!is.na(index)]
  index <- index[!is.na(index)]
  node_prz[index] <- terminal_values

  if(missing(dummies)||is.null(dummies)||is.na(dummies)){
    dummies <- terminal_names #re-assign this to allow for input
  }

  ## Prepare input file for MST-PCSF implementation in C++

  # Calculate the hub penalization scores
  node_degrees <- igraph::degree(ppi)
  hub_penalization <- - mu*node_degrees

  # Update the node prizes
  node_prizes <- b*node_prz
  index <- which(node_prizes==0)
  node_prizes[index] <- hub_penalization[index]

  # Construct the list of edges
  edges <- ends(ppi,es = E(ppi))
  from <- c(rep("DUMMY", length(dummies)), edges[,1])
  to <- c(dummies, edges[,2])

  cost <- c(rep(w, length(dummies)), E(ppi)$weight)

  #PCSF will faill if there are NAs in weights, this will check and fail gracefully
  if(any(is.na(E(ppi)$weight))){
    print("NAs found in the weight vector!");
    return (NULL);
  }

  ## Feed the input into the PCSF algorithm
  output <- XiaLabCppLib::call_sr(from,to,cost,node_names,node_prizes)

  # Check the size of output subnetwork and print a warning if it is 0
  if(length(output[[1]]) != 0){

    # Contruct an igraph object from the MST-PCSF output
    e <- data.frame(output[[1]], output[[2]], output[[3]])
    e <- e[which(e[,2]!="DUMMY"), ]
    names(e) <- c("from", "to", "weight")

    # Differentiate the type of nodes
    type <- rep("Steiner", length(output[[4]]))
    index <- match(terminal_names, output[[4]])
    index <- index[!is.na(index)]
    type[index] <- "Terminal"

    v <- data.frame(output[[4]], output[[5]], type)
    names(v) <- c("terminals", "prize", "type")
    subnet <- graph_from_data_frame(e,vertices=v,directed=F)
    #E(subnet)$weight <- as.numeric(output[[3]])
    subnet <- delete_vertices(subnet, "DUMMY")
    subnet <- delete_vertices(subnet, names(which(degree(subnet)==0)));
    return(subnet)

  } else{
    print("Subnetwork can not be identified for a given parameter set")
    return(NULL);
  }
}
