# Compute backbone layout for compartment networks and attach coordinates to JSON
# Expected inputs: elements (nodes/edges) list read from JSON; optional focus/algo

UpdateCompartmentLayout <- function(jsonNm = "proteoanalyst_1.json",
                                    algo = c("backbone", "fr", "kk")) {
  algo <- match.arg(algo)
  if (!file.exists(jsonNm)) {
    stop("JSON file not found: ", jsonNm)
  }
  graph <- jsonlite::read_json(jsonNm, simplifyVector = TRUE)
  if (is.null(graph$elements)) {
    stop("No elements found in JSON.")
  }

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph package is required for compartment layout.")
  }

  # Build igraph
  nodes <- graph$elements[sapply(graph$elements, function(e) e$group == "nodes")]
  edges <- graph$elements[sapply(graph$elements, function(e) e$group == "edges")]
  node_ids <- vapply(nodes, function(n) n$data$id, character(1))
  g <- make_empty_graph(n = length(nodes), directed = FALSE)
  V(g)$name <- node_ids

  if (length(edges)) {
    edge_pairs <- unlist(lapply(edges, function(e) c(e$data$source, e$data$target)))
    g <- add_edges(g, edge_pairs)
  }

  # Build membership by compartment for features; compartments themselves form their own group
  membership <- rep("ungrouped", length(nodes))
  names(membership) <- node_ids
  for (i in seq_along(nodes)) {
    if (!is.null(nodes[[i]]$data$type) && nodes[[i]]$data$type == "compartment") {
      membership[i] <- nodes[[i]]$data$id
    } else if (!is.null(nodes[[i]]$data$parent)) {
      membership[i] <- nodes[[i]]$data$parent
    }
  }
  comm <- igraph::make_clusters(g, membership = membership)

  edge.weights <- function(community, network, weight.within = 40, weight.between = 1) {
    bridges <- crossing(communities = community, graph = network)
    weights <- ifelse(test = bridges, yes = weight.between, no = weight.within)
    return(weights)
  }

  # Prefer backbone layout when requested; fall back to Fruchterman-Reingold with weights favoring intra-compartment cohesion
  coords <- switch(
    algo,
    backbone = {
      if (!requireNamespace("backbone", quietly = TRUE)) {
        layout_with_fr(g, weights = edge.weights(comm, g, 40, 1))
      } else {
        # use a sparse backbone approximation then layout_fr
        el <- igraph::as_data_frame(g, what = "edges")
        if (nrow(el) > 0) {
          mat <- Matrix::sparseMatrix(
            i = match(el$from, V(g)$name),
            j = match(el$to, V(g)$name),
            x = 1,
            dims = c(vcount(g), vcount(g))
          )
          bb <- try(backbone::sdsm(mat, signed = FALSE), silent = TRUE)
          if (inherits(bb, "try-error")) {
            layout_with_fr(g, weights = edge.weights(comm, g, 40, 1))
          } else {
            g2 <- graph_from_adjacency_matrix(as.matrix(bb), mode = "undirected", weighted = TRUE)
            layout_with_fr(g2, weights = igraph::E(g2)$weight)
          }
        } else {
          layout_with_fr(g, weights = edge.weights(comm, g, 40, 1))
        }
      }
    },
    fr = layout_with_fr(g, weights = edge.weights(comm, g, 40, 1)),
    kk = layout_with_kk(g)
  )

  # Normalize coordinates
  coords <- as.matrix(coords)
  coords <- scale(coords, center = TRUE, scale = FALSE)
  coords <- coords / max(abs(coords))

  # Attach coords back to nodes
  for (i in seq_along(nodes)) {
    nodes[[i]]$data$x <- coords[i, 1]
    nodes[[i]]$data$y <- coords[i, 2]
  }
  graph$elements[sapply(graph$elements, function(e) e$group == "nodes")] <- nodes

  outnm <- jsonNm
  jsonlite::write_json(graph, outnm, auto_unbox = TRUE, pretty = TRUE)
  return(outnm)
}
