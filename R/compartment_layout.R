# Compute backbone layout for compartment networks and attach coordinates to JSON
# Expected inputs: elements (nodes/edges) list read from JSON; optional focus/algo

UpdateCompartmentLayout <- function(jsonNm = "proteoanalyst_1.json",
                                    algo = c("backbone", "fr", "kk")) {
  algo <- match.arg(algo)
  if (!file.exists(jsonNm)) {
    AddErrMsg(paste0("JSON file not found: ", jsonNm));
    return(0);
  }
  graph <- jsonlite::read_json(jsonNm, simplifyVector = TRUE)
  if (is.null(graph$elements)) {
    AddErrMsg("No elements found in JSON.");
    return(0);
  }

  if (!requireNamespace("igraph", quietly = TRUE)) {
    AddErrMsg("igraph package is required for compartment layout.");
    return(0);
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


# Compute a deterministic 3D compartment-aware layout for localization/PPI JSON.
# The 2D fields posx/posy are preserved for the sigma viewer; the 3D viewer uses
# x/y/z and fixed coordinates fx/fy/fz.
UpdateCompartmentLayout3D <- function(jsonNm = "localization_network.json",
                                      algo = c("compartment", "sphere")) {
  algo <- match.arg(algo)
  if (!file.exists(jsonNm)) {
    AddErrMsg(paste0("JSON file not found: ", jsonNm))
    return("NA")
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    AddErrMsg("igraph package is required for 3D compartment layout.")
    return("NA")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    AddErrMsg("jsonlite package is required for 3D compartment layout.")
    return("NA")
  }

  graph <- jsonlite::read_json(jsonNm, simplifyVector = FALSE)
  if (is.null(graph$nodes) || length(graph$nodes) == 0) {
    AddErrMsg("No nodes found in network JSON.")
    return("NA")
  }

  `%||%` <- function(a, b) {
    if (is.null(a) || length(a) == 0 || is.na(a) || !nzchar(as.character(a))) b else a
  }

  sanitize_compartment <- function(x) {
    x <- as.character(x %||% "Unknown")
    x <- gsub("[^A-Za-z0-9_]", "_", x)
    if (!nzchar(x)) "Unknown" else x
  }

  safe_num <- function(x, default = 0) {
    val <- suppressWarnings(as.numeric(x))
    ifelse(is.finite(val), val, default)
  }

  normalize3d <- function(coords, radius = 1) {
    coords <- as.matrix(coords)
    if (ncol(coords) < 3) {
      coords <- cbind(coords, matrix(0, nrow = nrow(coords), ncol = 3 - ncol(coords)))
    }
    coords <- coords[, 1:3, drop = FALSE]
    coords[!is.finite(coords)] <- 0
    coords <- scale(coords, center = TRUE, scale = FALSE)
    max.dist <- max(sqrt(rowSums(coords^2)))
    if (!is.finite(max.dist) || max.dist == 0) {
      return(matrix(0, nrow = nrow(coords), ncol = 3))
    }
    coords / max.dist * radius
  }

  fibonacci_sphere <- function(n, radius = 1) {
    if (n <= 0) return(matrix(numeric(0), ncol = 3))
    if (n == 1) return(matrix(c(0, 0, 0), ncol = 3))
    i <- seq_len(n) - 1
    golden <- pi * (3 - sqrt(5))
    y <- 1 - (i / (n - 1)) * 2
    r <- sqrt(pmax(0, 1 - y * y))
    theta <- i * golden
    cbind(radius * cos(theta) * r, radius * y, radius * sin(theta) * r)
  }

  spacing_stats <- function(coords, min_gap) {
    n <- nrow(coords)
    nearest <- rep(Inf, n)
    close_pairs <- 0L
    if (n < 2) {
      return(list(close_pairs = 0L, tight_share = 0))
    }
    for (i in seq_len(n - 1)) {
      idx <- (i + 1):n
      diffs <- sweep(coords[idx, , drop = FALSE], 2, coords[i, ], "-")
      dists <- sqrt(rowSums(diffs * diffs))
      close_pairs <- close_pairs + sum(dists < min_gap, na.rm = TRUE)
      nearest[i] <- min(nearest[i], dists, na.rm = TRUE)
      nearest[idx] <- pmin(nearest[idx], dists)
    }
    nearest <- nearest[is.finite(nearest)]
    list(
      close_pairs = close_pairs,
      tight_share = if (length(nearest)) mean(nearest < min_gap) else 0
    )
  }

  deterministic_direction <- function(i, j, iter) {
    angle <- (i * 12.9898 + j * 78.233 + iter * 37.719) %% (2 * pi)
    z <- sin(i * 0.73 + j * 1.27 + iter * 0.19)
    r <- sqrt(max(0, 1 - z * z))
    c(cos(angle) * r, z, sin(angle) * r)
  }

  disperse_close_nodes <- function(coords, radius, min_gap, max_iter = 70) {
    coords <- as.matrix(coords)
    rn <- rownames(coords)
    n <- nrow(coords)
    if (n < 2 || n > 750 || !is.finite(min_gap) || min_gap <= 0) {
      return(list(coords = coords, radius = radius, adjusted = FALSE))
    }
    stats <- spacing_stats(coords, min_gap)
    if (stats$close_pairs == 0) {
      return(list(coords = coords, radius = radius, adjusted = FALSE))
    }

    work_radius <- max(radius, min(420, min_gap * (n ^ (1 / 3)) * 1.5))
    coords <- normalize3d(coords, radius = work_radius)

    for (iter in seq_len(max_iter)) {
      disp <- matrix(0, nrow = n, ncol = 3)
      close_pairs <- 0L

      for (i in seq_len(n - 1)) {
        idx <- (i + 1):n
        diffs <- sweep(coords[idx, , drop = FALSE], 2, coords[i, ], "-")
        dists <- sqrt(rowSums(diffs * diffs))
        close <- which(dists < min_gap)
        if (!length(close)) next

        close_pairs <- close_pairs + length(close)
        for (k in close) {
          j <- i + k
          d <- dists[[k]]
          direction <- if (is.finite(d) && d > 1e-6) diffs[k, ] / d else deterministic_direction(i, j, iter)
          push <- (min_gap - max(d, 1e-6)) * 0.5
          disp[i, ] <- disp[i, ] - direction * push
          disp[j, ] <- disp[j, ] + direction * push
        }
      }

      if (close_pairs == 0) break

      step_len <- sqrt(rowSums(disp * disp))
      max_step <- min_gap * 0.45
      step_scale <- ifelse(step_len > max_step & step_len > 0, max_step / step_len, 1)
      coords <- coords + disp * step_scale
      coords <- scale(coords, center = TRUE, scale = FALSE)

      dist_from_origin <- sqrt(rowSums(coords * coords))
      outside <- dist_from_origin > work_radius
      if (any(outside)) {
        coords[outside, ] <- coords[outside, , drop = FALSE] / dist_from_origin[outside] * work_radius
      }

      if (iter %% 15 == 0 && close_pairs > n) {
        work_radius <- min(420, work_radius * 1.08)
        coords <- normalize3d(coords, radius = work_radius)
      }
    }

    rownames(coords) <- rn
    actual_radius <- max(sqrt(rowSums(coords * coords)))
    list(coords = coords, radius = max(radius, actual_radius), adjusted = TRUE)
  }


  node_ids <- vapply(graph$nodes, function(n) as.character(n$id %||% ""), character(1))
  is_comp_node <- vapply(graph$nodes, function(n) {
    identical(as.character(n$type %||% ""), "compartment") ||
      identical(as.character(n$molType %||% ""), "compartment")
  }, logical(1))
  gene_idx <- which(!is_comp_node & nzchar(node_ids))
  if (length(gene_idx) == 0) {
    AddErrMsg("No protein/gene nodes found for 3D compartment layout.")
    return("NA")
  }

  gene_ids <- node_ids[gene_idx]

  # Derive a normalized category label for each node.
  # broad_category is already normalized by the R export; fall back to
  # normalizing the raw compartment string so that different raw IDs that
  # represent the same standard category (e.g. "Cell_membrane__Multi_pass"
  # and "Unknown") are merged into one compartment group.
  norm_category <- function(n) {
    bc <- as.character(n$broad_category %||% "")
    if (nzchar(bc) && bc != "Unknown") return(bc)
    raw <- as.character(n$compartment %||% n$location %||% "")
    raw <- gsub("_", " ", raw)
    raw <- trimws(raw)
    if (!nzchar(raw) || raw == "Unknown") return("Unknown")
    norm <- .normalizeBroadCategory(raw)
    if (length(norm) == 1 && nzchar(norm) && norm != "Unknown") norm else "Unknown"
  }

  comp_labels_raw <- vapply(graph$nodes[gene_idx], norm_category, character(1))
  # Use the normalized label as the grouping key (already human-readable)
  comp_ids <- comp_labels_raw
  names(comp_ids) <- gene_ids

  comp_labels <- split(comp_labels_raw, comp_ids)
  comp_labels <- vapply(comp_labels, function(x) x[[1]], character(1))

  comp_members <- split(gene_ids, comp_ids)
  comp_levels <- names(sort(vapply(comp_members, length, integer(1)), decreasing = TRUE))

  edge_df <- data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE)
  if (!is.null(graph$edges) && length(graph$edges) > 0) {
    edge_df <- do.call(rbind, lapply(graph$edges, function(e) {
      data.frame(
        from = as.character(e$source %||% ""),
        to = as.character(e$target %||% ""),
        stringsAsFactors = FALSE
      )
    }))
    edge_df <- edge_df[edge_df$from %in% gene_ids & edge_df$to %in% gene_ids & edge_df$from != edge_df$to, , drop = FALSE]
  }

  if (nrow(edge_df) > 0) {
    gene_graph <- igraph::graph_from_data_frame(edge_df, directed = FALSE, vertices = data.frame(name = gene_ids))
    gene_graph <- igraph::simplify(gene_graph, remove.multiple = TRUE, remove.loops = TRUE)
  } else {
    gene_graph <- igraph::make_empty_graph(n = length(gene_ids), directed = FALSE)
    igraph::V(gene_graph)$name <- gene_ids
  }

  comp_centers <- matrix(0, nrow = length(comp_levels), ncol = 3, dimnames = list(comp_levels, c("x", "y", "z")))
  center_radius <- max(390, min(960, 205 * sqrt(length(comp_levels))))
  if (length(comp_levels) > 1) {
    inter_edges <- data.frame(from = character(0), to = character(0), weight = numeric(0), stringsAsFactors = FALSE)
    if (nrow(edge_df) > 0) {
      edge_comp <- data.frame(
        from = comp_ids[edge_df$from],
        to = comp_ids[edge_df$to],
        stringsAsFactors = FALSE
      )
      edge_comp <- edge_comp[edge_comp$from != edge_comp$to, , drop = FALSE]
      if (nrow(edge_comp) > 0) {
        edge_comp$key <- apply(edge_comp, 1, function(x) paste(sort(x), collapse = "||"))
        tab <- sort(table(edge_comp$key), decreasing = TRUE)
        inter_edges <- do.call(rbind, lapply(names(tab), function(key) {
          parts <- strsplit(key, "\\|\\|")[[1]]
          data.frame(from = parts[1], to = parts[2], weight = as.numeric(tab[[key]]), stringsAsFactors = FALSE)
        }))
      }
    }

    if (algo == "compartment" && nrow(inter_edges) > 0) {
      comp_graph <- igraph::graph_from_data_frame(
        inter_edges,
        directed = FALSE,
        vertices = data.frame(name = comp_levels, stringsAsFactors = FALSE)
      )
      set.seed(42)
      center_coords <- tryCatch({
        igraph::layout_with_fr(
          comp_graph,
          weights = 1 + log1p(igraph::E(comp_graph)$weight),
          niter = 1000,
          dim = 3,
          grid = "nogrid"
        )
      }, error = function(e) {
        fibonacci_sphere(length(comp_levels), radius = 1)
      })
      comp_centers[igraph::V(comp_graph)$name, ] <- normalize3d(center_coords, radius = center_radius)
    } else {
      comp_centers[, ] <- fibonacci_sphere(length(comp_levels), radius = center_radius)
    }
  }

  coords_by_id <- matrix(0, nrow = length(gene_ids), ncol = 3, dimnames = list(gene_ids, c("x", "y", "z")))
  comp_radius <- setNames(numeric(length(comp_levels)), comp_levels)

  set.seed(42)
  for (comp_id in comp_levels) {
    ids <- comp_members[[comp_id]]
    m <- length(ids)
    local_radius <- max(28, min(150, 13.5 * sqrt(m)))

    if (m == 1) {
      local_coords <- matrix(c(0, 0, 0), ncol = 3)
    } else {
      sub_graph <- igraph::induced_subgraph(gene_graph, vids = ids)
      if (igraph::ecount(sub_graph) > 0 && m <= 600) {
        local_coords <- tryCatch({
          igraph::layout_with_fr(
            sub_graph,
            niter = ifelse(m > 250, 250, 600),
            dim = 3,
            start.temp = sqrt(m) * 1.5,
            grid = "nogrid"
          )
        }, error = function(e) {
          fibonacci_sphere(m, radius = 1)
        })
        local_coords <- normalize3d(local_coords, radius = local_radius)
        rownames(local_coords) <- igraph::V(sub_graph)$name
        local_coords <- local_coords[ids, , drop = FALSE]
      } else {
        local_coords <- fibonacci_sphere(m, radius = local_radius)
        rownames(local_coords) <- ids
      }
    }

    if (m > 2) {
      min_gap <- max(18, min(34, 9 + 2.7 * log1p(m)))
      spacing <- disperse_close_nodes(local_coords, local_radius, min_gap)
      local_coords <- spacing$coords
      local_radius <- max(local_radius, spacing$radius)
    }
    comp_radius[[comp_id]] <- max(28, local_radius * 1.08)

    center <- comp_centers[comp_id, ]
    coords_by_id[ids, ] <- sweep(local_coords, 2, center, "+")
  }

  for (idx in gene_idx) {
    id <- node_ids[[idx]]
    coord <- coords_by_id[id, ]
    graph$nodes[[idx]]$x <- round(coord[[1]], 2)
    graph$nodes[[idx]]$y <- round(coord[[2]], 2)
    graph$nodes[[idx]]$z <- round(coord[[3]], 2)
    graph$nodes[[idx]]$fx <- round(coord[[1]], 2)
    graph$nodes[[idx]]$fy <- round(coord[[2]], 2)
    graph$nodes[[idx]]$fz <- round(coord[[3]], 2)
    graph$nodes[[idx]]$posz <- round(coord[[3]], 2)
  }

  comp_volumes <- lapply(comp_levels, function(comp_id) {
    center <- comp_centers[comp_id, ]
    list(
      id = comp_id,
      label = as.character(comp_labels[[comp_id]] %||% comp_id),
      color = as.character(graph$nodes[[gene_idx[match(comp_members[[comp_id]][1], gene_ids)]]]$colorb %||% "#999999"),
      x = round(center[[1]], 2),
      y = round(center[[2]], 2),
      z = round(center[[3]], 2),
      radius = round(comp_radius[[comp_id]], 2),
      nodes = as.list(comp_members[[comp_id]])
    )
  })
  names(comp_volumes) <- comp_levels

  comp_node_idx <- which(is_comp_node & node_ids %in% comp_levels)
  for (idx in comp_node_idx) {
    comp_id <- node_ids[[idx]]
    center <- comp_centers[comp_id, ]
    graph$nodes[[idx]]$x <- round(center[[1]], 2)
    graph$nodes[[idx]]$y <- round(center[[2]], 2)
    graph$nodes[[idx]]$z <- round(center[[3]], 2)
    graph$nodes[[idx]]$fx <- round(center[[1]], 2)
    graph$nodes[[idx]]$fy <- round(center[[2]], 2)
    graph$nodes[[idx]]$fz <- round(center[[3]], 2)
    graph$nodes[[idx]]$radius3d <- round(comp_radius[[comp_id]], 2)
  }

  if (!is.null(graph$compartments) && length(graph$compartments) > 0) {
    compartment_names <- names(graph$compartments)
    if (is.null(compartment_names)) {
      compartment_names <- rep("", length(graph$compartments))
    }
    for (i in seq_along(graph$compartments)) {
      fallback_id <- if (length(compartment_names) >= i) compartment_names[[i]] else ""
      comp_id <- sanitize_compartment(graph$compartments[[i]]$id %||% fallback_id %||% "")
      if (comp_id %in% names(comp_volumes)) {
        graph$compartments[[i]]$x <- comp_volumes[[comp_id]]$x
        graph$compartments[[i]]$y <- comp_volumes[[comp_id]]$y
        graph$compartments[[i]]$z <- comp_volumes[[comp_id]]$z
        graph$compartments[[i]]$radius <- comp_volumes[[comp_id]]$radius
        graph$compartments[[i]]$nodes <- comp_volumes[[comp_id]]$nodes
      }
    }
  } else {
    graph$compartments <- comp_volumes
  }

  if (is.null(graph$metadata)) graph$metadata <- list()
  graph$metadata$layout3d <- list(
    type = "compartment3d",
    algorithm = algo,
    fixed = TRUE,
    compartmentCount = length(comp_levels),
    compartments = comp_volumes
  )
  graph$layout3d <- list(type = "compartment3d", fixed = TRUE)
  graph$threed <- TRUE

  jsonlite::write_json(graph, jsonNm, auto_unbox = TRUE, pretty = TRUE)
  return(jsonNm)
}
