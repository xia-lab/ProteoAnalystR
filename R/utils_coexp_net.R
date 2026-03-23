
# ====================================================================
# BuildIgraphFromCEM  —  make a data-derived gene–gene network
# --------------------------------------------------------------------
# file        : path to the qs-saved CEMiTool object (“cem.qs”)
# thresh      : numeric, keep edges with weight > thresh
# return      : igraph object with vertex/edge attributes
# ====================================================================
BuildIgraphFromCEM <- function(thresh    = 0.05,
                               layoutFun = igraph::layout_nicely) {

  ## ── 0 · packages ─────────────────────────────────────────────────
  library(CEMiTool)
  library(igraph)
  library(reshape2)

  ## ── 1 · read the CEMiTool object ────────────────────────────────
  if (!file.exists("cem.qs")) {
    stop("cem.qs file not found. Run CEMiTool analysis first.");
  }

  cem <- tryCatch({
    qs::qread("cem.qs")
  }, error = function(e) {
    stop(paste("Failed to read cem.qs:", conditionMessage(e)));
  })
  
  ## ── 2 · make sure we have an adjacency matrix -------------------
  ##     (CEMiTool stores β only if you explicitly asked for it.)
  get_beta <- function(cem) {
    # 1) try stored β
    if (!is.null(cem@parameters$beta))
      return(as.numeric(cem@parameters$beta))
    
    # 2) try default scale-free heuristic (≥ v1.29)
    beta <- tryCatch({
      get_cemitool_r2_beta(cem)[2]   # returns c(R2, β)
    }, error = function(e) NA)
    
    # 3) last-resort: WGCNA pickSoftThreshold on the expression matrix
    if (is.na(beta)) {
      expr <- CEMiTool:::get_expression(cem)   # matrix features × samples
      sft  <- WGCNA::pickSoftThreshold(t(expr), verbose = 0)
      beta <- sft$powerEstimate
      if (is.na(beta)) beta <- 6              # fallback default
    }
    beta
  }
  
  adj_missing <- is.null(cem@adjacency) ||
                 (is.matrix(cem@adjacency) && nrow(cem@adjacency) == 0)

  if (adj_missing) {
    beta <- get_beta(cem)
    cem  <- get_adj(cem, beta = beta)
  }
  adj <- adj_data(cem)                        # square features × features

  cat(sprintf("[BuildIgraphFromCEM] Adjacency matrix: %d x %d\n", nrow(adj), ncol(adj)))
  cat(sprintf("[BuildIgraphFromCEM] Sample adj rownames: %s\n", paste(head(rownames(adj), 5), collapse=", ")))
  cat(sprintf("[BuildIgraphFromCEM] Sample adj colnames: %s\n", paste(head(colnames(adj), 5), collapse=", ")))

  # Ensure adjacency matrix has proper feature names
  cat(sprintf("[BuildIgraphFromCEM] Adjacency has rownames: %s, colnames: %s\n",
              !is.null(rownames(adj)), !is.null(colnames(adj))))

  if(is.null(rownames(adj)) || is.null(colnames(adj))) {
    cat("[BuildIgraphFromCEM] WARNING: Adjacency matrix missing row/col names, setting from cem@module features\n")
    # Use the features from cem@module which should match the adjacency matrix
    rownames(adj) <- colnames(adj) <- cem@module$features
    cat(sprintf("[BuildIgraphFromCEM] After setting, sample adj rownames: %s\n", paste(head(rownames(adj), 5), collapse=", ")))
  }

  ## ── 3 · build edge list above threshold -------------------------
  edge.df <- melt(adj)
  # Convert factors to characters to preserve names correctly
  edge.df$Var1 <- as.character(edge.df$Var1)
  edge.df$Var2 <- as.character(edge.df$Var2)

  cat(sprintf("[BuildIgraphFromCEM] After melt, sample edge names: %s <-> %s\n",
              edge.df$Var1[1], edge.df$Var2[1]))

  edge.df <- subset(edge.df, value > thresh & Var1 != Var2)

  ## Check if we have any edges
  if (nrow(edge.df) == 0) {
    stop(paste0("No edges found with threshold > ", thresh, ". Try lowering the threshold or check if CEMiTool found modules."));
  }

  ## keep at most 2 000 heaviest edges
  if (nrow(edge.df) > 2000) {
    edge.df <- edge.df[order(edge.df$value, decreasing = TRUE), ]
    edge.df <- edge.df[1:2000, ]
  }

  g       <- graph_from_data_frame(edge.df, directed = FALSE)
  E(g)$weight <- edge.df$value
  
  
  
  ## ── 4 · add vertex-level annotations ----------------------------
  mod.df <- cem@module                       # cols: features, modules

  # Debug: check module assignments
  cat(sprintf("[BuildIgraphFromCEM] Total features in cem@module: %d\n", nrow(mod.df)))
  cat(sprintf("[BuildIgraphFromCEM] cem@module column names: %s\n", paste(colnames(mod.df), collapse=", ")))
  cat(sprintf("[BuildIgraphFromCEM] Unique modules: %s\n", paste(unique(mod.df$modules), collapse=", ")))
  cat(sprintf("[BuildIgraphFromCEM] Graph has %d vertices\n", vcount(g)))
  cat(sprintf("[BuildIgraphFromCEM] Sample graph vertex names: %s\n", paste(head(V(g)$name, 5), collapse=", ")))

  # Check which column has the feature names
  feature_col <- if("genes" %in% colnames(mod.df)) {
    "genes"
  } else if("features" %in% colnames(mod.df)) {
    "features"
  } else {
    # Use rownames if no explicit column
    NULL
  }

  if(is.null(feature_col)) {
    cat("[BuildIgraphFromCEM] Using rownames of cem@module as features\n")
    features <- rownames(mod.df)
  } else {
    cat(sprintf("[BuildIgraphFromCEM] Using column '%s' as features\n", feature_col))
    features <- mod.df[[feature_col]]
  }

  cat(sprintf("[BuildIgraphFromCEM] Sample cem@module features: %s\n", paste(head(features, 5), collapse=", ")))

  idx    <- match(V(g)$name, features)
  cat(sprintf("[BuildIgraphFromCEM] Matched vertices: %d out of %d\n", sum(!is.na(idx)), vcount(g)))
  V(g)$module <- mod.df$modules[idx]

  # Check for NA modules
  na_count <- sum(is.na(V(g)$module))
  if(na_count > 0) {
    cat(sprintf("[BuildIgraphFromCEM] WARNING: %d vertices have NA module assignments\n", na_count))
    # Assign unmatched vertices to a default module
    V(g)$module[is.na(V(g)$module)] <- "Not.Assigned"
  }

  cat(sprintf("[BuildIgraphFromCEM] Final module distribution:\n"))
  print(table(V(g)$module))

deg   <- igraph::degree(g)
ldeg  <- log10(deg + 1)                       # stabilise high degrees

# rescale helper
resc <- function(x) (x - min(x)) / (max(x) - min(x))

# colour ramp: yellow → dark red (works on white & black)
pal <- colorRampPalette(c("#FFD54F", "#FFA726", "#EF5350", "#B71C1C"))(10)

V(g)$color  <- pal[ ceiling( resc(ldeg) * 9 ) + 1 ]   # for light bg
V(g)$colorw <- V(g)$color                             # same for dark bg
  
  ## size by degree
  rescale <- function(x, from = 8, to = 20)
    (x - min(x)) / (max(x) - min(x)) * (to - from) + from
  V(g)$size <- rescale(log10(degree(g) + 1))
  #V(g)$size <- 8;

  ## ── 5 · 2-D layout coordinates ----------------------------------
  xy <- layoutFun(g)
  V(g)$posx <- xy[, 1]
  V(g)$posy <- xy[, 2]
  analSet <- qs::qread("analSet.qs");
  analSet$overall.graph <- g;
  analSet$overall.graph.orig <- g;

  saveSet(analSet, "analSet");

  return(1)
}

CorrIgraph2SigmaJS <- function(g,
                               netNm       = "coexp_net",
                               paramSet,
                               analSet,
                               projectPPI  = FALSE,
                               ppiDatabase = "") {

  nms <- V(g)$name

  # Detect ID type and apply appropriate mapping
  is.entrez <- all(grepl("^[0-9]+$", nms))
  is.uniprot <- any(grepl("^[A-Z][0-9][A-Z0-9]{3,4}[0-9]$", nms))

  # Get gene symbols for labels
  # First try to get from vertex attributes (if already set by ExtractModule)
  symVec <- V(g)$gene_symbol

  # If not available in attributes, try mapping
  if (is.null(symVec) || all(is.na(symVec))) {
    msg("[CorrIgraph2SigmaJS] No gene_symbol attribute, attempting ID mapping...")
    symVec <- tryCatch({
      if (is.uniprot) {
        msg("[CorrIgraph2SigmaJS] Detected UniProt IDs, using doUniprot2SymbolMapping")
        doUniprot2SymbolMapping(nms, paramSet$data.org, paramSet$data.idType)
      } else {
        msg("[CorrIgraph2SigmaJS] Detected Entrez IDs, using doEntrez2SymbolMapping")
        doEntrez2SymbolMapping(nms, paramSet$data.org, paramSet$data.idType)
      }
    }, error = function(e) {
      msg("[CorrIgraph2SigmaJS] WARNING: Symbol mapping failed: ", e$message)
      return(nms)  # Use IDs as labels if mapping fails
    })
  } else {
    msg("[CorrIgraph2SigmaJS] Using gene_symbol from vertex attributes")
  }

  # Ensure symVec has same length as nms
  if (is.null(symVec) || length(symVec) == 0 || length(symVec) != length(nms)) {
    msg("[CorrIgraph2SigmaJS] WARNING: Symbol mapping returned incorrect length, using IDs as labels")
    symVec <- nms
  }

  # Log first few symbols for debugging
  msg("[CorrIgraph2SigmaJS] First 5 symbols: ", paste(head(symVec, 5), collapse=", "))

  # Load UniProt and Entrez mappings from symbol.map.qs
  uniprot.vec <- rep(NA_character_, length(nms))
  entrez.vec <- rep(NA_character_, length(nms))
  symbol.vec.from.map <- rep(NA_character_, length(nms))
  #cat(sprintf("[CorrIgraph2SigmaJS] DEBUG: Checking for symbol.map.qs file...\n"))
  #cat(sprintf("[CorrIgraph2SigmaJS] DEBUG: Current directory: %s\n", getwd()))
  #cat(sprintf("[CorrIgraph2SigmaJS] DEBUG: symbol.map.qs exists: %s\n", file.exists("symbol.map.qs")))

  if (file.exists("symbol.map.qs")) {
    tryCatch({
      gene.map <- readDataQs("symbol.map.qs", paramSet$anal.type, paramSet$dataName)
      msg("[CorrIgraph2SigmaJS] Loaded symbol.map.qs with columns: ", paste(colnames(gene.map), collapse=", "))

      if (!is.null(gene.map) && is.data.frame(gene.map) &&
          "gene_id" %in% colnames(gene.map) && "uniprot" %in% colnames(gene.map)) {

        if (is.uniprot) {
          # Input IDs are UniProt - use them directly, and get Entrez and Symbol from mapping
          uniprot.vec <- nms
          hit.inx <- match(as.character(nms), as.character(gene.map$uniprot))
          entrez.vec <- gene.map$gene_id[hit.inx]

          # Get symbols from mapping if available
          if ("symbol" %in% colnames(gene.map)) {
            symbol.vec.from.map <- gene.map$symbol[hit.inx]
            msg("[CorrIgraph2SigmaJS] Mapped ", sum(!is.na(symbol.vec.from.map)), "/", length(nms), " UniProt to symbols from symbol.map.qs")
          }

          msg("[CorrIgraph2SigmaJS] Using UniProt IDs directly, mapped ", sum(!is.na(entrez.vec)), "/", length(nms), " to Entrez")
        } else {
          # Input IDs are Entrez - use them directly, and get UniProt and Symbol from mapping
          entrez.vec <- nms
          hit.inx <- match(as.character(nms), as.character(gene.map$gene_id))
          uniprot.vec <- gene.map$uniprot[hit.inx]

          # Get symbols from mapping if available
          if ("symbol" %in% colnames(gene.map)) {
            symbol.vec.from.map <- gene.map$symbol[hit.inx]
            msg("[CorrIgraph2SigmaJS] Mapped ", sum(!is.na(symbol.vec.from.map)), "/", length(nms), " Entrez to symbols from symbol.map.qs")
          }

          msg("[CorrIgraph2SigmaJS] Using Entrez IDs directly, mapped ", sum(!is.na(uniprot.vec)), "/", length(nms), " to UniProt")
        }
      }
    }, error = function(e) {
      msg("[CorrIgraph2SigmaJS] ERROR: Could not load ID mapping: ", e$message)
    })
  } else {
    # No mapping file - use what we have
    if (is.uniprot) {
      uniprot.vec <- nms
      msg("[CorrIgraph2SigmaJS] No symbol.map.qs found, using UniProt IDs as-is")
    } else {
      entrez.vec <- nms
      msg("[CorrIgraph2SigmaJS] No symbol.map.qs found, using Entrez IDs as-is")
    }
  }

  # If we got symbols from symbol.map.qs and the original symVec mapping failed or is incomplete,
  # use the symbols from symbol.map.qs instead
  if (!is.null(symbol.vec.from.map) && sum(!is.na(symbol.vec.from.map)) > 0) {
    # Count how many valid symbols we have from each source
    valid.from.map <- sum(!is.na(symbol.vec.from.map) & nzchar(symbol.vec.from.map))
    valid.from.mapping <- sum(!is.na(symVec) & nzchar(symVec) & symVec != nms)

    msg("[CorrIgraph2SigmaJS] Symbol comparison: symbol.map.qs=", valid.from.map,
        " vs mapping functions=", valid.from.mapping)

    # If symbol.map.qs has more valid symbols, use it
    if (valid.from.map > valid.from.mapping) {
      msg("[CorrIgraph2SigmaJS] Using symbols from symbol.map.qs (better coverage)")
      symVec <- symbol.vec.from.map
      # Fill in any NAs with original IDs
      na.inx <- is.na(symVec) | !nzchar(symVec)
      if (any(na.inx)) {
        symVec[na.inx] <- nms[na.inx]
      }
    }
  }

  # Final check: ensure no NAs in symVec
  na.inx <- is.na(symVec) | !nzchar(symVec)
  if (any(na.inx)) {
    msg("[CorrIgraph2SigmaJS] Replacing ", sum(na.inx), " NA/empty symbols with IDs")
    symVec[na.inx] <- nms[na.inx]
  }

  msg("[CorrIgraph2SigmaJS] Final symbols (first 5): ", paste(head(symVec, 5), collapse=", "))

  # Get phosphosite mapping if available (for phosphoproteomics data)
  phosphosite.map <- list()
  if (!is.null(paramSet$phospho.mapping) && !is.null(paramSet$phospho.mapping$entrez.to.phospho)) {
    phosphosite.map <- paramSet$phospho.mapping$entrez.to.phospho
    msg("[CorrIgraph2SigmaJS] Found phosphosite mapping for ", length(phosphosite.map), " proteins")
  }

  # Define compartment colors (matching PPI network export)
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

  nodes <- lapply(seq_len(vcount(g)), function(i) {
    v   <- V(g)[i]
    lbl <- if (!is.na(symVec[i]) && nzchar(symVec[i])) symVec[i] else v$name

    # Log label assignment for debugging
    if (i <= 3) {
      msg(sprintf("[CorrIgraph2SigmaJS] Node %d: v$name=%s, symVec[i]=%s, lbl=%s",
                  i, v$name, symVec[i], lbl))
    }

    # Determine primary ID based on input type
    # Use UniProt as primary ID if available (stable, unique identifier)
    # Gene symbol goes in label field for display
    primary.id <- if (!is.na(uniprot.vec[i]) && !is.null(uniprot.vec[i]) && nzchar(uniprot.vec[i])) {
      as.character(uniprot.vec[i])
    } else {
      as.character(v$name)
    }

    # Populate both ID types
    uniprot.id <- if (!is.na(uniprot.vec[i]) && !is.null(uniprot.vec[i]) && nzchar(uniprot.vec[i])) {
      as.character(uniprot.vec[i])
    } else {
      ""
    }

    entrez.id <- if (!is.na(entrez.vec[i]) && !is.null(entrez.vec[i]) && nzchar(entrez.vec[i])) {
      as.character(entrez.vec[i])
    } else {
      ""
    }

    # Get compartment info from vertex attributes (if available)
    broad_cat <- if (!is.null(v$broad_category)) as.character(v$broad_category) else "Unknown"
    main_loc <- if (!is.null(v$main_location)) as.character(v$main_location) else "Unknown"
    comp.id <- gsub("[^A-Za-z0-9_]", "_", broad_cat)
    comp.color <- category.colors[[broad_cat]]
    if (is.null(comp.color)) comp.color <- "#999999"

    # Use vertex attributes if they exist (set by ExtractModule), otherwise use compartment color as default
    node.size <- if (!is.null(v$size)) unclass(v$size)[1] else 2
    node.colorb <- if (!is.null(v$color)) as.character(v$color) else comp.color
    node.colorw <- if (!is.null(v$colorw)) as.character(v$colorw) else comp.color
    node.expr <- if (!is.null(v$expr)) unclass(v$expr)[1] else 0

    # Get expression colors (if available from ExtractModule)
    exp.colb <- if (!is.null(v$expcolb)) as.character(v$expcolb) else comp.color
    exp.colw <- if (!is.null(v$expcolw)) as.character(v$expcolw) else comp.color
    exp.colc <- if (!is.null(v$expcolc)) as.character(v$expcolc) else comp.color

    # Get topology colors (if available from ExtractModule)
    topo.colb <- if (!is.null(v$topocolb)) as.character(v$topocolb) else comp.color
    topo.colw <- if (!is.null(v$topocolw)) as.character(v$topocolw) else comp.color
    topo.colc <- if (!is.null(v$topocolc)) as.character(v$topocolc) else comp.color

    # Get phosphosite IDs for this protein (if phospho data)
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
      id        = primary.id,             # UniProt ID (stable identifier)
      label     = lbl,                    # Gene symbol (for display)
      uniprot   = uniprot.id,             # UniProt accession
      entrez    = entrez.id,              # Entrez gene ID
      size      = node.size,
      true_size = node.size,
      molType   = "gene",
      type      = "gene",
      colorb    = node.colorb,            # Default: compartment color
      colorw    = node.colorw,            # Default: compartment color
      expcolb   = exp.colb,               # Alternative: expression color
      expcolw   = exp.colw,               # Alternative: expression color
      expcolc   = exp.colc,               # Alternative: expression color
      topocolb  = topo.colb,              # Alternative: topology color
      topocolw  = topo.colw,              # Alternative: topology color
      topocolc  = topo.colc,              # Alternative: topology color
      exp       = node.expr,
      posx      = unclass(v$posx)[1],
      posy      = unclass(v$posy)[1],
      compartment = comp.id,
      location = main_loc,
      main_location = main_loc,
      broad_category = broad_cat
    )

    # Add phosphosites field if available
    if (!is.null(phosphosites)) {
      node.data$phosphosites <- phosphosites
    }

    node.data
  })
  
  
  ## ── edges with rescaled size 0.5–2.5  ----------------------------
  el <- igraph::as_data_frame(g, what = "edges")       # from, to, weight

  wMin <- min(el$weight)
  wMax <- max(el$weight)
  rescale <- function(x, from = 0.5, to = 2.5) {
    if (wMax == wMin) return((from + to) / 2)          # avoid 0/0
    (x - wMin) / (wMax - wMin) * (to - from) + from
  }

  # Query PPI database if projection is enabled
  ppi.edges <- NULL
  if (projectPPI && nzchar(ppiDatabase)) {
    msg("[CorrIgraph2SigmaJS] Querying PPI database: ", ppiDatabase)
    tryCatch({
      # Get SQLite path and organism from paramSet
      sqlite.path <- paramSet$sqlite.path
      org <- if (!is.null(paramSet$data.org)) paramSet$data.org else "hsa"

      if (is.null(sqlite.path) || !nzchar(sqlite.path)) {
        msg("[CorrIgraph2SigmaJS] ERROR: sqlite.path not set in paramSet")
        ppi.edges <- NULL
      } else {
        # Construct table name: organism_database (e.g., "hsa_string")
        table.nm <- paste0(org, "_", ppiDatabase)
        msg("[CorrIgraph2SigmaJS] Table name: ", table.nm)

        # PPI database uses Entrez IDs
        # Co-expression network vertex names are the feature IDs from the uploaded data
        # Need to ensure they're in Entrez format for PPI matching
        node.ids <- V(g)$name

        # Check if IDs are already Entrez-like (all numeric)
        is.entrez <- all(grepl("^[0-9]+$", node.ids))
        msg("[CorrIgraph2SigmaJS] Node IDs are Entrez-like: ", is.entrez)

        if (!is.entrez) {
          # IDs might be UniProt accessions, gene symbols, or other types
          # Try to convert to Entrez for PPI database matching
          msg("[CorrIgraph2SigmaJS] Converting node IDs to Entrez for PPI matching...")
          msg("[CorrIgraph2SigmaJS] Sample node IDs: ", paste(head(node.ids, 5), collapse=", "))

          # Check if IDs look like UniProt accessions (pattern: letter + digits + letters/digits)
          is.uniprot <- any(grepl("^[A-Z][0-9][A-Z0-9]{3,4}[0-9]$", node.ids))
          msg("[CorrIgraph2SigmaJS] IDs appear to be UniProt format: ", is.uniprot)

          entrez.ids <- NULL

          if (is.uniprot) {
            # Convert UniProt -> Entrez using database
            msg("[CorrIgraph2SigmaJS] Querying UniProt->Entrez mapping from database...")
            uniprot.map <- tryCatch({
              queryGeneDB("entrez_uniprot", org)
            }, error = function(e) {
              msg("[CorrIgraph2SigmaJS] queryGeneDB failed: ", e$message)
              NULL
            })

            if (!is.null(uniprot.map) && is.data.frame(uniprot.map)) {
              msg("[CorrIgraph2SigmaJS] UniProt map columns: ", paste(colnames(uniprot.map), collapse=", "))
              if ("accession" %in% colnames(uniprot.map) && "gene_id" %in% colnames(uniprot.map)) {
                # Match UniProt accessions to get Entrez IDs
                hit.inx <- match(node.ids, as.character(uniprot.map$accession))
                entrez.ids <- uniprot.map$gene_id[hit.inx]
                msg("[CorrIgraph2SigmaJS] Matched ", sum(!is.na(entrez.ids)), " UniProt accessions to Entrez")
              }
            }
          } else {
            # Try symbol -> entrez mapping
            msg("[CorrIgraph2SigmaJS] Trying symbol->Entrez mapping...")
            gene.map <- tryCatch({
              queryGeneDB("symbol", org)
            }, error = function(e) {
              msg("[CorrIgraph2SigmaJS] queryGeneDB failed: ", e$message)
              NULL
            })

            if (!is.null(gene.map) && is.data.frame(gene.map)) {
              msg("[CorrIgraph2SigmaJS] Symbol map columns: ", paste(colnames(gene.map), collapse=", "))
              # Try to match by symbol
              if ("symbol" %in% colnames(gene.map) && "gene_id" %in% colnames(gene.map)) {
                hit.inx <- match(toupper(node.ids), toupper(gene.map$symbol))
                entrez.ids <- gene.map$gene_id[hit.inx]
                msg("[CorrIgraph2SigmaJS] Matched ", sum(!is.na(entrez.ids)), " symbols to Entrez")
              }
            }
          }

          # Check results
          if (!is.null(entrez.ids)) {
            valid.inx <- !is.na(entrez.ids)
            if (sum(valid.inx) == 0) {
              msg("[CorrIgraph2SigmaJS] WARNING: Could not convert any node IDs to Entrez, PPI projection will not work")
              node.ids <- as.character(node.ids)
            } else {
              msg("[CorrIgraph2SigmaJS] Successfully converted ", sum(valid.inx), "/", length(node.ids), " node IDs to Entrez")
              # Keep only valid Entrez IDs for PPI query
              node.ids <- as.character(entrez.ids[valid.inx])
            }
          } else {
            msg("[CorrIgraph2SigmaJS] WARNING: Could not convert any node IDs to Entrez, PPI projection will not work")
            node.ids <- as.character(node.ids)
          }
        } else {
          node.ids <- as.character(node.ids)
        }

        # Query PPI database with Entrez IDs
        # Parameters: sqlite.path, table.nm (db name), q.vec (node IDs), requireExp, min.score
        requireExp <- FALSE  # Don't require experimental evidence
        min.score <- 0       # No minimum score filter (accept all interactions)

        ppi.edges <- QueryPpiSQLite(sqlite.path, table.nm, node.ids, requireExp, min.score)

        if (!is.null(ppi.edges) && nrow(ppi.edges) > 0) {
          msg("[CorrIgraph2SigmaJS] Found ", nrow(ppi.edges), " PPI interactions from ", table.nm)

          # Show sample PPI edges
          sample.size <- min(10, nrow(ppi.edges))
          msg("[CorrIgraph2SigmaJS] Sample PPI edges (first ", sample.size, "):")
          for (j in 1:sample.size) {
            msg("  ", ppi.edges$id1[j], " <-> ", ppi.edges$id2[j])
          }

          # Show PPI edge ID types
          msg("[CorrIgraph2SigmaJS] PPI database ID1 sample: ", paste(head(ppi.edges$id1, 5), collapse=", "))
          msg("[CorrIgraph2SigmaJS] PPI database ID2 sample: ", paste(head(ppi.edges$id2, 5), collapse=", "))

          # Create mapping from original IDs to Entrez IDs for edge matching
          # Store as named vector: original_id -> entrez_id
          if (!is.entrez && exists("entrez.ids") && !is.null(entrez.ids) && length(entrez.ids) > 0) {
            # Build full mapping including NAs
            original.ids <- V(g)$name
            id.to.entrez <- setNames(entrez.ids, original.ids)

            # Log mapping info
            valid.map.count <- sum(!is.na(id.to.entrez))
            msg("[CorrIgraph2SigmaJS] ID mapping: ", valid.map.count, " valid mappings out of ", length(id.to.entrez))
            msg("[CorrIgraph2SigmaJS] Sample ID mapping (first 5):")
            for (j in 1:min(5, length(id.to.entrez))) {
              msg("  ", names(id.to.entrez)[j], " -> ", if(!is.na(id.to.entrez[j])) id.to.entrez[j] else "NA")
            }

            ppi.edges$id.mapping <- list(id.to.entrez)  # Store for later use (includes NAs)
          } else {
            # IDs are already Entrez, no mapping needed
            msg("[CorrIgraph2SigmaJS] IDs are already Entrez, no mapping needed")
            msg("[CorrIgraph2SigmaJS] Sample node IDs for PPI matching: ", paste(head(node.ids, 5), collapse=", "))
            ppi.edges$id.mapping <- list(NULL)
          }
        } else {
          msg("[CorrIgraph2SigmaJS] No PPI interactions found in ", table.nm)
          ppi.edges <- NULL
        }
      }
    }, error = function(e) {
      msg("[CorrIgraph2SigmaJS] ERROR querying PPI database: ", e$message)
      ppi.edges <- NULL
    })
  }

  # Create mapping from original vertex names to primary IDs (UniProt preferred)
  vertex.name.to.primary <- setNames(
    sapply(seq_along(uniprot.vec), function(i) {
      if (!is.na(uniprot.vec[i]) && !is.null(uniprot.vec[i]) && nzchar(uniprot.vec[i])) {
        as.character(uniprot.vec[i])
      } else {
        as.character(nms[i])
      }
    }),
    nms
  )

  # Track PPI overlap statistics
  ppi.overlap.count <- 0
  total.edges <- nrow(el)

  # Log edge matching details for first few edges (debugging)
  debug.edge.count <- min(5, total.edges)

  edges <- lapply(seq_len(nrow(el)), function(i) {
    w <- as.numeric(el$weight[i])
    src.orig <- as.character(el$from[i])
    tgt.orig <- as.character(el$to[i])

    # Map to primary IDs (UniProt if available)
    src.primary <- vertex.name.to.primary[src.orig]
    tgt.primary <- vertex.name.to.primary[tgt.orig]

    # Check if this edge is in PPI database
    has.ppi <- FALSE
    if (!is.null(ppi.edges) && nrow(ppi.edges) > 0) {
      # PPI database uses id1 and id2 columns (Entrez IDs)
      # Co-expression edge uses src/tgt from original data

      # Get Entrez IDs for src and tgt (use original vertex names for lookup)
      if (!is.null(ppi.edges$id.mapping[[1]])) {
        # Use mapping: original_id -> entrez_id
        id.map <- ppi.edges$id.mapping[[1]]
        src.entrez <- id.map[src.orig]
        tgt.entrez <- id.map[tgt.orig]

        # Debug logging for first few edges
        if (i <= debug.edge.count) {
          msg("[CorrIgraph2SigmaJS] Edge ", i, ": ", src.orig, " -> ", tgt.orig)
          msg("  src.entrez=", if(!is.na(src.entrez)) src.entrez else "NA",
              ", tgt.entrez=", if(!is.na(tgt.entrez)) tgt.entrez else "NA")
        }

        # Skip if either ID couldn't be mapped
        if (!is.na(src.entrez) && !is.na(tgt.entrez)) {
          # Check both directions (undirected network)
          has.ppi <- any((ppi.edges$id1 == src.entrez & ppi.edges$id2 == tgt.entrez) |
                         (ppi.edges$id1 == tgt.entrez & ppi.edges$id2 == src.entrez))

          if (i <= debug.edge.count) {
            msg("  Checking PPI: ", has.ppi)
            if (has.ppi) {
              msg("  MATCHED in PPI database!")
            }
          }
        } else {
          if (i <= debug.edge.count) {
            msg("  Skipped (unmapped Entrez IDs)")
          }
        }
      } else {
        # IDs are already Entrez, direct comparison
        has.ppi <- any((ppi.edges$id1 == src.orig & ppi.edges$id2 == tgt.orig) |
                       (ppi.edges$id1 == tgt.orig & ppi.edges$id2 == src.orig))

        if (i <= debug.edge.count) {
          msg("[CorrIgraph2SigmaJS] Edge ", i, " (Entrez): ", src.orig, " -> ", tgt.orig, " PPI=", has.ppi)
        }
      }
    }

    # Count PPI overlaps
    if (has.ppi) {
      ppi.overlap.count <<- ppi.overlap.count + 1
    }

    # Set edge color and size based on PPI support
    if (has.ppi) {
      edge.color <- "#ff7f00"  # Orange for PPI-supported
      edge.size <- rescale(w) * 1.8  # 1.8x thicker
    } else {
      edge.color <- "#cccccc"  # Light gray for co-expression only
      edge.size <- rescale(w)
    }

    list(
      id     = paste0("e", i),
      source = src.primary,           # Use primary ID (UniProt preferred)
      target = tgt.primary,           # Use primary ID (UniProt preferred)
      weight = w,                     # keeps the raw weight
      size   = edge.size,             # stroke width (adjusted for PPI)
      color  = edge.color,            # edge color (based on PPI support)
      ppi    = has.ppi,               # boolean flag
      ppiDb  = if(has.ppi) ppiDatabase else ""  # which database
    )
  })

  # Log final PPI overlap statistics
  msg("[CorrIgraph2SigmaJS] ===== PPI OVERLAP SUMMARY =====")
  msg("[CorrIgraph2SigmaJS] Total co-expression edges: ", total.edges)
  msg("[CorrIgraph2SigmaJS] Edges with PPI support: ", ppi.overlap.count)
  msg("[CorrIgraph2SigmaJS] Overlap percentage: ",
      round(100 * ppi.overlap.count / total.edges, 2), "%")
  msg("[CorrIgraph2SigmaJS] ===============================")
  
  
  dataSet <- readDataset(paramSet$dataName);
  nodeTable <- BuildNodeTable(g,paramSet,dataSet,analSet);

  # Detect ID type and map to symbols for node labels
  node.ids <- V(g)$name
  is.entrez.ids <- all(grepl("^[0-9]+$", node.ids))
  is.uniprot.ids <- any(grepl("^[A-Z][0-9][A-Z0-9]{3,4}[0-9]$", node.ids))

  initsbls <- tryCatch({
    if (is.uniprot.ids) {
      doUniprot2SymbolMapping(node.ids, paramSet$data.org, paramSet$data.idType)
    } else {
      doEntrez2SymbolMapping(node.ids, paramSet$data.org, paramSet$data.idType)
    }
  }, error = function(e) {
    msg("[CorrIgraph2SigmaJS] WARNING: Symbol mapping failed: ", e$message)
    node.ids  # Use IDs as labels if mapping fails
  })

  # Ensure same length
  if (is.null(initsbls) || length(initsbls) != length(node.ids)) {
    initsbls <- node.ids
  }

  names(initsbls) <- V(g)$name
  ppi.net <- list();

  ppi.net[["node.data"]] <- data.frame(Id=V(g)$name, Label=unname(initsbls));
  ppi.net <<- ppi.net;
  ## ── 3 · assemble JSON payload -----------------------------------
  netData <- list(nodes            = nodes,
                  edges            = edges,
                  backgroundColor  = list("#f5f5f5", "#0066CC"),
                  naviString       = "Correlation Network",
                  org              = paramSet$data.org,
                  proteinlist=initsbls, 
                  nodeTable = nodeTable);
  
  fileNm <- paste0(netNm, ".json")
  jsonlite::write_json(netData, fileNm, auto_unbox = TRUE)
  
  ## track for later download
  paramSet$partialToBeSaved <- c(paramSet$partialToBeSaved, fileNm)
  paramSet$jsonNms$coexpNet <- basename(fileNm)
  saveSet(paramSet, "paramSet")
  analSet$corNet <- netData
  saveSet(analSet, "analSet")
  
  invisible(netData)
}

SplitIgraphByModule <- function(g, keepXTalk = FALSE) {

  stopifnot("module" %in% vertex_attr_names(g))

  mods <- sort(unique(V(g)$module))
  cat(sprintf("[SplitIgraphByModule] Found %d unique modules: %s\n", length(mods), paste(mods, collapse=", ")))

  subG <- setNames(vector("list", length(mods)), mods)
  
  for (m in mods) {
    features.m <- V(g)[module == m]
    cat(sprintf("[SplitIgraphByModule] Module '%s' has %d vertices\n", m, length(features.m)))

    if (keepXTalk) {
      # keep all edges touching those features
      subG[[m]] <- induced_subgraph(g, vids = features.m)
    } else {
      # keep *only* edges whose BOTH endpoints are in module m
      eKeep <- E(g)[inc(V(g)[module == m])]
      eKeep <- eKeep[ which(ends(g, eKeep)[,1] %in% features.m$name &
                              ends(g, eKeep)[,2] %in% features.m$name) ]
      subG[[m]] <- subgraph_from_edges(g, eKeep)
    }

    cat(sprintf("[SplitIgraphByModule] Subgraph for module '%s': %d nodes, %d edges\n",
                m, vcount(subG[[m]]), ecount(subG[[m]])))
  }

  cat(sprintf("[SplitIgraphByModule] Returning %d subgraphs\n", length(subG)))
  subG
}

GenerateCEMModuleNetworks <- function(fileName  = "coexp_network",
                                      thresh    = 0.05,
                                      keepXTalk = FALSE,
                                      minNodeNum = 3) {

  paramSet <- readSet(paramSet, "paramSet")
  analSet  <- readSet(analSet,  "analSet")
  library(igraph);

  print(names(analSet));
  g.all   <- analSet$overall.graph.orig

  # Check if graph exists
  if (is.null(g.all)) {
    AddErrMsg("overall.graph.orig is NULL. BuildIgraphFromCEM may have failed.");
    return(0);
  }

  g.byMod <- SplitIgraphByModule(g.all, keepXTalk = keepXTalk)

  comps <- g.byMod

  # Check if we have any modules
  if(length(comps) == 0) {
    AddErrMsg("No modules found in the CEMiTool network. CEMiTool may not have identified any co-expression modules.");
    return(0);
  }

  # first compute subnet stats
  net.stats <- ComputeSubnetStats(comps);
  ord.inx <- order(net.stats[,1], decreasing=TRUE);
  net.stats <- net.stats[ord.inx,];
  comps <- comps[ord.inx];
  names(comps) <- rownames(net.stats) <- c(paste("module", 1:length(comps), sep=""));

  # note, we report stats for all nets (at least 3 nodes);
  hit.inx <- net.stats$Node >= minNodeNum;
  comps <- comps[hit.inx];

  # Check if filtering removed all modules
  if(length(comps) == 0) {
    AddErrMsg(sprintf("All modules were filtered out (minNodeNum=%d). Try lowering the minimum node threshold.", minNodeNum));
    return(0);
  }


  # now record
  net.stats <<- net.stats;
  sub.stats <- unlist(lapply(comps, vcount));

  
  ## ── 3 · write JSON for *only the first* module in the list ─────
  firstMod <- names(g.byMod)[1]              # e.g. "M1"
  netNm <- fileName;
  analSet$ppi.comps <- comps;
  saveSet(analSet, "analSet")
  return(c(vcount(g.all), ecount(g.all), length(comps), sub.stats));
}


BuildNodeTable <- function(g,
                           paramSet,
                           dataSet = NULL,
                           analSet = NULL) {

  ids    <- V(g)$name

  # Detect ID type and apply appropriate mapping
  # Co-expression networks can use Entrez IDs or UniProt accessions
  is.entrez <- all(grepl("^[0-9]+$", ids))
  is.uniprot <- any(grepl("^[A-Z][0-9][A-Z0-9]{3,4}[0-9]$", ids))

  labels <- tryCatch({
    if (is.uniprot) {
      msg("[BuildNodeTable] Detected UniProt IDs, using doUniprot2SymbolMapping")
      doUniprot2SymbolMapping(ids, paramSet$data.org, paramSet$data.idType)
    } else {
      msg("[BuildNodeTable] Detected Entrez IDs, using doEntrez2SymbolMapping")
      doEntrez2SymbolMapping(ids, paramSet$data.org, paramSet$data.idType)
    }
  }, error = function(e) {
    msg("[BuildNodeTable] WARNING: ID mapping failed: ", e$message)
    # Return IDs as labels if mapping fails
    return(ids)
  })

  # Ensure labels has same length as ids
  if (is.null(labels) || length(labels) == 0) {
    msg("[BuildNodeTable] WARNING: ID mapping returned empty vector, using IDs as labels")
    labels <- ids
  } else if (length(labels) != length(ids)) {
    msg("[BuildNodeTable] WARNING: Label length mismatch (", length(labels), " vs ", length(ids), "), using IDs as labels")
    labels <- ids
  }

  anal.type <- paramSet$anal.type
  ## --- topological features --------------------------------------
  deg  <- igraph::degree(g)

  # PERFORMANCE FIX (Issue #8): Adaptive betweenness calculation
  # For large networks (>1000 nodes), use cutoff to limit path length
  # Reduces complexity from O(V³) to approximately O(V²) for dense graphs
  n_nodes <- vcount(g)
  if (n_nodes < 1000) {
    # Small networks: use exact betweenness (fast enough)
    btw <- igraph::betweenness(g)
  } else {
    # Large networks: use path length cutoff for approximation
    # cutoff=3 considers paths up to length 3 (captures local centrality)
    # 10-100x faster for networks with 1000+ nodes
    btw <- igraph::betweenness(g, cutoff = 3)
  }

  ## --- expression values (log2FC) --------------------------------
  expr <- rep(0, length(ids))         # default NA

  if (anal.type == "onedata") {

    tbl <- dataSet$comp.res
    inx <- match(ids, rownames(tbl))
    expr <- tbl[inx, paramSet$selectedFactorInx]

  } else if (anal.type == "metadata") {

    if (paramSet$selDataNm == "meta_default") {
      tbl  <- analSet$meta.mat.all
      sy   <- doEntrez2SymbolMapping(rownames(tbl),
                                     paramSet$data.org,
                                     paramSet$data.idType)
      inx  <- match(ids, sy)
      expr <- analSet$meta.avgFC[rownames(tbl)][inx]

    } else {                     # metadata but user-selected dataset
      ds   <- readDataset(paramSet$selDataNm)
      tbl  <- ds$comp.res
      sy   <- doEntrez2SymbolMapping(rownames(tbl),
                                     paramSet$data.org,
                                     paramSet$data.idType)
      inx  <- match(ids, sy)
      expr <- tbl[inx, "logFC"]
    }
  }
  ## for "proteinlist" expr stays NA  (per requirement)

  ## --- load uniprot mapping --------------------------------------
  uniprot.vec <- rep("", length(ids))
  if (file.exists("symbol.map.qs")) {
    tryCatch({
      gene.map <- readDataQs("symbol.map.qs", paramSet$anal.type, paramSet$dataName)
      if (!is.null(gene.map) && is.data.frame(gene.map) &&
          "gene_id" %in% colnames(gene.map) && "uniprot" %in% colnames(gene.map)) {
        # Match node IDs (Entrez) to gene.map
        hit.inx <- match(as.character(ids), as.character(gene.map$gene_id))
        uniprot.vec.temp <- gene.map$uniprot[hit.inx]
        # Convert NA to empty string
        uniprot.vec <- ifelse(is.na(uniprot.vec.temp), "", as.character(uniprot.vec.temp))
      }
    }, error = function(e) {
      # Silently ignore errors
    })
  }

  ## --- assemble data-frame ---------------------------------------
  node.df <- data.frame(
    id          = ids,
    label       = ifelse(is.na(labels) | labels == "", ids, labels),
    uniprot     = uniprot.vec,
    degree      = deg,
    betweenness = btw,
    expr        = expr,
    stringsAsFactors = FALSE
  )

  node.df <- node.df[order(node.df$degree, decreasing = TRUE), ]

  node.df
}



ComputeSubnetStats <- function(comps){
  library(igraph);

  # Check if comps is empty or has issues
  if(length(comps) == 0) {
    warning("ComputeSubnetStats: comps is empty");
    return(data.frame(Node = integer(0), Edge = integer(0), Query = integer(0)));
  }

  net.stats <- as.data.frame(matrix(0, ncol = 3, nrow = length(comps)));
  colnames(net.stats) <- c("Node", "Edge", "Query");

  for(i in 1:length(comps)){
    # Check if this component exists and is valid
    if(is.null(comps[[i]])) {
      warning(sprintf("ComputeSubnetStats: comps[[%d]] is NULL", i));
      net.stats[i,] <- c(0, 0, 0);
      next;
    }

    g <- comps[[i]];
    num_nodes <- vcount(g);
    num_edges <- ecount(g);
    # Count query nodes (from original upload)
    num_query <- if (!is.null(V(g)$is_query)) {
      sum(V(g)$is_query, na.rm = TRUE)
    } else {
      0  # If is_query attribute doesn't exist, default to 0
    }
    net.stats[i,] <- c(num_nodes, num_edges, num_query);
  }
  return(net.stats);
}
# ====================================================================
# filterNetByThresh  —  edge-weight filtering for an igraph object
# --------------------------------------------------------------------
# g          : igraph object that already has a numeric edge attribute
#              called 'weight'
# thresh     : keep edges with weight > thresh
# maxEdges   : cap the network at this many heaviest edges (NULL = no cap)
# rmIsolated : TRUE → delete nodes that become isolated after filtering
# layoutFun  : (optional) layout recalculation if you need updated coords
# return     : list(graph = <filtered igraph>,
#                   stats = c(nodes, edges, n.components))
# ====================================================================
FilterNetByThresh <- function(thresh      = 0.05,
                                 maxEdges    = 2000,
                                 rmIsolated  = TRUE) {
  # save.image("filter.RDAta");
  analSet  <- readSet(analSet,  "analSet")
  overall.graph <- analSet$overall.graph;
  g <- overall.graph;
  if (!"weight" %in% edge_attr_names(g))
    stop("edge attribute 'weight' not found")

  # ── 1 · keep only edges above threshold ───────────────────────────
  g <- subgraph_from_edges(g, E(g)[weight > thresh], delete.vertices = FALSE)

  # ── 2 · cap total edges if requested ──────────────────────────────
  if (!is.null(maxEdges) && ecount(g) > maxEdges) {
    el <- igraph::as_data_frame(g, what = "edges")
    el <- el[order(el$weight, decreasing = TRUE), ][seq_len(maxEdges), ]
    g  <- graph_from_data_frame(el,
                                directed = FALSE,
                                vertices = igraph::as_data_frame(g, what = "vertices"))
    E(g)$weight <- el$weight                    # restore weights
  }

  # ── 3 · optionally drop newly-isolated nodes ──────────────────────
  if (rmIsolated) {
    iso <- which(igraph::degree(g) == 0)
    if (length(iso) > 0)
      g <- delete_vertices(g, iso)
  }

  # ── 4 · final tidy-up (loops / multi-edges) ───────────────────────
  g <- simplify(g, edge.attr.comb = list("first"))

  # ── 6 · book-keeping & return  ────────────────────────────────────
  current.msg <<-
    paste("FilterNetByThreshold:",
          vcount(g), "nodes and", ecount(g), "edges retained at thresh >", thresh);
  

    analSet <- DecomposeGraph(overall.graph,analSet);
    substats <- analSet$substats;
  outStats <- c(vcount(g), ecount(g), length(substats))
  analSet$overall.graph <- g;
    return(saveSet(analSet, "analSet", outStats));
}


FilterBipartiNet <- function(nd.type, min.dgr, min.btw){
    paramSet <- readSet(paramSet, "paramSet");
    analSet <- readSet(analSet, "analSet");
    overall.graph <- analSet$overall.graph
    all.nms <- V(overall.graph)$name;
    edge.mat <- as_edgelist(overall.graph);
    dgrs <- degree(overall.graph);
    nodes2rm.dgr <- nodes2rm.btw <- NULL;

    if(nd.type == "gene"){
        hit.inx <- all.nms %in% edge.mat[,1];
    }else if(nd.type=="other"){
        hit.inx <- all.nms %in% edge.mat[,2];
    }else{ # all
        hit.inx <- rep(TRUE, length(all.nms));
    }

    if(min.dgr > 0){
        rm.inx <- dgrs <= min.dgr & hit.inx;
        nodes2rm.dgr <- V(overall.graph)$name[rm.inx];
    }
    if(min.btw > 0){
        # PERFORMANCE FIX (Issue #8): Adaptive betweenness calculation
        # Use cutoff for large networks to reduce O(V³) complexity
        n_nodes <- vcount(overall.graph)
        if (n_nodes < 1000) {
            btws <- betweenness(overall.graph)
        } else {
            btws <- betweenness(overall.graph, cutoff = 3)
        }
        rm.inx <- btws <= min.btw & hit.inx;
        nodes2rm.btw <- V(overall.graph)$name[rm.inx];
    }

    nodes2rm <- unique(c(nodes2rm.dgr, nodes2rm.btw));
    overall.graph <- simplify(delete_vertices(overall.graph, nodes2rm));
    current.msg <<- paste("A total of", length(nodes2rm) , "was reduced.");
    analSet <- DecomposeGraph(overall.graph,analSet);
    substats <- analSet$substats;
    if(!is.null(substats)){
        output <- c(vcount(overall.graph), ecount(overall.graph), length(analSet$ppi.comps), substats);
    }else{
        output <- 0;
    }
    analSet$overall.graph <- overall.graph;

    return(saveSet(analSet, "analSet", output));
}

PrepareCoexpNetwork <- function(net.nm, jsonNm, projectPPI = FALSE, ppiDatabase = ""){
   analSet <- readSet(analSet, "analSet");
   paramSet <- readSet(paramSet, "paramSet");
   #print(analSet$ppi.comps);

   my.ppi <- analSet$ppi.comps[[net.nm]];
   nd.nms <- V(my.ppi)$name;
  CorrIgraph2SigmaJS(my.ppi,
                     netNm       = jsonNm,
                     paramSet    = paramSet,
                     analSet     = analSet,
                     projectPPI  = projectPPI,
                     ppiDatabase = ppiDatabase)

   current.net.nm <<- net.nm;
   paramSet$current.net.nm <- net.nm;
   saveSet(paramSet, "paramSet");

   return(saveSet(analSet, "analSet", 1));

}
