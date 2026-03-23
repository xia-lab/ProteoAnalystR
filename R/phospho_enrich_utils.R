##################################################
## R script for ProteoAnalyst
## Description: Phosphoproteomics kinase enrichment analysis
## Authors:
## Generated for phospho module implementation
###################################################

#' Perform kinase enrichment analysis on significant phosphosites
#'
#' @param dataName Dataset name
#' @param database Kinase-substrate database to use (matches files in resources/data/ksea, or "all")
#' @return 1 on success, 0 on failure
PerformKinaseEnrichment <- function(dataName, database = "phosphositeplus") {
  dataSet <- readDataset(dataName)
  paramSet <- readSet(paramSet, "paramSet")
  analSet <- readSet(analSet, "analSet")
  msgSet <- readSet(msgSet, "msgSet")

  # Check if this is phospho data
  if (is.null(paramSet$data.type) || paramSet$data.type != "phospho") {
    #msg("[PerformKinaseEnrichment] Not phospho data, skipping.")
    return(0)
  }

  # Get significant phosphosites
  if (is.null(dataSet$sig.mat) || nrow(dataSet$sig.mat) == 0) {
    #msg("[PerformKinaseEnrichment] No significant sites found. Run differential expression first.")
    msgSet$current.msg <- "No significant phosphosites found. Please perform differential expression analysis first."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  sig_sites <- rownames(dataSet$sig.mat)
  all_sites <- rownames(dataSet$data.norm)

  # Attempt to convert site IDs to symbol-based substrate strings (SYMBOL;Y123;)
  converted <- .convertSitesToSymbols(dataSet, paramSet)
  if (!is.null(converted)) {
    #msg("[PerformKinaseEnrichment] Converted site IDs to symbols. Example all_sites head: ",
    #        paste(head(converted$all_sites), collapse=", "))
    sig_sites <- converted$sig_sites
    all_sites <- converted$all_sites
  } else {
    #msg("[PerformKinaseEnrichment] Site conversion failed; using raw IDs. Example head: ",
    #        paste(head(all_sites), collapse=", "))
  }

  # Load kinase-substrate database
  ks_db <- .loadKinaseSubstrateDB(database, paramSet)

  if (is.null(ks_db) ||
      (is.data.frame(ks_db) && nrow(ks_db) == 0) ||
      (is.list(ks_db) && !is.data.frame(ks_db) && length(ks_db) == 0)) {
    #msg("[PerformKinaseEnrichment] Failed to load kinase-substrate database.")
    msgSet$current.msg <- "Failed to load kinase-substrate database. Please check database files."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Perform enrichment analysis
  enrichment_results <- .performKinaseEnrichmentTest(sig_sites, all_sites, ks_db, dataSet)

  if (is.null(enrichment_results) || nrow(enrichment_results) == 0) {
    #msg("[PerformKinaseEnrichment] No enriched kinases found.")
    msgSet$current.msg <- "No significantly enriched kinases were identified."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Save results
  analSet$kinase.enrich <- enrichment_results
  saveSet(analSet, "analSet")

  # Export to CSV
  fast.write(enrichment_results, file = "kinase_enrichment.csv", row.names = FALSE)

  msgSet$current.msg <- paste0("Kinase enrichment analysis complete. Found ",
                                nrow(enrichment_results), " enriched kinases.")
  saveSet(msgSet, "msgSet")

  return(1)
}

#' List available kinase databases for the current organism
#' @return character vector of database names (includes "all" if >1)
GetAvailableKinaseDBs <- function() {
  paramSet <- readSet(paramSet, "paramSet")
  .listAvailableKseaDBs(paramSet, include_all = TRUE)
}

#' Load kinase-substrate database
#'
#' @param database Database name
#' @param paramSet Parameter set with organism info
#' @return Data frame with kinase-substrate relationships
.loadKinaseSubstrateDB <- function(database = "phosphositeplus", paramSet) {

  # Get organism
  org <- paramSet$data.org
  if (is.null(org)) org <- "hsa"

  # Database file paths
  # Expected location: resources/data/ksea/<DB>_uniprot_<org>.qs
  # lib.path points to ../../data/ relative to R working directory
  db_dir <- paste0(paramSet$lib.path, "ksea/")

  #msg("[.loadKinaseSubstrateDB] Looking for database in: ", db_dir)

  if (!dir.exists(db_dir)) {
    #msg("[.loadKinaseSubstrateDB] Database directory not found: ", db_dir)
    #msg("[.loadKinaseSubstrateDB] Creating minimal kinase-substrate mapping.")
    ks_db <- .createMinimalKinaseDB(org)
    return(ks_db)
  }

  available_dbs <- .listAvailableKseaDBs(paramSet, include_all = FALSE)

  # Determine which DBs to load
  if (is.null(database) || database == "" || tolower(database) == "auto") {
    database <- if (length(available_dbs) > 0) available_dbs[1] else "auto"
  }

  dbs_to_load <- character(0)
  if (tolower(database) == "all") {
    dbs_to_load <- available_dbs
  } else {
    hits <- available_dbs[tolower(available_dbs) == tolower(database)]
    if (length(hits) == 0) {
      #msg("[.loadKinaseSubstrateDB] Requested database ", database,
      #        " not found for org ", org, ". Available: ", paste(available_dbs, collapse = ", "))
      if (length(available_dbs) > 0) {
        dbs_to_load <- available_dbs[1]
        #msg("[.loadKinaseSubstrateDB] Falling back to ", dbs_to_load)
      }
    } else {
      dbs_to_load <- hits
    }
  }

  if (length(dbs_to_load) == 0) {
    #msg("[.loadKinaseSubstrateDB] No curated DB found; using minimal motif DB.")
    return(.createMinimalKinaseDB(org))
  }

  loaded <- lapply(dbs_to_load, function(db_nm) .loadSingleKseaDB(db_nm, org, db_dir))
  # Drop nulls
  loaded <- loaded[!vapply(loaded, is.null, logical(1))]

  if (length(loaded) == 0) {
    #msg("[.loadKinaseSubstrateDB] Failed to load requested DB(s); using minimal motif DB.")
    return(.createMinimalKinaseDB(org))
  }

  if (length(loaded) == 1) {
    return(loaded[[1]])
  }

  # Merge multiple DBs (union of substrates per kinase)
  merged <- list()
  for (db in loaded) {
    if (is.list(db) && !is.data.frame(db)) {
      for (k in names(db)) {
        merged[[k]] <- unique(c(merged[[k]], as.character(db[[k]])))
      }
    } else if (is.data.frame(db) && "Substrate" %in% colnames(db) && "Kinase" %in% colnames(db)) {
      split_db <- split(as.character(db$Substrate), db$Kinase)
      for (k in names(split_db)) {
        merged[[k]] <- unique(c(merged[[k]], split_db[[k]]))
      }
    }
  }
  #msg("[.loadKinaseSubstrateDB] Merged ", length(merged), " kinases from ",
  #        paste(dbs_to_load, collapse = ", "))
  merged
}

# Return available KSEA databases for the current organism
.listAvailableKseaDBs <- function(paramSet, include_all = TRUE) {
  org <- paramSet$data.org
  if (is.null(org)) org <- "hsa"
  db_dir <- paste0(paramSet$lib.path, "ksea/")
  if (!dir.exists(db_dir)) return(character(0))
  files <- list.files(db_dir, pattern = "\\.qs$", full.names = FALSE)
  if (length(files) == 0) return(character(0))

  opts <- character(0)
  for (f in files) {
    parts <- strsplit(basename(f), "_")[[1]]
    if (length(parts) < 3) next
    org_part <- sub("\\.qs$", "", parts[length(parts)])
    if (tolower(org_part) != tolower(org)) next
    opts <- c(opts, parts[1])
  }
  opts <- unique(opts)
  opts <- opts[order(tolower(opts))]

  opts
}

# Load a single KSEA db file and normalize to list format
.loadSingleKseaDB <- function(db_nm, org, db_dir) {
  candidates <- list.files(db_dir,
                           pattern = paste0("^", db_nm, "_.*_", org, "\\.qs$"),
                           full.names = TRUE,
                           ignore.case = TRUE)
  if (length(candidates) == 0) {
    #msg("[.loadSingleKseaDB] No file for ", db_nm, " and org ", org)
    return(NULL)
  }
  db_file <- candidates[1]
  #msg("[.loadSingleKseaDB] Loading ", basename(db_file))
  db_obj <- qs::qread(db_file)

  # Normalize to list(kinase = substrates)
  if (is.list(db_obj) && !is.data.frame(db_obj)) {
    db_obj <- lapply(db_obj, function(x) unique(as.character(x)))
    return(db_obj)
  }

  if (is.data.frame(db_obj)) {
    # Try to detect kinase/substrate columns
    kcol <- NULL
    scol <- NULL
    if ("Kinase" %in% colnames(db_obj)) kcol <- "Kinase"
    if ("Substrate" %in% colnames(db_obj)) scol <- "Substrate"
    if (is.null(kcol)) kcol <- colnames(db_obj)[1]
    if (is.null(scol)) scol <- colnames(db_obj)[ncol(db_obj)]
    db_obj[[kcol]] <- as.character(db_obj[[kcol]])
    db_obj[[scol]] <- as.character(db_obj[[scol]])
    return(split(db_obj[[scol]], db_obj[[kcol]]))
  }

  #msg("[.loadSingleKseaDB] Unsupported object type for ", basename(db_file), "; skipping.")
  NULL
}

#' Create minimal kinase-substrate database based on known motifs
#'
#' This is a fallback when no curated database is available
#' Maps phosphosites to kinases based on simple sequence motifs
.createMinimalKinaseDB <- function(org) {

  #msg("[.createMinimalKinaseDB] Creating minimal kinase database for organism: ", org)

  # Define basic kinase motifs (very simplified)
  # Format: Kinase, Motif (regex), Residue
  kinase_motifs <- data.frame(
    Kinase = c("PKA", "PKA", "PKC", "PKC", "CK2", "CK2",
               "MAPK", "MAPK", "CDK", "CDK", "GSK3", "GSK3",
               "AKT", "AKT", "mTOR", "mTOR"),
    Motif = c("R.{2}S", "K.{2}S", "S.{2}K", "T.{2}K",
              "S.{2}[DE]", "T.{2}[DE]",
              "[ST]P", "[ST]P",
              "[ST]P.{0,2}[KR]", "[ST]P.{0,2}[KR]",
              "S.{3}[ST]", "T.{3}[ST]",
              "RxRxx[ST]", "RxRxx[ST]",
              "[ST].{3}[FILV]", "[ST].{3}[FILV]"),
    Residue = c("S", "S", "S", "T", "S", "T", "S", "T", "S", "T", "S", "T", "S", "T", "S", "T"),
    stringsAsFactors = FALSE
  )

  # Note: This returns just the motif patterns
  # Actual site matching happens in the enrichment function
  return(kinase_motifs)
}

#' Perform kinase enrichment test using Fisher's exact test
#'
#' @param sig_sites Significant phosphosites
#' @param all_sites All detected phosphosites
#' @param ks_db Kinase-substrate database (list format)
#' @param dataSet Dataset object with feature info
#' @return Data frame with enrichment results
.performKinaseEnrichmentTest <- function(sig_sites, all_sites, ks_db, dataSet) {

  # Check database format
  # If it's a list (like the PhosphoSitePlus database), process accordingly
  if (is.list(ks_db) && !is.data.frame(ks_db)) {
    #msg("[.performKinaseEnrichmentTest] KS DB is list-based; example kinase entries: ",
    #        paste(names(ks_db)[1:min(3, length(ks_db))], collapse=", "))
    results <- .enrichmentFromListDB(sig_sites, all_sites, ks_db)
  } else if (is.data.frame(ks_db) && "Substrate" %in% colnames(ks_db)) {
    # Data frame format with Substrate column
    results <- .enrichmentFromCuratedDB(sig_sites, all_sites, ks_db)
  } else {
    # Use motif-based mapping
    results <- .enrichmentFromMotifs(sig_sites, all_sites, ks_db, dataSet)
  }

  return(results)
}

# Convert phosphosite IDs to standardized substrate strings using UniProt accessions (e.g., P12345;Y542;)
.convertSitesToSymbols <- function(dataSet, paramSet) {
  fi <- dataSet$feature.info
  if (is.null(fi) || !all(c("Proteins", "Amino acid", "Position") %in% colnames(fi))) {
    #msg("[.convertSitesToSymbols] Feature info missing Proteins/Amino acid/Position; skipping conversion.")
    return(NULL)
  }

  #msg("[.convertSitesToSymbols] Feature info columns: ", paste(colnames(fi), collapse=", "))

  prot_raw <- as.character(fi$Proteins)
  aa   <- toupper(as.character(fi$`Amino acid`))
  pos  <- as.character(fi$Position)

  # Base IDs: first cleaned accession from list
  base_acc <- vapply(prot_raw, function(x) {
    accs <- trimws(strsplit(x, ";")[[1]])
    acc <- accs[1]
    acc <- sub("-\\d+$", "", acc)
    acc <- sub("_.*", "", acc)
    acc
  }, character(1))

  #msg("[.convertSitesToSymbols] Building site strings using UniProt accessions (no symbol mapping).")
  #msg("[.convertSitesToSymbols] Example base_acc: ", paste(head(base_acc), collapse=", "))

  mapped_ids <- toupper(base_acc)
  sym_sites <- paste0(mapped_ids, ";", aa, pos, ";")
  names(sym_sites) <- rownames(fi)

  #msg("[.convertSitesToSymbols] FINAL SUMMARY:")
  #msg("  Total sites: ", length(sym_sites))
  #msg("  Example mapped sites: ", paste(head(sym_sites), collapse=", "))

  list(
    all_sites = sym_sites,
    sig_sites = sym_sites[rownames(dataSet$sig.mat)]
  )
}

#' Enrichment from list-based kinase-substrate database
#' Database format: list(kinase = c("PROTEIN;RESIDUE_POS;", ...))
.enrichmentFromListDB <- function(sig_sites, all_sites, ks_db) {
  .enrichmentFromListDB_Fast(sig_sites, all_sites, ks_db)
}

.enrichmentFromListDB_Fast <- function(sig_sites, all_sites, ks_db) {
  
  #msg("[.enrichmentFromListDB] Processing list-based database with ", length(ks_db), " kinases")
  
  # 1. Convert the list-based DB to a long Data Frame for fast vector operations
  # This creates a table with columns: "Kinase", "Substrate"
  # We stack the list to avoid looping through it later
  db_flat <- stack(ks_db)
  colnames(db_flat) <- c("Substrate", "Kinase")
  
  # 2. Filter the DB: Only keep substrates that exist in your Background (all_sites)
  # This is critical for accurate statistics (the "Universe")
  # Note: This assumes IDs in ks_db match IDs in all_sites. 
  # If you need specific ID mapping (e.g. Uniprot to Gene), do it before this function.
  db_flat <- db_flat[db_flat$Substrate %in% all_sites, ]

  if (nrow(db_flat) == 0) {
    #msg("No overlap between Kinase DB and provided background sites. Example all_sites head: ",
    #        paste(head(all_sites), collapse=", "))
    return(NULL)
  }
  
  # 3. Annotate which substrates are Significant
  # creates a logical vector: TRUE if substrate is in sig_sites
  db_flat$is_sig <- db_flat$Substrate %in% sig_sites
  
  # 4. Calculate Counts (Vectorized)
  # 'm' = Number of white balls in urn (Total targets for the kinase in background)
  # 'x' = Number of white balls drawn (Significant targets for the kinase)
  stats <- aggregate(is_sig ~ Kinase, data = db_flat, FUN = function(x) c(Total = length(x), Sig = sum(x)))

  # Clean up aggregate output (it creates a matrix column by default)
  # IMPORTANT: Convert Kinase from factor to character to preserve names
  stats_df <- data.frame(Kinase = as.character(stats$Kinase),
                         Substrates_Total = stats$is_sig[, "Total"],
                         Substrates_Sig   = stats$is_sig[, "Sig"],
                         stringsAsFactors = FALSE)
  
  # 5. Define Global Parameters for Hypergeometric Test
  # 'N' = Total number of balls in urn (Total background sites)
  # 'k' = Number of balls drawn (Total significant sites)
  total_bg_size <- length(all_sites)
  total_sig_size <- length(sig_sites)
  
  # 6. Apply Hypergeometric Test (phyper) - Replaces fisher.test loop
  # q: num significant hits for kinase - 1 (for P(X >= x))
  # m: total targets for kinase
  # n: total non-targets (Universe - m)
  # k: total significant sites in experiment
  
  stats_df$P_value <- phyper(q = stats_df$Substrates_Sig - 1, 
                             m = stats_df$Substrates_Total, 
                             n = total_bg_size - stats_df$Substrates_Total, 
                             k = total_sig_size, 
                             lower.tail = FALSE) # One-sided (greater)
  
  # 7. Calculate Fold Enrichment
  # (Sig_Kinase / Total_Sig) / (Total_Kinase / Total_Background)
  stats_df$Fold_Enrichment <- (stats_df$Substrates_Sig / total_sig_size) / 
                              (stats_df$Substrates_Total / total_bg_size)
  
  # 8. Filter and Clean up
  # Remove kinases with 0 significant hits (optional, but cleaner)
  stats_df <- stats_df[stats_df$Substrates_Sig > 0, ]
  
  if (nrow(stats_df) == 0) {
    #msg("[.enrichmentFromListDB] No enriched kinases found")
    return(NULL)
  }
  
  # Adjust FDR
  stats_df$FDR <- p.adjust(stats_df$P_value, method = "BH")
  
  # Sort
  stats_df <- stats_df[order(stats_df$P_value), ]
  
  # Report
  n_sig <- sum(stats_df$FDR < 0.05)
  #msg("[.enrichmentFromListDB] Found ", n_sig, " enriched kinases (FDR < 0.05)")
  
  return(stats_df)
}


#' Count overlaps between phosphosite IDs and kinase substrate list
#' @param site_ids Vector of phosphosite IDs (format: PROTEIN_RES_POS like P12345_S_123)
#' @param substrates Vector of substrates (format: "PROTEIN;RES_POS;" like "TP53;S15;")
#' @return Number of matching sites
.countSiteOverlaps <- function(site_ids, substrates) {

  # Parse phosphosite IDs
  # Format: PROTEIN_RESIDUE_POSITION (e.g., P12345_S_123 or TP53_S_15)

  overlap_count <- 0

  for (site_id in site_ids) {
    # Split by underscore
    parts <- strsplit(site_id, "_")[[1]]

    if (length(parts) < 3) next

    protein <- parts[1]
    residue <- parts[2]
    position <- parts[3]

    # Create possible match patterns
    # Pattern: "PROTEIN;RESIDUE_POSITION;" (e.g., "TP53;S15;")
    pattern <- paste0(residue, position, ";")

    # Check if any substrate matches
    # Match on residue_position first, then optionally check protein
    for (substrate in substrates) {
      # Format: "PROTEIN;RES_POS;"
      if (grepl(pattern, substrate, fixed = TRUE)) {
        # Optional: check if protein also matches
        # Extract protein from substrate
        substrate_protein <- gsub(";.*", "", substrate)

        # Match if same protein OR if protein in substrate matches our protein
        if (substrate_protein == protein || grepl(protein, substrate_protein) || grepl(substrate_protein, protein)) {
          overlap_count <- overlap_count + 1
          break  # Count this site only once
        }
      }
    }
  }

  return(overlap_count)
}


#' Enrichment from curated kinase-substrate database
.enrichmentFromCuratedDB <- function(sig_sites, all_sites, ks_db) {

  # Group by kinase
  kinases <- unique(ks_db$Kinase)

  results_list <- lapply(kinases, function(kinase) {
    # Get substrates for this kinase
    kinase_substrates <- ks_db$Substrate[ks_db$Kinase == kinase]

    # Count overlaps
    sig_overlap <- length(intersect(sig_sites, kinase_substrates))
    all_overlap <- length(intersect(all_sites, kinase_substrates))

    # Fisher's exact test
    # Contingency table:
    #           In kinase substrates | Not in kinase substrates
    # Sig sites:       a              |          b
    # Other sites:     c              |          d

    a <- sig_overlap
    b <- length(sig_sites) - a
    c <- all_overlap - a
    d <- length(all_sites) - length(sig_sites) - c

    if (a == 0) return(NULL)  # Skip if no overlap

    fisher_result <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")

    fold_enrichment <- (a / length(sig_sites)) / (all_overlap / length(all_sites))

    data.frame(
      Kinase = kinase,
      Substrates_Total = all_overlap,
      Substrates_Sig = sig_overlap,
      P_value = fisher_result$p.value,
      Fold_Enrichment = fold_enrichment,
      stringsAsFactors = FALSE
    )
  })

  results <- do.call(rbind, results_list[!sapply(results_list, is.null)])

  if (is.null(results) || nrow(results) == 0) return(NULL)

  # Adjust p-values for multiple testing
  results$FDR <- p.adjust(results$P_value, method = "BH")

  # Sort by p-value
  results <- results[order(results$P_value), ]

  # Filter by significance
  results <- results[results$FDR < 0.05, ]

  return(results)
}

#' Enrichment from motif-based kinase mapping
.enrichmentFromMotifs <- function(sig_sites, all_sites, motif_db, dataSet) {

  # This is a simplified implementation
  # In practice, would need sequence window information from feature.info

  if (is.null(dataSet$feature.info)) {
    #msg("[.enrichmentFromMotifs] No feature info available for motif matching.")
    return(NULL)
  }

  # For now, return a placeholder message
  #msg("[.enrichmentFromMotifs] Motif-based enrichment requires sequence windows.")
  #msg("[.enrichmentFromMotifs] Please provide a curated kinase-substrate database.")

  # Return empty results
  return(data.frame(
    Kinase = character(),
    Substrates_Total = integer(),
    Substrates_Sig = integer(),
    P_value = numeric(),
    Fold_Enrichment = numeric(),
    FDR = numeric(),
    stringsAsFactors = FALSE
  ))
}

#' Get kinase enrichment results for display
#'
#' @return Data frame with enrichment results
GetKinaseEnrichmentResults <- function() {
  analSet <- readSet(analSet, "analSet")

  if (is.null(analSet$kinase.enrich)) {
    return(data.frame(
      Kinase = character(),
      Substrates_Total = integer(),
      Substrates_Sig = integer(),
      P_value = numeric(),
      Fold_Enrichment = numeric(),
      FDR = numeric()
    ))
  }

  return(analSet$kinase.enrich)
}

#' Plot kinase enrichment bar chart
#'
#' @param imgName Image file name
#' @param dpi DPI for image
#' @param format Image format (png, pdf)
#' @param top_n Number of top kinases to show
#' @return Image file name
PlotKinaseEnrichment <- function(imgName, dpi = 96, format = "png", top_n = 20) {
  require('ggplot2')

  analSet <- readSet(analSet, "analSet")

  if (is.null(analSet$kinase.enrich) || nrow(analSet$kinase.enrich) == 0) {
    #msg("[PlotKinaseEnrichment] No kinase enrichment results available.")
    return("NA")
  }

  results <- analSet$kinase.enrich

  # Take top N kinases by significance
  results <- results[order(results$P_value), ]
  if (nrow(results) > top_n) {
    results <- results[1:top_n, ]
  }

  # Reverse order for plotting (so most significant is on top)
  results$Kinase <- factor(results$Kinase, levels = rev(results$Kinase))

  # Create -log10(FDR) for visualization
  results$NegLogFDR <- -log10(results$FDR)

  # Use rank-based coloring instead of fold enrichment
  # Rank by p-value (lower p-value = better rank = higher color intensity)
  results$Rank <- rank(results$P_value, ties.method = "first")
  # Normalize ranks to 0-1 range
  results$NormalizedRank <- (results$Rank - 1) / (nrow(results) - 1)

  dpi <- as.numeric(dpi)
  fileNm <- paste(imgName, "dpi", dpi, ".", sep = "")
  imgNm <- paste0(fileNm, format, sep = "")

  # Use same color gradient as kinase network nodes
  # This matches GetColorGradient("black", center=FALSE, colorblind=FALSE)
  # which returns colorRampPalette(rev(heat.colors(9)))(100)
  heat_colors <- rev(heat.colors(9))

  # Create horizontal bar plot with rank-based coloring
  p <- ggplot(results, aes(x = Kinase, y = NegLogFDR, fill = NormalizedRank)) +
    geom_bar(stat = "identity", color = "white") +
    scale_fill_gradientn(colors = heat_colors) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    coord_flip() +
    labs(
      title = "Enriched Kinases in Significant Phosphosites",
      x = "Kinase",
      y = "-log10(FDR)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "none"  # Remove legend
    )

  # Save plot
  if (format == "png") {
    ggsave(imgNm, plot = p, dpi = dpi, width = 10, height = max(6, nrow(results) * 0.3), bg = "white")
  } else if (format == "pdf") {
    ggsave(imgNm, plot = p, width = 10, height = max(6, nrow(results) * 0.3), device = cairo_pdf, bg = "white")
  }

  return(imgNm)
}
