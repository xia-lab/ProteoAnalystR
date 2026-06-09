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

  if ("Substrates_Total" %in% colnames(enrichment_results)) {
    enrichment_results <- enrichment_results[as.numeric(enrichment_results$Substrates_Total) >= 5, ]
  }

  if (nrow(enrichment_results) == 0) {
    msgSet$current.msg <- "No kinases with sufficient substrate coverage found (minimum 5 substrates required)."
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
  unsupported <- c("PhosphoELM", "PhosphoPICK")
  opts <- opts[!tolower(opts) %in% tolower(unsupported)]
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
  db_obj <- ov_qs_read(db_file)

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

.phosphoCleanSequence <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  gsub("[^A-Za-z]", "", x)
}

.phosphoContextFromWindow <- function(window, residue = NA_character_) {
  raw <- as.character(window)
  if (is.na(raw) || raw == "") return(NULL)
  clean <- .phosphoCleanSequence(raw)
  if (nchar(clean) == 0) return(NULL)

  chars.raw <- strsplit(raw, "")[[1]]
  letters.raw <- chars.raw[grepl("[A-Za-z]", chars.raw)]
  lower.idx <- which(letters.raw %in% c("s", "t", "y"))
  center <- if (length(lower.idx) > 0) lower.idx[1] else ceiling(nchar(clean) / 2)

  clean <- toupper(clean)
  if (center < 1 || center > nchar(clean)) return(NULL)
  center.res <- substr(clean, center, center)
  residue <- toupper(as.character(residue))
  if (!is.na(residue) && residue %in% c("S", "T", "Y") && !(center.res %in% c("S", "T", "Y"))) {
    hits <- gregexpr(residue, clean, fixed = TRUE)[[1]]
    if (hits[1] != -1) {
      center <- hits[which.min(abs(hits - ceiling(nchar(clean) / 2)))]
      center.res <- substr(clean, center, center)
    }
  }
  if (!(center.res %in% c("S", "T", "Y"))) return(NULL)

  list(sequence = clean, center = center, residue = center.res)
}

.phosphoContextFromModifiedPeptide <- function(mod.seq, stripped.seq = NA_character_) {
  mod.seq <- as.character(mod.seq)
  if (is.na(mod.seq) || mod.seq == "") return(NULL)
  phospho.patterns <- c("([STY])\\(Phospho[^)]*\\)", "([STY])\\(UniMod:21\\)")
  hit.starts <- integer(0)
  for (pat in phospho.patterns) {
    hits <- gregexpr(pat, mod.seq, perl = TRUE)[[1]]
    if (hits[1] != -1) hit.starts <- c(hit.starts, hits)
  }
  hit.starts <- sort(unique(hit.starts))
  if (length(hit.starts) == 0) return(NULL)

  center.start <- hit.starts[1]
  residue <- substr(mod.seq, center.start, center.start)
  before <- substr(mod.seq, 1, center.start)
  clean.before <- gsub("\\([^)]+\\)", "", before)
  center <- nchar(.phosphoCleanSequence(clean.before))

  clean.seq <- if (!is.na(stripped.seq) && as.character(stripped.seq) != "") {
    .phosphoCleanSequence(stripped.seq)
  } else {
    .phosphoCleanSequence(gsub("\\([^)]+\\)", "", mod.seq))
  }
  clean.seq <- toupper(clean.seq)
  if (nchar(clean.seq) == 0 || center < 1 || center > nchar(clean.seq)) return(NULL)
  list(sequence = clean.seq, center = center, residue = toupper(residue))
}

.getRelativeAA <- function(sequence, center, offset) {
  pos <- center + offset
  if (pos < 1 || pos > nchar(sequence)) return(NA_character_)
  substr(sequence, pos, pos)
}

.extractPhosphoMotifContexts <- function(dataSet) {
  fi <- dataSet$feature.info
  if (is.null(fi) || nrow(fi) == 0) return(data.frame())
  if (is.null(rownames(fi))) {
    if ("id" %in% colnames(fi)) rownames(fi) <- as.character(fi$id)
  }

  all.ids <- rownames(dataSet$data.norm)
  fi <- fi[rownames(fi) %in% all.ids, , drop = FALSE]
  if (nrow(fi) == 0) return(data.frame())

  first_col <- function(cands) {
    hit <- cands[cands %in% colnames(fi)]
    if (length(hit) > 0) hit[1] else NA_character_
  }

  seq.col <- first_col(c("Sequence window", "SequenceWindow", "Sequence.Window", "sequence_window", "Window"))
  mod.col <- first_col(c("Modified.Sequence", "Modified sequence", "Modified.Sequence.Unique", "Modified Peptide"))
  strip.col <- first_col(c("Stripped.Sequence", "Stripped sequence", "Peptide", "Sequence"))
  residue.col <- first_col(c("Amino acid", "Residue", "residue"))
  position.col <- first_col(c("Position", "position"))
  gene.col <- first_col(c("Gene names", "Gene", "Genes", "Symbol"))
  protein.col <- first_col(c("Proteins", "Protein", "ProteinID", "UniProt"))

  rows <- lapply(seq_len(nrow(fi)), function(i) {
    residue <- if (!is.na(residue.col)) fi[[residue.col]][i] else NA_character_
    ctx <- NULL
    source <- NA_character_
    if (!is.na(seq.col)) {
      ctx <- .phosphoContextFromWindow(fi[[seq.col]][i], residue)
      if (!is.null(ctx)) source <- seq.col
    }
    if (is.null(ctx) && !is.na(mod.col)) {
      stripped <- if (!is.na(strip.col)) fi[[strip.col]][i] else NA_character_
      ctx <- .phosphoContextFromModifiedPeptide(fi[[mod.col]][i], stripped)
      if (!is.null(ctx)) source <- mod.col
    }
    if (is.null(ctx)) return(NULL)

    data.frame(
      Site = rownames(fi)[i],
      Protein = if (!is.na(protein.col)) as.character(fi[[protein.col]][i]) else NA_character_,
      Gene = if (!is.na(gene.col)) as.character(fi[[gene.col]][i]) else NA_character_,
      Residue = ctx$residue,
      Position = if (!is.na(position.col)) as.character(fi[[position.col]][i]) else NA_character_,
      SequenceContext = ctx$sequence,
      Center = ctx$center,
      ContextSource = source,
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

.classifyPhosphoMotifs <- function(contexts) {
  if (is.null(contexts) || nrow(contexts) == 0) return(data.frame())
  motif.defs <- list(
    list(Name = "Proline-directed S/T-P", Pattern = "pS/pT followed by P", Family = "Proline-directed",
         Match = function(seq, cen, res) res %in% c("S", "T") && identical(.getRelativeAA(seq, cen, 1), "P")),
    list(Name = "SP motif", Pattern = "pS-P", Family = "Proline-directed",
         Match = function(seq, cen, res) res == "S" && identical(.getRelativeAA(seq, cen, 1), "P")),
    list(Name = "TP motif", Pattern = "pT-P", Family = "Proline-directed",
         Match = function(seq, cen, res) res == "T" && identical(.getRelativeAA(seq, cen, 1), "P")),
    list(Name = "Basophilic upstream", Pattern = "R/K at -1 to -3", Family = "Basophilic",
         Match = function(seq, cen, res) res %in% c("S", "T") && any(vapply(-3:-1, function(k) .getRelativeAA(seq, cen, k) %in% c("R", "K"), logical(1)), na.rm = TRUE)),
    list(Name = "RxxS/T", Pattern = "R at -3 before pS/pT", Family = "Basophilic",
         Match = function(seq, cen, res) res %in% c("S", "T") && identical(.getRelativeAA(seq, cen, -3), "R")),
    list(Name = "Acidic downstream", Pattern = "D/E at +1 to +3", Family = "Acidic",
         Match = function(seq, cen, res) res %in% c("S", "T") && any(vapply(1:3, function(k) .getRelativeAA(seq, cen, k) %in% c("D", "E"), logical(1)), na.rm = TRUE)),
    list(Name = "CK2-like S/TxxD/E", Pattern = "D/E at +3 after pS/pT", Family = "Acidic",
         Match = function(seq, cen, res) res %in% c("S", "T") && .getRelativeAA(seq, cen, 3) %in% c("D", "E")),
    list(Name = "Tyrosine phosphorylation", Pattern = "pY", Family = "Tyrosine",
         Match = function(seq, cen, res) res == "Y")
  )

  motif.mat <- do.call(cbind, lapply(motif.defs, function(def) {
    vapply(seq_len(nrow(contexts)), function(i) {
      isTRUE(def$Match(contexts$SequenceContext[i], contexts$Center[i], contexts$Residue[i]))
    }, logical(1))
  }))
  colnames(motif.mat) <- vapply(motif.defs, `[[`, character(1), "Name")

  long <- do.call(rbind, lapply(seq_along(motif.defs), function(j) {
    hit <- motif.mat[, j]
    if (!any(hit)) return(NULL)
    data.frame(
      Site = contexts$Site[hit],
      Motif = motif.defs[[j]]$Name,
      MotifFamily = motif.defs[[j]]$Family,
      Pattern = motif.defs[[j]]$Pattern,
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(long)) long <- data.frame(Site = character(), Motif = character(), MotifFamily = character(), Pattern = character())

  contexts$Motifs <- vapply(seq_len(nrow(contexts)), function(i) {
    hits <- colnames(motif.mat)[motif.mat[i, ]]
    if (length(hits) == 0) "Unclassified" else paste(hits, collapse = "; ")
  }, character(1))

  list(contexts = contexts, long = long, motif.defs = motif.defs)
}

CanPerformMotifEnrichment <- function(dataName) {
  dataSet <- readDataset(dataName)
  paramSet <- readSet(paramSet, "paramSet")

  if (is.null(paramSet$data.type) || paramSet$data.type != "phospho") return(0L)
  if (is.null(dataSet$sig.mat) || nrow(dataSet$sig.mat) == 0) return(0L)
  if (is.null(dataSet$data.norm) || nrow(dataSet$data.norm) == 0) return(0L)

  contexts <- .extractPhosphoMotifContexts(dataSet)
  if (is.null(contexts) || nrow(contexts) == 0) return(0L)

  sig.sites <- intersect(rownames(dataSet$sig.mat), contexts$Site)
  all.sites <- intersect(rownames(dataSet$data.norm), contexts$Site)
  if (length(sig.sites) == 0 || length(all.sites) < 10) return(0L)

  classified <- .classifyPhosphoMotifs(contexts[contexts$Site %in% all.sites, , drop = FALSE])
  motif.long <- classified$long
  if (is.null(motif.long) || nrow(motif.long) == 0) return(0L)

  motif.names <- unique(motif.long$Motif)
  usable <- any(vapply(motif.names, function(motif) {
    sub <- motif.long[motif.long$Motif == motif, , drop = FALSE]
    bg.hit <- length(unique(intersect(sub$Site, all.sites)))
    sig.hit <- length(unique(intersect(sub$Site, sig.sites)))
    bg.hit >= 5 && sig.hit > 0
  }, logical(1)))
  if (!usable) return(0L)

  return(1L)
}

PerformMotifEnrichment <- function(dataName) {
  dataSet <- readDataset(dataName)
  paramSet <- readSet(paramSet, "paramSet")
  analSet <- readSet(analSet, "analSet")
  msgSet <- readSet(msgSet, "msgSet")

  if (is.null(paramSet$data.type) || paramSet$data.type != "phospho") {
    msgSet$current.msg <- "Motif analysis is only available for phosphoproteomics data."
    saveSet(msgSet, "msgSet")
    return(0)
  }
  if (is.null(dataSet$sig.mat) || nrow(dataSet$sig.mat) == 0) {
    msgSet$current.msg <- "No significant phosphosites found. Please perform differential expression analysis first."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  contexts <- .extractPhosphoMotifContexts(dataSet)
  if (is.null(contexts) || nrow(contexts) == 0) {
    msgSet$current.msg <- "Motif analysis requires sequence windows or modified peptide sequences in feature information. The current dataset does not contain usable sequence context."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  sig.sites <- intersect(rownames(dataSet$sig.mat), contexts$Site)
  all.sites <- intersect(rownames(dataSet$data.norm), contexts$Site)
  if (length(sig.sites) == 0 || length(all.sites) < 10) {
    msgSet$current.msg <- "Insufficient significant phosphosites with usable sequence context for motif enrichment."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  classified <- .classifyPhosphoMotifs(contexts[contexts$Site %in% all.sites, , drop = FALSE])
  motif.long <- classified$long
  if (nrow(motif.long) == 0) {
    msgSet$current.msg <- "No recognized phosphosite motif classes were found in the available sequence contexts."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  total.bg <- length(unique(all.sites))
  total.sig <- length(unique(sig.sites))
  motif.names <- unique(motif.long$Motif)

  results <- do.call(rbind, lapply(motif.names, function(motif) {
    sub <- motif.long[motif.long$Motif == motif, , drop = FALSE]
    bg.hit.sites <- unique(intersect(sub$Site, all.sites))
    sig.hit.sites <- unique(intersect(sub$Site, sig.sites))
    bg.hit <- length(bg.hit.sites)
    sig.hit <- length(sig.hit.sites)
    if (bg.hit < 5 || sig.hit == 0) return(NULL)
    p.val <- phyper(sig.hit - 1, bg.hit, total.bg - bg.hit, total.sig, lower.tail = FALSE)
    fold <- (sig.hit / total.sig) / (bg.hit / total.bg)
    data.frame(
      Motif = motif,
      Family = sub$MotifFamily[1],
      Pattern = sub$Pattern[1],
      Background_Sites = bg.hit,
      Significant_Sites = sig.hit,
      P_value = p.val,
      Fold_Enrichment = fold,
      Example_Sites = paste(head(sig.hit.sites, 8), collapse = "; "),
      stringsAsFactors = FALSE
    )
  }))

  if (is.null(results) || nrow(results) == 0) {
    msgSet$current.msg <- "No motif classes had enough significant and background sites for enrichment testing."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  results$FDR <- p.adjust(results$P_value, method = "BH")
  results <- results[order(results$P_value), ]

  context.out <- classified$contexts
  context.out$Significant <- context.out$Site %in% sig.sites

  # Annotate each site with its cellular compartment (best-effort)
  context.out$Compartment <- .mapSitesToCompartment(context.out$Protein, paramSet, analSet)

  analSet$motif.enrich   <- results
  analSet$motif.context  <- context.out
  analSet$motif.long     <- classified$long   # one row per site-motif assignment
  saveSet(analSet, "analSet")

  fast.write(results, file = "motif_enrichment.csv", row.names = FALSE)
  fast.write(context.out, file = "phosphosite_motif_context.csv", row.names = FALSE)

  msgSet$current.msg <- paste0(
    "Motif analysis complete. Tested ", nrow(results),
    " motif class(es) using ", total.sig, " significant phosphosite(s) and ",
    total.bg, " background phosphosite(s) with sequence context."
  )
  saveSet(msgSet, "msgSet")
  return(1)
}

# Write motif enrichment data as Plotly-ready JSON for the interactive browser chart.
# Includes both flat enrichment results and (when possible) per-compartment enrichment.
GetMotifEnrichmentJSON <- function(jsonNm) {
  require(jsonlite)
  analSet  <- readSet(analSet, "analSet")
  paramSet <- readSet(paramSet, "paramSet")
  .dbg <- function(...) {
    txt <- sprintf(...)
    message(txt)
    cat(paste0(Sys.time(), " ", txt, "\n"), file = "motif_enrich_debug.log", append = TRUE)
  }
  .dbg("[GetMotifEnrichmentJSON] called with jsonNm=%s; motif.enrich rows=%s; wd=%s",
       jsonNm, if (!is.null(analSet$motif.enrich)) nrow(analSet$motif.enrich) else "NULL",
       getwd())
  if (is.null(analSet$motif.enrich) || nrow(analSet$motif.enrich) == 0) return("NA")

  res <- analSet$motif.enrich
  res$NegLogFDR         <- round(-log10(pmax(as.numeric(res$FDR), .Machine$double.xmin)), 4)
  res$Fold_Enrichment   <- round(as.numeric(res$Fold_Enrichment), 3)
  res$P_value           <- signif(as.numeric(res$P_value), 3)
  res$FDR               <- signif(as.numeric(res$FDR), 3)

  motifs.list <- lapply(seq_len(nrow(res)), function(i) {
    r <- res[i, , drop = FALSE]
    list(
      Motif             = r$Motif,
      Family            = r$Family,
      Pattern           = if ("Pattern" %in% names(r)) r$Pattern else "",
      Background_Sites  = as.integer(r$Background_Sites),
      Significant_Sites = as.integer(r$Significant_Sites),
      P_value           = r$P_value,
      FDR               = r$FDR,
      Fold_Enrichment   = r$Fold_Enrichment,
      NegLogFDR         = r$NegLogFDR,
      Example_Sites     = if ("Example_Sites" %in% names(r) && !is.na(r$Example_Sites)) r$Example_Sites else ""
    )
  })

  # Per-compartment enrichment (recomputed from motif.long + motif.context compartment column)
  comp.rows <- list()
  has.compartment <- FALSE

  ml  <- analSet$motif.long
  ctx <- analSet$motif.context
  if (!is.null(ml) && !is.null(ctx) && "Compartment" %in% colnames(ctx) &&
      any(!is.na(ctx$Compartment))) {

    site.comp <- setNames(as.character(ctx$Compartment), as.character(ctx$Site))
    sig.mask  <- setNames(as.logical(ctx$Significant),   as.character(ctx$Site))
    ml$Compartment <- site.comp[as.character(ml$Site)]

    comps <- sort(unique(ml$Compartment[!is.na(ml$Compartment)]))
    if (length(comps) > 0) {
      has.compartment <- TRUE
      for (comp in comps) {
        sub.ml     <- ml[!is.na(ml$Compartment) & ml$Compartment == comp, , drop = FALSE]
        all.sites  <- unique(sub.ml$Site)
        sig.sites  <- unique(sub.ml$Site[!is.na(sig.mask[sub.ml$Site]) & sig.mask[sub.ml$Site]])
        total.bg   <- length(all.sites)
        total.sig  <- length(sig.sites)
        if (total.sig == 0) next

        for (motif in unique(sub.ml$Motif)) {
          m.sub   <- sub.ml[sub.ml$Motif == motif, , drop = FALSE]
          bg.hit  <- length(unique(intersect(m.sub$Site, all.sites)))
          sig.hit <- length(unique(intersect(m.sub$Site, sig.sites)))
          if (bg.hit < 3 || sig.hit == 0) next
          p.val <- phyper(sig.hit - 1, bg.hit, total.bg - bg.hit, total.sig, lower.tail = FALSE)
          fold  <- (sig.hit / total.sig) / (bg.hit / total.bg)
          fdr   <- NA  # FDR computed after
          comp.rows[[length(comp.rows) + 1]] <- list(
            Compartment       = comp,
            Motif             = motif,
            Family            = m.sub$MotifFamily[1],
            Background_Sites  = bg.hit,
            Significant_Sites = sig.hit,
            P_value           = signif(p.val, 3),
            Fold_Enrichment   = round(fold, 3),
            NegLogFDR         = NA_real_  # filled below
          )
        }
      }
      # FDR per compartment, fill NegLogFDR
      if (length(comp.rows) > 0) {
        all.p <- vapply(comp.rows, `[[`, numeric(1), "P_value")
        all.f <- p.adjust(all.p, method = "BH")
        for (k in seq_along(comp.rows)) {
          comp.rows[[k]]$FDR      <- signif(all.f[k], 3)
          comp.rows[[k]]$NegLogFDR <- round(-log10(max(all.f[k], .Machine$double.xmin)), 4)
        }
      }
    }
  }

  # ── Site-level data for Manhattan plot ──────────────────────────────────────
  # Each phosphosite gets: motif category (x), -log10(adj.P.Val) from DE (y),
  # significance status, gene label, compartment (if available).
  # Sites assigned to multiple motifs appear once each (deduplicated to primary motif).
  sites.list <- list()
  ml  <- if (is.null(ml))  analSet$motif.long    else ml
  ctx <- if (is.null(ctx)) analSet$motif.context else ctx

  if (!is.null(ml) && nrow(ml) > 0) {
    # Deduplicate: one motif per site (first assignment when ordered by Family priority)
    fam.order <- c("Proline-directed", "Basophilic", "Acidophilic", "Hydrophobic", "Other")
    ml$FamRank <- match(ml$MotifFamily, fam.order)
    ml$FamRank[is.na(ml$FamRank)] <- 99L
    ml <- ml[order(ml$Site, ml$FamRank), , drop = FALSE]
    ml.primary <- ml[!duplicated(ml$Site), , drop = FALSE]

    # Pull per-site DE stats from comp.res
    ds <- tryCatch(readDataset(paramSet$dataName, quiet = TRUE), error = function(e) NULL)
    cr <- if (!is.null(ds)) ds$comp.res else NULL

    site.fdr  <- rep(NA_real_, nrow(ml.primary))
    site.fc   <- rep(NA_real_, nrow(ml.primary))
    if (!is.null(cr) && "adj.P.Val" %in% colnames(cr)) {
      idx <- match(as.character(ml.primary$Site), rownames(cr))
      site.fdr <- suppressWarnings(as.numeric(cr$adj.P.Val[idx]))
      if ("logFC" %in% colnames(cr))
        site.fc <- suppressWarnings(as.numeric(cr$logFC[idx]))
    }

    # Attach compartment
    site.comp <- if (!is.null(ctx) && "Compartment" %in% colnames(ctx))
      setNames(as.character(ctx$Compartment), as.character(ctx$Site)) else character(0)
    site.sig <- if (!is.null(ctx) && "Significant" %in% colnames(ctx))
      setNames(as.logical(ctx$Significant), as.character(ctx$Site)) else logical(0)

    # Gene lookup from motif context
    site.gene <- if (!is.null(ctx) && "Gene" %in% colnames(ctx))
      setNames(as.character(ctx$Gene), as.character(ctx$Site)) else character(0)

    sites.list <- lapply(seq_len(nrow(ml.primary)), function(i) {
      s    <- as.character(ml.primary$Site[i])
      fdr  <- site.fdr[i]
      nlp  <- if (!is.na(fdr) && fdr > 0) round(-log10(fdr), 4) else NA_real_
      list(
        Site        = s,
        Motif       = ml.primary$Motif[i],
        Family      = ml.primary$MotifFamily[i],
        NegLogFDR   = nlp,
        logFC       = if (!is.na(site.fc[i])) round(site.fc[i], 3) else NA_real_,
        FDR         = if (!is.na(fdr)) signif(fdr, 3) else NA_real_,
        Significant = isTRUE(site.sig[s]),
        Gene        = if (!is.na(site.gene[s])) site.gene[s] else "",
        Compartment = if (length(site.comp) > 0 && !is.na(site.comp[s])) site.comp[s] else ""
      )
    })
  }

  # ── Output ───────────────────────────────────────────────────────────────
  out <- list(
    motifs            = motifs.list,       # per-motif enrichment summary (bubble view)
    sites             = sites.list,        # per-site DE stats (Manhattan view)
    compartmentMotifs = comp.rows,         # per-compartment enrichment (compartment view)
    hasCompartment    = has.compartment,
    hasSites          = length(sites.list) > 0
  )

  json.path <- paste0(jsonNm, ".json")
  .dbg("[GetMotifEnrichmentJSON] writing to: %s", json.path)
  tryCatch(
    jsonlite::write_json(out, json.path, auto_unbox = TRUE, null = "null"),
    error = function(e) .dbg("[GetMotifEnrichmentJSON] write_json ERROR: %s", e$message)
  )
  .dbg("[GetMotifEnrichmentJSON] file.exists=%s; motifs=%d; sites=%d",
       file.exists(json.path), length(out$motifs), length(out$sites))
  return(json.path)
}

# Map a vector of protein IDs (UniProt or mixed format) to Broad.category compartment labels.
# Returns a character vector of equal length, NA for unmapped proteins.
.mapSitesToCompartment <- function(protein.ids, paramSet, analSet) {
  result <- rep(NA_character_, length(protein.ids))
  if (is.null(protein.ids) || length(protein.ids) == 0) return(result)

  org <- if (!is.null(paramSet$data.org) && nzchar(paramSet$data.org)) paramSet$data.org else "hsa"
  loc.path <- paste0(paramSet$lib.path, org, "/", org, "_localization.qs")
  if (!file.exists(loc.path)) return(result)
  loc <- tryCatch(ov_qs_read(loc.path), error = function(e) NULL)
  if (is.null(loc) || !all(c("EntrezID", "Broad.category") %in% colnames(loc))) return(result)

  loc$Broad.category[is.na(loc$Broad.category) | loc$Broad.category == ""] <- "Unknown"
  if (!"Main.location" %in% colnames(loc)) {
    loc$Main.location <- "Unknown"
  }
  comp.res <- .paResolveCompartmentAnnotations(loc$Broad.category, loc$Main.location)
  comp.map <- setNames(comp.res$primary, as.character(loc$EntrezID))

  # Build UniProt -> Entrez map: prefer dataset map, fall back to gene DB
  up2ent <- analSet$uniprot_to_entrez_map
  if (is.null(up2ent) || length(up2ent) == 0) {
    gdb <- tryCatch(queryGeneDB("entrez_uniprot", org), error = function(e) NULL)
    if (!is.null(gdb) && is.data.frame(gdb) && "gene_id" %in% colnames(gdb) && "accession" %in% colnames(gdb))
      up2ent <- setNames(as.character(gdb$gene_id), as.character(gdb$accession))
  }

  clean_id <- function(x) {
    # strip sp|ACC|GENE or ACC_SPECIES formats to primary accession
    x <- sub("^[a-z]{2}\\|([A-Z0-9_-]+)\\|.*", "\\1", x)
    sub("_[A-Z]+$", "", x)
  }

  for (i in seq_along(protein.ids)) {
    pid <- as.character(protein.ids[i])
    if (is.na(pid) || !nzchar(pid)) next
    ent <- NA_character_
    if (!is.null(up2ent)) {
      ent <- up2ent[pid]
      if (is.na(ent)) ent <- up2ent[clean_id(pid)]
    }
    if (!is.na(ent) && nzchar(ent)) {
      comp <- comp.map[ent]
      if (!is.na(comp)) result[i] <- comp
    }
  }
  result
}

GetMotifEnrichmentResults <- function() {
  analSet <- readSet(analSet, "analSet")
  if (is.null(analSet$motif.enrich)) {
    return(data.frame(
      Motif = character(),
      Family = character(),
      Pattern = character(),
      Background_Sites = integer(),
      Significant_Sites = integer(),
      P_value = numeric(),
      Fold_Enrichment = numeric(),
      FDR = numeric(),
      Example_Sites = character()
    ))
  }
  analSet$motif.enrich
}

PlotMotifEnrichment <- function(imgName, dpi = 96, format = "png", top_n = 20, byCompartment = FALSE) {
  require(ggplot2)
  require(Cairo)
  analSet <- readSet(analSet, "analSet")
  if (is.null(analSet$motif.enrich) || nrow(analSet$motif.enrich) == 0) return("NA")
  dpi <- suppressWarnings(as.numeric(dpi))
  if (is.na(dpi) || dpi <= 0) dpi <- 96

  # ── Palette shared between both modes ──────────────────────────────────────
  family.cols <- c(
    "Proline-directed"  = "#4e79a7",
    "Basophilic"        = "#f28e2b",
    "Acidophilic"       = "#e15759",
    "Hydrophobic"       = "#76b7b2",
    "Other"             = "#b07aa1"
  )

  # ── Helper: build enrichment rows for one set of sites ─────────────────────
  .enrich_one <- function(motif.long, sig.mask, top_n) {
    motif.names <- unique(motif.long$Motif)
    all.sites  <- unique(motif.long$Site)
    sig.sites  <- unique(motif.long$Site[motif.long$Site %in% names(sig.mask)[sig.mask]])
    total.bg   <- length(all.sites)
    total.sig  <- length(sig.sites)
    if (total.sig == 0) return(NULL)
    rows <- lapply(motif.names, function(motif) {
      sub      <- motif.long[motif.long$Motif == motif, , drop = FALSE]
      bg.hit   <- length(unique(intersect(sub$Site, all.sites)))
      sig.hit  <- length(unique(intersect(sub$Site, sig.sites)))
      if (bg.hit < 3 || sig.hit == 0) return(NULL)
      p.val <- phyper(sig.hit - 1, bg.hit, total.bg - bg.hit, total.sig, lower.tail = FALSE)
      fold  <- (sig.hit / total.sig) / (bg.hit / total.bg)
      data.frame(
        Motif              = motif,
        Family             = sub$MotifFamily[1],
        Background_Sites   = bg.hit,
        Significant_Sites  = sig.hit,
        P_value            = p.val,
        Fold_Enrichment    = fold,
        stringsAsFactors   = FALSE
      )
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) == 0) return(NULL)
    out <- do.call(rbind, rows)
    out$FDR <- p.adjust(out$P_value, method = "BH")
    out <- out[order(out$P_value), , drop = FALSE]
    if (nrow(out) > top_n) out <- out[seq_len(top_n), , drop = FALSE]
    out
  }

  # ── Standard (no compartment) mode ─────────────────────────────────────────
  if (!isTRUE(byCompartment)) {
    res <- analSet$motif.enrich
    res <- res[order(res$P_value), , drop = FALSE]
    if (nrow(res) > top_n) res <- res[seq_len(top_n), , drop = FALSE]
    res$NegLogFDR       <- -log10(pmax(as.numeric(res$FDR), .Machine$double.xmin))
    res$Significant_Sites <- as.numeric(res$Significant_Sites)
    res$Background_Sites  <- as.numeric(res$Background_Sites)
    res$Motif <- factor(res$Motif, levels = res$Motif[order(res$NegLogFDR)])

    p <- ggplot(res, aes(x = Motif, y = NegLogFDR,
                         size = Significant_Sites, fill = Family)) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                 color = "#999999", linewidth = 0.4) +
      geom_point(shape = 21, color = "#333333", stroke = 0.4, alpha = 0.88) +
      geom_text(aes(label = paste0(Significant_Sites, "/", Background_Sites)),
                vjust = -1.2, size = 2.8, color = "#444444") +
      scale_size_continuous(range = c(4, 14), name = "Sig. sites") +
      scale_fill_manual(values = family.cols, name = "Motif family",
                        na.value = "#aaaaaa") +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      theme_bw(base_size = 11) +
      theme(axis.text.x    = element_text(angle = 35, hjust = 1),
            panel.grid.major.x = element_blank(),
            panel.grid.minor   = element_blank(),
            legend.position    = "right",
            plot.title         = element_text(hjust = 0.5, face = "bold")) +
      labs(title = "Phosphosite Motif Enrichment",
           x = NULL, y = "-log₁₀(FDR)")

    n.motifs <- nrow(res)
    w <- max(5.5, 0.65 * n.motifs + 3.5)
    h <- 4.8

  } else {
    # ── Compartment mode ─────────────────────────────────────────────────────
    ctx     <- analSet$motif.context
    ml      <- analSet$motif.long
    if (is.null(ctx) || is.null(ml) || !"Compartment" %in% colnames(ctx)) {
      # Fall back to standard mode silently
      return(PlotMotifEnrichment(imgName, dpi = dpi, format = format,
                                 top_n = top_n, byCompartment = FALSE))
    }

    # Attach compartment to motif.long via site-level lookup
    site.comp <- setNames(as.character(ctx$Compartment), as.character(ctx$Site))
    ml$Compartment <- site.comp[as.character(ml$Site)]

    sig.mask <- setNames(as.logical(ctx$Significant), as.character(ctx$Site))

    comps <- sort(unique(ml$Compartment[!is.na(ml$Compartment)]))
    if (length(comps) == 0) {
      return(PlotMotifEnrichment(imgName, dpi = dpi, format = format,
                                 top_n = top_n, byCompartment = FALSE))
    }

    all.rows <- lapply(comps, function(comp) {
      sub.ml <- ml[!is.na(ml$Compartment) & ml$Compartment == comp, , drop = FALSE]
      rows <- .enrich_one(sub.ml, sig.mask, top_n)
      if (!is.null(rows)) { rows$Compartment <- comp; rows }
    })
    all.rows <- all.rows[!vapply(all.rows, is.null, logical(1))]
    if (length(all.rows) == 0) {
      return(PlotMotifEnrichment(imgName, dpi = dpi, format = format,
                                 top_n = top_n, byCompartment = FALSE))
    }
    res <- do.call(rbind, all.rows)
    res$NegLogFDR        <- -log10(pmax(as.numeric(res$FDR), .Machine$double.xmin))
    res$Significant_Sites <- as.numeric(res$Significant_Sites)
    res$Background_Sites  <- as.numeric(res$Background_Sites)

    # Order motifs within each compartment by FDR
    motif.order <- unique(res$Motif[order(res$NegLogFDR, decreasing = TRUE)])
    res$Motif <- factor(res$Motif, levels = motif.order)

    p <- ggplot(res, aes(x = Motif, y = NegLogFDR,
                         size = Significant_Sites, fill = Family)) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                 color = "#999999", linewidth = 0.4) +
      geom_point(shape = 21, color = "#333333", stroke = 0.4, alpha = 0.88) +
      scale_size_continuous(range = c(3, 12), name = "Sig. sites") +
      scale_fill_manual(values = family.cols, name = "Motif family",
                        na.value = "#aaaaaa") +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
      facet_wrap(~ Compartment, scales = "free_x", ncol = 2) +
      theme_bw(base_size = 10) +
      theme(axis.text.x    = element_text(angle = 35, hjust = 1, size = 8),
            panel.grid.major.x = element_blank(),
            panel.grid.minor   = element_blank(),
            legend.position    = "right",
            strip.text         = element_text(face = "bold"),
            plot.title         = element_text(hjust = 0.5, face = "bold")) +
      labs(title = "Motif Enrichment by Cellular Compartment",
           x = NULL, y = "-log₁₀(FDR)")

    n.comp <- length(comps)
    w <- 9.5
    h <- max(4.5, ceiling(n.comp / 2) * 3.5 + 1.2)
  }

  imgNm <- paste0(imgName, "dpi", dpi, ".", format)
  Cairo(file = imgNm, width = w, height = h, unit = "in",
        dpi = dpi, bg = "white", type = format)
  print(p)
  invisible(dev.off())
  imgNm
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

#' Plot kinase enrichment bubble chart
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

  # Create -log10(P_value) for visualization. Higher values are more significant.
  pval_for_plot <- pmax(as.numeric(results$P_value), .Machine$double.xmin)
  results$NegLogFDR <- -log10(pval_for_plot)
  results$FoldEnrichmentPlot <- as.numeric(results$Fold_Enrichment)
  results$FoldEnrichmentPlot[!is.finite(results$FoldEnrichmentPlot)] <- NA
  if (all(is.na(results$FoldEnrichmentPlot)) &&
      all(c("Substrates_Sig", "Substrates_Total") %in% colnames(results))) {
    results$FoldEnrichmentPlot <- as.numeric(results$Substrates_Sig) /
      pmax(as.numeric(results$Substrates_Total), 1)
  }
  results$FoldEnrichmentPlot[is.na(results$FoldEnrichmentPlot)] <- 0
  results$FoldEnrichmentPlot <- pmax(results$FoldEnrichmentPlot, 1e-3)
  results$HitCount <- if ("Substrates_Sig" %in% colnames(results)) {
    as.numeric(results$Substrates_Sig)
  } else {
    rep(1, nrow(results))
  }
  results$HitCount[!is.finite(results$HitCount) | results$HitCount < 1] <- 1

  dpi <- as.numeric(dpi)
  fileNm <- paste(imgName, "dpi", dpi, ".", sep = "")
  imgNm <- paste0(fileNm, format, sep = "")

  # Map significance directly to color: pale = less significant, red = more significant.
  significance_colors <- c("#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026")

  x_breaks <- pretty(c(0, max(results$FoldEnrichmentPlot, na.rm = TRUE)))
  x_upper <- max(x_breaks)

  # Create bubble plot: enrichment magnitude, support, and significance in one view.
  p <- ggplot(results, aes(x = FoldEnrichmentPlot, y = Kinase)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray70") +
    geom_point(aes(size = HitCount, color = NegLogFDR), alpha = 0.9) +
    scale_color_gradientn(colors = significance_colors, name = "-log10(p-value)") +
    scale_size_area(max_size = 11, name = "Sig. substrates") +
    scale_x_continuous(
      limits = c(0, x_upper),
      breaks = x_breaks,
      labels = function(x) format(x, trim = TRUE, scientific = FALSE)
    ) +
    labs(
      x = "Fold Enrichment",
      y = "Kinase"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      legend.position = "right"
    )

  # Save plot
  if (format == "png") {
    ggsave(imgNm, plot = p, dpi = dpi, width = 7, height = 6, bg = "white")
  } else if (format == "pdf") {
    ggsave(imgNm, plot = p, width = 7, height = 6, device = cairo_pdf, bg = "white")
  }

  # Store in imgSet for report/slides inclusion
  imgSet <- readSet(imgSet, "imgSet")
  imgSet$kinase_enrichment <- imgNm
  saveSet(imgSet, "imgSet")

  return(imgNm)
}

# ── Cellular compartment enrichment ───────────────────────────────────────────
# Tests whether significant phosphosites are over-represented in any of the
# broad cellular compartments from the localization database (same data used
# for the PPI network compartment layout).

CanPerformCompartmentEnrichment <- function(dataName) {
  dataSet  <- readDataset(dataName)
  paramSet <- readSet(paramSet, "paramSet")
  if (is.null(paramSet$data.type) || paramSet$data.type != "phospho") return(0L)
  if (is.null(dataSet$sig.mat)   || nrow(dataSet$sig.mat)   == 0) return(0L)
  if (is.null(dataSet$data.norm) || nrow(dataSet$data.norm) == 0) return(0L)
  # Need feature.info with protein IDs to map to compartments
  fi <- dataSet$feature.info
  if (is.null(fi) || nrow(fi) == 0) return(0L)
  return(1L)
}

PerformCompartmentEnrichment <- function(dataName) {
  dataSet  <- readDataset(dataName)
  paramSet <- readSet(paramSet, "paramSet")
  analSet  <- readSet(analSet,  "analSet")
  msgSet   <- readSet(msgSet,   "msgSet")

  if (is.null(paramSet$data.type) || paramSet$data.type != "phospho") {
    msgSet$current.msg <- "Compartment enrichment is only available for phosphoproteomics data."
    saveSet(msgSet, "msgSet"); return(0L)
  }
  if (is.null(dataSet$sig.mat) || nrow(dataSet$sig.mat) == 0) {
    msgSet$current.msg <- "No significant phosphosites found. Run differential expression analysis first."
    saveSet(msgSet, "msgSet"); return(0L)
  }

  all.sites <- rownames(dataSet$data.norm)
  sig.sites <- rownames(dataSet$sig.mat)

  # Get protein IDs for each site from feature.info
  fi <- dataSet$feature.info
  if (is.null(fi) || nrow(fi) == 0) {
    msgSet$current.msg <- "Feature information is missing. Cannot map phosphosites to compartments."
    saveSet(msgSet, "msgSet"); return(0L)
  }
  prot.col <- {
    cands <- c("Proteins", "Protein", "ProteinID", "UniProt")
    hit <- cands[cands %in% colnames(fi)]
    if (length(hit) > 0) hit[1] else NA_character_
  }
  if (is.na(prot.col)) {
    msgSet$current.msg <- "No protein ID column found in feature information."
    saveSet(msgSet, "msgSet"); return(0L)
  }

  # Map all sites to compartments
  prot.ids <- as.character(fi[[prot.col]][match(all.sites, rownames(fi))])
  comp.vec  <- .mapSitesToCompartment(prot.ids, paramSet, analSet)
  names(comp.vec) <- all.sites

  # Keep only sites that could be mapped (ignore Unknown and NA)
  mappable.mask <- !is.na(comp.vec) & comp.vec != "Unknown"
  all.mapped    <- all.sites[mappable.mask]
  sig.mapped    <- intersect(sig.sites, all.mapped)

  if (length(sig.mapped) == 0) {
    msgSet$current.msg <- "No significant phosphosites could be mapped to a cellular compartment."
    saveSet(msgSet, "msgSet"); return(0L)
  }

  total.bg  <- length(all.mapped)
  total.sig <- length(sig.mapped)
  comps     <- sort(unique(comp.vec[all.mapped]))

  results <- do.call(rbind, lapply(comps, function(comp) {
    bg.sites  <- all.mapped[comp.vec[all.mapped] == comp]
    sig.in    <- intersect(sig.mapped, bg.sites)
    bg.hit    <- length(bg.sites)
    sig.hit   <- length(sig.in)
    if (bg.hit == 0) return(NULL)
    p.val <- phyper(sig.hit - 1, bg.hit, total.bg - bg.hit, total.sig, lower.tail = FALSE)
    fold  <- (sig.hit / total.sig) / (bg.hit / total.bg)
    data.frame(
      Compartment      = comp,
      Total_Sites      = bg.hit,
      Sig_Sites        = sig.hit,
      P_value          = p.val,
      Fold_Enrichment  = fold,
      stringsAsFactors = FALSE
    )
  }))

  if (is.null(results) || nrow(results) == 0) {
    msgSet$current.msg <- "Compartment enrichment test produced no results."
    saveSet(msgSet, "msgSet"); return(0L)
  }

  results$FDR <- p.adjust(results$P_value, method = "BH")
  results <- results[order(results$P_value), ]

  analSet$comp.enrich <- results
  saveSet(analSet, "analSet")
  fast.write(results, file = "compartment_enrichment.csv", row.names = FALSE)

  msgSet$current.msg <- paste0(
    "Compartment enrichment complete. Tested ", nrow(results), " compartment(s) using ",
    total.sig, " significant and ", total.bg, " background phosphosite(s)."
  )
  saveSet(msgSet, "msgSet")
  return(1L)
}

GetCompartmentEnrichmentResults <- function() {
  analSet <- readSet(analSet, "analSet")
  if (is.null(analSet$comp.enrich)) {
    return(data.frame(
      Compartment     = character(),
      Total_Sites     = integer(),
      Sig_Sites       = integer(),
      P_value         = numeric(),
      Fold_Enrichment = numeric(),
      FDR             = numeric()
    ))
  }
  analSet$comp.enrich
}

GetCompartmentEnrichmentJSON <- function(jsonNm) {
  require(jsonlite)
  analSet <- readSet(analSet, "analSet")
  if (is.null(analSet$comp.enrich) || nrow(analSet$comp.enrich) == 0) return("NA")

  res <- analSet$comp.enrich
  res$NegLogFDR       <- round(-log10(pmax(as.numeric(res$FDR), .Machine$double.xmin)), 4)
  res$Fold_Enrichment <- round(as.numeric(res$Fold_Enrichment), 3)
  res$P_value         <- signif(as.numeric(res$P_value), 3)
  res$FDR             <- signif(as.numeric(res$FDR), 3)

  rows <- lapply(seq_len(nrow(res)), function(i) {
    r <- res[i, , drop = FALSE]
    list(
      Compartment     = r$Compartment,
      Total_Sites     = as.integer(r$Total_Sites),
      Sig_Sites       = as.integer(r$Sig_Sites),
      P_value         = r$P_value,
      FDR             = r$FDR,
      Fold_Enrichment = r$Fold_Enrichment,
      NegLogFDR       = r$NegLogFDR
    )
  })

  json.path <- paste0(jsonNm, ".json")
  jsonlite::write_json(list(compartments = rows), json.path, auto_unbox = TRUE, null = "null")
  return(json.path)
}

PlotCompartmentEnrichment <- function(imgName, dpi = 96, format = "png") {
  require(ggplot2)
  require(Cairo)
  analSet <- readSet(analSet, "analSet")
  if (is.null(analSet$comp.enrich) || nrow(analSet$comp.enrich) == 0) return("NA")
  dpi <- suppressWarnings(as.numeric(dpi)); if (is.na(dpi) || dpi <= 0) dpi <- 96

  res <- analSet$comp.enrich
  res$NegLogFDR <- -log10(pmax(as.numeric(res$FDR), .Machine$double.xmin))
  res$Compartment <- factor(res$Compartment, levels = res$Compartment[order(res$NegLogFDR)])

  p <- ggplot(res, aes(x = NegLogFDR, y = Compartment, fill = Fold_Enrichment)) +
    geom_col(width = 0.65) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = "#888888", linewidth = 0.4) +
    scale_fill_gradient2(low = "#4393c3", mid = "#f7f7f7", high = "#d6604d",
                         midpoint = 1, name = "Fold\nenrichment") +
    labs(x = expression(-log[10](FDR)), y = NULL,
         title = "Cellular Compartment Enrichment") +
    theme_bw(base_size = 11) +
    theme(panel.grid.major.y = element_blank(),
          plot.title = element_text(size = 12, colour = "#333333"))

  imgNm <- paste0(imgName, "dpi", dpi, ".", format)
  if (format == "png") {
    ggsave(imgNm, plot = p, dpi = dpi, width = 7, height = 4, bg = "white")
  } else {
    ggsave(imgNm, plot = p, width = 7, height = 4, device = cairo_pdf, bg = "white")
  }

  imgSet <- readSet(imgSet, "imgSet")
  imgSet$comp_enrichment <- imgNm
  saveSet(imgSet, "imgSet")
  return(imgNm)
}

# ── Phospho-site occupancy via protein reference ──────────────────────────────
#
# Computes per-site phosphorylation occupancy as log2(phosphosite / protein).
# Requires paramSet$protein.ref (log2 protein abundances from proteinGroups.txt
# or a DIA-NN pg_matrix).  Both MaxQuant-phospho and DIA-NN-phospho formats are
# supported as long as the protein reference has been loaded.
#
DetectPhosphoOccupancyBySite <- function(dataName) {
  msgSet   <- readSet(msgSet,   "msgSet")
  dataSet  <- readDataset(dataName)
  paramSet <- readSet(paramSet, "paramSet")

  fail <- function(msg) {
    message("[PhosphoOccupancy] FAIL: ", msg)
    msgSet$current.msg <- msg
    saveSet(msgSet, "msgSet")
    return(0L)
  }

  if (!isTRUE(paramSet$has.protein.ref) || is.null(paramSet$protein.ref))
    return(fail(paste0(
      "No protein reference available. ",
      "Please provide a total proteome file (MaxQuant proteinGroups.txt or DIA-NN pg_matrix.tsv) ",
      "so phosphorylation occupancy can be estimated."
    )))

  prot.ref <- paramSet$protein.ref   # log2 protein abundances: proteins × samples

  phospho.mat <- dataSet$data.norm
  if (is.null(phospho.mat) || nrow(phospho.mat) == 0)
    return(fail("No normalized phospho data found. Run normalization first."))

  # --- conditions ---
  meta.info <- dataSet$meta.info
  disc.inx  <- dataSet$disc.inx
  if (is.null(meta.info) || is.list(meta.info) && !is.data.frame(meta.info)) {
    disc.inx  <- meta.info$disc.inx
    meta.info <- meta.info$meta.info
  }
  if (is.null(disc.inx) || sum(disc.inx) == 0)
    return(fail("No categorical condition found in metadata."))

  cond.col <- names(which(disc.inx))[1]
  groups   <- as.character(meta.info[[cond.col]])
  names(groups) <- rownames(meta.info)

  # align to common samples present in both matrices
  common.smp <- intersect(intersect(colnames(phospho.mat), names(groups)),
                          colnames(prot.ref))
  if (length(common.smp) < 4)
    return(fail(paste0(
      "Only ", length(common.smp),
      " sample(s) overlap between the phospho data and the protein reference. ",
      "Check that the protein reference file uses the same sample names."
    )))

  phospho.mat <- phospho.mat[, common.smp, drop = FALSE]
  prot.ref    <- prot.ref[,    common.smp, drop = FALSE]
  groups      <- groups[common.smp]
  uniq.grp    <- unique(groups)
  if (length(uniq.grp) < 2) return(fail("At least two conditions required."))

  g1.idx <- groups == uniq.grp[1]
  g2.idx <- groups == uniq.grp[2]
  if (sum(g1.idx) < 2 || sum(g2.idx) < 2)
    return(fail("At least 2 replicates per condition required."))

  # --- site-to-protein mapping ---
  site.ids   <- rownames(phospho.mat)
  prot.ids   <- sub("_[STY]_[0-9]+$", "", site.ids)
  residues   <- regmatches(site.ids, regexpr("[STY](?=_[0-9]+$)", site.ids, perl = TRUE))
  residues[lengths(regmatches(site.ids, gregexpr("[STY](?=_[0-9]+$)", site.ids, perl = TRUE))) == 0] <- ""

  # gene name lookup from feature.info
  gene.lookup <- setNames(prot.ids, site.ids)
  if (!is.null(dataSet$feature.info)) {
    fi       <- dataSet$feature.info
    gene.col <- intersect(c("Gene names", "Gene", "GN", "Mapped Genes"), colnames(fi))[1]
    if (!is.na(gene.col)) {
      gv <- fi[[gene.col]]
      gv <- sapply(strsplit(as.character(gv), ";"), function(x) x[1])
      idx <- match(site.ids, rownames(fi))
      valid <- !is.na(idx)
      gene.lookup[site.ids[valid]] <- ifelse(
        is.na(gv[idx[valid]]) | gv[idx[valid]] == "",
        prot.ids[valid], gv[idx[valid]]
      )
    }
  }

  # --- MSstatsPTM-style adjustment ---
  # Fit limma on site matrix and protein matrix separately, then adjust
  # site-level fold changes for protein-level changes (delta-method SE propagation).
  require(limma)

  grp.fac  <- factor(groups)
  design   <- model.matrix(~ 0 + grp.fac)
  colnames(design) <- levels(grp.fac)

  g2.nm <- make.names(uniq.grp[2])
  g1.nm <- make.names(uniq.grp[1])
  colnames(design) <- make.names(colnames(design))
  contrast.vec <- setNames(c(-1, 1), c(g1.nm, g2.nm))
  contrast.mat <- limma::makeContrasts(
    contrasts = paste0("`", g2.nm, "` - `", g1.nm, "`"),
    levels = design
  )

  .fitLimma <- function(mat, des, cmat) {
    keep <- apply(mat, 1, function(x) sum(!is.na(x)) >= 2)
    mat  <- mat[keep, , drop = FALSE]
    fit  <- tryCatch(limma::lmFit(mat, des), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    fit2 <- tryCatch({
      f2 <- limma::contrasts.fit(fit, cmat)
      limma::eBayes(f2)
    }, error = function(e) NULL)
    fit2
  }

  fit.site <- .fitLimma(phospho.mat, design, contrast.mat)
  if (is.null(fit.site)) return(fail("limma fitting failed on site matrix."))

  # restrict protein matrix to proteins that appear in site list
  matched.prot.ids <- unique(prot.ids[prot.ids %in% rownames(prot.ref)])
  if (length(matched.prot.ids) == 0)
    return(fail("No phosphosites matched to the protein reference. Check that protein IDs align."))
  prot.mat.sub <- prot.ref[matched.prot.ids, , drop = FALSE]

  fit.prot <- .fitLimma(prot.mat.sub, design, contrast.mat)
  if (is.null(fit.prot)) return(fail("limma fitting failed on protein reference matrix."))

  # extract per-feature logFC, posterior SE, and total df from each fit
  .extractStats <- function(fit) {
    lfc  <- fit$coefficients[, 1]
    se   <- sqrt(fit$s2.post) * fit$stdev.unscaled[, 1]
    # df.total may be scalar (global eBayes) or vector (robust/trend)
    dft  <- if (length(fit$df.total) == 1) rep(fit$df.total, nrow(fit)) else fit$df.total
    names(dft) <- rownames(fit)
    list(lfc = lfc, se = se, df = dft)
  }

  ss <- .extractStats(fit.site)
  ps <- .extractStats(fit.prot)

  prot.ref.ids <- rownames(prot.ref)

  results.list <- lapply(seq_len(nrow(phospho.mat)), function(i) {
    pid <- prot.ids[i]
    sid <- site.ids[i]
    res <- residues[i]

    # skip sites not present in the limma fit (filtered for low coverage)
    if (!(sid %in% names(ss$lfc))) return(NULL)

    # find matching protein
    pidx <- which(prot.ref.ids == pid)
    if (length(pidx) == 0) pidx <- which(grepl(pid, prot.ref.ids, fixed = TRUE))
    if (length(pidx) == 0) return(NULL)
    pname <- prot.ref.ids[pidx[1]]
    if (!(pname %in% names(ps$lfc))) return(NULL)

    site.lfc <- ss$lfc[sid];  site.se <- ss$se[sid];  site.df <- ss$df[sid]
    prot.lfc <- ps$lfc[pname]; prot.se <- ps$se[pname]; prot.df <- ps$df[pname]
    if (any(is.na(c(site.lfc, site.se, prot.lfc, prot.se)))) return(NULL)

    # MSstatsPTM delta-method correction
    adj.lfc <- site.lfc - prot.lfc
    adj.se  <- sqrt(site.se^2 + prot.se^2)
    t.stat  <- adj.lfc / adj.se
    # Satterthwaite df approximation
    adj.df  <- (site.se^2 + prot.se^2)^2 /
               (site.se^4 / site.df + prot.se^4 / prot.df)
    p.occ   <- 2 * pt(-abs(t.stat), df = adj.df)
    p.site  <- 2 * pt(-abs(site.lfc / site.se), df = site.df)

    # mean log2(phospho/protein) per condition for display
    phos.vec <- as.numeric(phospho.mat[i, ])
    prot.vec <- as.numeric(prot.ref[pidx[1], ])
    occ.vec  <- phos.vec - prot.vec
    occ.g1   <- mean(occ.vec[g1.idx], na.rm = TRUE)
    occ.g2   <- mean(occ.vec[g2.idx], na.rm = TRUE)

    mod.name <- if (nchar(res) > 0) paste0("Phospho(", res, ")") else "Phospho"

    data.frame(
      Gene             = gene.lookup[[sid]],
      Peptide          = sid,
      Modification     = mod.name,
      Mod.Sig          = tolower(paste0("phospho_", res)),
      Precursors.Unmod = 0L,
      Precursors.Mod   = as.integer(sum(!is.na(phos.vec))),
      Occupancy.Cond1  = round(occ.g1, 3),
      Occupancy.Cond2  = round(occ.g2, 3),
      Delta.Occupancy  = round(adj.lfc, 3),
      Occ.Pvalue       = signif(p.occ, 3),
      Total.Pvalue     = signif(p.site, 3),
      stringsAsFactors = FALSE
    )
  })

  results <- do.call(rbind, Filter(Negate(is.null), results.list))
  if (is.null(results) || nrow(results) == 0)
    return(fail(paste0(
      "No phosphosites could be matched to the protein reference. ",
      "Verify that the protein IDs in the reference match those in the phospho data."
    )))

  results$Occ.FDR      <- signif(p.adjust(results$Occ.Pvalue, method = "BH"), 3)
  results              <- results[order(results$Occ.Pvalue), ]
  results$Cond1.Label  <- uniq.grp[1]
  results$Cond2.Label  <- uniq.grp[2]

  # Use UniProt protein IDs (extracted from site IDs) for compartment annotation.
  # This is more reliable than gene symbol lookup because:
  # (a) some organisms (e.g. yeast) use systematic ORF names in the gene DB that
  #     differ from the common names in DIA-NN Genes columns, and
  # (b) the entrez_uniprot table provides a direct accession → Entrez mapping.
  analSet <- readSet(analSet, "analSet")
  result.prot.ids <- sub("_[STY]_[0-9]+$", "", results$Peptide)
  comp.vec <- .mapSitesToCompartment(result.prot.ids, paramSet, analSet)
  comp.vec[is.na(comp.vec)] <- "Unknown"
  results$Parent.Compartment      <- comp.vec
  results$Parent.All.Compartments <- comp.vec
  results$Landscape.Group         <- comp.vec

  message("[PhosphoOccupancy] ", nrow(results), " sites with matched protein reference")
  ov_qs_save(results, "ptm_occupancy_results.qs")
  fast.write(results,  "ptm_occupancy_results.csv")

  msgSet$current.msg <- paste0(
    "Phosphorylation occupancy analysis complete. ", nrow(results), " site(s) analyzed."
  )
  saveSet(msgSet,  "msgSet")
  saveSet(paramSet, "paramSet")
  invisible(RegisterData(dataSet, output = 1L))
  return(nrow(results))
}

PlotPhosphoOccupancyProfile <- function(dataName, imageName, siteId,
                                        format = "png", dpi = 96,
                                        paletteOpt = "default") {
  require(ggplot2)
  require(Cairo)

  dataSet  <- readDataset(dataName)
  paramSet <- readSet(paramSet, "paramSet")

  if (!isTRUE(paramSet$has.protein.ref) || is.null(paramSet$protein.ref))
    return(invisible(0))

  phospho.mat <- dataSet$data.norm
  prot.ref    <- paramSet$protein.ref
  if (is.null(phospho.mat) || !(siteId %in% rownames(phospho.mat)))
    return(invisible(0))

  # condition labels
  meta.info <- dataSet$meta.info
  disc.inx  <- dataSet$disc.inx
  if (is.null(disc.inx)) {
    disc.inx  <- meta.info$disc.inx
    meta.info <- meta.info$meta.info
  }
  cond.col <- names(which(disc.inx))[1]
  groups   <- as.character(meta.info[[cond.col]])
  names(groups) <- rownames(meta.info)

  # find matching protein
  prot.id <- sub("_[STY]_[0-9]+$", "", siteId)
  pidx    <- which(rownames(prot.ref) == prot.id)
  if (length(pidx) == 0) pidx <- which(grepl(prot.id, rownames(prot.ref), fixed = TRUE))
  if (length(pidx) == 0) return(invisible(0))

  common.smp <- intersect(intersect(colnames(phospho.mat), names(groups)), colnames(prot.ref))
  if (length(common.smp) < 2) return(invisible(0))

  phos.vec <- as.numeric(phospho.mat[siteId,   common.smp])
  prot.vec <- as.numeric(prot.ref[pidx[1], common.smp])
  occ.vec  <- phos.vec - prot.vec   # log2 ratio

  df <- data.frame(
    Sample    = common.smp,
    Condition = groups[common.smp],
    Occupancy = occ.vec,
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$Occupancy), ]
  if (nrow(df) == 0) return(invisible(0))

  col <- GetGroupPalette(df$Condition, paletteOpt)
  df$Condition <- factor(df$Condition, levels = names(col))

  trunc.id <- if (nchar(siteId) > 36) paste0(substr(siteId, 1, 36), "…") else siteId

  myplot <- ggplot(df, aes(x = Condition, y = Occupancy, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 1.8, alpha = 0.7, color = "#555555") +
    scale_fill_manual(values = col, guide = "none") +
    theme_bw(base_size = 10) +
    theme(
      axis.title.x       = element_blank(),
      axis.text.x        = element_text(size = 9),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(size = 11, hjust = 0.5),
      plot.subtitle      = element_text(size = 9,  hjust = 0.5, color = "#666666")
    ) +
    ylab("log₂(phosphosite / protein)") +
    ggtitle(trunc.id, subtitle = paste("Protein:", prot.id))

  imgName <- paste0(imageName, "dpi", dpi, ".", format)
  Cairo(file = imgName, width = 5, height = 5, unit = "in", dpi = dpi, bg = "white", type = format)
  invisible(print(myplot))
  invisible(dev.off())
  return(invisible(1))
}
