##################################################
## R script for ProteoAnalyst
## Description: Pathway impact analysis using ORA for proteins/enzymes
## Mirrors MetaboAnalyst CalculateOraScore but for protein (Entrez) inputs
##################################################

CalculateEnzymePathwayOra <- function(dataName, topoCode = "rbc") {
  paramSet <- readSet(paramSet, "paramSet")
  msgSet   <- readSet(msgSet, "msgSet")
  dataSet  <- readDataset(dataName)

  # For protein-list uploads sig.mat is absent — use the full uploaded list as the hit set
  if (is.null(dataSet$sig.mat) || nrow(dataSet$sig.mat) == 0) {
    if (is.null(dataSet$norm.mat) || nrow(dataSet$norm.mat) == 0) {
      msgSet$current.msg <- "No proteins found. Please upload data first."
      saveSet(msgSet, "msgSet")
      return(0)
    }
    sig.vec <- rownames(dataSet$norm.mat)
  } else {
    # Get significant protein IDs (UniProt IDs as row names of sig.mat)
    sig.vec <- rownames(dataSet$sig.mat)
  }
  sig.vec <- sub("_[A-Z]_\\d+$", "", sig.vec)   # strip phosphosite annotations
  sig.vec <- sub("-\\d+$", "", sig.vec)           # strip isoform suffixes
  sig.vec <- unique(trimws(sig.vec))

  # Convert UniProt → Entrez IDs
  org <- paramSet$data.org
  uniprot.map <- queryGeneDB("entrez_uniprot", org)
  hit.inx       <- match(sig.vec, uniprot.map[, "accession"])
  mapped.entrez <- uniprot.map[hit.inx, "gene_id"]  # parallel to sig.vec
  valid.inx     <- !is.na(mapped.entrez)
  entrez.paired <- mapped.entrez[valid.inx]          # Entrez IDs with valid mapping
  uniprot.paired <- sig.vec[valid.inx]               # corresponding UniProt accessions
  sig.entrez    <- unique(entrez.paired)

  if (length(sig.entrez) == 0) {
    msgSet$current.msg <- "No significant proteins could be mapped to KEGG gene IDs."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Build symbol → UniProt map (for network viewer highlighting)
  sym.paired <- doEntrez2SymbolMapping(entrez.paired, org, "entrez")
  valid.sym  <- !is.na(sym.paired) & nchar(trimws(sym.paired)) > 0
  if (any(valid.sym)) {
    sym.uni.raw <- uniprot.paired[valid.sym]
    names(sym.uni.raw) <- sym.paired[valid.sym]
    # For genes with multiple UniProt accessions, collapse with semicolon
    ora.enzyme.symbol.to.uniprot <<- tapply(sym.uni.raw, names(sym.uni.raw),
                                             function(x) paste(unique(x), collapse = ";"))
  } else {
    ora.enzyme.symbol.to.uniprot <<- character(0)
  }

  # Load KEGG gene sets for organism
  setres <- .loadEnrichLib("kegg", paramSet)
  if (!is.list(setres)) {
    msgSet$current.msg <- "Failed to load KEGG pathway library."
    saveSet(msgSet, "msgSet")
    return(0)
  }
  current.featureset <- setres$current.featureset
  current.setids     <- setres$current.setids    # named vector: pathway name → KEGG pathway ID

  # Determine background universe
  universeOpt <- if (!is.null(paramSet$universeOpt)) paramSet$universeOpt else "library"
  if (universeOpt == "uploaded") {
    all.vec    <- rownames(dataSet$norm.mat)
    all.vec    <- sub("_[A-Z]_\\d+$", "", all.vec)
    all.vec    <- sub("-\\d+$", "", all.vec)
    all.vec    <- unique(trimws(all.vec))
    hit.inx2   <- match(all.vec, uniprot.map[, "accession"])
    all.entrez <- unique(uniprot.map[hit.inx2, "gene_id"])
    all.entrez <- all.entrez[!is.na(all.entrez)]
  } else {
    all.entrez <- unique(unlist(current.featureset))
  }

  # Filter feature sets to universe members with >= 2 genes
  current.featureset <- lapply(current.featureset, function(x) x[x %in% all.entrez])
  keep.inx           <- sapply(current.featureset, length) >= 2
  current.featureset <- current.featureset[keep.inx]

  if (length(current.featureset) == 0) {
    msgSet$current.msg <- "No KEGG pathways contain enough mapped proteins."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  q.size     <- length(sig.entrez)
  uniq.count <- length(all.entrez)

  set.size   <- sapply(current.featureset, length)
  hits.query <- lapply(current.featureset, function(x) sig.entrez[sig.entrez %in% x])
  hit.num    <- sapply(hits.query, length)

  # Keep only pathways with at least one hit
  keep          <- hit.num > 0
  set.size      <- set.size[keep]
  hit.num       <- hit.num[keep]
  hits.query    <- hits.query[keep]
  ora.enzyme.hits.list <<- hits.query   # pathway name → Entrez IDs vector

  if (length(hit.num) == 0) {
    msgSet$current.msg <- "No significant pathways found."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Hypergeometric test
  raw.pvals            <- phyper(hit.num - 1, set.size, uniq.count - set.size, q.size, lower.tail = FALSE)
  raw.pvals[is.nan(raw.pvals)] <- 1
  holm.pvals           <- p.adjust(raw.pvals, "holm")
  fdr.pvals            <- p.adjust(raw.pvals, "fdr")
  expected             <- q.size * (set.size / uniq.count)

  # Impact score: fraction of pathway covered by significant proteins
  impact <- hit.num / set.size

  # Restrict to metabolic pathways only (reactions > 0 in KEGG JSON) immediately
  # after computation, before any sorting or cutoffs.
  json.dir           <- paste0(paramSet$lib.path, "kegg_pws/json/")
  metabolic.ids.file <- paste0(json.dir, "metabolic_ko_ids.txt")
  metabolic.ko.ids   <- if (file.exists(metabolic.ids.file)) {
    trimws(readLines(metabolic.ids.file, warn = FALSE))
  } else {
    character(0)
  }
  all.kegg.ids <- as.vector(current.setids[names(set.size)])
  ko.norm.all  <- sub("^[a-z]{2,4}(\\d{5})$", "ko\\1", all.kegg.ids)
  has.json.all <- file.exists(paste0(json.dir, ko.norm.all, ".json"))
  is.metabolic <- (length(metabolic.ko.ids) == 0 | ko.norm.all %in% metabolic.ko.ids) & has.json.all

  set.size   <- set.size[is.metabolic]
  hit.num    <- hit.num[is.metabolic]
  raw.pvals  <- raw.pvals[is.metabolic]
  holm.pvals <- holm.pvals[is.metabolic]
  fdr.pvals  <- fdr.pvals[is.metabolic]
  expected   <- expected[is.metabolic]
  impact     <- impact[is.metabolic]

  if (length(set.size) == 0) {
    msgSet$current.msg <- "No metabolic pathways found with mapped proteins."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Sort by raw p-value and take top 10 (show top results even if none reach p < 0.05)
  ord        <- order(raw.pvals)
  path.names <- names(set.size)[ord]
  kegg.ids   <- as.vector(current.setids[path.names])
  ord        <- ord[seq_len(min(10, length(ord)))]
  path.names <- path.names[seq_len(min(10, length(path.names)))]
  kegg.ids   <- kegg.ids[seq_len(min(10, length(kegg.ids)))]

  # Result matrix: Total, Expected, Hits, RawP, logP, HolmP, FDR, Impact
  res.mat <- cbind(
    set.size[ord],
    round(expected[ord], 3),
    hit.num[ord],
    raw.pvals[ord],
    -log10(raw.pvals[ord]),
    holm.pvals[ord],
    fdr.pvals[ord],
    impact[ord]
  )
  rownames(res.mat) <- path.names

  ora.enzyme.mat     <<- res.mat
  ora.enzyme.paths   <<- path.names
  ora.enzyme.kegg.ids <<- kegg.ids

  # Persist to Arrow for table display
  result.df <- data.frame(
    row_names_id = path.names,
    total        = set.size[ord],
    expected     = round(expected[ord], 3),
    hits         = as.integer(hit.num[ord]),
    rawP         = raw.pvals[ord],
    holmP        = holm.pvals[ord],
    fdr          = fdr.pvals[ord],
    impact       = impact[ord],
    stringsAsFactors = FALSE
  )
  arrow::write_feather(result.df, "pathway_enzyme_ora.arrow")

  msgSet$current.msg <- paste0("Pathway impact analysis complete. ", nrow(result.df), " enriched pathways found.")
  saveSet(msgSet, "msgSet")
  return(1)
}

GetEnzymeOra.pathNames <- function() {
  return(ora.enzyme.paths)
}

GetEnzymeOra.keggIDs <- function() {
  return(ora.enzyme.kegg.ids)
}

GetEnzymeOra.mat <- function() {
  return(ora.enzyme.mat)
}

GetEnzymeHitsSymbols <- function(pathwayQuery) {
  if (!exists("ora.enzyme.hits.list") || length(ora.enzyme.hits.list) == 0) return(character(0))

  if (pathwayQuery %in% names(ora.enzyme.hits.list)) {
    entrez.ids <- ora.enzyme.hits.list[[pathwayQuery]]
  } else if (exists("ora.enzyme.kegg.ids") && exists("ora.enzyme.paths") &&
             pathwayQuery %in% ora.enzyme.kegg.ids) {
    path.name <- ora.enzyme.paths[match(pathwayQuery, ora.enzyme.kegg.ids)]
    if (!is.na(path.name) && path.name %in% names(ora.enzyme.hits.list)) {
      entrez.ids <- ora.enzyme.hits.list[[path.name]]
    } else {
      return(character(0))
    }
  } else {
    return(character(0))
  }

  if (length(entrez.ids) == 0) return(character(0))
  paramSet <- readSet(paramSet, "paramSet")
  sym.vec <- doEntrez2SymbolMapping(entrez.ids, paramSet$data.org, "entrez")
  return(unique(sym.vec[!is.na(sym.vec) & nchar(trimws(sym.vec)) > 0]))
}

GetEnzymeSymbolUniprotJson <- function() {
  if (!exists("ora.enzyme.symbol.to.uniprot") || length(ora.enzyme.symbol.to.uniprot) == 0) {
    return("{}")
  }
  m <- ora.enzyme.symbol.to.uniprot
  pairs <- paste0('"', gsub('"', '\\\\"', names(m)), '":"', gsub('"', '\\\\"', as.character(m)), '"')
  return(paste0("{", paste(pairs, collapse = ","), "}"))
}
