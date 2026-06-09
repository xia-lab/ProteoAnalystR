##################################################
## R script for ProteoAnalyst
## Description: Pathway impact analysis using ORA for proteins/enzymes
## Mirrors MetaboAnalyst CalculateOraScore but for protein (Entrez) inputs
##################################################

.paStandardCompartmentLevels <- function() {
  c(
    "Nucleus",
    "Cell surface & adhesion",
    "Cytoskeleton",
    "Endomembrane",
    "Mitochondria & metabolic organelles",
    "Cytosol",
    "Extracellular",
    "Unknown"
  )
}

.paFormatPathwayCompartment <- function(tab, clear.primary.cutoff = 0.5, clear.margin.cutoff = 0.2) {
  tab <- tab[tab > 0]
  if (length(tab) == 0) {
    return(list(primary = "Unknown", score = 0, display = "Unknown"))
  }

  priority <- .paStandardCompartmentLevels()
  tab <- tab[priority[priority %in% names(tab)]]
  ord <- order(as.numeric(tab), decreasing = TRUE)
  tab <- tab[ord]

  top.count <- as.numeric(tab[1])
  total.count <- sum(as.numeric(tab))
  top.frac <- top.count / total.count
  second.frac <- if (length(tab) > 1) as.numeric(tab[2]) / total.count else 0
  top.name <- names(tab)[1]

  if (length(tab) == 1 || top.frac >= clear.primary.cutoff || (top.frac - second.frac) >= clear.margin.cutoff) {
    return(list(primary = top.name, score = top.frac, display = top.name))
  }

  comp.fracs <- as.numeric(tab) / total.count
  main.comps <- names(tab)[(top.frac - comp.fracs) < clear.margin.cutoff]
  if (length(main.comps) < 2 && length(tab) >= 2) main.comps <- names(tab)[1:2]
  if (length(main.comps) == 2) {
    return(list(primary = paste(main.comps, collapse = " + "), score = top.frac, display = paste(main.comps, collapse = " + ")))
  }

  list(primary = "Mixed compartments", score = top.frac, display = "Mixed compartments")
}

.curatedKeggPathwayCompartment <- function(path.name, kegg.id = "") {
  pname <- tolower(paste(path.name, kegg.id, sep = " "))
  primary <- "Unknown"

  if (grepl(paste0("oxidative phosphorylation|citrate cycle|tca cycle|fatty acid degradation|",
                   "valine, leucine and isoleucine degradation|",
                   "propanoate metabolism|butanoate metabolism|ketone bod|mitophagy|",
                   "mitochondrial"), pname)) {
    primary <- "Mitochondria & metabolic organelles"
  } else if (grepl(paste0("protein processing in endoplasmic reticulum|protein export|",
                          "n-glycan biosynthesis|glycosphingolipid|glycosaminoglycan|",
                          "lysosome|peroxisome|endocytosis|phagosome|autophagy|",
                          "snare|vesicle|endoplasmic reticulum|golgi"), pname)) {
    primary <- "Endomembrane"
  } else if (grepl(paste0("ribosome|proteasome|aminoacyl-trna biosynthesis|glycolysis|",
                          "pentose phosphate|purine metabolism|pyrimidine metabolism|",
                          "amino sugar|glutathione metabolism|rna degradation"), pname)) {
    primary <- "Cytosol"
  } else if (grepl(paste0("dna replication|base excision repair|nucleotide excision repair|",
                          "mismatch repair|homologous recombination|non-homologous|",
                          "rna polymerase|basal transcription|spliceosome|cell cycle|",
                          "ribosome biogenesis|chromatin|nuclear"), pname)) {
    primary <- "Nucleus"
  } else if (grepl(paste0("focal adhesion|ecm-receptor|cell adhesion molecules|adherens junction|",
                          "tight junction|gap junction|neuroactive ligand|cytokine-cytokine|",
                          "receptor interaction|synaptic"), pname)) {
    primary <- "Cell surface & adhesion"
  } else if (grepl("regulation of actin cytoskeleton|actin cytoskeleton|microtubule|axon guidance", pname)) {
    primary <- "Cytoskeleton"
  } else if (grepl("complement and coagulation|extracellular matrix|collagen|secreted", pname)) {
    primary <- "Extracellular"
  }

  list(
    primary = primary,
    all = primary,
    method = if (primary == "Unknown") "inferred" else "curated",
    score = if (primary == "Unknown") 0 else 1,
    mapped = 0L,
    distribution = primary
  )
}

.loadPathwayLocalization <- function(org, lib.path) {
  loc.path <- paste0(lib.path, org, "/", org, "_localization.qs")
  if (!file.exists(loc.path)) loc.path <- paste0(lib.path, org, "/localization.qs")
  if (!file.exists(loc.path)) return(NULL)

  loc <- try(ov_qs_read(loc.path), silent = TRUE)
  if (inherits(loc, "try-error") || is.null(loc)) return(NULL)
  if (!all(c("EntrezID", "Broad.category") %in% colnames(loc))) return(NULL)
  if (!"Main.location" %in% colnames(loc)) loc$Main.location <- "Unknown"

  loc$EntrezID <- as.character(loc$EntrezID)
  loc$Broad.category[is.na(loc$Broad.category) | !nzchar(as.character(loc$Broad.category))] <- "Unknown"
  loc$Main.location[is.na(loc$Main.location) | !nzchar(as.character(loc$Main.location))] <- "Unknown"
  loc$Primary.category <- .paResolveCompartmentAnnotations(loc$Broad.category, loc$Main.location)$primary
  loc
}

.inferKeggPathwayCompartment <- function(path.genes, loc.data) {
  if (is.null(loc.data) || length(path.genes) == 0) {
    return(list(primary = "Unknown", all = "Unknown", method = "inferred",
                score = 0, mapped = 0L, distribution = "Unknown"))
  }

  genes <- unique(as.character(path.genes))
  loc.rows <- loc.data[loc.data$EntrezID %in% genes, , drop = FALSE]
  if (nrow(loc.rows) == 0) {
    return(list(primary = "Unknown", all = "Unknown", method = "inferred",
                score = 0, mapped = 0L, distribution = "Unknown"))
  }

  comps <- .normalizeBroadCategory(loc.rows$Primary.category)
  comps <- comps[!is.na(comps) & nzchar(comps)]
  known.comps <- comps[comps != "Unknown"]
  count.comps <- if (length(known.comps) > 0) known.comps else comps
  if (length(count.comps) == 0) count.comps <- "Unknown"

  tab <- table(factor(count.comps, levels = .paStandardCompartmentLevels()))
  tab <- tab[tab > 0]
  if (length(tab) == 0) {
    primary <- "Unknown"
    score <- 0
    distribution <- "Unknown"
  } else {
    formatted <- .paFormatPathwayCompartment(tab)
    primary <- formatted$primary
    score <- formatted$score
    distribution <- paste(paste0(names(tab), " (", as.integer(tab), ")"), collapse = "; ")
  }

  all.comps <- names(tab)
  if (length(all.comps) == 0) all.comps <- "Unknown"
  list(
    primary = primary,
    all = paste(all.comps, collapse = "; "),
    method = "inferred",
    score = score,
    mapped = length(unique(loc.rows$EntrezID)),
    distribution = distribution
  )
}

.annotateKeggPathwayCompartments <- function(path.names, kegg.ids, featureset, org, lib.path) {
  loc.data <- .loadPathwayLocalization(org, lib.path)
  rows <- lapply(seq_along(path.names), function(i) {
    curated <- .curatedKeggPathwayCompartment(path.names[i], if (length(kegg.ids) >= i) kegg.ids[i] else "")
    if (!is.null(curated$primary) && curated$primary != "Unknown") {
      ann <- curated
      genes <- if (path.names[i] %in% names(featureset)) featureset[[path.names[i]]] else character(0)
      inferred <- .inferKeggPathwayCompartment(genes, loc.data)
      ann$mapped <- inferred$mapped
      ann$distribution <- inferred$distribution
    } else {
      genes <- if (path.names[i] %in% names(featureset)) featureset[[path.names[i]]] else character(0)
      ann <- .inferKeggPathwayCompartment(genes, loc.data)
    }
    data.frame(
      Primary.Compartment = ann$primary,
      All.Compartments = ann$all,
      Compartment.Method = ann$method,
      Compartment.Score = round(as.numeric(ann$score), 3),
      Compartment.Mapped = as.integer(ann$mapped),
      Compartment.Distribution = ann$distribution,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

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

  # Keep pathways that have a KEGG JSON map available for the embedded viewer.
  # Do not restrict to metabolic pathways; signaling and disease maps should
  # remain visible when they are enriched and have pathway JSON data.
  json.dir      <- paste0(paramSet$lib.path, "kegg_pws/json/")
  all.kegg.ids <- as.vector(current.setids[names(set.size)])
  ko.norm.all  <- sub("^[a-z]{2,4}(\\d{5})$", "ko\\1", all.kegg.ids)
  has.json.all <- file.exists(paste0(json.dir, ko.norm.all, ".json"))

  set.size   <- set.size[has.json.all]
  hit.num    <- hit.num[has.json.all]
  raw.pvals  <- raw.pvals[has.json.all]
  holm.pvals <- holm.pvals[has.json.all]
  fdr.pvals  <- fdr.pvals[has.json.all]
  expected   <- expected[has.json.all]
  impact     <- impact[has.json.all]

  if (length(set.size) == 0) {
    msgSet$current.msg <- "No KEGG pathways with available pathway maps found."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Sort by raw p-value and take top 20 (show top results even if none reach p < 0.05)
  ord        <- order(raw.pvals)
  path.names <- names(set.size)[ord]
  kegg.ids   <- as.vector(current.setids[path.names])
  ord        <- ord[seq_len(min(20, length(ord)))]
  path.names <- path.names[seq_len(min(20, length(path.names)))]
  kegg.ids   <- kegg.ids[seq_len(min(20, length(kegg.ids)))]

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
  ora.enzyme.pathway.compartments <<- .annotateKeggPathwayCompartments(
    path.names,
    kegg.ids,
    current.featureset,
    org,
    if (exists("api.lib.path")) api.lib.path else paramSet$lib.path
  )

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
    primary_compartment = ora.enzyme.pathway.compartments$Primary.Compartment,
    all_compartments = ora.enzyme.pathway.compartments$All.Compartments,
    compartment_method = ora.enzyme.pathway.compartments$Compartment.Method,
    compartment_score = ora.enzyme.pathway.compartments$Compartment.Score,
    compartment_mapped = ora.enzyme.pathway.compartments$Compartment.Mapped,
    compartment_distribution = ora.enzyme.pathway.compartments$Compartment.Distribution,
    stringsAsFactors = FALSE
  )
  arrow::write_feather(result.df, "pathway_enzyme_ora.arrow")

  msgSet$current.msg <- paste0("Pathway impact analysis complete. ", nrow(result.df), " enriched pathways found.")
  saveSet(msgSet, "msgSet")
  return(1)
}

CalculateEnzymePathwayGsea <- function(dataName) {
  require(fgsea)
  paramSet <- readSet(paramSet, "paramSet")
  msgSet   <- readSet(msgSet, "msgSet")
  dataSet  <- readDataset(dataName)

  comp.res <- dataSet$comp.res
  if (is.null(comp.res) || nrow(comp.res) == 0) {
    msgSet$current.msg <- "No differential analysis results found. Please run differential analysis first."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Build ranked vector: UniProt IDs → logFC
  lfc <- comp.res[, "logFC"]
  names(lfc) <- rownames(comp.res)
  lfc <- lfc[is.finite(lfc)]
  if (length(lfc) == 0) {
    msgSet$current.msg <- "No finite fold-change values available for GSEA ranking."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Strip phosphosite and isoform suffixes from UniProt IDs
  clean.ids <- sub("_[A-Z]_\\d+$", "", names(lfc))
  clean.ids <- sub("-\\d+$", "", clean.ids)

  # Collapse multiple phosphosites per protein: keep largest absolute logFC
  lfc.df <- data.frame(uniprot = clean.ids, lfc = as.numeric(lfc), stringsAsFactors = FALSE)
  lfc.df <- lfc.df[order(-abs(lfc.df$lfc)), ]
  lfc.df <- lfc.df[!duplicated(lfc.df$uniprot), ]

  org <- paramSet$data.org
  uniprot.map <- queryGeneDB("entrez_uniprot", org)
  hit.inx <- match(lfc.df$uniprot, uniprot.map[, "accession"])
  entrez.ids <- uniprot.map[hit.inx, "gene_id"]
  valid <- !is.na(entrez.ids)
  ranked.vec <- setNames(lfc.df$lfc[valid], entrez.ids[valid])
  # Deduplicate Entrez (keep max absolute)
  ranked.vec <- ranked.vec[order(-abs(ranked.vec))]
  ranked.vec <- ranked.vec[!duplicated(names(ranked.vec))]

  # Build symbol → UniProt map (for network viewer highlighting)
  sym.paired <- doEntrez2SymbolMapping(names(ranked.vec), org, "entrez")
  uniprot.paired <- lfc.df$uniprot[valid][!duplicated(names(ranked.vec))]
  valid.sym <- !is.na(sym.paired) & nchar(trimws(sym.paired)) > 0
  if (any(valid.sym)) {
    sym.uni.raw <- uniprot.paired[valid.sym]
    names(sym.uni.raw) <- sym.paired[valid.sym]
    ora.enzyme.symbol.to.uniprot <<- tapply(sym.uni.raw, names(sym.uni.raw),
                                             function(x) paste(unique(x), collapse = ";"))
  } else {
    ora.enzyme.symbol.to.uniprot <<- character(0)
  }

  setres <- .loadEnrichLib("kegg", paramSet)
  if (!is.list(setres)) {
    msgSet$current.msg <- "Failed to load KEGG pathway library."
    saveSet(msgSet, "msgSet")
    return(0)
  }
  current.featureset <- setres$current.featureset
  current.setids     <- setres$current.setids

  # Run fgsea
  gsea.res <- tryCatch(
    fgsea::fgsea(pathways = current.featureset, stats = ranked.vec,
                 minSize = 3, maxSize = 500, scoreType = "std"),
    error = function(e) { message("[GSEA] fgsea error: ", e$message); NULL }
  )
  if (is.null(gsea.res) || nrow(gsea.res) == 0) {
    msgSet$current.msg <- "GSEA returned no results. Try ORA instead."
    saveSet(msgSet, "msgSet")
    return(0)
  }
  gsea.res <- gsea.res[order(gsea.res$pval), ]

  # Keep only pathways with available KEGG JSON maps
  json.dir    <- paste0(paramSet$lib.path, "kegg_pws/json/")
  path.nms    <- as.character(gsea.res$pathway)
  kegg.ids.all <- as.vector(current.setids[path.nms])
  ko.norm     <- sub("^[a-z]{2,4}(\\d{5})$", "ko\\1", kegg.ids.all)
  has.json    <- file.exists(paste0(json.dir, ko.norm, ".json")) & !is.na(kegg.ids.all)
  gsea.res    <- gsea.res[has.json, ]

  if (nrow(gsea.res) == 0) {
    msgSet$current.msg <- "No KEGG pathways with available pathway maps found."
    saveSet(msgSet, "msgSet")
    return(0)
  }

  # Build hits.query for network highlighting: leading edge → Entrez → UniProt
  le.list <- gsea.res$leadingEdge
  names(le.list) <- gsea.res$pathway
  ora.enzyme.hits.list <<- le.list

  # Top 20 by p-value
  top.n       <- min(20, nrow(gsea.res))
  gsea.top    <- gsea.res[seq_len(top.n), ]
  path.names  <- as.character(gsea.top$pathway)
  kegg.ids    <- as.vector(current.setids[path.names])

  # Result matrix matching ORA column layout:
  # Total, Expected(=0 for GSEA), Hits(=leadingEdge size), RawP, logP, HolmP, FDR, NES
  leading.sizes <- sapply(gsea.top$leadingEdge, length)
  set.sizes     <- sapply(current.featureset[path.names], length)
  nes.score     <- as.numeric(gsea.top$NES)

  res.mat <- cbind(
    set.sizes,
    rep(0, top.n),
    leading.sizes,
    as.numeric(gsea.top$pval),
    -log10(pmax(as.numeric(gsea.top$pval), 1e-300)),
    p.adjust(as.numeric(gsea.top$pval), "holm"),
    as.numeric(gsea.top$padj),
    nes.score
  )
  rownames(res.mat) <- path.names

  ora.enzyme.mat      <<- res.mat
  ora.enzyme.paths    <<- path.names
  ora.enzyme.kegg.ids <<- kegg.ids
  ora.enzyme.pathway.compartments <<- .annotateKeggPathwayCompartments(
    path.names, kegg.ids, current.featureset, org,
    if (exists("api.lib.path")) api.lib.path else paramSet$lib.path
  )

  result.df <- data.frame(
    row_names_id = path.names,
    total        = set.sizes,
    expected     = rep(0, top.n),
    hits         = as.integer(leading.sizes),
    rawP         = as.numeric(gsea.top$pval),
    holmP        = p.adjust(as.numeric(gsea.top$pval), "holm"),
    fdr          = as.numeric(gsea.top$padj),
    impact       = nes.score,
    kegg_id      = kegg.ids,
    stringsAsFactors = FALSE
  )

  msgSet$current.msg <- paste0("Pathway GSEA complete. ", nrow(result.df), " enriched pathways found.")
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

GetEnzymeOra.primaryCompartments <- function() {
  if (!exists("ora.enzyme.pathway.compartments") || is.null(ora.enzyme.pathway.compartments)) {
    return(rep("Unknown", length(ora.enzyme.paths)))
  }
  return(as.character(ora.enzyme.pathway.compartments$Primary.Compartment))
}

GetEnzymeOra.allCompartments <- function() {
  if (!exists("ora.enzyme.pathway.compartments") || is.null(ora.enzyme.pathway.compartments)) {
    return(rep("Unknown", length(ora.enzyme.paths)))
  }
  return(as.character(ora.enzyme.pathway.compartments$All.Compartments))
}

GetEnzymeOra.compartmentMethods <- function() {
  if (!exists("ora.enzyme.pathway.compartments") || is.null(ora.enzyme.pathway.compartments)) {
    return(rep("inferred", length(ora.enzyme.paths)))
  }
  return(as.character(ora.enzyme.pathway.compartments$Compartment.Method))
}

GetEnzymeOra.compartmentScores <- function() {
  if (!exists("ora.enzyme.pathway.compartments") || is.null(ora.enzyme.pathway.compartments)) {
    return(rep(0, length(ora.enzyme.paths)))
  }
  return(as.numeric(ora.enzyme.pathway.compartments$Compartment.Score))
}

GetEnzymeOra.compartmentDistributions <- function() {
  if (!exists("ora.enzyme.pathway.compartments") || is.null(ora.enzyme.pathway.compartments)) {
    return(rep("Unknown", length(ora.enzyme.paths)))
  }
  return(as.character(ora.enzyme.pathway.compartments$Compartment.Distribution))
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

GetEnzymeHitsCompartmentsJson <- function(pathwayQuery) {
  if (!exists("ora.enzyme.hits.list") || length(ora.enzyme.hits.list) == 0) return("{}")

  if (pathwayQuery %in% names(ora.enzyme.hits.list)) {
    entrez.ids <- ora.enzyme.hits.list[[pathwayQuery]]
  } else if (exists("ora.enzyme.kegg.ids") && exists("ora.enzyme.paths") &&
             pathwayQuery %in% ora.enzyme.kegg.ids) {
    path.name <- ora.enzyme.paths[match(pathwayQuery, ora.enzyme.kegg.ids)]
    if (!is.na(path.name) && path.name %in% names(ora.enzyme.hits.list)) {
      entrez.ids <- ora.enzyme.hits.list[[path.name]]
    } else return("{}")
  } else return("{}")

  if (length(entrez.ids) == 0) return("{}")
  paramSet <- readSet(paramSet, "paramSet")
  org <- paramSet$data.org

  sym.vec <- doEntrez2SymbolMapping(as.character(entrez.ids), org, "entrez")
  valid <- !is.na(sym.vec) & nchar(trimws(sym.vec)) > 0
  if (!any(valid)) return("{}")
  sym.vec <- sym.vec[valid]
  entrez.ids <- as.character(entrez.ids)[valid]

  lib.path <- if (exists("api.lib.path")) api.lib.path else paramSet$lib.path
  loc.path <- paste0(lib.path, org, "/", org, "_localization.qs")
  if (!file.exists(loc.path)) loc.path <- paste0(lib.path, org, "/localization.qs")
  if (!file.exists(loc.path)) return("{}")

  loc.data <- try(ov_qs_read(loc.path), silent = TRUE)
  if (inherits(loc.data, "try-error") || is.null(loc.data)) return("{}")
  if (!all(c("EntrezID", "Broad.category") %in% colnames(loc.data))) return("{}")

  loc.entrez <- as.character(loc.data$EntrezID)
  comps <- vapply(seq_along(sym.vec), function(i) {
    eid <- entrez.ids[i]
    if (eid %in% loc.entrez) {
      loc.rows <- loc.data[loc.entrez == eid, , drop = FALSE]
      broad <- paste(unique(as.character(loc.rows$Broad.category)), collapse = "; ")
      main <- if ("Main.location" %in% colnames(loc.rows)) paste(unique(as.character(loc.rows$Main.location)), collapse = "; ") else NA_character_
      .paPrimaryCompartment(broad, main)$primary
    } else "Unknown"
  }, character(1))

  pairs <- paste0('"', gsub('"', '\\\\"', sym.vec), '":"',
                  gsub('"', '\\\\"', comps), '"')
  paste0("{", paste(pairs, collapse = ","), "}")
}

GetEnzymeSymbolUniprotJson <- function() {
  if (!exists("ora.enzyme.symbol.to.uniprot") || length(ora.enzyme.symbol.to.uniprot) == 0) {
    return("{}")
  }
  m <- ora.enzyme.symbol.to.uniprot
  pairs <- paste0('"', gsub('"', '\\\\"', names(m)), '":"', gsub('"', '\\\\"', as.character(m)), '"')
  return(paste0("{", paste(pairs, collapse = ","), "}"))
}
