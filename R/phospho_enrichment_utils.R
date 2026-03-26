##################################################
## R script for ProteoAnalyst
## Description: Phosphoproteomics enrichment analysis utilities
## Handles phosphosite IDs (uniprot+site) for enrichment analysis
##
## Authors:
## Jeff Xia, jeff.xia@mcgill.ca
## Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
###################################################

#' Strip phosphosite IDs to get UniProt part only
#' Handles both ID formats: "P12345_S123" and "Q80XQ2_S_578"
#' @param phosphosite.ids Character vector of phosphosite IDs
#' @return Character vector of UniProt IDs (e.g., "P12345", "Q80XQ2")
.stripPhosphositeToUniprot <- function(phosphosite.ids) {
  if (is.null(phosphosite.ids) || length(phosphosite.ids) == 0) {
    return(character(0))
  }

  # Split on underscore and take first part (UniProt ID with possible isoform)
  uniprot.ids <- sapply(strsplit(as.character(phosphosite.ids), "_"), function(x) x[1])

  # Strip isoform suffixes (e.g., "A2AJT9-3" -> "A2AJT9")
  # This is essential for matching against annotation databases
  uniprot.ids <- sub("-\\d+$", "", uniprot.ids)

  # Preserve names if they exist
  if (!is.null(names(phosphosite.ids))) {
    names(uniprot.ids) <- names(phosphosite.ids)
  }

  return(uniprot.ids)
}

.mapUniprotToEntrezForEnrichment <- function(uniprot.ids, paramSet) {
  org <- paramSet$data.org
  normalized.ids <- trimws(as.character(uniprot.ids))
  db.map <- queryGeneDB("entrez_uniprot", org)

  if (is.null(db.map) || !is.data.frame(db.map) || nrow(db.map) == 0 ||
      !all(c("accession", "gene_id") %in% colnames(db.map))) {
    warning("[Phospho Enrichment] entrez_uniprot mapping table unavailable for ", org)
    return(rep(NA_character_, length(normalized.ids)))
  }

  hit.inx <- match(normalized.ids, as.character(db.map$accession))
  entrez.ids <- as.character(db.map$gene_id[hit.inx])
  names(entrez.ids) <- names(uniprot.ids)
  entrez.ids
}

#' Check if data type is phosphoproteomics
#' @return Logical TRUE if phospho data
.isPhosphoData <- function() {
  paramSet <- readSet(paramSet, "paramSet")
  is_phospho <- (!is.null(paramSet$data.type) && paramSet$data.type == "phospho")
  return(is_phospho)
}

#' Collapse phosphosite IDs to gene-level
#' @param phosphosite.ids Character vector of phosphosite IDs
#' @param paramSet Parameter set containing organism info
#' @return List with entrez IDs and mapping from entrez to original phosphosites
.collapsePhosphositesToGenes <- function(phosphosite.ids, paramSet) {
  org <- paramSet$data.org

  # Strip to UniProt IDs
  uniprot.ids <- .stripPhosphositeToUniprot(phosphosite.ids)

  # Convert UniProt to Entrez using existing mapping
  entrez.ids <- .mapUniprotToEntrezForEnrichment(uniprot.ids, paramSet)

  mapped.inx <- !is.na(entrez.ids) & nzchar(as.character(entrez.ids))
  cat(sprintf("[Phospho Enrichment] UniProt->Entrez mapped %d/%d phosphosite entries (%.1f%%)\n",
              sum(mapped.inx), length(entrez.ids), 100 * sum(mapped.inx) / max(1, length(entrez.ids))))


  # Build gene-level representation from existing mappings
  # Create a named vector: phosphosite -> entrez
  phospho.to.entrez <- entrez.ids
  names(phospho.to.entrez) <- phosphosite.ids

  # Get unique entrez IDs
  unique.entrez <- unique(entrez.ids[mapped.inx])
  cat(sprintf("[Phospho Enrichment] Collapsed to %d unique Entrez IDs\n", length(unique.entrez)))

  # Create reverse mapping: entrez -> list of phosphosites
  entrez.to.phospho <- lapply(unique.entrez, function(eid) {
    phosphosite.ids[which(phospho.to.entrez == eid & !is.na(phospho.to.entrez))]
  })
  names(entrez.to.phospho) <- unique.entrez

  return(list(
    entrez.ids = unique.entrez,
    phospho.to.entrez = phospho.to.entrez,
    entrez.to.phospho = entrez.to.phospho,
    uniprot.ids = uniprot.ids,
    mapped.inx = mapped.inx
  ))
}

#' Wrapper for doEntrez2SymbolMapping that handles phosphosites
#' @param id.vec Character vector - phosphosite IDs if phospho data, else entrez IDs
#' @param data.org Organism
#' @param data.idType ID type
#' @return Character vector of gene symbols
doEntrez2SymbolMappingPhospho <- function(id.vec, data.org = "NA", data.idType = "NA") {
  if (.isPhosphoData()) {
    # For phospho data, strip to UniProt, convert to Entrez, then to Symbol
    paramSet <- readSet(paramSet, "paramSet")
    uniprot.ids <- .stripPhosphositeToUniprot(id.vec)

    # Convert UniProt to Entrez
    entrez.ids <- .mapUniprotToEntrezForEnrichment(uniprot.ids, paramSet)

    # Convert Entrez to Symbol
    symbols <- doEntrez2SymbolMapping(entrez.ids, data.org, "entrez")

    # Preserve original phosphosite IDs as names
    names(symbols) <- id.vec

    return(symbols)
  } else {
    # Regular data - use original function
    return(doEntrez2SymbolMapping(id.vec, data.org, data.idType))
  }
}

#' Wrapper for doEntrezIDAnot that handles phosphosites
#' @param id.vec Character vector - phosphosite IDs if phospho data, else entrez IDs
#' @param data.org Organism
#' @param data.idType ID type
#' @return Data frame with gene_id, symbol, name
doEntrezIDAnotPhospho <- function(id.vec, data.org = "NA", data.idType = "NA") {
  if (.isPhosphoData()) {
    # For phospho data, strip to UniProt, convert to Entrez, then annotate
    paramSet <- readSet(paramSet, "paramSet")
    uniprot.ids <- .stripPhosphositeToUniprot(id.vec)

    # Convert UniProt to Entrez
    entrez.ids <- .mapUniprotToEntrezForEnrichment(uniprot.ids, paramSet)

    # Get annotation for Entrez IDs
    anot.df <- doEntrezIDAnot(entrez.ids, data.org, "entrez")

    # Add original phosphosite IDs as rownames
    rownames(anot.df) <- id.vec

    return(anot.df)
  } else {
    # Regular data - use original function
    return(doEntrezIDAnot(id.vec, data.org, data.idType))
  }
}

#' Perform enrichment analysis for phosphosite data
#' @param dataSet Dataset object
#' @param file.nm Output file name
#' @param fun.type Function library type
#' @param phosphosite.vec Named vector of phosphosite IDs (names are symbols)
#' @param vis.type Visualization type
#' @return Result code
.performEnrichAnalysisPhospho <- function(dataSet, file.nm, fun.type, phosphosite.vec, vis.type) {
  paramSet <- readSet(paramSet, "paramSet")

  cat(sprintf("[Phospho Enrichment] Input: %d phosphosite IDs\n", length(phosphosite.vec)))

  # Collapse phosphosites to genes
  collapse.res <- .collapsePhosphositesToGenes(phosphosite.vec, paramSet)
  entrez.vec <- collapse.res$entrez.ids
  entrez.to.phospho <- collapse.res$entrez.to.phospho
  phospho.to.entrez <- collapse.res$phospho.to.entrez

  # Get gene symbols for entrez IDs
  sym.vec <- doEntrez2SymbolMapping(entrez.vec, paramSet$data.org, "entrez")
  names(entrez.vec) <- sym.vec
  cat(sprintf("[Phospho Enrichment] Unique Entrez IDs entering ORA: %d\n", length(entrez.vec)))


  # Perform enrichment analysis with entrez IDs
  # Temporarily store phosphosite mapping in paramSet
  paramSet$phospho.mapping <- list(
    entrez.to.phospho = entrez.to.phospho,
    phospho.to.entrez = phospho.to.entrez,
    original.phosphosites = phosphosite.vec
  )
  saveSet(paramSet, "paramSet")

  # Call the standard enrichment analysis
  res <- .performEnrichAnalysisPhosphoInternal(dataSet, file.nm, fun.type, entrez.vec, vis.type)

  return(res)
}

#' Modified performEnrichAnalysis that returns phosphosite IDs in fun.anot
#' This is called internally by .performEnrichAnalysisPhospho
.performEnrichAnalysisPhosphoInternal <- function(dataSet, file.nm, fun.type, ora.vec, vis.type) {
  dataSet <<- dataSet

  msgSet <- readSet(msgSet, "msgSet")
  paramSet <- readSet(paramSet, "paramSet")
  require(dplyr)

  # prepare lib
  setres <- .loadEnrichLib(fun.type, paramSet)
  current.featureset <- setres$current.featureset
  current.setids <<- setres$current.setids

  # prepare query
  ora.nms <- names(ora.vec)

  if (is.null(ora.nms)) {
    ora.nms <- ora.vec
    names(ora.vec) <- ora.vec
  }

  # Get phospho mapping if it exists
  phospho.mapping <- paramSet$phospho.mapping
  has.phospho <- !is.null(phospho.mapping)


  if (has.phospho) {
    cat(sprintf("  entrez.to.phospho: class=%s, length=%d\n",
                class(phospho.mapping$entrez.to.phospho), length(phospho.mapping$entrez.to.phospho)))
    cat(sprintf("  Sample entrez IDs: %s\n", paste(head(names(phospho.mapping$entrez.to.phospho), 5), collapse=", ")))
  }

  if (paramSet$universe.opt == "library") {
    # OPTION 1: Use all genes from pathway libraries as universe
    # This is already in Entrez ID space, no conversion needed
    current.universe <- unique(unlist(current.featureset))
    #cat(sprintf("[Phospho Enrichment DEBUG] Using library as universe (Entrez IDs)\n"))
    #cat(sprintf("  Universe size: %d\n", length(current.universe)))
  } else {
    # OPTION 2: Use uploaded features as universe
    # For phospho data, these are phosphosite IDs that need conversion to Entrez
    if (paramSet$anal.type == "onedata") {
      data.anot <- .get.annotated.data()
      current.universe <- rownames(data.anot)
    } else if (paramSet$anal.type == "metadata") {
      inmex <- qs::qread("inmex_meta.qs")
      current.universe <- rownames(inmex$data)
    } else {
      if (!is.null(paramSet$backgroundUniverse)) {
        current.universe <- paramSet$backgroundUniverse
      } else {
        current.universe <- unique(unlist(current.featureset))
      }
    }

    # For phospho data, convert universe from phosphosite IDs to Entrez IDs
    if (has.phospho) {

      # The current.universe contains all phosphosite IDs from the data
      # We need to convert ALL of them to Entrez IDs, not just the significant ones
      universe.phosphosites <- current.universe

      # Strip to UniProt IDs
      universe.uniprots <- .stripPhosphositeToUniprot(universe.phosphosites)

      # Convert UniProt to Entrez
      universe.entrez <- .mapUniprotToEntrezForEnrichment(universe.uniprots, paramSet)

      # Remove NAs and get unique Entrez IDs
      universe.entrez <- unique(universe.entrez[!is.na(universe.entrez)])

      current.universe <- universe.entrez
    }
  }

  # also make sure pathways only contain features measured in experiment
  if (file.exists("data.anot.qs")) {
    current.featureset <- lapply(current.featureset, function(x) {
      x[x %in% current.universe]
    })
    inds <- lapply(current.featureset, length) > 0
    current.featureset <- current.featureset[inds]
  }

  # prepare for the result table
  set.size <- length(current.featureset)
  res.mat <- matrix(0, nrow = set.size, ncol = 5)
  rownames(res.mat) <- names(current.featureset)
  colnames(res.mat) <- c("Total", "Expected", "Hits", "Pval", "FDR")

  q.size <- length(ora.vec)

  # get the matched query for each pathway
  hits.query <- lapply(current.featureset, function(x) {
    ora.nms[ora.vec %in% unlist(x)]
  })

  names(hits.query) <- names(current.featureset)
  hit.num <- unlist(lapply(hits.query, function(x) {
    length(unique(x))
  }), use.names = FALSE)

  gene.vec <- current.universe
  sym.vec <- doEntrez2SymbolMapping(gene.vec, paramSet$data.org, "entrez")
  gene.nms <- sym.vec

  current.featureset.symb <- lapply(current.featureset, function(x) {
    gene.nms[gene.vec %in% unlist(x)]
  })

  # total unique gene number
  uniq.count <- length(current.universe)

  # unique gene count in each pathway
  set.size <- unlist(lapply(current.featureset, length))

  res.mat[, 1] <- set.size
  res.mat[, 2] <- q.size * (set.size / uniq.count)
  res.mat[, 3] <- hit.num

  # use lower.tail = F for P(X>x)
  raw.pvals <- phyper(hit.num - 1, set.size, uniq.count - set.size, q.size, lower.tail = F)
  # Replace NaN values with 1
  raw.pvals[is.nan(raw.pvals)] <- 1

  res.mat[, 4] <- raw.pvals
  res.mat[, 5] <- p.adjust(raw.pvals, "fdr")

  # now, clean up result, synchronize with hit.query
  res.mat <- res.mat[hit.num > 0, , drop = F]
  hits.query <- hits.query[hit.num > 0]

  if (nrow(res.mat) > 1) {
    # order by p value
    ord.inx <- order(res.mat[, 4])
    res.mat <- signif(res.mat[ord.inx, ], 3)
    hits.query <- hits.query[ord.inx]

    res.mat.all <- as.data.frame(res.mat)
    res.mat.all$Pathway <- rownames(res.mat)
    res.mat.all$features <- rep("NA", nrow(res.mat))

    # For phospho data, convert hits.query back to phosphosite IDs
    if (has.phospho) {
      entrez.to.phospho <- phospho.mapping$entrez.to.phospho

      # Convert gene symbols back to phosphosite IDs
      hits.query.phospho <- lapply(hits.query, function(gene.symbols) {
        # gene.symbols are the names, ora.vec values are the entrez IDs
        # Get entrez IDs for these symbols
        entrez.ids <- ora.vec[ora.nms %in% gene.symbols]

        # Map each entrez ID to its phosphosites
        phosphosites <- unlist(lapply(entrez.ids, function(eid) {
          if (as.character(eid) %in% names(entrez.to.phospho)) {
            entrez.to.phospho[[as.character(eid)]]
          } else {
            character(0)
          }
        }))

        return(unique(phosphosites))
      })

      # Keep only phosphosites that exist in the user's input list to avoid missing nodes downstream
      # NOTE: Only apply this filtering for "genelist" visualization type
      # For "volcano" type, the phosphosites come directly from differential analysis, not from list.features
      if (vis.type == "genelist") {
        allowed_sites <- analSet$list.features
        if (!is.null(allowed_sites)) {
          hits.query.phospho <- lapply(hits.query.phospho, function(sites) {
            sites[sites %in% allowed_sites]
          })
        }
      }

      # Replace hits.query with phosphosite version
      hits.query.original <- hits.query  # Keep gene symbol version for features column
      hits.query <- hits.query.phospho   # Use phosphosite version for fun.anot
    } else {
      hits.query.original <- hits.query
    }

    # Iterate through the list and add comma-separated values to the data frame
    for (name in names(hits.query.original)) {
      if (name %in% res.mat.all$Pathway) {
        res.mat.all[which(res.mat.all$Pathway == name), "features"] <-
          paste(hits.query.original[[name]], collapse = ",")
      }
    }

    res.mat.all <- res.mat.all[which(res.mat.all$features != "NA"), ]
    res.mat.all$Pathway <- NULL
    pws <- rownames(res.mat[which(res.mat.all$features != "NA"), ])
    fun.ids2 <- as.vector(setres$current.setids[pws])
    resTable.all <- data.frame(Pathway = pws, ID = fun.ids2, res.mat.all)

    csv.nm <- paste(file.nm, ".csv", sep = "")
    write.csv(resTable.all, file = csv.nm, row.names = F)

    imp.inx <- res.mat[, 4] <= 0.05
    imp.inx[is.na(imp.inx)] <- F
    if (sum(imp.inx) < 10) {
      # too little left, give the top ones
      topn <- ifelse(nrow(res.mat) > 10, 10, nrow(res.mat))
      res.mat <- res.mat[1:topn, ]
      hits.query <- hits.query[1:topn]
    } else {
      res.mat <- res.mat[imp.inx, ]
      hits.query <- hits.query[imp.inx]
    if (sum(imp.inx) > 120) {
        # now, clean up result, synchronize with hit.query
        res.mat <- res.mat[1:120, ]
        hits.query <- hits.query[1:120]
    }
  }
  } else {
    msgSet$current.msg <- "No overlap between queried features and pathway library!"
    return(0)
  }

  # Check for and handle duplicate row names in enr.mat
  if (any(duplicated(rownames(res.mat)))) {
    res.mat <- res.mat[!duplicated(rownames(res.mat)), ]
    hits.query <- hits.query[match(rownames(res.mat), names(hits.query))]
    print("Duplicates in enr.mat were removed.")
  } else {
    res.mat <- res.mat
  }

  resTable <- data.frame(Pathway = rownames(res.mat), res.mat)

  # Persist finalized hits.query for downstream visualizations (ridgeline/network).
  # For phospho data this must be phosphosite IDs, not intermediate gene symbols.
  qs::qsave(hits.query, "hits_query.qs")

  # msg("[PhosphoEnrichAnalysis] Saving enr.mat.qs to: ", getwd())
  # msg("[PhosphoEnrichAnalysis] res.mat dimensions: ", nrow(res.mat), " x ", ncol(res.mat))
  qs:::qsave(res.mat, "enr.mat.qs")

  # Verify file was written
  if (file.exists("enr.mat.qs")) {
    # msg("[PhosphoEnrichAnalysis] Successfully saved enr.mat.qs (", file.size("enr.mat.qs"), " bytes)")
  } else {
    # warning("[PhosphoEnrichAnalysis] WARNING: enr.mat.qs was not created!")
  }

  msgSet$current.msg <- "Functional enrichment analysis was completed"

  # write json
  # fun.anot now contains phosphosite IDs if phospho data, gene symbols otherwise
  fun.anot <- hits.query
  total <- resTable$Total
  if (length(total) == 1) {
    total <- matrix(total)
  }
  fun.pval <- resTable$Pval
  if (length(fun.pval) == 1) {
    fun.pval <- matrix(fun.pval)
  }
  fun.padj <- resTable$FDR
  if (length(fun.padj) == 1) {
    fun.padj <- matrix(fun.padj)
  }
  hit.num <- paste0(resTable$Hits, "/", resTable$Total)
  if (length(hit.num) == 1) {
    hit.num <- matrix(hit.num)
  }
  fun.ids <- as.vector(setres$current.setids[resTable$Pathway])

  resTable$IDs <- fun.ids
  if (length(fun.ids) == 1) {
    fun.ids <- matrix(fun.ids)
  }

  json.res <- list(
    fun.link = setres$current.setlink[1],
    fun.anot = fun.anot,  # Contains phosphosite IDs if phospho data
    fun.ids = fun.ids,
    fun.pval = fun.pval,
    fun.padj = fun.padj,
    hit.num = hit.num,
    total = total
  )
  json.mat <- rjson::toJSON(json.res)
  json.nm <- paste(file.nm, ".json", sep = "")

  sink(json.nm)
  cat(json.mat)
  sink()

  # write csv
  fun.hits <<- hits.query
  fun.pval <<- fun.pval
  hit.num <<- resTable$Hits

  paramSet$partialToBeSaved <- c(paramSet$partialToBeSaved, c(json.nm))

  imgSet <- readSet(imgSet, "imgSet")
  rownames(resTable) <- NULL
  imgSet$enrTables[[vis.type]] <- list()
  imgSet$enrTables[[vis.type]]$table <- resTable
  imgSet$enrTables[[vis.type]]$library <- fun.type
  imgSet$enrTables[[vis.type]]$algo <- "Overrepresentation Analysis"
  imgSet$enrTables[[vis.type]]$current.featureset <- current.featureset

  saveSet(msgSet, "msgSet")
  saveSet(imgSet, "imgSet")
  saveSet(paramSet, "paramSet")

  return(1)
}
