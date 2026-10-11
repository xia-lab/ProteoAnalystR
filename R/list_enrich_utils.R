##################################################
## R script for ProteoAnalyst
## Description: what a protein LIST is tested as in enrichment analysis.
##
## One approach, the same in every Xia Lab tool (MetaboAnalystR SetRankedListCutoff,
## MicrobiomeAnalystR SetTseaListInput):
##   1. A COMPLETE RANKED LIST (every protein tested, with a score; p and FDR columns optional):
##      (a) its top hits are tested for over-representation against the WHOLE LIST as background;
##      (b) the whole ranking is tested by rank (preranked fgsea).
##   2. A SIGNIFICANT LIST plus the proteins TESTED in the study: over-representation against the
##      tested proteins that are in the library; a set counts only its tested members.
##   3. A SIGNIFICANT LIST alone: over-representation against the library, said as the weaker claim.
## What was tested is recorded in paramSet (enrich.mode, enrich.background, enrich.cutoff) and read
## back, in words, by GetListEnrichInputMsg.
###################################################

# ---- the score of a list ---------------------------------------------------------------------

# A header cell that names a p-value or an adjusted p-value (FDR, q-value).
.paListIsPLabel <- function(label) {
  grepl("^(p|pval|pvalue|p[._ -]?value|fdr|padj|p[._ -]?adj.*|adj[._ -]?p.*|q|qval|qvalue|q[._ -]?value|bh)$",
        tolower(trimws(as.character(label))))
}

# A header cell that names an ADJUSTED p-value (FDR, q-value).
.paListIsFdrLabel <- function(label) {
  grepl("^(fdr|padj|p[._ -]?adj.*|adj[._ -]?p.*|q|qval|qvalue|q[._ -]?value|bh)$",
        tolower(trimws(as.character(label))))
}

# What the user's score column is called, for the sentences ("|logFC| >= 1.21"): its header when it has
# a short one, else "score".
.paListScoreLabel <- function(hdr) {
  if (is.null(hdr) || length(hdr) < 2) return("score")
  lab <- trimws(as.character(hdr[2]))
  if (is.na(lab) || !nzchar(lab) || nchar(lab) > 24) return("score")
  lab
}

# What kind of score a list carries: "p" (small is strong: a p-value or FDR, named so by its header or
# every value in [0, 1]) or "effect" (a fold change or a statistic: large in size is strong).
.paListScoreKind <- function(score, label = "") {
  if (.paListIsPLabel(label)) return("p")
  sc <- suppressWarnings(as.numeric(score))
  sc <- sc[is.finite(sc)]
  if (length(sc) == 0) return(NA_character_)
  if (all(sc >= 0 & sc <= 1)) "p" else "effect"
}

# How strong each protein's score is, larger = stronger, in either direction for an effect.
.paListStrength <- function(score, kind) {
  score <- suppressWarnings(as.numeric(score))
  if (identical(kind, "p")) {
    p <- pmax(score, .Machine$double.xmin)
    -log10(p)
  } else {
    abs(score)
  }
}

# The score, in words: "|logFC|" for an effect (the cutoff is on its size), "P.Value" for a p-like score.
.paListScoreText <- function(label, kind) {
  if (identical(kind, "p")) label else paste0("|", label, "|")
}

.paListN <- function(x) format(x, big.mark = ",", trim = TRUE, scientific = FALSE)

.paListCut <- function(x) format(signif(as.numeric(x), 3), trim = TRUE, scientific = FALSE)

.paListDirText <- function(dir) if (identical(dir, "<=")) "≤" else "≥"

# The smallest value per group, NA when a group has none.
.paMinBy <- function(v, g) {
  if (length(v) == 0) return(stats::setNames(numeric(0), character(0)))
  out <- tapply(v, g, function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE))
  stats::setNames(as.numeric(out), names(out))
}

# The p-value and FDR columns of a list pasted with more than two columns, by protein id: the third
# column is the p-value and the fourth the FDR unless the header names them otherwise. NULL for a
# column that holds no number.
.paListExtraColumns <- function(gene.mat, hdr) {
  ids <- trimws(gene.mat[, 1])
  out <- list(p = NULL, fdr = NULL)
  for (j in 3:min(4, ncol(gene.mat))) {
    role <- if (j == 3) "p" else "fdr"
    if (!is.null(hdr) && length(hdr) >= j) {
      if (.paListIsFdrLabel(hdr[j])) role <- "fdr" else if (.paListIsPLabel(hdr[j])) role <- "p"
    }
    if (!is.null(out[[role]])) next
    v <- suppressWarnings(as.numeric(gene.mat[, j]))
    if (any(is.finite(v))) out[[role]] <- .paMinBy(v, ids)
  }
  out
}

# gene_id (Entrez), score, p and fdr of the proteins MapListIds kept; p and fdr are NA for a list without them.
.paListBuildStats <- function(prot.mat, mapped, extra) {
  ids <- rownames(prot.mat)
  if (is.null(ids) || length(ids) == 0) return(NULL)
  st <- data.frame(gene_id = ids, score = as.numeric(prot.mat[, 1]), p = NA_real_, fdr = NA_real_,
                   stringsAsFactors = FALSE)
  for (col in c("p", "fdr")) {
    v <- extra[[col]]
    if (is.null(v)) next
    by.gene <- .paMinBy(v[as.character(mapped$orig)], as.character(mapped$gene_id))
    st[[col]] <- as.numeric(by.gene[st$gene_id])
  }
  st
}

# ---- reading a list ------------------------------------------------------------------------

# The proteins of one list, as MapListIds left them: gene_id (Entrez), score, p, fdr.
# Returns list(stats, label, kind, name), or NULL when the list is not there.
.paListStats <- function(paramSet, dataName = "") {
  if (is.null(dataName) || length(dataName) != 1 || is.na(dataName) || !nzchar(dataName)) {
    dataName <- if (!is.null(paramSet$selDataNm) && nzchar(paramSet$selDataNm)) paramSet$selDataNm else "datalist1"
  }
  ds <- readDataset(dataName, quiet = TRUE)
  if (is.null(ds)) return(NULL)
  st <- ds$list.stats
  if (is.null(st) && !is.null(ds$prot.mat) && !is.null(rownames(ds$prot.mat))) {
    # a list mapped before the p and FDR columns were kept
    st <- data.frame(gene_id = rownames(ds$prot.mat), score = as.numeric(ds$prot.mat[, 1]),
                     p = NA_real_, fdr = NA_real_, stringsAsFactors = FALSE)
  }
  if (is.null(st) || nrow(st) == 0) return(NULL)
  label <- if (!is.null(ds$score.label) && nzchar(ds$score.label)) ds$score.label else "score"
  list(stats = st, label = label, kind = .paListScoreKind(st$score, label), name = dataName)
}

# The Entrez ids of a list's proteins. A list is kept in the ID space the user gave (MapListIds keeps UniProt
# accessions as they are), the gene-set libraries are in Entrez: UniProt accessions are mapped the way
# .performEnrichAnalysis maps its query (phosphosite and isoform suffixes dropped). NA where unmapped.
.paListEntrez <- function(ids, paramSet) {
  ids <- trimws(as.character(ids))
  non.na <- ids[!is.na(ids) & nzchar(ids)]
  if (length(non.na) > 0 && mean(grepl("^[0-9]+$", non.na)) >= 0.8) return(ids)
  norm <- sub("_[A-Z]_\\d+$", "", ids)
  norm <- sub("-\\d+$", "", norm)
  uniprot.map <- queryGeneDB("entrez_uniprot", paramSet$data.org)
  as.character(uniprot.map[match(norm, uniprot.map[, "accession"]), "gene_id"])
}

# ---- which proteins are the top hits -------------------------------------------------------

# The hits of a complete ranked list. `spec` is list(by, value, dir):
#   by = "top"   the top `value` share (a fraction, 0 < value <= 1) by the size of the score;
#   by = "score" the proteins whose |score| is >= or <= `value` (dir), a p-like score as the p itself;
#   by = "fdr"   the proteins whose FDR column is below `value`;
#   by = "p"     the proteins whose p-value column is below `value`.
# Returns list(ok = FALSE, msg) when nothing can be compared (no column, nobody passes, everybody
# passes), else list(ok = TRUE, hits, by, value, dir, threshold, n.hits, n.list, label, kind).
.paListHits <- function(stats, kind, label, spec, default.frac = 0.10) {
  n <- nrow(stats)
  nn <- .paListN(n)
  by <- if (is.null(spec) || is.null(spec$by)) "top" else spec$by
  value <- if (is.null(spec) || is.null(spec$value)) default.frac else spec$value
  ids <- as.character(stats$gene_id)
  strength <- .paListStrength(stats$score, kind)
  dir <- NA_character_

  if (by == "top") {
    n.hit <- floor(n * value + 1e-9)
    n.hit <- max(1L, min(n.hit, n))
    ord <- order(-strength, seq_len(n), na.last = TRUE)
    pick <- ord[seq_len(n.hit)]
    pass <- rep(FALSE, n); pass[pick] <- TRUE
    if (identical(kind, "p")) {
      dir <- "<="; threshold <- max(abs(stats$score[pick]))
    } else {
      dir <- ">="; threshold <- min(abs(stats$score[pick]))
    }
  } else if (by == "score") {
    dir <- spec$dir
    if (is.null(dir) || length(dir) != 1 || is.na(dir)) dir <- if (identical(kind, "p")) "<=" else ">="
    val <- abs(suppressWarnings(as.numeric(stats$score)))
    pass <- is.finite(val) & (if (dir == ">=") val >= value * (1 - 1e-9) else val <= value * (1 + 1e-9))
    threshold <- value
    if (!any(pass)) {
      return(list(ok = FALSE, msg = paste0("No protein has ", .paListScoreText(label, kind), " ", dir, " ",
                                           .paListCut(value), " among the ", nn, " proteins of your list. ",
                                           "Change the cutoff or its direction.")))
    }
  } else {
    col <- if (by == "fdr") stats$fdr else stats$p
    word <- if (by == "fdr") "FDR" else "p-value"
    if (is.null(col) || !any(is.finite(col))) {
      return(list(ok = FALSE, msg = paste0("Your list has no ", word, " column, so a ", word,
                                           " cutoff cannot be used. Choose the top percentage or a score cutoff instead.")))
    }
    pass <- is.finite(col) & col < value
    threshold <- value
    if (!any(pass)) {
      return(list(ok = FALSE, msg = paste0("No protein of your list has ", if (by == "fdr") "an " else "a ", word,
                                           " below ", .paListCut(value), ". Change the cutoff.")))
    }
  }

  if (all(pass)) {
    return(list(ok = FALSE, msg = paste0("Every one of the ", nn, " proteins of your list passes this cutoff, so there ",
                                         "is nothing to compare the top hits with. Choose a stricter cutoff.")))
  }
  list(ok = TRUE, hits = unique(ids[pass]), by = by, value = value, dir = dir, threshold = threshold,
       n.hits = sum(pass), n.list = n, label = label, kind = kind)
}

# ---- the input -----------------------------------------------------------------------------

.paListRefuse <- function(msg) {
  if (!exists("current.msg", envir = .GlobalEnv, inherits = FALSE)) current.msg <<- character(0)
  AddErrMsg(msg)
  0
}

#' Say what kind of protein list was entered, and how a complete ranked list is cut.
#' @description A protein list is either a set of significant proteins (\code{list.type = "sig"}: tested
#' against the proteins tested in the study when they were given, see \code{ReadUniverseFile}, else
#' against the library) or a COMPLETE ranked list (\code{"ranked"}: every protein tested, with a score and
#' optionally its p-value and FDR). A complete list is its own background: the existing enrichment
#' entry points (\code{PerformDefaultEnrichment}, \code{InitEnrichmentNetwork},
#' \code{PerformHeatmapEnrichment}, \code{PerformGSEA}) then test its top hits against the whole list, and
#' the whole ranking by rank. What was tested is read back, in words, by \code{GetListEnrichInputMsg}.
#' @param dataName The list (default: the first, "datalist1").
#' @param list.type "sig" or "ranked".
#' @param cutoff.by For a ranked list: "top" (the top share of the list by the size of its score),
#' "score" (the proteins whose |score| passes the cutoff in \code{cutoff.dir}), "fdr" or "p" (the proteins
#' below the cutoff in the list's own FDR or p-value column; refused when the list has none).
#' @param cutoff The share for "top" (above 0, at most 1; 0.10 is the top 10 percent), the threshold otherwise.
#' @param cutoff.dir For "score" only: ">=" (|score| at least the cutoff: a fold change or a statistic) or
#' "<=" (at most: a p-value-like score). When not given it follows the kind of score the list carries.
#' @return 1 on success; 0 with the reason in the error message when the choice cannot be used (no score,
#' no such column, a cutoff nobody passes or everybody passes).
#' @export
SetListEnrichInput <- function(dataName = "", list.type = "sig", cutoff.by = "top", cutoff = 0.10, cutoff.dir = NULL) {
  paramSet <- readSet(paramSet, "paramSet")
  if (!identical(paramSet$anal.type, "proteinlist")) {
    return(.paListRefuse("The input type of a list applies to a protein list."))
  }
  if (!(is.character(list.type) && length(list.type) == 1 && list.type %in% c("sig", "ranked"))) {
    return(.paListRefuse("A protein list is either a significant list (sig) or a complete ranked list (ranked)."))
  }
  if (!is.character(dataName) || length(dataName) != 1 || is.na(dataName)) dataName <- ""

  if (list.type == "sig") {
    paramSet$list.enrich <- list(dataName = dataName, list.type = "sig")
    paramSet$enrich.used <- NULL
    saveSet(paramSet, "paramSet")
    return(1)
  }

  if (!is.null(paramSet$numOfLists) && paramSet$numOfLists > 1) {
    return(.paListRefuse("A complete ranked list is one list; paste only one list."))
  }
  if (!(is.character(cutoff.by) && length(cutoff.by) == 1 && cutoff.by %in% c("top", "score", "fdr", "p"))) {
    return(.paListRefuse("The top hits are chosen by the top percentage (top), a score cutoff (score), an FDR (fdr) or a p-value (p)."))
  }
  cutoff <- suppressWarnings(as.numeric(cutoff))
  if (length(cutoff) != 1 || !is.finite(cutoff)) {
    return(.paListRefuse("The cutoff for the top hits must be a number."))
  }
  if (cutoff.by == "top" && (cutoff <= 0 || cutoff > 1)) {
    return(.paListRefuse("The top percentage must be above 0 and at most 100 percent."))
  }
  if (cutoff.by %in% c("fdr", "p") && (cutoff <= 0 || cutoff > 1)) {
    return(.paListRefuse("An FDR or p-value cutoff must be above 0 and at most 1."))
  }
  if (cutoff.by == "score" && cutoff < 0) {
    return(.paListRefuse("The score cutoff is compared with |score|, so it must be zero or larger."))
  }
  if (!is.null(cutoff.dir) && !(length(cutoff.dir) == 1 && !is.na(cutoff.dir) && cutoff.dir %in% c(">=", "<="))) {
    return(.paListRefuse("The direction of a score cutoff is >= (|score| at least the cutoff) or <= (|score| at most the cutoff)."))
  }

  lst <- .paListStats(paramSet, dataName)
  if (is.null(lst)) {
    return(.paListRefuse("The protein list is not there; upload it again."))
  }
  if (!any(is.finite(lst$stats$score) & lst$stats$score != 0)) {
    return(.paListRefuse(paste0("Your list has no score, so it cannot be used as a complete ranked list. Give each protein a score ",
                                "(a fold change or a statistic) in the second column, or use it as a significant list.")))
  }

  dir <- if (cutoff.by == "score") {
    if (is.null(cutoff.dir)) { if (identical(lst$kind, "p")) "<=" else ">=" } else cutoff.dir
  } else NA_character_
  spec <- list(by = cutoff.by, value = cutoff, dir = dir)
  hits <- .paListHits(lst$stats, lst$kind, lst$label, spec)
  if (!isTRUE(hits$ok)) return(.paListRefuse(hits$msg))

  paramSet$list.enrich <- list(dataName = lst$name, list.type = "ranked", cutoff.by = cutoff.by,
                               cutoff = cutoff, cutoff.dir = dir)
  paramSet$enrich.used <- NULL
  saveSet(paramSet, "paramSet")
  1
}

# The saved input for this list, or NULL when the list is not a ranked one.
.paListRankedSpec <- function(paramSet) {
  le <- paramSet$list.enrich
  if (!is.null(le) && identical(le$list.type, "ranked")) le else NULL
}

# ---- what an enrichment run tests ---------------------------------------------------------

# The background and the query of an over-representation run on a protein list.
# `query` holds the Entrez ids the caller tests. Returns list(mode, universe, label, ...) where
# universe is NULL when the library is the background, or list(refused = TRUE, msg) when the saved
# cutoff cannot be applied to the list now.
.paListBackground <- function(paramSet, query) {
  le <- .paListRankedSpec(paramSet)
  if (!is.null(le)) {
    lst <- .paListStats(paramSet, le$dataName)
    if (is.null(lst)) return(list(refused = TRUE, msg = "The protein list is not there; upload it again."))
    spec <- list(by = le$cutoff.by, value = le$cutoff, dir = le$cutoff.dir)
    hits <- .paListHits(lst$stats, lst$kind, lst$label, spec)
    if (!isTRUE(hits$ok)) return(list(refused = TRUE, msg = hits$msg))
    cut <- list(kind = hits$by, value = hits$value, dir = hits$dir, threshold = hits$threshold)
    # a query that is not the top hits (nodes the user selected, say) is said to be a selection
    hit.entrez <- unique(stats::na.omit(.paListEntrez(hits$hits, paramSet)))
    if (!setequal(unique(stats::na.omit(as.character(query))), hit.entrez)) {
      cut <- list(kind = "selection", value = NA_real_, dir = NA_character_, threshold = NA_real_)
    }
    universe <- unique(stats::na.omit(.paListEntrez(lst$stats$gene_id, paramSet)))
    return(list(mode = "ranked", universe = universe,
                label = "All proteins of your complete list", n.list = nrow(lst$stats),
                tested = nrow(lst$stats), score.label = lst$label, score.kind = lst$kind,
                cutoff = cut, n.hits = hits$n.hits))
  }
  uni <- paramSet$backgroundUniverse
  if (!identical(paramSet$universe.opt, "library") && !is.null(uni) && length(uni) > 1) {
    return(list(mode = "ora_measured", universe = unique(stats::na.omit(.paListEntrez(uni, paramSet))),
                label = "Proteins tested in your study (uploaded)", n.list = NA_integer_, tested = length(unique(uni)),
                cutoff = list(kind = "none", value = NA_real_, dir = NA_character_, threshold = NA_real_)))
  }
  list(mode = "ora_catalog", universe = NULL, label = "Gene Set Library", n.list = NA_integer_,
       tested = NA_integer_,
       cutoff = list(kind = "none", value = NA_real_, dir = NA_character_, threshold = NA_real_))
}

# The Entrez ids an over-representation run tests for this list: `ids` (the whole list) for a significant
# list, the top hits for a complete ranked one. NULL, with the reason in the error message, when the saved
# cutoff cannot be applied.
.paListQueryIds <- function(paramSet, ids) {
  le <- .paListRankedSpec(paramSet)
  if (is.null(le)) return(ids)
  lst <- .paListStats(paramSet, le$dataName)
  if (is.null(lst)) { .paListRefuse("The protein list is not there; upload it again."); return(NULL) }
  hits <- .paListHits(lst$stats, lst$kind, lst$label, list(by = le$cutoff.by, value = le$cutoff, dir = le$cutoff.dir))
  if (!isTRUE(hits$ok)) { .paListRefuse(hits$msg); return(NULL) }
  hits$hits
}

# The ranking of a complete ranked list for preranked GSEA: a named vector (Entrez id -> score) with the
# strongest up first, and what the run will record. A p-like score ranks as -log10(p). NULL when the list
# is not a complete ranked one.
.paListRankedVec <- function(paramSet) {
  le <- .paListRankedSpec(paramSet)
  if (is.null(le)) return(NULL)
  lst <- .paListStats(paramSet, le$dataName)
  if (is.null(lst)) return(NULL)
  st <- lst$stats
  score <- if (identical(lst$kind, "p")) .paListStrength(st$score, "p") else as.numeric(st$score)
  vec <- stats::setNames(score, as.character(st$gene_id))
  vec <- vec[is.finite(vec)]
  list(vec = vec[order(-vec)],
       bg = list(mode = "ranked", label = "All proteins of your complete list", n.list = nrow(st),
                 tested = nrow(st), score.label = lst$label, score.kind = lst$kind,
                 cutoff = list(kind = "none", value = NA_real_, dir = NA_character_, threshold = NA_real_)))
}

# Record what a run tested, in paramSet, so every surface can say it.
.paListRecordRun <- function(paramSet, bg, method, query.input, query.size, bg.size, sets.tested = NA_integer_) {
  paramSet$enrich.mode <- bg$mode
  paramSet$enrich.background <- list(label = bg$label, size = bg.size)
  paramSet$enrich.cutoff <- list(kind = bg$cutoff$kind, value = bg$cutoff$value, dir = bg$cutoff$dir)
  paramSet$enrich.used <- list(mode = bg$mode, method = method, label = bg$label,
                               n.list = bg$n.list, tested = bg$tested, n.hits = bg$n.hits,
                               query.input = query.input, query.size = query.size, bg.size = bg.size,
                               sets.tested = sets.tested, cutoff = bg$cutoff,
                               score.label = bg$score.label, score.kind = bg$score.kind,
                               library = paramSet$init.lib)
  paramSet$universe.opt.readable <- bg$label
  paramSet
}

# ---- saying it -----------------------------------------------------------------------------

# One plain sentence from a record of what was (or is about to be) tested.
.paListSentence <- function(u) {
  if (is.null(u) || is.null(u$mode)) return("")
  lab <- if (is.null(u$score.label)) "score" else u$score.label
  kind <- if (is.null(u$score.kind)) "effect" else u$score.kind
  stext <- .paListScoreText(lab, kind)
  inlib <- function(n.bg, n.all) {
    if (is.null(n.bg) || is.na(n.bg) || is.null(n.all) || is.na(n.all) || n.bg >= n.all) ""
    else paste0(" (", .paListN(n.bg), " of them in the gene-set library)")
  }

  if (identical(u$mode, "ranked") && identical(u$method, "gsea")) {
    return(paste0("All ", .paListN(u$n.list), " proteins of your complete list tested by rank (preranked GSEA on ", lab, ")"))
  }
  if (identical(u$mode, "ranked")) {
    cut <- u$cutoff
    n <- .paListN(u$n.list)
    tail <- paste0(" tested against all ", n, " as background", inlib(u$bg.size, u$n.list))
    if (identical(cut$kind, "top")) {
      thr <- if (is.null(cut$threshold) || is.na(cut$threshold)) "" else
        paste0(", ", stext, " ", .paListDirText(cut$dir), " ", .paListCut(cut$threshold))
      return(paste0("Top ", format(round(cut$value * 100, 2), trim = TRUE), "% of the ", n,
                    " proteins in your complete list (", .paListN(u$n.hits), thr, ")", tail))
    }
    if (identical(cut$kind, "score")) {
      return(paste0("The ", .paListN(u$n.hits), " proteins with ", stext, " ", .paListDirText(cut$dir), " ",
                    .paListCut(cut$value), " of the ", n, " in your complete list", tail))
    }
    if (cut$kind %in% c("fdr", "p")) {
      return(paste0("The ", .paListN(u$n.hits), " proteins with ", if (cut$kind == "fdr") "FDR" else "p-value", " below ",
                    .paListCut(cut$value), " of the ", n, " in your complete list", tail))
    }
    return(paste0("The ", .paListN(u$query.input), " proteins you selected from your complete list of ", n, tail))
  }
  if (identical(u$mode, "ora_measured")) {
    q <- if (!is.null(u$query.input) && !is.na(u$query.input)) paste0(.paListN(u$query.input), " significant proteins") else "significant proteins"
    left <- if (!is.null(u$query.input) && !is.null(u$query.size) && !is.na(u$query.input) && !is.na(u$query.size) &&
                u$query.size < u$query.input)
      paste0(" (", .paListN(u$query.input - u$query.size), " of them not among the tested proteins, left out)") else ""
    return(paste0("Your ", q, left, " tested against the ", .paListN(u$tested), " proteins tested in your study as background",
                  inlib(u$bg.size, u$tested)))
  }
  q <- if (!is.null(u$query.input) && !is.na(u$query.input)) paste0(.paListN(u$query.input), " significant proteins") else "significant proteins"
  size <- if (!is.null(u$bg.size) && !is.na(u$bg.size)) paste0(.paListN(u$bg.size), " ") else ""
  paste0("Your ", q, " tested against all ", size, "proteins in the gene-set library as background, ",
         "a weaker claim than testing against the proteins tested in your study")
}

#' What the last enrichment on a protein list tested, in words
#' @description One plain sentence, for example "Top 10% of the 1,532 proteins in your complete list
#' (153, |logFC| >= 1.21) tested against all 1,532 as background", the same sentence for the proteins
#' tested in the study, or for the library as background (said as the weaker claim). Before any
#' enrichment has run it says what the saved input will test. Empty when the analysis is not a protein
#' list. There is no final full stop; callers add their own.
#' @param dataName The list (default: the first).
#' @return A character string.
#' @export
GetListEnrichInputMsg <- function(dataName = "") {
  paramSet <- readSet(paramSet, "paramSet")
  if (!identical(paramSet$anal.type, "proteinlist")) return("")
  if (!is.null(paramSet$enrich.used)) return(.paListSentence(paramSet$enrich.used))

  # nothing has run yet: what the saved input will test
  le <- .paListRankedSpec(paramSet)
  if (!is.null(le)) {
    lst <- .paListStats(paramSet, le$dataName)
    if (is.null(lst)) return("")
    hits <- .paListHits(lst$stats, lst$kind, lst$label, list(by = le$cutoff.by, value = le$cutoff, dir = le$cutoff.dir))
    if (!isTRUE(hits$ok)) return("")
    return(.paListSentence(list(mode = "ranked", method = "ora", n.list = nrow(lst$stats), n.hits = hits$n.hits,
                                score.label = lst$label, score.kind = lst$kind,
                                cutoff = list(kind = hits$by, value = hits$value, dir = hits$dir, threshold = hits$threshold))))
  }
  lst <- .paListStats(paramSet, dataName)
  if (is.null(lst)) return("")
  bg <- .paListBackground(paramSet, character(0))
  .paListSentence(list(mode = bg$mode, method = "ora", query.input = nrow(lst$stats), tested = bg$tested))
}
