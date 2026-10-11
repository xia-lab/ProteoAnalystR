# Protein-LIST enrichment: what a list is tested as (SetListEnrichInput / GetListEnrichInputMsg and the
# existing ORA / GSEA entry points that honour it).
#
# The approach is the same in every Xia Lab tool:
#   1. a COMPLETE ranked list: its top hits (a top share, a score cutoff with an explicit direction, an FDR or
#      p-value cutoff) against the WHOLE list as background, and the whole ranking by rank (preranked fgsea);
#   2. a SIGNIFICANT list plus the proteins TESTED: against the tested proteins that are in the library;
#   3. a SIGNIFICANT list alone: against the library, said as the weaker claim.
# The tests run the package's real storage (readSet / saveSet / readDataset) and the real enrichment code in a
# scratch directory; only the reference databases (gene-set library, UniProt -> Entrez, Entrez -> symbol) are
# small stand-ins, so every expected number below can be worked out by hand.

R_DIR <- normalizePath(file.path(testthat::test_path(), "..", "..", "R"), mustWork = FALSE)
R_FILES <- c("misc_utils.R", "data_utils_general.R", "list_enrich_utils.R", "enrich_utils.R",
             "data_maplist_utils.R", "utils_gsea.R", "report_summary.R", "enrich_path_impact.R")

# ---- the world of the tests -------------------------------------------------------------------
# Proteins P01..P60 are Entrez genes 1..60. The gene-set library:
#   A  1:10                  A2  1-4 and 11-16        B  11:20      C  21:30      H  31:40
#   D  51:60 (none of them is in the lists below)     E  1,2 (two members)        F  3 (one member)
#   G  31:40 and 51:60 (half of it is outside a list of P01..P40)
LIB <- list(A = as.character(1:10), A2 = as.character(c(1:4, 11:16)), B = as.character(11:20),
            C = as.character(21:30), H = as.character(31:40), D = as.character(51:60),
            E = c("1", "2"), F = "3", G = as.character(c(31:40, 51:60)))
LIB2 <- list(X = as.character(1:6), Y = as.character(5:12), Z = as.character(c(1, 2, 20:25)))

.load_list_env <- function() {
  e <- new.env()
  for (f in R_FILES) {
    skip_if_not(file.exists(file.path(R_DIR, f)), paste(f, "not found"))
    sys.source(file.path(R_DIR, f), envir = e)
  }
  # state the package keeps in globals: give it a home in this environment, not the user's session
  e$paramSet <- list(objName = "paramSet"); e$msgSet <- list(objName = "msgSet")
  e$cmdSet <- list(objName = "cmdSet"); e$current.msg <- character(0)
  e$dataSet <- list(annotated = FALSE)     # a global in a session; PerformDefaultEnrichment hands it on unread for a list
  # what the pathway-impact analysis leaves in globals for its viewer
  for (g in c("ora.enzyme.symbol.to.uniprot", "ora.enzyme.hits.list", "ora.enzyme.mat", "ora.enzyme.paths",
              "ora.enzyme.kegg.ids", "ora.enzyme.pathway.compartments")) assign(g, NULL, envir = e)
  e$lib <- LIB
  e$queryGeneDB <- function(db.nm, org) {
    ids <- 1:60
    if (identical(db.nm, "entrez_uniprot")) data.frame(accession = sprintf("P%02d", ids), gene_id = as.character(ids), stringsAsFactors = FALSE)
    else data.frame(gene_id = as.character(ids), symbol = sprintf("G%02d", ids), stringsAsFactors = FALSE)
  }
  e$doEntrez2SymbolMapping <- function(entrez.vec, data.org = "NA", data.idType = "NA") ifelse(is.na(entrez.vec), NA_character_, sprintf("G%02d", suppressWarnings(as.integer(as.character(entrez.vec)))))
  e$.loadEnrichLib <- function(fun.type, paramSet) {
    lib <- get("lib", envir = e)
    list(current.featureset = lib, current.setids = stats::setNames(sprintf("hsa%05d", seq_along(lib)), names(lib)), current.setlink = "http://x/")
  }
  # UniProt accessions are kept as they are by MapListIds
  e$.doGeneIDMapping <- function(feature.vec, idType, paramSet, outputType = "vec", keepNA = FALSE)
    data.frame(accession = feature.vec, gene_id = feature.vec, orig = feature.vec, unmapped = FALSE, stringsAsFactors = FALSE)
  e$.addUniprotColumn <- function(symbol.map, idType, feature.vec, mapped.ids, paramSet) symbol.map
  e$saveDataQs <- function(data, name, module.nm, dataName) invisible(NULL)
  e$getEntrezTableName <- function(data.org, data.idType) "entrez"
  e$.annotateKeggPathwayCompartments <- function(path.names, kegg.ids, featureset, org, lib.path)
    data.frame(Primary.Compartment = rep("Cytosol", length(path.names)), All.Compartments = "Cytosol",
               Compartment.Method = "test", Compartment.Score = 1, Compartment.Mapped = 1, Compartment.Distribution = "",
               stringsAsFactors = FALSE)
  e$RecordRCommand <- function(...) invisible(NULL)
  e
}

# a fresh analysis folder with a protein-list session in it
.new_session <- function(e) {
  dir <- tempfile("list_enrich_"); dir.create(dir)
  withr::local_dir(dir, .local_envir = parent.frame())
  e$paramSet <- list(objName = "paramSet", anal.type = "proteinlist", data.org = "hsa", data.idType = "uniprot",
                     universe.opt = "uploaded", universe.opt.readable = "Uploaded data", numOfLists = 1,
                     on.public.web = TRUE, init.lib = "kegg", selDataNm = "datalist1", mdata.all = list(),
                     partialToBeSaved = character(0), lib.path = paste0(dir, "/lib/"), api.bool = FALSE)
  e$msgSet <- list(objName = "msgSet", current.msg = character(0))
  e$current.msg <- character(0)
  # a pathway map for every set of the library, which the pathway-impact analysis requires
  dir.create(file.path(dir, "lib", "kegg_pws", "json"), recursive = TRUE)
  for (i in seq_along(LIB)) writeLines("{}", file.path(dir, "lib", "kegg_pws", "json", sprintf("ko%05d.json", i)))
  e$saveSet(e$paramSet, "paramSet"); e$saveSet(e$msgSet, "msgSet")
  e$saveSet(list(objName = "analSet"), "analSet"); e$saveSet(list(objName = "cmdSet"), "cmdSet")
  e$saveSet(list(objName = "imgSet", enrTables = list(), featureList = list()), "imgSet")
  e$paramSet <- NULL
  invisible(dir)
}

# a mapped list as MapListIds leaves it: proteins P01..Pn, with a score (and p, fdr)
.save_list <- function(e, n = 40, score = NULL, p = NULL, fdr = NULL, label = "logFC", name = "datalist1", sig.ids = NULL) {
  ids <- sprintf("P%02d", seq_len(n))
  if (is.null(score)) {
    # P01..P04 stand out, then a smooth fall; signs alternate so that both directions occur
    score <- c(3.0, 2.8, 2.6, 2.4, seq(1.2, 0.05, length.out = n - 4)) * rep(c(1, -1), length.out = n)
    score[1:4] <- abs(score[1:4])
  }
  st <- data.frame(gene_id = ids, score = score, p = if (is.null(p)) NA_real_ else p,
                   fdr = if (is.null(fdr)) NA_real_ else fdr, stringsAsFactors = FALSE)
  if (is.null(sig.ids)) sig.ids <- ids
  ds <- list(name = name, list.stats = st, score.label = label, prot.mat = matrix(score, ncol = 1, dimnames = list(ids, "x")),
             sig.mat = matrix(1, ncol = 1, nrow = length(sig.ids), dimnames = list(sig.ids, "x")))
  e$RegisterData(ds)
  ps <- e$readSet(paramSet, "paramSet"); ps$all.ent.mat <- ds$prot.mat; e$saveSet(ps, "paramSet")
  st
}

.msgs <- function(e) e$current.msg
.ps <- function(e) e$readSet(paramSet, "paramSet")

# ---- which proteins are the top hits ------------------------------------------------------------

test_that("top share: the strongest |score| in either direction, with the realised threshold", {
  e <- .load_list_env()
  st <- data.frame(gene_id = sprintf("P%02d", 1:10), score = c(-5, 4, 3, -2, 1, 0.5, 0.4, 0.3, 0.2, 0.1), p = NA_real_, fdr = NA_real_)
  h <- e$.paListHits(st, "effect", "logFC", list(by = "top", value = 0.2))
  expect_true(h$ok)
  expect_equal(h$hits, c("P01", "P02"))
  expect_equal(h$n.hits, 2L); expect_equal(h$dir, ">="); expect_equal(h$threshold, 4)
  # a share that rounds to no protein still tests one
  expect_equal(e$.paListHits(st, "effect", "logFC", list(by = "top", value = 0.01))$n.hits, 1L)
})

test_that("direction is defaulted by the kind of score, and changeable", {
  e <- .load_list_env()
  eff <- data.frame(gene_id = sprintf("P%02d", 1:6), score = c(2, -1.5, 1, 0.5, 0.2, -0.1), p = NA_real_, fdr = NA_real_)
  pl <- data.frame(gene_id = sprintf("P%02d", 1:6), score = c(0.001, 0.01, 0.04, 0.2, 0.6, 0.9), p = NA_real_, fdr = NA_real_)
  expect_equal(e$.paListScoreKind(eff$score, "logFC"), "effect")
  expect_equal(e$.paListScoreKind(pl$score, "score"), "p")
  expect_equal(e$.paListScoreKind(eff$score, "adj.P.Val"), "p")   # a header that names a p-value decides
  # default: |score| >= x for an effect, score <= x for a p-like score
  expect_equal(e$.paListHits(eff, "effect", "logFC", list(by = "score", value = 1))$hits, c("P01", "P02", "P03"))
  expect_equal(e$.paListHits(pl, "p", "score", list(by = "score", value = 0.05))$hits, c("P01", "P02", "P03"))
  # explicit direction wins: |score| <= 1 for the same effect list
  h <- e$.paListHits(eff, "effect", "logFC", list(by = "score", value = 1, dir = "<="))
  expect_equal(h$hits, c("P03", "P04", "P05", "P06")); expect_equal(h$dir, "<=")
  # and the top share of a p-like score is its smallest p-values
  expect_equal(e$.paListHits(pl, "p", "score", list(by = "top", value = 0.5))$hits, c("P01", "P02", "P03"))
})

test_that("FDR and p cutoffs use the list's own columns", {
  e <- .load_list_env()
  st <- data.frame(gene_id = sprintf("P%02d", 1:5), score = c(3, 2, 1, -1, -2), p = c(0.001, 0.01, 0.04, 0.2, 0.5),
                   fdr = c(0.01, 0.04, 0.2, 0.5, 0.9))
  expect_equal(e$.paListHits(st, "effect", "logFC", list(by = "fdr", value = 0.05))$hits, c("P01", "P02"))
  expect_equal(e$.paListHits(st, "effect", "logFC", list(by = "p", value = 0.05))$hits, c("P01", "P02", "P03"))
})

test_that("a cutoff nobody passes, everybody passes, or with no such column is refused in words", {
  e <- .load_list_env()
  st <- data.frame(gene_id = sprintf("P%02d", 1:5), score = c(3, 2, 1, -1, -2), p = NA_real_, fdr = NA_real_)
  none <- e$.paListHits(st, "effect", "logFC", list(by = "score", value = 50, dir = ">="))
  expect_false(none$ok); expect_match(none$msg, "No protein has \\|logFC\\| >= 50 among the 5 proteins")
  all <- e$.paListHits(st, "effect", "logFC", list(by = "score", value = 0))
  expect_false(all$ok); expect_match(all$msg, "Every one of the 5 proteins")
  nofdr <- e$.paListHits(st, "effect", "logFC", list(by = "fdr", value = 0.05))
  expect_false(nofdr$ok); expect_match(nofdr$msg, "no FDR column")
  nop <- e$.paListHits(st, "effect", "logFC", list(by = "p", value = 0.05))
  expect_false(nop$ok); expect_match(nop$msg, "no p-value column")
  st$fdr <- c(0.4, 0.5, 0.6, 0.7, 0.8)
  low <- e$.paListHits(st, "effect", "logFC", list(by = "fdr", value = 0.05))
  expect_false(low$ok); expect_match(low$msg, "below 0.05")
})

# ---- SetListEnrichInput: the contract ----------------------------------------------------------

test_that("SetListEnrichInput returns 1 on success and 0 with a reason on refusal", {
  e <- .load_list_env(); .new_session(e); .save_list(e)
  expect_equal(e$SetListEnrichInput("datalist1", "ranked", "top", 0.10), 1)
  expect_equal(.ps(e)$list.enrich$list.type, "ranked")
  expect_equal(e$SetListEnrichInput("", "sig"), 1)
  expect_equal(.ps(e)$list.enrich$list.type, "sig")

  expect_equal(e$SetListEnrichInput("datalist1", "bogus"), 0)
  expect_match(tail(.msgs(e), 1), "significant list")
  expect_equal(e$SetListEnrichInput("datalist1", "ranked", "nope", 0.1), 0)
  expect_equal(e$SetListEnrichInput("datalist1", "ranked", "top", 0), 0)
  expect_match(tail(.msgs(e), 1), "above 0 and at most 100 percent")
  expect_equal(e$SetListEnrichInput("datalist1", "ranked", "top", "ten"), 0)
  expect_equal(e$SetListEnrichInput("datalist1", "ranked", "score", 1, "~"), 0)
  expect_match(tail(.msgs(e), 1), "direction of a score cutoff")
  expect_equal(e$SetListEnrichInput("datalist1", "ranked", "fdr", 0.05), 0)
  expect_match(tail(.msgs(e), 1), "no FDR column")
  expect_equal(e$SetListEnrichInput("datalist1", "ranked", "score", 99), 0)
  expect_match(tail(.msgs(e), 1), "No protein has")
})

test_that("a ranked list needs a score and one list", {
  e <- .load_list_env(); .new_session(e)
  .save_list(e, score = rep(0, 40))                     # a one-column list is padded with zeros
  expect_equal(e$SetListEnrichInput("datalist1", "ranked"), 0)
  expect_match(tail(.msgs(e), 1), "has no score")
  expect_equal(e$SetListEnrichInput("datalist1", "sig"), 1)    # still usable as a significant list
  .save_list(e)
  ps <- .ps(e); ps$numOfLists <- 2; e$saveSet(ps, "paramSet")
  expect_equal(e$SetListEnrichInput("datalist1", "ranked"), 0)
  expect_match(tail(.msgs(e), 1), "one list")
})

test_that("the direction of a score cutoff is recorded: defaulted by the kind of score, or the user's", {
  e <- .load_list_env(); .new_session(e)
  .save_list(e)                                         # an effect (a fold change)
  expect_equal(e$SetListEnrichInput("datalist1", "ranked", "score", 1), 1)
  expect_equal(.ps(e)$list.enrich$cutoff.dir, ">=")
  expect_equal(e$SetListEnrichInput("datalist1", "ranked", "score", 1, "<="), 1)
  expect_equal(.ps(e)$list.enrich$cutoff.dir, "<=")
  .save_list(e, score = seq(0.001, 0.9, length.out = 40), label = "P.Value")    # a p-like score
  expect_equal(e$SetListEnrichInput("datalist1", "ranked", "score", 0.05), 1)
  expect_equal(.ps(e)$list.enrich$cutoff.dir, "<=")
  expect_equal(e$SetListEnrichInput("datalist1", "ranked", "score", 0.5, ">="), 1)
  expect_equal(.ps(e)$list.enrich$cutoff.dir, ">=")
})

# ---- saying it -------------------------------------------------------------------------------------

test_that("GetListEnrichInputMsg says the contract's sentence for a complete list", {
  e <- .load_list_env(); .new_session(e)
  # 1,532 proteins; the 153rd strongest has |logFC| 1.2149
  score <- c(seq(5, 1.2149, length.out = 153), seq(1.2, 0.01, length.out = 1379)) * rep(c(1, -1), length.out = 1532)
  .save_list(e, n = 1532, score = score)
  expect_equal(e$SetListEnrichInput("datalist1", "ranked", "top", 0.10), 1)
  expect_equal(e$GetListEnrichInputMsg(),
               "Top 10% of the 1,532 proteins in your complete list (153, |logFC| ≥ 1.21) tested against all 1,532 as background")
  expect_equal(e$SetListEnrichInput("datalist1", "ranked", "score", 1, ">="), 1)
  expect_match(e$GetListEnrichInputMsg(), "^The [0-9,]+ proteins with \\|logFC\\| ≥ 1 of the 1,532 in your complete list tested against all 1,532 as background$")
})

test_that("the sentence for a significant list names its background, and the library one says it is weaker", {
  e <- .load_list_env(); .new_session(e); .save_list(e, n = 8)
  expect_equal(e$SetListEnrichInput("datalist1", "sig"), 1)
  expect_match(e$GetListEnrichInputMsg(), "tested against all proteins in the gene-set library as background, a weaker claim")
  ps <- .ps(e); ps$backgroundUniverse <- sprintf("P%02d", 1:30); e$saveSet(ps, "paramSet")
  expect_equal(e$GetListEnrichInputMsg(), "Your 8 significant proteins tested against the 30 proteins tested in your study as background")
  expect_equal(.ps(e)$list.enrich$list.type, "sig")
})

# ---- the three modes, through the real over-representation code ---------------------------------

test_that("a complete ranked list: top hits against the whole list; sets and hits restricted to it", {
  e <- .load_list_env(); .new_session(e); .save_list(e, n = 40)
  expect_equal(e$SetListEnrichInput("datalist1", "ranked", "top", 0.10), 1)       # the top 4: P01..P04
  expect_equal(e$PerformDefaultEnrichment("datalist1", "enr_ranked", "kegg"), 1)
  u <- .ps(e)$enrich.used
  expect_equal(.ps(e)$enrich.mode, "ranked")
  expect_equal(c(u$n.list, u$n.hits, u$query.size, u$bg.size), c(40, 4, 4, 40))
  # D has no member in the list and F one: neither is tested. A, A2, B, C, H, E and G (10 of its 20) are.
  expect_equal(u$sets.tested, 7L)
  res <- utils::read.csv("enr_ranked.csv", stringsAsFactors = FALSE)
  expect_false(any(c("D", "F") %in% res$Pathway))
  # a set counts its members in the list: G is 10, not 20 (the 10 in 51:60 are not in the list)
  a <- res[res$Pathway == "A", ]
  expect_equal(a$Total, 10); expect_equal(a$Hits, 4)
  # N is the 40 proteins of the list that are in the library; the test is on that, adjusted over the 7 tested sets
  raw <- vapply(list(A = c(10, 4), A2 = c(10, 4), B = c(10, 0), C = c(10, 0), H = c(10, 0), E = c(2, 2), G = c(10, 0)),
                function(x) stats::phyper(x[2] - 1, x[1], 40 - x[1], 4, lower.tail = FALSE), 0)
  raw[is.nan(raw)] <- 1
  expect_equal(a$Pval, signif(raw[["A"]], 3)); expect_equal(a$FDR, signif(stats::p.adjust(raw, "fdr")[["A"]], 3))
  expect_equal(u$cutoff$kind, "top")
  expect_equal(.ps(e)$enrich.cutoff, list(kind = "top", value = 0.1, dir = ">="))
  expect_equal(.ps(e)$enrich.background$label, "All proteins of your complete list")
  expect_equal(.ps(e)$enrich.background$size, 40L)
  expect_equal(.ps(e)$universe.opt.readable, "All proteins of your complete list")
  expect_match(e$GetListEnrichInputMsg(), "^Top 10% of the 40 proteins in your complete list \\(4, \\|logFC\\| ≥ 2.4\\) tested against all 40 as background$")
})

test_that("a significant list with the proteins tested: a set counts its tested members, hits outside them are left out", {
  e <- .load_list_env(); .new_session(e)
  .save_list(e, n = 8)                                                    # P01..P08 are significant
  ps <- .ps(e); ps$all.ent.mat <- matrix(1, nrow = 8, dimnames = list(c(sprintf("P%02d", 1:6), "P35", "P36"), "x"))
  ps$backgroundUniverse <- sprintf("P%02d", 1:30); e$saveSet(ps, "paramSet")   # what ReadUniverseFile leaves
  expect_equal(e$SetListEnrichInput("datalist1", "sig"), 1)
  expect_equal(e$PerformDefaultEnrichment("datalist1", "enr_measured", "kegg"), 1)
  u <- .ps(e)$enrich.used
  expect_equal(.ps(e)$enrich.mode, "ora_measured")
  # P35 and P36 are in the library (set H) but not among the proteins tested: left out, not counted in N
  expect_equal(c(u$query.input, u$query.size, u$bg.size, u$tested), c(8, 6, 30, 30))
  # H and G have no tested member, D none, F one: not tested. A, A2, B, C and E are.
  expect_equal(u$sets.tested, 5L)
  res <- utils::read.csv("enr_measured.csv", stringsAsFactors = FALSE)
  expect_false(any(c("H", "G", "D", "F") %in% res$Pathway))
  expect_equal(res$Total[res$Pathway == "A"], 10)
  raw <- vapply(list(A = c(10, 6), A2 = c(10, 4), B = c(10, 0), C = c(10, 0), E = c(2, 2)),
                function(x) stats::phyper(x[2] - 1, x[1], 30 - x[1], 6, lower.tail = FALSE), 0)
  raw[is.nan(raw)] <- 1
  expect_equal(res$FDR[res$Pathway == "A"], signif(stats::p.adjust(raw, "fdr")[["A"]], 3))
  expect_equal(.ps(e)$backgroundUniverse, sprintf("P%02d", 1:30))        # what the user uploaded is left as it was
  expect_equal(e$GetListEnrichInputMsg(),
               "Your 8 significant proteins (2 of them not among the tested proteins, left out) tested against the 30 proteins tested in your study as background")
})

test_that("a significant list alone is tested against the library and says so: the label never claims uploaded data", {
  e <- .load_list_env(); .new_session(e); .save_list(e, n = 8)
  ps <- .ps(e); ps$all.ent.mat <- matrix(1, nrow = 6, dimnames = list(sprintf("P%02d", 1:6), "x")); e$saveSet(ps, "paramSet")
  expect_equal(.ps(e)$universe.opt, "uploaded")                           # the default, with nothing uploaded
  expect_equal(e$PerformDefaultEnrichment("datalist1", "enr_cat", "kegg"), 1)
  expect_equal(.ps(e)$enrich.mode, "ora_catalog")
  expect_equal(.ps(e)$universe.opt.readable, "Gene Set Library")
  expect_equal(.ps(e)$enrich.background, list(label = "Gene Set Library", size = 50L))   # 1:40 and 51:60 are in some set
  expect_null(.ps(e)$backgroundUniverse)
  expect_match(e$GetListEnrichInputMsg(), "tested against all 50 proteins in the gene-set library as background, a weaker claim")
  # a second library is tested against ITS proteins, not the first run's library left behind as a "universe"
  e$lib <- LIB2
  expect_equal(e$PerformDefaultEnrichment("datalist1", "enr_cat2", "kegg"), 1)
  expect_null(.ps(e)$backgroundUniverse)
  expect_equal(.ps(e)$enrich.background$size, 18L)                        # the union of 1:6, 5:12 and 1,2,20:25
})

test_that("only sets with at least two tested members are tested, and the adjustment runs over all of them", {
  e <- .load_list_env(); .new_session(e); .save_list(e, n = 40)
  expect_equal(e$SetListEnrichInput("datalist1", "ranked", "top", 0.10), 1)
  e$PerformDefaultEnrichment("datalist1", "enr_two", "kegg")
  # F (one member) is not in the family the p-values are adjusted over
  expect_equal(.ps(e)$enrich.used$sets.tested, 7L)
})

test_that("a ranked list tests its whole ranking by rank, over the proteins of the list only", {
  skip_if_not_installed("fgsea")
  e <- .load_list_env(); .new_session(e); .save_list(e, n = 40)
  expect_equal(e$SetListEnrichInput("datalist1", "ranked"), 1)
  res <- suppressWarnings(e$PerformGSEA("datalist1", "gsea_ranked", "kegg", "gseanet", "gsea", 1, "multi", ""))
  expect_equal(res, 1)
  expect_equal(.ps(e)$enrich.mode, "ranked")
  expect_equal(.ps(e)$enrich.used$method, "gsea")
  tab <- utils::read.csv("gsea_ranked.csv", stringsAsFactors = FALSE)
  # G has 20 members in the library and 10 in the list: its size is 10; D (none in the list) is not there
  expect_equal(tab$Total[tab$Name == "G"], 10)
  expect_false("D" %in% tab$Name)
  expect_equal(e$GetListEnrichInputMsg(), "All 40 proteins of your complete list tested by rank (preranked GSEA on logFC)")
})

test_that("a significant list has nothing to rank", {
  skip_if_not_installed("fgsea")
  e <- .load_list_env(); .new_session(e); .save_list(e, n = 40)
  expect_equal(e$SetListEnrichInput("datalist1", "sig"), 1)
  expect_equal(suppressWarnings(e$PerformGSEA("datalist1", "gsea_sig", "kegg", "gseanet", "gsea", 1, "multi", "")), 0)
})

test_that("a cutoff that no longer fits the list is refused at run time, not tested as empty", {
  e <- .load_list_env(); .new_session(e); .save_list(e, n = 40)
  expect_equal(e$SetListEnrichInput("datalist1", "ranked", "score", 1.5), 1)
  .save_list(e, n = 40, score = rep(c(0.5, -0.4), 20))                   # the list was replaced: nobody reaches 1.5
  expect_equal(e$PerformDefaultEnrichment("datalist1", "enr_refused", "kegg"), 0)
  expect_match(tail(.msgs(e), 1), "No protein has")
  expect_false(file.exists("enr_refused.csv"))
})

# ---- reading a pasted list ------------------------------------------------------------------------

test_that("a pasted ranked list keeps its score name, p-value and FDR through MapListIds", {
  e <- .load_list_env(); .new_session(e)
  txt <- paste(c("Protein\tlogFC\tP.Value\tadj.P.Val",
                 sprintf("P%02d\t%s\t%s\t%s", 1:6, c(2.5, -1.2, 0.9, 0.4, -0.2, 0.1),
                         c(0.001, 0.01, 0.03, 0.2, 0.5, 0.9), c(0.01, 0.04, 0.08, 0.4, 0.7, 0.95))), collapse = "\n")
  expect_equal(e$MapListIds("datalist1", txt, "hsa", "uniprot"), 1)
  ds <- e$readDataset("datalist1")
  expect_equal(ds$score.label, "logFC")
  expect_equal(ds$list.stats$gene_id, sprintf("P%02d", 1:6))
  expect_equal(ds$list.stats$p[1:3], c(0.001, 0.01, 0.03))
  expect_equal(ds$list.stats$fdr[1:3], c(0.01, 0.04, 0.08))
  expect_equal(e$SetListEnrichInput("datalist1", "ranked", "fdr", 0.05), 1)
  expect_equal(e$.paListHits(ds$list.stats, "effect", "logFC", list(by = "fdr", value = 0.05))$hits, c("P01", "P02"))
  # a new list starts with no chosen input type and nothing tested
  expect_equal(e$MapListIds("datalist1", txt, "hsa", "uniprot"), 1)
  expect_null(.ps(e)$list.enrich)
})

test_that("a third and fourth column are read by their header, else as p-value and FDR", {
  e <- .load_list_env()
  m <- cbind(c("P01", "P02", "P03"), c("1", "2", "3"), c("0.1", "0.2", "0.3"))
  out <- e$.paListExtraColumns(m, c("Protein", "logFC", "adj.P.Val"))      # a third column named as an FDR is the FDR
  expect_null(out$p); expect_equal(unname(out$fdr), c(0.1, 0.2, 0.3))
  out <- e$.paListExtraColumns(m, NULL)                                      # no header: the third is the p-value
  expect_equal(unname(out$p), c(0.1, 0.2, 0.3)); expect_null(out$fdr)
})

test_that("the analysis summary names the background a list was tested against", {
  e <- .load_list_env(); .new_session(e); .save_list(e, n = 40)
  e$SetListEnrichInput("datalist1", "ranked", "top", 0.10)
  e$PerformDefaultEnrichment("datalist1", "enr_report", "kegg")
  ps <- .ps(e)
  an <- list(enr.mat = utils::read.csv("enr_report.csv"))
  txt <- paste(e$.pa_enrich_section(ps, an), collapse = "\n")
  expect_match(txt, "Background / universe: +All proteins of your complete list")
  expect_match(txt, "Background feature count: +40")
  expect_match(txt, "What was tested: +Top 10% of the 40 proteins")
})

# ---- the pathway-impact analysis honours the same input ---------------------------------------

test_that("pathway impact (ORA): a ranked list tests its top hits against the whole list", {
  e <- .load_list_env(); .new_session(e); .save_list(e, n = 40)
  expect_equal(e$SetListEnrichInput("datalist1", "ranked", "top", 0.10), 1)
  expect_equal(e$CalculateEnzymePathwayOra("datalist1"), 1)
  u <- .ps(e)$enrich.used
  expect_equal(.ps(e)$enrich.mode, "ranked")
  expect_equal(c(u$n.hits, u$query.size, u$bg.size, u$sets.tested), c(4, 4, 40, 7))
  tab <- as.data.frame(arrow::read_feather("pathway_enzyme_ora.arrow"))
  a <- tab[tab$row_names_id == "A", ]
  expect_equal(c(a$total, a$hits), c(10, 4))
  expect_equal(a$rawP, stats::phyper(3, 10, 30, 4, lower.tail = FALSE))
  expect_false(any(c("D", "F", "H", "G") %in% tab$row_names_id))     # no hit, or not tested at all
})

test_that("pathway impact (ORA): a significant list is tested against the proteins tested, not the library", {
  e <- .load_list_env(); .new_session(e)
  .save_list(e, n = 8, sig.ids = c(sprintf("P%02d", 1:6), "P35", "P36"))
  ps <- .ps(e); ps$backgroundUniverse <- sprintf("P%02d", 1:30); e$saveSet(ps, "paramSet")
  expect_equal(e$SetListEnrichInput("datalist1", "sig"), 1)
  expect_equal(e$CalculateEnzymePathwayOra("datalist1"), 1)
  u <- .ps(e)$enrich.used
  expect_equal(.ps(e)$enrich.mode, "ora_measured")
  expect_equal(c(u$query.input, u$query.size, u$bg.size, u$sets.tested), c(8, 6, 30, 5))
  tab <- as.data.frame(arrow::read_feather("pathway_enzyme_ora.arrow"))
  expect_equal(tab$total[tab$row_names_id == "A"], 10)
  expect_equal(tab$rawP[tab$row_names_id == "A"], stats::phyper(5, 10, 20, 6, lower.tail = FALSE))
})

test_that("pathway impact (ORA): a significant list alone falls back to the library without claiming uploaded data", {
  e <- .load_list_env(); .new_session(e)
  .save_list(e, n = 8, sig.ids = sprintf("P%02d", 1:6))
  expect_equal(e$SetListEnrichInput("datalist1", "sig"), 1)
  expect_equal(e$CalculateEnzymePathwayOra("datalist1"), 1)
  expect_equal(.ps(e)$enrich.mode, "ora_catalog")
  expect_equal(.ps(e)$universe.opt.readable, "Gene Set Library")
  expect_false(any(grepl("No measured proteome", .ps(e)$enrich.used$label)))
})

test_that("pathway impact (GSEA) ranks a complete ranked list by its own score", {
  skip_if_not_installed("fgsea")
  e <- .load_list_env(); .new_session(e); .save_list(e, n = 40)
  expect_equal(e$SetListEnrichInput("datalist1", "ranked"), 1)
  expect_equal(suppressWarnings(e$CalculateEnzymePathwayGsea("datalist1")), 1)
  expect_equal(.ps(e)$enrich.mode, "ranked"); expect_equal(.ps(e)$enrich.used$method, "gsea")
  expect_equal(e$GetListEnrichInputMsg(), "All 40 proteins of your complete list tested by rank (preranked GSEA on logFC)")
})
