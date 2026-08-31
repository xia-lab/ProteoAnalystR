# Controlled examples for identifier mapping and duplicate-collapse policies.

R_DIR <- normalizePath(file.path(testthat::test_path(), "..", "..", "R"), mustWork = FALSE)

.source_in_clean_env <- function(file) {
  skip_if_not(file.exists(file.path(R_DIR, file)), paste(file, "not found"))
  e <- new.env(parent = globalenv())
  sys.source(file.path(R_DIR, file), envir = e)
  e
}

test_that("one-to-many mappings select the first database row", {
  e <- .source_in_clean_env("data_idanot_utils.R")
  map_ids <- get("doIdMappingGeneric", envir = e)
  gene.map <- data.frame(
    accession = c("A", "A", "B"),
    gene_id = c("101", "102", "201"),
    symbol = c("FIRST", "SECOND", "ONLY"),
    name = c("first", "second", "only"),
    stringsAsFactors = FALSE
  )

  expect_identical(map_ids(c("A", "B", "X"), gene.map,
                           "accession", "gene_id", "vec"),
                   c("101", "201", "X"))
  expect_identical(map_ids("A", gene.map, "accession", "symbol", "vec"),
                   "FIRST")
})

test_that("UniProt suffix cleanup is confined to lookup normalization", {
  e <- .source_in_clean_env("data_idanot_utils.R")
  normalize_ids <- get(".normalizeUniProtAccessionsForLookup", envir = e)

  expect_identical(
    normalize_ids(c("sp|P12345-2|PROT_HUMAN", "Q9TEST_S_17", "A0ABC1;Q00001")),
    c("P12345", "Q9TEST", "A0ABC1")
  )
})

test_that("ranked duplicate IDs retain max absolute statistic deterministically", {
  e <- .source_in_clean_env("misc_utils.R")
  collapse <- get(".paCollapseRankedStats", envir = e)

  out <- collapse(c(2, -5, 5, 1, -3), c("G1", "G1", "G1", "G2", "G2"))
  expect_identical(names(out), c("G1", "G2"))
  expect_equal(unname(out), c(-5, -3))
})

test_that("quantitative duplicate rows use the selected sample-wise aggregator", {
  e <- .source_in_clean_env("misc_utils.R")
  assign("readSet", function(x, name) x, envir = e)
  assign("saveSet", function(x, name) invisible(NULL), envir = e)
  remove_duplicates <- get("RemoveDuplicates", envir = e)
  mat <- matrix(c(1, 3, 5, 7, 10, 20), nrow = 3, byrow = TRUE,
                dimnames = list(c("P1", "P1", "P2"), c("S1", "S2")))

  out <- remove_duplicates(mat, "mean", quiet = TRUE,
                           paramSet = list(numOfLists = 1), msgSet = list())[[1]]
  expect_equal(unname(out["P1", ]), c(3, 5))
  expect_equal(unname(out["P2", ]), c(10, 20))
})

test_that("phosphosite representative uses statistic magnitude, not sign", {
  e <- .source_in_clean_env("phospho_data_utils.R")
  collapse <- get("phosphoCollapseLocal", envir = e)
  mat <- matrix(c(1, 2, 10, 20, 100, 200), nrow = 3, byrow = TRUE,
                dimnames = list(c("site_a", "site_b", "site_c"), c("S1", "S2")))

  out <- collapse(mat, id = c("P1", "P1", "P2"), stat = c(-8, 3, 1))
  expect_identical(rownames(out), c("P1", "P2"))
  expect_equal(unname(out["P1", ]), c(1, 2))
})

test_that("phosphosite ORA uses unique mapped genes and retains provenance", {
  e <- .source_in_clean_env("phospho_enrichment_utils.R")
  assign(".isEntrezLikeForPhosphoEnrichment", function(ids) FALSE, envir = e)
  assign(".mapUniprotToEntrezForEnrichment", function(ids, paramSet) {
    c("101", "101", "202", NA_character_)
  }, envir = e)
  collapse <- get(".collapsePhosphositesToGenes", envir = e)
  sites <- c("P1_S_1", "P1_T_2", "P2_Y_3", "UNMAPPED_S_4")

  out <- collapse(sites, list(data.org = "hsa"))
  expect_identical(out$entrez.ids, c("101", "202"))
  expect_identical(out$entrez.to.phospho[["101"]], sites[1:2])
  expect_false(out$mapped.inx[4])
})
