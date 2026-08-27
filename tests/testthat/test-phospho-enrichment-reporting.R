# Regression tests for the self-documenting phosphosite enrichment tables.

R_DIR <- normalizePath(file.path(testthat::test_path(), "..", "..", "R"), mustWork = FALSE)
PHOSPHO_ENRICH <- file.path(R_DIR, "phospho_enrich_utils.R")

.load_phospho_reporting_env <- function() {
  skip_if_not(file.exists(PHOSPHO_ENRICH), "phospho_enrich_utils.R not found")
  e <- new.env()
  sys.source(PHOSPHO_ENRICH, envir = e)
  e
}

test_that("reported DE rule preserves raw-p and custom thresholds", {
  e <- .load_phospho_reporting_env()
  thresholds <- get(".phosphoSigThresholds", envir = e)
  annotate <- get(".annotateEnrichReport", envir = e)
  param <- list(BHth = 0.01, fc.thresh = 1.25, pval.selection = "raw")

  expect_equal(thresholds(param), list(fdr = 0.01, fc = 1.25, selection = "raw"))

  out <- annotate(
    data.frame(P_value = c(0.001, 0.2), FDR = c(0.02, 0.2)),
    param, n_tested = 17, n_background = 250, enrich_fdr = 0.05
  )
  expect_identical(out$DE_Sig_Type, c("raw", "raw"))
  expect_equal(out$DE_Sig_Cutoff, c(0.01, 0.01))
  expect_equal(out$DE_log2FC_Cutoff, c(1.25, 1.25))
  expect_identical(out$Significant, c("Yes", "No"))
})

test_that("adjusted-p selection is labeled explicitly", {
  e <- .load_phospho_reporting_env()
  thresholds <- get(".phosphoSigThresholds", envir = e)
  expect_identical(
    thresholds(list(BHth = 0.05, fc.thresh = 0, pval.selection = "fdr"))$selection,
    "FDR-adjusted"
  )
})
