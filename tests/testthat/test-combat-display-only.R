# ComBat may be requested for the exploratory PCA display, but the matrix used
# by every biomarker model must remain the linear train-fit/test-apply space.

R_DIR <- normalizePath(file.path(testthat::test_path(), "..", "..", "R"), mustWork = FALSE)
BIOMARKER <- file.path(R_DIR, "biomarker_utils.R")

test_that("ComBat output is display-only and predictive data are linear", {
  skip_if_not_installed("limma")
  skip_if_not_installed("sva")
  skip_if_not(file.exists(BIOMARKER), "biomarker_utils.R not found")

  e <- new.env(parent = globalenv())
  sys.source(BIOMARKER, envir = e)
  set.seed(91)
  n <- 24; p <- 30
  cls <- factor(rep(c("ctrl", "case"), each = n / 2))
  batch <- factor(rep(c("B1", "B2", "B3"), times = n / 3))
  x <- matrix(rnorm(n * p), nrow = n,
              dimnames = list(paste0("S", seq_len(n)), paste0("F", seq_len(p))))
  x <- x + outer(as.numeric(batch) - 1, rnorm(p, 1.5, 0.2))
  meta <- data.frame(Class = cls, Batch = batch, row.names = rownames(x))
  ds <- list(data.norm.transposed = x, data.norm = t(x), meta.info = meta, cls = cls)

  assign("msgSet", list(), envir = e)
  assign("paramSet", list(), envir = e)
  assign("adj.vec", "Batch", envir = e)
  assign("readDataset", function(name) ds, envir = e)
  assign("readSet", function(x, name) list(), envir = e)
  assign("saveSet", function(x, name) invisible(NULL), envir = e)
  assign("msg", function(...) invisible(NULL), envir = e)
  assign("RegisterData", function(x) {
    assign("captured.dataset", x, envir = e)
    1
  }, envir = e)

  run <- get("PerformCovariateAdjustmentForROC", envir = e)
  expect_identical(run("toy", "Class", use.combat = TRUE), 1)
  out <- get("captured.dataset", envir = e)
  linear <- get(".LinearAdjustmentForBiomarker", envir = e)(t(x), meta, "Class", "Batch")

  expect_identical(out$combat.scope, "exploratory_display_only")
  expect_equal(out$data.norm.transposed, t(linear), tolerance = 1e-8)
  expect_equal(out$biomarker.model.data, t(linear), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(out$data.adjusted.transposed,
                                out$biomarker.model.data, tolerance = 1e-8)))
})
