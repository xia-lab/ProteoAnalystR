# Regression tests for covariate/batch adjustment of hold-out / new-sample
# prediction inputs (.AdjustNewSamplesForPrediction in biomarker_utils.R).
#
# Guards against the train/serve mismatch where the final model is trained on the
# covariate-adjusted matrix but new/unlabelled samples are scored un-adjusted.

R_DIR <- normalizePath(file.path(testthat::test_path(), "..", "..", "R"), mustWork = FALSE)
BIOMARKER <- file.path(R_DIR, "biomarker_utils.R")

.load_biomarker_env <- function() {
  skip_if_not(file.exists(BIOMARKER), "biomarker_utils.R not found")
  e <- new.env()
  sys.source(BIOMARKER, envir = e)
  e
}

# Shared synthetic dataset: real class effect + a large per-feature batch shift.
.make_fixture <- function(seed = 42) {
  set.seed(seed)
  p <- 60; n <- 40
  batch <- factor(rep(c("A", "B"), each = n / 2))
  cls   <- factor(rep(c("ctrl", "case"), times = n / 2))
  mu        <- matrix(rnorm(p * n), nrow = n)
  class.eff <- outer(ifelse(cls == "case", 1, 0), rnorm(p, 0, 0.8))
  batch.eff <- outer(ifelse(batch == "B", 1, 0), rnorm(p, 3, 0.5))
  X <- mu + class.eff + batch.eff
  rownames(X) <- paste0("s", 1:n); colnames(X) <- paste0("f", 1:p)

  new.idx <- c(sample(which(batch == "A"), 5), sample(which(batch == "B"), 5))
  tr.idx  <- setdiff(1:n, new.idx)
  meta.all <- data.frame(Class = as.character(cls), Batch = as.character(batch),
                         row.names = rownames(X), stringsAsFactors = FALSE)
  meta.tr  <- meta.all[tr.idx, ]; meta.tr$Class <- factor(meta.tr$Class)

  dataSet <- list(
    covariate.adjusted = TRUE,
    meta.info.original = meta.all,
    meta.info          = meta.all,
    cv.adjust = list(primary.condition = "Class", adj.vec = "Batch",
                     base.data = X[tr.idx, ], meta.info = meta.tr)
  )
  list(X = X, batch = batch, new.idx = new.idx, tr.idx = tr.idx,
       meta.all = meta.all, meta.tr = meta.tr, dataSet = dataSet)
}

test_that("projecting labelled rows reproduces limma::removeBatchEffect exactly", {
  skip_if_not_installed("limma")
  e <- .load_biomarker_env()
  f <- get(".AdjustNewSamplesForPrediction", envir = e)
  fx <- .make_fixture()

  ref <- t(limma::removeBatchEffect(
    t(fx$X[fx$tr.idx, ]),
    covariates = model.matrix(~Batch, fx$meta.tr)[, -1, drop = FALSE],
    design     = model.matrix(~Class, fx$meta.tr)))

  lab <- f(fx$dataSet, fx$X[fx$tr.idx, ])
  expect_equal(max(abs(lab - ref)), 0, tolerance = 1e-8)
})

test_that("full-training predictive adjustment is the same linear feature space", {
  skip_if_not_installed("limma")
  e <- .load_biomarker_env()
  adjust.full <- get(".LinearAdjustmentForBiomarker", envir = e)
  fx <- .make_fixture()

  got <- adjust.full(t(fx$X[fx$tr.idx, ]), fx$meta.tr, "Class", "Batch")
  expected <- limma::removeBatchEffect(
    t(fx$X[fx$tr.idx, ]),
    covariates = model.matrix(~Batch, fx$meta.tr)[, -1, drop = FALSE],
    design = model.matrix(~Class, fx$meta.tr))

  expect_equal(got, expected, tolerance = 1e-10)
})

test_that("batch shift is removed from new/held-out samples (the bug)", {
  e <- .load_biomarker_env()
  f <- get(".AdjustNewSamplesForPrediction", envir = e)
  fx <- .make_fixture()

  new.raw <- fx$X[fx$new.idx, ]
  new.adj <- f(fx$dataSet, new.raw)
  bnew <- fx$batch[fx$new.idx]

  raw.gap <- mean(abs(colMeans(new.raw[bnew == "B", ]) - colMeans(new.raw[bnew == "A", ])))
  adj.gap <- mean(abs(colMeans(new.adj[bnew == "B", ]) - colMeans(new.adj[bnew == "A", ])))
  expect_gt(raw.gap, 1)              # raw new rows carry the batch shift
  expect_lt(adj.gap, raw.gap / 3)    # adjustment collapses it
})

test_that("unseen batch level rejects prediction with a warning", {
  e <- .load_biomarker_env()
  f <- get(".AdjustNewSamplesForPrediction", envir = e)
  fx <- .make_fixture()
  new.raw <- fx$X[fx$new.idx, ]

  bad <- fx$dataSet
  meta.bad <- fx$meta.all; meta.bad[fx$new.idx[1], "Batch"] <- "C"  # level not in training
  bad$meta.info.original <- meta.bad
  expect_warning(out <- f(bad, new.raw))
  expect_null(out)
})

test_that("no covariate adjustment active is a no-op", {
  e <- .load_biomarker_env()
  f <- get(".AdjustNewSamplesForPrediction", envir = e)
  fx <- .make_fixture()
  new.raw <- fx$X[fx$new.idx, ]
  expect_identical(f(list(covariate.adjusted = FALSE), new.raw), new.raw)
})
