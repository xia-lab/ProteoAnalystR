# Regression tests for the leakage-free property of the fold-wise covariate
# adjustment (.RemoveCovariateEffectFold in biomarker_utils.R). The correction
# applied to a held-out fold must be estimated on the TRAINING rows only, so:
#   (a) the training-fold adjusted values do not depend on held-out expression, and
#   (b) perturbing a held-out sample's expression passes straight through the
#       correction (the subtracted covariate term is unchanged) — i.e. held-out
#       values never enter the fitted adjustment.

R_DIR <- normalizePath(file.path(testthat::test_path(), "..", "..", "R"), mustWork = FALSE)
BIOMARKER <- file.path(R_DIR, "biomarker_utils.R")

.load_env <- function() {
  skip_if_not(file.exists(BIOMARKER), "biomarker_utils.R not found")
  e <- new.env(); sys.source(BIOMARKER, envir = e); e
}

.fixture <- function(seed = 7) {
  set.seed(seed)
  p <- 40; n <- 30
  batch <- factor(rep(c("A", "B"), length.out = n))
  cls   <- factor(rep(c("ctrl", "case"), each = n / 2))
  X <- matrix(rnorm(p * n), nrow = n)
  X <- X + outer(ifelse(batch == "B", 1, 0), rnorm(p, 2, 0.5))   # batch effect
  meta <- data.frame(Class = cls, Batch = batch)
  d <- get(".BuildCovariateDesign", envir = .env)(meta, "Class", "Batch")
  tr <- 1:20; te <- 21:30
  list(X = X, prim = d$primary.design, cov = d$covariate.design, tr = tr, te = te)
}

test_that("held-out expression does not affect the TRAINING-fold correction", {
  .env <<- .load_env()
  fold <- get(".RemoveCovariateEffectFold", envir = .env)
  fx <- .fixture()
  x.tr <- fx$X[fx$tr, ]; x.te <- fx$X[fx$te, ]
  a1 <- fold(x.tr, x.te, fx$prim[fx$tr, , drop = FALSE],
             fx$cov[fx$tr, , drop = FALSE], fx$cov[fx$te, , drop = FALSE])
  # perturb ONLY the held-out expression, refit
  x.te2 <- x.te + matrix(rnorm(length(x.te), 0, 5), nrow = nrow(x.te))
  a2 <- fold(x.tr, x.te2, fx$prim[fx$tr, , drop = FALSE],
             fx$cov[fx$tr, , drop = FALSE], fx$cov[fx$te, , drop = FALSE])
  expect_equal(a1$train, a2$train, tolerance = 1e-10)   # training output invariant
})

test_that("perturbing held-out expression passes straight through the correction", {
  .env <<- .load_env()
  fold <- get(".RemoveCovariateEffectFold", envir = .env)
  fx <- .fixture()
  x.tr <- fx$X[fx$tr, ]; x.te <- fx$X[fx$te, ]
  delta <- matrix(rnorm(length(x.te), 0, 3), nrow = nrow(x.te))
  a1 <- fold(x.tr, x.te,          fx$prim[fx$tr, , drop=FALSE], fx$cov[fx$tr, , drop=FALSE], fx$cov[fx$te, , drop=FALSE])
  a2 <- fold(x.tr, x.te + delta,  fx$prim[fx$tr, , drop=FALSE], fx$cov[fx$tr, , drop=FALSE], fx$cov[fx$te, , drop=FALSE])
  # the correction term is identical (fit on training), so the difference == delta
  expect_equal(unname(a2$test - a1$test), unname(delta), tolerance = 1e-10)
})

test_that("no covariate design is a safe no-op", {
  .env <<- .load_env()
  fold <- get(".RemoveCovariateEffectFold", envir = .env)
  fx <- .fixture()
  out <- fold(fx$X[fx$tr, ], fx$X[fx$te, ], fx$prim[fx$tr, , drop = FALSE], NULL, NULL)
  expect_equal(out$train, fx$X[fx$tr, ])
  expect_equal(out$test,  fx$X[fx$te, ])
})

test_that("a non-estimable training-fold adjustment fails explicitly", {
  .env <<- .load_env()
  fold <- get(".RemoveCovariateEffectFold", envir = .env)
  fx <- .fixture()
  # The requested covariate is identical to the protected class term in this
  # training split, so its effect cannot be estimated separately.
  aliased <- fx$prim[, 2, drop = FALSE]
  expect_error(
    fold(fx$X[fx$tr, ], fx$X[fx$te, ],
         fx$prim[fx$tr, , drop = FALSE],
         aliased[fx$tr, , drop = FALSE], aliased[fx$te, , drop = FALSE]),
    "not estimable"
  )
})

test_that("permutation adjustment protects the permuted, not original, outcome", {
  .env <<- .load_env()
  adjust.perm <- get(".AdjustPermutationFold", envir = .env)
  fold <- get(".RemoveCovariateEffectFold", envir = .env)
  fx <- .fixture()
  perm <- factor(sample(factor(rep(c("ctrl", "case"), each = 15))))
  ctx <- list(covariate.design = fx$cov,
              primary.design = matrix(999, nrow = nrow(fx$X), ncol = 2))

  got <- adjust.perm(ctx, fx$X[fx$tr, ], fx$X[fx$te, ], fx$tr, fx$te, perm)
  perm.design <- model.matrix(~ factor(perm, levels = levels(factor(perm))))
  expected <- fold(fx$X[fx$tr, ], fx$X[fx$te, ],
                   perm.design[fx$tr, , drop = FALSE],
                   fx$cov[fx$tr, , drop = FALSE], fx$cov[fx$te, , drop = FALSE])

  expect_equal(got$train, expected$train, tolerance = 1e-10)
  expect_equal(got$test, expected$test, tolerance = 1e-10)
})
