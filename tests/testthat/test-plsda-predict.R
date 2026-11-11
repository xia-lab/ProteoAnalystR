# Regression tests for the PLS-DA branch of Predict.class (biomarker_utils.R).
#
# Guards against (1) unnamed PLS scores, which left every sample at 0.5 in the
# Prob. View and collapsed the confusion matrix to one row, and (2) per-batch
# min-max rescaling, which gave NaN for one new sample and forced two new
# samples of the same class to opposite classes.

R_DIR <- normalizePath(file.path(testthat::test_path(), "..", "..", "R"), mustWork = FALSE)
BIOMARKER <- file.path(R_DIR, "biomarker_utils.R")

test_that("PLS-DA predictions are named, in [0,1] and independent of batch size", {
  skip_if_not(file.exists(BIOMARKER), "biomarker_utils.R not found")
  skip_if_not_installed("pls")
  e <- new.env()
  sys.source(BIOMARKER, envir = e)

  set.seed(1)
  n <- 40; p <- 5
  cls <- factor(rep(c("Tumour", "Others"), each = n / 2), levels = c("Tumour", "Others"))
  X <- matrix(rnorm(n * p), n) + outer(ifelse(cls == "Others", 2, 0), rep(1, p))
  rownames(X) <- paste0("s", 1:n); colnames(X) <- paste0("f", 1:p)
  tr <- c(1:15, 21:35); te <- setdiff(1:n, tr)

  prob <- e$Predict.class(X[tr, ], cls[tr], X[te, ], "pls", 2)
  expect_identical(names(prob), rownames(X)[te])
  expect_true(all(prob >= 0 & prob <= 1))
  # second level ("Others") is the positive class, as for SVM / RF
  expect_true(mean((prob > 0.5) == (cls[te] == "Others")) >= 0.9)

  one <- e$Predict.class(X[tr, ], cls[tr], X[te[1], , drop = FALSE], "pls", 2)
  expect_false(is.nan(one))
  expect_equal(unname(one), unname(prob[1]))

  two <- e$Predict.class(X[tr, ], cls[tr], X[te[1:2], , drop = FALSE], "pls", 2)
  expect_equal(two, prob[1:2])

  res <- e$Predict.class(X[tr, ], cls[tr], X[te, ], "pls", 2, imp.out = TRUE)
  expect_identical(names(res$prob.out), rownames(X)[te])
})
