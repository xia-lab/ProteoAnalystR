# Regression tests for phosphoproteomics group-aware hybrid imputation.
R_DIR <- normalizePath(file.path(testthat::test_path(), "..", "..", "R"), mustWork = FALSE)

load_group_imputer <- function() {
  e <- new.env(parent = globalenv())
  sys.source(file.path(R_DIR, "meta_proc.R"), envir = e)
  sys.source(file.path(R_DIR, "data_impute.R"), envir = e)
  get(".imputeMinGroupAware", envir = e)
}

test_that("sporadic phosphosite gaps use reproducible group mean plus noise", {
  f <- load_group_imputer()
  mat <- rbind(
    site1 = c(20, 22, NA, 30, 32, 31),
    site2 = c(10, 14, NA, 40, 44, NA)
  )
  colnames(mat) <- paste0("s", seq_len(ncol(mat)))
  grp <- c("A", "A", "A", "B", "B", "B")

  x1 <- f(mat, grp, seed = 91L)
  x2 <- f(mat, grp, seed = 91L)
  x3 <- f(mat, grp, seed = 92L)

  expect_equal(x1, x2)
  expect_false(isTRUE(all.equal(x1[is.na(mat)], x3[is.na(mat)])))
  expect_equal(x1[!is.na(mat)], mat[!is.na(mat)])
  expect_false(isTRUE(all.equal(x1["site1", "s3"], mean(c(20, 22)))))
  expect_false(isTRUE(all.equal(x1["site2", "s6"], mean(c(40, 44)))))
})

test_that("censored groups use the canonical scale-aware LoD", {
  f <- load_group_imputer()
  mat <- rbind(site = c(20, 21, 22, NA, NA, NA))
  grp <- c("A", "A", "A", "B", "B", "B")

  out <- f(mat, grp, seed = 7L)

  expect_equal(out[1, 4:6], rep(20 - log2(5), 3), tolerance = 1e-12)
})

test_that("stochastic phosphosite imputation does not alter caller RNG state", {
  f <- load_group_imputer()
  mat <- rbind(site = c(20, 22, NA, 30, 31, 32))
  grp <- c("A", "A", "A", "B", "B", "B")

  set.seed(314)
  before <- .Random.seed
  invisible(f(mat, grp, seed = 123L))
  expect_equal(.Random.seed, before)
})
