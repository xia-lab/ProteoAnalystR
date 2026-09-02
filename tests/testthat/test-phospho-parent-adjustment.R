# Regression tests for strict sample-wise phosphosite/parent-protein adjustment.

R_DIR <- normalizePath(file.path(testthat::test_path(), "..", "..", "R"), mustWork = FALSE)

.load_norm_utils <- function() {
  skip_if_not(file.exists(file.path(R_DIR, "norm_utils.R")), "norm_utils.R not found")
  e <- new.env(parent = globalenv())
  sys.source(file.path(R_DIR, "norm_utils.R"), envir = e)
  get(".correctPhosphoByProteinAbundance", envir = e)
}

test_that("parent adjustment never retains an unadjusted site value", {
  adjust <- .load_norm_utils()
  phospho <- matrix(
    c(10, 11, 12,
      20, 21, 22,
      30, 31, 32),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("P1234_S_10", "P9999_T_2", "PAMB_Y_3"),
                    c("S1", "S2", "S3"))
  )
  protein <- matrix(
    c(4, NA,
      8, 9,
      5, 6,
      7, 8),
    nrow = 4, byrow = TRUE,
    dimnames = list(c("P1234", "P12345", "PAMB-2", "PAMB-3"),
                    c("S1", "S2"))
  )

  out <- adjust(phospho, protein)

  # P1234 is an exact, unique match. Missing parent data and the sample absent
  # from the reference remain NA rather than falling back to site intensities.
  expect_identical(rownames(out), "P1234_S_10")
  expect_equal(unname(out[1, ]), c(6, NA, NA))
  # P9999 is unmatched; PAMB maps to two isoform rows and is ambiguous. Both
  # are removed from the protein-adjusted matrix.
  expect_false(any(c("P9999_T_2", "PAMB_Y_3") %in% rownames(out)))
})

test_that("canonical matching does not confuse accession prefixes", {
  adjust <- .load_norm_utils()
  phospho <- matrix(10, nrow = 1,
                    dimnames = list("P1234_S_10", "S1"))
  protein <- matrix(3, nrow = 1,
                    dimnames = list("P12345", "S1"))

  out <- adjust(phospho, protein)
  expect_equal(nrow(out), 0)
})
