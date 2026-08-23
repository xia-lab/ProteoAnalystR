# Smoke test: every declared dependency is installed and loadable.
# Reads Imports (hard requirements) + Suggests (optional) from DESCRIPTION so the
# check stays in sync with the package metadata.

.read_desc_field <- function(field) {
  desc <- normalizePath(file.path(testthat::test_path(), "..", "..", "DESCRIPTION"),
                        mustWork = FALSE)
  if (!file.exists(desc)) return(character(0))
  d <- read.dcf(desc)
  if (!(field %in% colnames(d))) return(character(0))
  raw <- d[1, field]
  pk <- trimws(unlist(strsplit(raw, "[,\n]")))
  pk <- sub("\\s*\\(.*\\)$", "", pk)             # drop version constraints
  pk <- pk[nzchar(pk) & pk != "R"]
  unique(pk)
}

test_that("all Imports (hard dependencies) are installed and loadable", {
  imports <- .read_desc_field("Imports")
  skip_if(length(imports) == 0, "no Imports field found")
  missing <- imports[!vapply(imports, requireNamespace, logical(1), quietly = TRUE)]
  expect_true(length(missing) == 0,
              info = paste("Missing required packages:", paste(missing, collapse = ", ")))
})

test_that("key analysis backends load and report a version", {
  # The statistical engines the pipelines depend on directly.
  core <- c("MSstats", "limma", "DEqMS", "edgeR", "DESeq2", "impute", "imputeLCMD",
            "preprocessCore", "data.table", "qs", "matrixStats", "igraph", "fgsea",
            "WGCNA", "pROC")
  for (p in core) {
    ok <- requireNamespace(p, quietly = TRUE)
    expect_true(ok, info = paste0(p, " is not installed"))
    if (ok) expect_silent(v <- as.character(utils::packageVersion(p)))
  }
})

test_that("optional (Suggests) packages, if declared, are reported", {
  suggests <- .read_desc_field("Suggests")
  skip_if(length(suggests) == 0, "no Suggests field")
  present <- suggests[vapply(suggests, requireNamespace, logical(1), quietly = TRUE)]
  # Not a hard failure — record which optional method packages are available.
  testthat::expect_type(present, "character")
  message("Optional packages present: ",
          if (length(present)) paste(present, collapse = ", ") else "(none)")
})
