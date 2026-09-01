# Regression tests for the native feature-level MSstatsPTM workflow.

R_DIR <- normalizePath(file.path(testthat::test_path(), "..", "..", "R"), mustWork = FALSE)
PHOSPHO_ENRICH <- file.path(R_DIR, "phospho_enrich_utils.R")
PHOSPHO_DATA <- file.path(R_DIR, "phospho_data_utils.R")

.load_msstatsptm_helpers <- function() {
  skip_if_not(file.exists(PHOSPHO_ENRICH), "phospho_enrich_utils.R not found")
  e <- new.env(parent = globalenv())
  sys.source(PHOSPHO_ENRICH, envir = e)
  e
}

test_that("already-summarized site-by-run tables do not masquerade as raw features", {
  e <- .load_msstatsptm_helpers()
  ready <- get(".msstatsptmRawReady", envir = e)
  summarized <- data.frame(
    ProteinName = c("P1_K10", "P1_K10"),
    Run = c("r1", "r2"), Condition = c("A", "B"),
    BioReplicate = c("1", "1"), Intensity = c(20, 21)
  )

  expect_false(ready(list(PTM = summarized, PROTEIN = summarized), c("r1", "r2")))
})

test_that("canonical upload retains both converter-level branches", {
  skip_if_not(file.exists(PHOSPHO_DATA), "phospho_data_utils.R not found")
  skip_if_not_installed("MSstatsPTM")
  data_env <- new.env(parent = emptyenv())
  utils::data("raw.input", package = "MSstatsPTM", envir = data_env)
  raw <- get("raw.input", envir = data_env)

  td <- tempfile("msstatsptm-upload-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  ptm.file <- file.path(td, "ptm.csv")
  protein.file <- file.path(td, "protein.csv")
  data.table::fwrite(raw$PTM, ptm.file)
  data.table::fwrite(raw$PROTEIN, protein.file)

  e <- new.env(parent = globalenv())
  sys.source(PHOSPHO_DATA, envir = e)
  parsed <- get(".readMSstatsPTMLong", envir = e)(ptm.file, protein.file)
  expect_s3_class(parsed$msstatsptm.raw$PTM, "data.frame")
  expect_s3_class(parsed$msstatsptm.raw$PROTEIN, "data.frame")
  expect_equal(nrow(parsed$msstatsptm.raw$PTM), nrow(raw$PTM))
  expect_equal(nrow(parsed$msstatsptm.raw$PROTEIN), nrow(raw$PROTEIN))
  expect_true(all(c("PeptideSequence", "Intensity") %in%
                    colnames(parsed$msstatsptm.raw$PTM)))
  expect_equal(as.character(parsed$msstatsptm.raw$PTM$ProteinName),
               as.character(parsed$msstatsptm.raw$PTM$SiteName))
})

test_that("converter-level rows run through dataSummarizationPTM with MBimpute", {
  skip_if_not_installed("MSstatsPTM")
  e <- .load_msstatsptm_helpers()
  ready <- get(".msstatsptmRawReady", envir = e)
  summarize.raw <- get(".summarizeMSstatsPTMRaw", envir = e)

  data_env <- new.env(parent = emptyenv())
  utils::data("raw.input", package = "MSstatsPTM", envir = data_env)
  raw <- get("raw.input", envir = data_env)

  # Keep the complete paired example. Arbitrary site subsets can legitimately
  # remove the protein features needed by TMP and make the test exercise an
  # invalid input rather than the application route.
  raw$PTM <- as.data.frame(raw$PTM)
  raw$PROTEIN <- as.data.frame(raw$PROTEIN)
  runs <- intersect(unique(as.character(raw$PTM$Run)), unique(as.character(raw$PROTEIN$Run)))
  run_meta <- unique(raw$PTM[, c("Run", "Condition", "BioReplicate")])
  run_meta <- run_meta[match(runs, run_meta$Run), ]
  groups <- setNames(as.character(run_meta$Condition), runs)
  subjects <- setNames(as.character(run_meta$BioReplicate), runs)

  expect_true(ready(raw, runs))
  summarized <- suppressWarnings(
    summarize.raw(raw, runs, groups, subjects, mbimpute = TRUE)
  )
  expect_type(summarized, "list")
  expect_named(summarized, c("PTM", "PROTEIN"))
  expect_gt(nrow(summarized$PTM$ProteinLevelData), 0)
  expect_gt(nrow(summarized$PROTEIN$ProteinLevelData), 0)
  expect_true(all(c("Protein", "LogIntensities", "NumImputedFeature") %in%
                    colnames(summarized$PTM$ProteinLevelData)))
  expect_identical(attr(summarized, "pa.log.transform"), "log2")
})
