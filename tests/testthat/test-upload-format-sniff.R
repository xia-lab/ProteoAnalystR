# Regression tests for wrong-tab detection at upload (.sniffToolFormat /
# .wrongTabMessage): a file exported by one tool but uploaded on another tool's
# tab must be reported by name instead of failing inside the reader.

R_DIR <- normalizePath(file.path(testthat::test_path(), "..", "..", "R"), mustWork = FALSE)
READTABLE <- file.path(R_DIR, "data_readtable_utils.R")

.load_readtable <- function() {
  skip_if_not(file.exists(READTABLE), "data_readtable_utils.R not found")
  e <- new.env(parent = globalenv())
  sys.source(READTABLE, envir = e)
  e
}

.write_header <- function(cols, eol = "\n") {
  f <- tempfile(fileext = ".tsv")
  row <- rep("1", length(cols))
  writeChar(paste0(paste(cols, collapse = "\t"), eol, paste(row, collapse = "\t"), eol),
            f, eos = NULL)
  f
}

HEADERS <- list(
  spectronaut = c("PG.ProteinAccessions", "PG.ProteinGroups", "PG.Genes",
                  "[1] run_A.htrms.PG.Quantity", "[2] run_B.htrms.PG.Quantity"),
  diann = c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
            "First.Protein.Description", "run_A.d", "run_B.d"),
  maxquant = c("Protein IDs", "Majority protein IDs", "LFQ intensity A", "LFQ intensity B"),
  fragpipe = c("Protein", "Protein ID", "Entry Name", "A MaxLFQ Intensity", "B MaxLFQ Intensity")
)

test_that("each tool's header is recognised as that tool only", {
  e <- .load_readtable()
  for (fmt in names(HEADERS)) {
    f <- .write_header(HEADERS[[fmt]])
    expect_identical(e$.sniffToolFormat(f)$formats, fmt, info = fmt)
    expect_null(e$.wrongTabMessage(f, fmt), info = fmt)
  }
})

test_that("a file on another tool's tab names the tool it came from", {
  e <- .load_readtable()
  f <- .write_header(HEADERS$spectronaut)
  msg <- e$.wrongTabMessage(f, "diann")
  expect_match(msg, "looks like a Spectronaut export")
  expect_match(msg, "uploaded on the DIA-NN tab")
  expect_match(msg, "Spectronaut Outputs tab")
})

test_that("CR-only line endings (old Mac exports) are still sniffed", {
  e <- .load_readtable()
  f <- .write_header(HEADERS$spectronaut, eol = "\r")
  expect_identical(e$.sniffToolFormat(f)$formats, "spectronaut")
})

test_that("unrecognised generic tables are not rejected", {
  e <- .load_readtable()
  f <- .write_header(c("#NAME", "sample1", "sample2"))
  expect_length(e$.sniffToolFormat(f)$formats, 0)
  expect_null(e$.wrongTabMessage(f, "diann"))
})
