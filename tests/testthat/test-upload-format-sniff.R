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

# Spectronaut pivot reports: per-run quantity columns with "Filtered" cells must
# all be kept as runs (with Filtered -> NA), and an empty PG.Genes is not a run.
test_that("Spectronaut pivot keeps every run column and reads Filtered as NA", {
  e <- .load_readtable()
  dat <- data.frame(
    PG.ProteinAccessions = c("P1", "P2", "P3"),
    PG.Genes = c(NA, NA, NA),
    "[1] run_A.htrms.PG.Quantity" = c("Filtered", "10.5", "20"),
    "[2] run_B.htrms.PG.Quantity" = c("1", "NaN", "3"),
    "[3] run_C.htrms.PG.Quantity" = c("4", "5", "6"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  m <- e$.specIntensityMatrix(dat, annot.cols = c("PG.ProteinAccessions"))
  expect_identical(colnames(m), c("run_A", "run_B", "run_C"))
  expect_true(is.na(m[1, "run_A"]))
  expect_true(is.na(m[2, "run_B"]))
  expect_equal(unname(m[2, "run_A"]), 10.5)
})

test_that("generic numeric detection ignores Filtered tokens and empty columns", {
  e <- .load_readtable()
  dat <- data.frame(ID = c("P1", "P2"), Empty = c(NA, NA),
                    s1 = c("Filtered", "2"), s2 = c("3", "4"),
                    Note = c("x", "y"), stringsAsFactors = FALSE)
  m <- e$.specIntensityMatrix(dat, annot.cols = "ID")
  expect_identical(colnames(m), c("s1", "s2"))
})

test_that("run-name cleaning keeps the original names when it would collide", {
  e <- .load_readtable()
  expect_identical(e$.cleanSpectronautRunName("[12] F1RS_TB_BR_1_10.htrms.PG.Quantity"), "F1RS_TB_BR_1_10")
  expect_identical(e$.cleanSpectronautRunName("[1] x.raw.PEP.Quantity"), "x")
  clash <- c("[1] a.PG.Quantity", "[2] a.PG.Quantity")
  expect_identical(e$.cleanSpectronautRunName(clash), clash)
})
