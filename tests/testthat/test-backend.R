# Lightweight, dependency-free regression tests for the ProteoAnalyst backend.
R_DIR <- normalizePath(file.path(testthat::test_path(), "..", "..", "R"), mustWork = FALSE)

test_that("every backend R file parses without a syntax error", {
  skip_if_not(dir.exists(R_DIR), "R/ source dir not found")
  files <- list.files(R_DIR, pattern = "\\.R$", full.names = TRUE)
  expect_gt(length(files), 0)
  for (f in files) expect_silent(parse(f))
})

test_that("replace_extension_with_qs maps tabular uploads (incl .xls) to .qs", {
  skip_if_not(file.exists(file.path(R_DIR, "data_utils_general.R")))
  e <- new.env(); sys.source(file.path(R_DIR, "data_utils_general.R"), envir = e)
  f <- get("replace_extension_with_qs", envir = e)
  expect_equal(f("evidence.txt"), "evidence.qs")
  expect_equal(f("data.csv"),     "data.qs")
  expect_equal(f("report.tsv"),   "report.qs")
  # regression: .xls / .xlsx must be remapped so ov_qs_save() never overwrites the upload
  expect_equal(f("Spectronaut_input.xls"),  "Spectronaut_input.qs")
  expect_equal(f("book.xlsx"),               "book.qs")
})

test_that("PTM site id parsing handles labelled / no-underscore forms (no crash)", {
  # Mirrors the parent/residue extraction fixed in phospho_enrich_utils.R.
  site.ids <- c("A0AVT1_K_490", "O43175_K351 (heavy)", "Protein_1_S_1")
  prot <- sub("_[A-Z]_?[0-9]+( *\\((?:heavy|light)\\))?$", "", site.ids,
              ignore.case = TRUE, perl = TRUE)
  m <- regexpr("[A-Z](?=_?[0-9]+( *\\((?:heavy|light)\\))?$)", site.ids,
               perl = TRUE, ignore.case = TRUE)
  res <- rep("", length(site.ids)); res[m > 0] <- toupper(regmatches(site.ids, m))
  expect_equal(prot, c("A0AVT1", "O43175", "Protein_1"))
  expect_equal(res,  c("K", "K", "S"))
  expect_length(res, length(site.ids))   # aligned: the crash was a length mismatch
})
