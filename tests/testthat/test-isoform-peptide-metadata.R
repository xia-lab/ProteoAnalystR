# Per-peptide isoform metadata (diann_metadata.qs) written by the MaxQuant
# peptides.txt and Spectronaut peptide readers: ids must equal the matrix row
# ids, accessions must equal the peptide-to-protein map, and no DIA-NN-only
# column may appear (they would switch on DIA-NN filters / modification fraction).

R_DIR <- normalizePath(file.path(testthat::test_path(), "..", "..", "R"), mustWork = FALSE)
READTABLE <- file.path(R_DIR, "data_readtable_utils.R")
ISO_COLS <- c("Protein.IDs", "Protein.Group", "Gene", "Proteotypic", "Stripped.Sequence")

.load_readtable <- function() {
  skip_if_not(file.exists(READTABLE), "data_readtable_utils.R not found")
  e <- new.env(parent = globalenv())
  sys.source(READTABLE, envir = e)
  e$ov_qs_save <- function(obj, file, ...) saveRDS(obj, file)
  e
}

.write_tsv <- function(df, sep = "\t") {
  f <- tempfile(fileext = ".txt")
  write.table(df, f, sep = sep, quote = FALSE, row.names = FALSE, na = "")
  f
}

.check_meta <- function(res) {
  meta <- readRDS("diann_metadata.qs")
  expect_identical(colnames(meta), ISO_COLS)
  expect_identical(meta$Protein.IDs, rownames(res$data))
  expect_identical(meta$Protein.Group, res$prot.map$Protein[match(meta$Protein.IDs, res$prot.map$Peptide)])
  meta
}

test_that("MaxQuant peptides.txt writes per-peptide isoform metadata", {
  e <- .load_readtable()
  withr::local_dir(withr::local_tempdir())
  pep <- data.frame(
    Sequence = c("AAAK", "CCCK", "DDDK", "EEEK", "FFFK", "FFFK", "GGGK"),
    Proteins = c("P1-2", "P1;P1-2", "CON__P9", "REV__P8", "Q2", "Q3", "Q4"),
    "Leading razor protein" = c("P1-2", "P1", "CON__P9", "REV__P8", "Q2", "Q3", "Q4"),
    "Gene names" = c("G1;G1B", "G1", "KRT1", "", "G2", "G3", "G4"),
    "Unique (Proteins)" = c("yes", "no", "yes", "yes", "yes", "yes", "yes"),
    PEP = 0.001,
    Reverse = c("", "", "", "+", "", "", ""),
    "Potential contaminant" = c("", "", "+", "", "", "", ""),
    "Intensity A1" = c(10, 20, 30, 40, 50, 60, 0),
    "Intensity B1" = c(11, 21, 31, 41, 51, 61, 0),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  res <- e$.readMaxQuantPeptides(.write_tsv(pep), list(quantType = "peptide", removeContaminants = TRUE, removeDecoys = TRUE))
  meta <- .check_meta(res)
  # contaminant, decoy and all-zero rows stay out; duplicate sequences get unique ids
  expect_identical(meta$Protein.IDs, c("AAAK", "CCCK", "FFFK", "FFFK.1"))
  expect_identical(meta$Proteotypic, c(1L, 0L, 1L, 1L))
  expect_identical(meta$Gene, c("G1", "G1", "G2", "G3"))
  expect_identical(meta$Stripped.Sequence, c("AAAK", "CCCK", "FFFK", "FFFK"))
  expect_false(any(grepl("CON__|REV__", meta$Protein.Group)))

  # older tables without "Unique (Proteins)": a single accession is proteotypic
  pep[["Unique (Proteins)"]] <- NULL
  res <- e$.readMaxQuantPeptides(.write_tsv(pep), list(quantType = "peptide", removeContaminants = TRUE, removeDecoys = TRUE))
  expect_identical(.check_meta(res)$Proteotypic, c(1L, 0L, 1L, 1L))
})

.spec_long <- function() {
  data.frame(
    R.FileName = rep(c("r1", "r2"), each = 5),
    PG.ProteinAccessions = rep(c("sp|P1-2|A_HUMAN", "P1;P1-2", "P1;P1-2", "Q5", "Q6"), 2),
    PG.Genes = rep(c("G1", "G1;G1B", "G1;G1B", "G5", "G6"), 2),
    PEP.StrippedSequence = rep(c("AAAK", "CCCK", "CCCK", "HHHK", "HHHK"), 2),
    PEP.IsProteotypic = rep(c("True", "False", "False", "False", "False"), 2),
    EG.ModifiedSequence = rep(c("_AAAK_", "_CCCK_", "_CCCK_", "_HHHK_", "_HHHK_"), 2),
    FG.Charge = rep(c(2, 2, 3, 2, 2), 2),
    FG.Quantity = c(10, 20, 22, 30, 31, 11, 21, 23, 32, 33),
    stringsAsFactors = FALSE
  )
}

test_that("Spectronaut long precursor report writes metadata keyed like the matrix rows", {
  e <- .load_readtable()
  withr::local_dir(withr::local_tempdir())
  res <- e$.readSpectronaut(.write_tsv(.spec_long()), list(inputType = "peptide", featureLevel = "precursor"))
  expect_identical(res$format, "spectronaut-peptide")
  meta <- .check_meta(res)
  # charge states collapse into one feature; a sequence under two protein groups gets two ids
  ord <- match(c("_AAAK_", "_CCCK_", "_HHHK_", "_HHHK_.1"), meta$Protein.IDs)
  expect_false(anyNA(ord))
  expect_identical(meta$Protein.Group[ord[1]], "P1-2")   # sp|..| form cleaned like prot.map
  expect_identical(meta$Proteotypic[ord], c(1L, 0L, 0L, 0L))
  expect_identical(meta$Gene[ord[2]], "G1")
  expect_identical(meta$Stripped.Sequence[ord], c("AAAK", "CCCK", "HHHK", "HHHK"))
})

test_that("Spectronaut peptide pivot report writes metadata; no PEP.IsProteotypic writes none", {
  e <- .load_readtable()
  withr::local_dir(withr::local_tempdir())
  piv <- data.frame(
    PG.ProteinAccessions = c("P1-2", "P1;P1-2", "Q7"),
    PG.Genes = c("G1", "G1", "G7"),
    PEP.StrippedSequence = c("AAAK", "CCCK", "MMMK"),
    # flagged, but a two-accession group cannot name the isoform: must come out 0
    PEP.IsProteotypic = c("True", "True", "True"),
    EG.ModifiedSequence = c("_AAAK_", "_CCCK_", "_M[Oxidation (M)]MMK_"),
    "[1] r1.PEP.Quantity" = c(10, 20, 30),
    "[2] r2.PEP.Quantity" = c(11, 21, 31),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  res <- e$.readSpectronaut(.write_tsv(piv), list(inputType = "peptide", featureLevel = "precursor"))
  meta <- .check_meta(res)
  expect_identical(meta$Proteotypic, c(1L, 0L, 1L))
  expect_identical(meta$Stripped.Sequence, c("AAAK", "CCCK", "MMMK"))

  unlink("diann_metadata.qs")
  piv$PEP.IsProteotypic <- NULL
  res <- e$.readSpectronaut(.write_tsv(piv), list(inputType = "peptide", featureLevel = "precursor"))
  expect_identical(res$format, "spectronaut-peptide")
  expect_false(file.exists("diann_metadata.qs"))
})

test_that("GetPeptideMetaFlags reads the metadata once and refreshes when it changes", {
  e <- .load_readtable()
  withr::local_dir(withr::local_tempdir())
  withr::defer(if (exists(".pa.meta.flags", envir = .GlobalEnv)) rm(".pa.meta.flags", envir = .GlobalEnv))
  reads <- 0L
  e$ov_qs_exists <- function(file) file.exists(file)
  e$ov_qs_read <- function(file, ...) { reads <<- reads + 1L; readRDS(file) }

  expect_identical(e$GetPeptideMetaFlags(), integer(5))   # nothing uploaded
  e$.savePeptideIsoformMeta(c("a", "b"), c("P1-2", "P1"), c("yes", "no"), stripped = c("AAK", "CCK"))
  saveRDS(data.frame(Peptide = c("a", "b"), Protein = c("P1-2", "P1")), "peptide_to_protein_map.qs")
  # flags: PEP, Peptide.Count, proteotypic filter, isoform usage, modification fraction
  expect_identical(e$GetPeptideMetaFlags(), c(0L, 0L, 1L, 1L, 0L))
  n <- reads
  e$GetPeptideMetaFlags(); e$GetPeptideMetaFlags()
  expect_identical(reads, n)   # cached until the files change

  saveRDS(data.frame(Proteotypic = 1L, Stripped.Sequence = "CK", Modified.Sequence = "C(UniMod:4)K", PEP = 0.01),
          "diann_metadata.qs")
  expect_identical(e$GetPeptideMetaFlags(), c(1L, 0L, 0L, 1L, 0L))   # fixed UniMod:4 alone: no modification fraction
  saveRDS(data.frame(Proteotypic = 1L, Stripped.Sequence = "MK", Modified.Sequence = "M(UniMod:35)K", PEP = 0.01),
          "diann_metadata.qs")
  expect_identical(e$GetPeptideMetaFlags(), c(1L, 0L, 0L, 1L, 1L))
})
