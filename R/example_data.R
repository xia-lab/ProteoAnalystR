# =============================================================================
# On-demand access to the benchmark example datasets (MaxQuant / Spectronaut /
# MSstatsPTM). The raw inputs are large (the Spectronaut report alone is ~1 GB),
# so they are NOT bundled in the package. Instead the package ships a manifest of
# their public MassIVE locations (checksummed) and downloads them on request.
# See inst/extdata/README.md.
# =============================================================================

.example_manifest_path <- function() {
  p <- system.file("extdata", "benchmark_manifest.tsv", package = "ProteoAnalystR")
  if (!nzchar(p)) {
    # dev fallback: source tree layout when the package is not installed
    p <- file.path("inst", "extdata", "benchmark_manifest.tsv")
  }
  p
}

.read_example_manifest <- function() {
  p <- .example_manifest_path()
  if (!file.exists(p)) stop("Example manifest not found: ", p)
  utils::read.delim(p, stringsAsFactors = FALSE)
}

.massive_endpoint <- "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile"

#' List the bundled benchmark example datasets
#'
#' Returns a data frame describing the MaxQuant, Spectronaut, and MSstatsPTM
#' example datasets used in the ProteoAnalyst benchmarks. The data themselves are
#' public MassIVE deposits and are fetched on demand with
#' \code{\link{fetch_example_data}} (not bundled, due to size).
#'
#' @return A data frame with one row per file (dataset, role, format, accession,
#'   expected_bytes, sha256, relative_path).
#' @export
list_example_datasets <- function() {
  mf <- .read_example_manifest()
  mf$expected_bytes <- suppressWarnings(as.numeric(mf$expected_bytes))
  mf
}

#' Fetch a benchmark example dataset from its public repository
#'
#' Downloads the requested files for one example dataset from MassIVE into
#' \code{dest_dir}, verifying each file's size (with tolerance for the servlet's
#' CRLF->LF canonicalization of text files). The MaxQuant \code{evidence.txt}
#' (~143 MB) and especially the Spectronaut report (~1 GB) are large; a message
#' warns before large downloads.
#'
#' @param dataset One of "maxquant", "spectronaut", "ptm".
#' @param dest_dir Directory to download into (created if needed).
#' @param roles Which file roles to fetch. Default: the analysis inputs +
#'   metadata. Use "reference" to also get the MSstats(PTM) reference result, or
#'   \code{NULL}/"all" for everything.
#' @param verify Logical; verify downloaded file sizes against the manifest.
#' @param quiet Logical; suppress progress messages.
#' @return A named character vector of local file paths (names = roles).
#' @examples
#' \dontrun{
#'   ## MaxQuant example (evidence.txt + annotation):
#'   f <- fetch_example_data("maxquant", tempdir())
#'   ## then run it through the ProteoAnalystR pipeline exactly as the app does:
#'   ## ReadTabExpressData(f[["input"]], f[["metadata"]], ...)
#' }
#' @export
fetch_example_data <- function(dataset = c("maxquant", "spectronaut", "ptm"),
                               dest_dir = file.path(tempdir(), "ProteoAnalystR_examples"),
                               roles = c("input", "input_protein", "metadata"),
                               verify = TRUE, quiet = FALSE) {
  dataset <- match.arg(dataset)
  mf <- .read_example_manifest()
  sub <- mf[mf$dataset == dataset, , drop = FALSE]
  if (!is.null(roles) && !identical(roles, "all")) {
    sub <- sub[sub$role %in% roles, , drop = FALSE]
  }
  if (nrow(sub) == 0) stop("No files for dataset '", dataset, "' with roles: ",
                           paste(roles, collapse = ", "))
  total_mb <- round(sum(suppressWarnings(as.numeric(sub$expected_bytes)), na.rm = TRUE) / 1e6)
  if (!quiet && total_mb > 100) {
    message(sprintf("Note: '%s' download is ~%d MB. This can take a while.", dataset, total_mb))
  }
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

  out <- character(nrow(sub)); names(out) <- sub$role
  for (i in seq_len(nrow(sub))) {
    descriptor <- sub$file_descriptor[i]
    dest <- file.path(dest_dir, basename(sub$relative_path[i]))
    expected <- suppressWarnings(as.numeric(sub$expected_bytes[i]))

    if (file.exists(dest) && (!verify || .size_ok(dest, expected))) {
      if (!quiet) message("cached    ", basename(dest))
    } else {
      url <- paste0(.massive_endpoint, "?file=", utils::URLencode(descriptor, reserved = TRUE),
                    "&forceDownload=true")
      if (!quiet) message("fetching  ", basename(dest), " (", sub$role[i], ")")
      ok <- .download(url, dest)
      if (!ok) stop("Download failed for ", basename(dest), " (", descriptor, ")")
      if (verify && !.size_ok(dest, expected)) {
        stop("Size mismatch for ", basename(dest), ": expected ", expected,
             ", got ", file.info(dest)$size,
             ". Re-run, or verify manually against sha256 ", sub$sha256[i])
      }
    }
    out[i] <- dest
  }
  out
}

# size check with the same CRLF-canonicalization tolerance the shell fetcher uses
.size_ok <- function(path, expected) {
  if (is.na(expected)) return(TRUE)
  actual <- file.info(path)$size
  if (isTRUE(actual == expected)) return(TRUE)
  if (grepl("\\.(csv|txt|tsv|xml|R|Rmd|log|xls)$", path, ignore.case = TRUE)) {
    nlines <- tryCatch(length(readLines(path, warn = FALSE)), error = function(e) 0L)
    diff <- expected - actual
    return(diff >= 0 && diff <= nlines)
  }
  FALSE
}

.download <- function(url, dest) {
  tmp <- paste0(dest, ".part")
  ok <- FALSE
  if (requireNamespace("curl", quietly = TRUE)) {
    ok <- tryCatch({ curl::curl_download(url, tmp, quiet = TRUE); TRUE },
                   error = function(e) FALSE)
  }
  if (!ok) {
    ok <- tryCatch(utils::download.file(url, tmp, mode = "wb", quiet = TRUE) == 0L,
                   error = function(e) FALSE)
  }
  if (ok && file.exists(tmp)) { file.rename(tmp, dest); return(TRUE) }
  if (file.exists(tmp)) unlink(tmp)
  FALSE
}
