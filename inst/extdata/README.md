# ProteoAnalystR benchmark example datasets

These are the public controlled-mixture / spike-in datasets used to benchmark
ProteoAnalyst (MaxQuant DDA, Spectronaut DIA, MSstatsPTM ubiquitination). The raw
inputs are large — the Spectronaut report alone is ~1 GB — so they are **not
bundled** in the package. Instead this directory ships a checksummed manifest
(`benchmark_manifest.tsv`) of their permanent public locations, and the package
downloads them on demand.

## Datasets

| dataset | mode / tool | accession | ground truth |
|---|---|---|---|
| `maxquant`    | DDA / MaxQuant     | RMSV000000249.2 (Choi 2017, iPRG-2015) | 6 spike-ins at known fmol amounts |
| `spectronaut` | DIA / Spectronaut  | RMSV000000250.2 (Navarro 2016, LFQbench HYE124) | organism ratios (H 1:1, Y 2:1, E 1:4) |
| `ptm`         | PTM / MSstatsPTM   | MSV000088971 / RMSV000000669.1 | heavy KGG spike-in ubiquitination sites |

Each dataset lists an `input` (and, for PTM, an `input_protein`), a `metadata`
annotation, and a `reference` (the published MSstats / MSstatsPTM result, for
head-to-head comparison).

## Usage

```r
library(ProteoAnalystR)

# what's available
list_example_datasets()

# fetch the MaxQuant example (evidence.txt + annotation) into a temp dir
files <- fetch_example_data("maxquant", dest_dir = tempdir())
files[["input"]]      # -> .../Choi2017_DDA_MaxQuant_evidence.txt
files[["metadata"]]   # -> .../Choi2017_DDA_MaxQuant_annotation.csv

# also get the MSstats reference result for comparison
fetch_example_data("maxquant", roles = "all")

# Spectronaut is ~1 GB; fetch_example_data() warns before large downloads
fetch_example_data("spectronaut")
```

Downloaded files are size-verified against the manifest (with tolerance for the
MassIVE servlet's CRLF→LF text canonicalization). SHA-256 hashes are in the
manifest for manual verification.

## Reproducible full benchmarks

The complete benchmark suite (all datasets, scoring runners, and reports that
reproduce the manuscript numbers) lives in the application repository under
`reviewer_revision/benchmarks/` (`fetch_benchmarks.sh`, `MANIFEST.tsv`,
`CHECKSUMS.sha256`, `scoring/`). This package-level manifest is a self-contained
subset for quick, in-R examples.
