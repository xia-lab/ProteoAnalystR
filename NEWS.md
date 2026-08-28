# ProteoAnalystR NEWS

## Version 1.0.0 — 2026-08 (manuscript release)

First versioned, installable release of the ProteoAnalyst analysis backend,
accompanying the *Protein Science* manuscript. Now a standard R package
with `DESCRIPTION`, `NAMESPACE`, roxygen docs under `man/`, and a `testthat`
suite. A citable Zenodo DOI will be minted for this release at publication.

### Readers & ingestion
- Native MaxQuant `evidence.txt` ingestion at peptide/precursor level;
  duplicate precursor/run matches resolved by maximum intensity.
- MSstats-long duplicate aggregation fixed: `reshape2::dcast()` had silently
  converted duplicate protein/run observations into feature counts; now uses
  median intensity.
- Spectronaut: peptide/protein groups no longer merged during reshaping;
  UniProt accessions such as `1/sp|P37108|SRP14_HUMAN` parsed correctly (pipes
  are not treated as protein-group delimiters).
- Paired canonical MSstatsPTM / PTM long-table ingestion (SiteName, parent
  ProteinName, run design, KGG sites, paired protein abundance); repeated
  site/run and protein/run values aggregated by median; protein references
  aligned and log2-converted when supplied on a linear scale.
- Parent-protein extraction preserves identifiers containing underscores.

### Statistics
- **Enrichment background defaults to the measured proteome/phosphoproteome**
  (`paramSet$universe.opt = "uploaded"`) for protein ORA, volcano ORA, and
  phospho pathway enrichment; the KEGG pathway-impact enrichment now reads the
  canonical `universe.opt` toggle instead of a stray `universeOpt`.
- Fold-wise covariate/batch removal for cross-validated biomarker analysis
  (`.RemoveCovariateEffectFold`), eliminating held-out leakage.

### Fixes
- `data_impute.R`: corrected KNN feature×sample orientation in `knn_var` and
  `knn_smp` (previously crashed with "invalid 'row.names' length").
- `data_utils_general.R` `replace_extension_with_qs()`: now maps `.xls`/`.xlsx`
  (as well as `.csv`/`.txt`/`.tsv`) to `.qs`, preventing a `.xls` upload from
  being overwritten by a `.qs` blob.
- `phospho_enrich_utils.R` `DetectPhosphoOccupancyBySite()`: builds a
  full-length, index-aligned residue vector (previously `regmatches()` dropped
  non-matching sites and misaligned residues, crashing on labelled KGG sites);
  handles MSstatsPTM-style ids without an underscore before the position and
  with a trailing `(heavy)`/`(light)` label.

### Packaging
- `Version: 1.0.0`, `License: MIT`, `Depends: R (>= 4.5)`.
- Pinned dependency environment documented in
  `reviewer_revision/reproducibility/DEPENDENCIES.md` (R 4.5.3 / Bioconductor 3.22).
