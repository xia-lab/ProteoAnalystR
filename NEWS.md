# ProteoAnalystR NEWS

## Unreleased

### Statistics
- `DetectPTMOccupancy()`: the per-pair Welch t-test on logit occupancies is
  replaced by a limma-moderated test (lmFit + eBayes robust across all
  mod/unmod pairs; Welch fallback if limma is unavailable). The occupancy
  estimator is unchanged. On the semi-synthetic occupancy mixture, scored
  through the production end-to-end run reported in Supporting Information
  Table S11, this raises BH-FDR power at |delta occupancy| = 0.1/0.2 from
  0.21/0.54 to 0.54/0.93 and discrimination AUC from 0.958 to 0.972, while
  keeping null FPR (BH 0.014) and abundance-confounded calls (BH 0.019) under
  the nominal 5%. The parallel total-abundance test is moderated the same way.

### Reproducibility & examples
- Bundled two small, self-contained example datasets in the package
  (`inst/extdata/bundled/`), accessible offline via `bundled_example()` and
  described by `list_bundled_datasets()`: `fragpipe_proteomics` (FragPipe LFQ
  protein matrix, human, IDHmut vs IDHwt) and `diann_phospho` (DIA-NN
  phosphosite report, human, CTL vs DRUG). Two modalities, two upstream tools.
- New vignette `replay-example-datasets` walks a user from a bundled input
  through reader → normalization → differential expression →
  `analysis_summary.txt`, so the package can be installed and the analyses
  replayed locally.

### Analysis summary report
- `WriteAnalysisSummary()` now emits dedicated **Cellular-Compartment
  Enrichment** and **Motif / Sequence-Context Enrichment** sections (previously
  only Kinase enrichment was reported), and the Kinase section now states the
  number of significant phosphosites tested, the DE-significance definition, the
  measured background size, and the enrichment FDR cutoff — matching what the
  enrichment tabs show in the UI.

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
- Paired converter-level PTM and protein rows are now retained for the optional
  native MSstatsPTM engine, which runs `dataSummarizationPTM` (log2 transform,
  equalize-medians normalization, TMP with MBimpute) before
  `groupComparisonPTM`. Already-summarized uploads are detected
  explicitly and use the existing matrix fallback; result exports identify the
  path and the applied summarization settings.
- Parent-protein extraction preserves identifiers containing underscores.

### Statistics
- **Enrichment background defaults to the measured proteome/phosphoproteome**
  (`paramSet$universe.opt = "uploaded"`) for protein ORA, volcano ORA, and
  phospho pathway enrichment; the KEGG pathway-impact enrichment now reads the
  canonical `universe.opt` toggle instead of a stray `universeOpt`.
- Fold-wise covariate/batch removal for cross-validated biomarker analysis
  (`.RemoveCovariateEffectFold`), eliminating held-out leakage.

### Fixes
- `data_impute.R`: phosphoproteomics `min` imputation now uses a reproducible,
  PhosR-style site/condition draw from `N(group mean, group SD)` for sporadic
  within-group gaps, while retaining scale-aware LoD fills for censored groups.
  This avoids the variance deflation of deterministic group-mean replacement.
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
