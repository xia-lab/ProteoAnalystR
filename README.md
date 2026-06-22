# ProteoAnalystR

**ProteoAnalystR** contains the R functions and libraries that power the
[ProteoAnalyst](https://www.proteoanalyst.ca) web server — a comprehensive
platform for the statistical, functional, and network-based analysis of
**proteomics** and **phosphoproteomics** data.

The package implements the complete analysis backend: data reading and
processing, normalization and imputation, differential and meta-analysis,
functional enrichment, biomarker exploration, and network-based
interpretation. Every analysis run on the web server is carried out by the
functions in this repository, which keeps the public tool and the underlying
R code fully reproducible.

> **Note:** ProteoAnalystR accompanies the ProteoAnalyst web application. It is
> not a standalone CRAN-style package — the code here is maintained as embedded
> runtime code that is synchronized with the web platform. APIs may change
> between releases.

## Installation

ProteoAnalystR is not on CRAN. To install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("xia-lab/ProteoAnalystR", build_vignettes = FALSE)
```

The package depends on a number of Bioconductor packages (e.g. `limma`,
`edgeR`, `DESeq2`, `DEqMS`, `sva`, `fgsea`, `WGCNA`). If `install_github`
fails to resolve these automatically, install them first with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "DESeq2", "DEqMS", "sva", "fgsea",
                       "preprocessCore", "WGCNA"))
```

## Functional modules

The R functions are organized by analysis stage. Key source files in `R/`:

| Area | Files | Purpose |
| --- | --- | --- |
| **Data import & processing** | `data_readtable_utils.R`, `data_idanot_utils.R`, `data_utils_general.R`, `data_maplist_utils.R`, `data_feattable_utils.R`, `utils_arrow.R` | Reading expression/abundance tables, ID annotation, feature handling |
| **Normalization & imputation** | `norm_utils.R`, `data_impute.R` | Sample/feature normalization and missing-value imputation |
| **Differential analysis** | `data_deanal_utils.R`, `covariate_utils.R`, `volcano_utils.R` | Differential expression (limma/edgeR/DESeq2/DEqMS), covariate adjustment, volcano plots |
| **Meta-analysis** | `data_utils_meta.R`, `metadata_utils.R`, `meta_methods.R`, `meta_proc.R` | Integration and statistical meta-analysis across multiple datasets |
| **Functional enrichment** | `enrich_utils.R`, `enrich_path_impact.R`, `utils_gsea.R`, `utils_enrichnet.R`, `_utils_regenrich.R` | ORA / GSEA, pathway impact, enrichment networks |
| **Phosphoproteomics** | `phospho_data_utils.R`, `phospho_qc_utils.R`, `phospho_enrich_utils.R`, `phospho_enrichment_utils.R`, `phospho_network_utils.R` | PTM-specific QC, enrichment, and kinase/network analysis |
| **Network analysis** | `net_lib_utils.R`, `ppi_network_utils.R`, `graph_utils_general.R`, `utils_coexp.R`, `utils_coexp_net.R`, `compartment_layout.R` | PPI queries, co-expression networks (WGCNA/CEMiTool), graph layout |
| **Dimensionality reduction** | `dim_red_utils.R`, `utils_scatter3d.R` | PCA / 2D & 3D ordination |
| **Biomarker analysis** | `biomarker_utils.R` | ROC-based biomarker and feature-selection workflows |
| **Visualization** | `heatmap_wrapper_utils.R`, `utils_heatmap_list.R`, `utils_heatmap_table.R`, `plot_violin_utils.R`, `qc_graphics.R`, `utils_ridgeline.R`, `upset_utils.R` | Heatmaps, violin/ridgeline plots, QC graphics, UpSet plots |
| **Helpers** | `helper_functions.R`, `misc_utils.R` | Shared utilities |

## Dependencies

ProteoAnalystR relies on a mix of statistical, visualization, dimensionality
reduction, network, and utility libraries:

- **Statistics & differential analysis:** `limma`, `edgeR`, `DESeq2`, `DEqMS`, `sva`, `fgsea`, `caret`, `preprocessCore`
- **Data manipulation:** `dplyr`, `tidyr`, `tibble`, `readr`, `reshape`, `reshape2`, `pryr`, `RSQLite`
- **Visualization:** `ggplot2`, `plotly`, `ggrepel`, `ggpubr`, `ggridges`, `Cairo`, `RColorBrewer`, `gridExtra`, `lattice`, `png`, `see`
- **Network & graph analysis:** `igraph`, `graphlayouts`, `WGCNA`, `CEMiTool`
- **3D & geometry / JSON:** `rgl`, `alphashape3d`, `ks`, `jsonlite`, `RJSONIO`, `rjson`

The exact packages used depend on the module and workflow being executed.

## Citation

If you use ProteoAnalyst or ProteoAnalystR in your work, please cite the
ProteoAnalyst web server (manuscript in preparation). For the most current
citation, see [https://www.proteoanalyst.ca](https://www.proteoanalyst.ca).

## Bugs & feedback

- Report issues on the [GitHub issue tracker](https://github.com/xia-lab/ProteoAnalystR/issues)
- Ask questions on the [OmicsForum](https://omicsforum.ca/c/proteoanalyst/22)

## License

ProteoAnalystR is released under the [MIT License](LICENSE).
Copyright (c) 2026 Jeff Xia.

---

> **ProteoAnalystR is under active development — we cannot guarantee full functionality.**
