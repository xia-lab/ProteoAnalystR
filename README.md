# ProteoAnalystR

### Note - ProteoAnalystR is under active development - we cannot guarantee full functionality ###

## Description

**ProteoAnalystR** contains the R functions and libraries underlying the ProteoAnalyst web server for proteomics and phosphoproteomics data analysis. It provides the statistical, functional, and visual analysis routines used by the platform for tasks such as data processing, differential analysis, enrichment analysis, biomarker exploration, and network-based interpretation.

ProteoAnalystR is intended to accompany the ProteoAnalyst web application. It is not a standalone CRAN-style package. The code in this directory is maintained as embedded runtime analysis code and synchronized with the web platform.

## R libraries used in ProteoAnalystR

ProteoAnalystR relies on a mix of statistical, visualization, dimensionality reduction, network, and utility libraries. Commonly used packages in the current codebase include:

- Core statistics and differential analysis: `limma`, `edgeR`, `DESeq2`, `DEqMS`, `sva`, `metagenomeSeq`
- Data manipulation and utilities: `dplyr`, `tidyr`, `tibble`, `readr`, `pryr`, `RSQLite`
- Visualization: `ggplot2`, `plotly`, `ggrepel`, `ggpubr`, `ggridges`, `Cairo`, `RColorBrewer`, `gridExtra`, `lattice`, `png`, `see`
- Network and graph analysis: `igraph`, `graphlayouts`, `WGCNA`, `CEMiTool`
- Enrichment and dimensionality reduction: `fgsea`, `uwot`, `Rtsne`
- Other supporting packages used in specific workflows: `caret`, `preprocessCore`, `reshape`, `reshape2`, `jsonlite`, `RJSONIO`, `rjson`, `rgl`, `alphashape3d`, `ks`

The exact packages used depend on the module and workflow being executed.

## Note for Developers

Keep this documentation synchronized with the corresponding JSF pages, Java calls, and R workflow behavior in the main ProteoAnalyst application.
