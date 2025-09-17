# Repository Catalogue
This repository contains code used for Nanostring GeoMx analysis. 


# DSP analysis and visualization

This Script contains the R code and workflow used for a manuscript (To be Updated).

---

## 📦 Requirements

- R version: ≥ 4.2.0
- Recommended environment management: [renv](https://rstudio.github.io/renv/)
- Required R packages:
  - tidyverse
  - readxl
  - openxlsx
  - DESeq2
  - clusterProfiler
  - org.Hs.eg.db
  - enrichplot
  - EnhancedVolcano
  - pheatmap
  - ggrepel
  - sessioninfo

To install missing packages:

```r
install.packages(c("tidyverse", "readxl", "openxlsx", "ggrepel"))
BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db",
                       "enrichplot", "EnhancedVolcano", "pheatmap"))
