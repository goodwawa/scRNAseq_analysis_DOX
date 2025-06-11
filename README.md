##scRNAseq_analysis_DOX Pipeline 
For manuscript: Endothelial transcription factor EB protects against doxorubicin-induced endothelial toxicity and cardiac dysfunction
This repository contains the R scripts used to perform the single-cell-RNA-seq (scRNA-seq) analysis in our study investigating the protective role of endothelial TFEB against doxorubicin-induced vascular and cardiac damage.

##Software & Dependencies
Core Requirements
- R version 4.5.0
- Seurat v5.3.0
- DoubletFinder v2.0.6
- msigdbr v24.1.0
- CellChat v2.2.0
- SCPA v1.6.2
- MAST v1.33.0
- tidyverse (includes ggplot2, dplyr, etc.)
- ggrepel
- future

## Analysis Workflow
1. **01_QC_and_filtering.R** – Load raw single-cell RNA-seq data (10X Genomics format); perform quality control based on mitochondrial content and gene counts; identifies and removes doublets using DoubletFinder.
2. **02_Normalization_clustering_integration.R** – Normalize, identifies highly variable features and performs dimensionality reduction (PCA, UMAP), cluster cells, and individual samples integration.
3. **03_Differential_expression.R** – FindMarkers() for DEGs.
4. **04_Cellchat_analysis.R** – CellChat analysis for endothelial cell-cardiomyocyte communications under the vehicle and DOX treated mouse cardiac tissue.
5. **05_Pathway_analysis.R** – SCPA pathway enrichment analysis using msigdbr hallmark and individual genesets.

## Notes
- Input files are expected in 'data/' folder (not included)
-  The raw and processed scRNA-seq data are currently under peer review. Upon acceptance, the dataset will be deposited and publicly available via the NCBI GEO database.
