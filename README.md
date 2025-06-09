# scRNAseq_analysis_DOX Pipeline for [Endothelial transcription factor EB protects against doxorubicin-induced endothelial toxicity and cardiac dysfunction ]
This repository contains the R scripts used to perform single-cell RNA-seq analysis in our study "[Endothelial transcription factor EB protects against doxorubicin-induced endothelial toxicity and cardiac dysfunction]".
##Software & Dependencies
-R## Software & Dependencies
- R version 4.2.0
- Seurat v4.3.0
- tidyverse, clusterProfiler, msigdbr, etc.

## Analysis Overview
1. **01_QC_and_filtering.R** – Load raw data, perform quality control, filter low-quality cells.
2. **02_Normalization_and_clustering.R** – Normalize, find variable features, cluster.
3. **03_Integration.R** – Integration across samples using Seurat anchors.
4. **04_Differential_expression.R** – FindMarkers() for DEGs.
5. **05_Pathway_analysis.R** – GSEA or enrichment analysis.

## Notes
- Input files are expected in `data/` folder (not included).
- Data sharing available upon request due to IRB restrictions.
## Script Description
**01_QC_and_filtering.R**
Performs the following:
-Loads filtered 10x Genomics matrix
-Filters cells based on mitochondrial content and feature count
-Normalizes and scales the data
-Runs PCA for dimensionality reduction
-Uses DoubletFinder to score and remove doublets

##Requirements
-R version =
-Seurat v
-DoubletFinder v
Install packages (if not already installed):
