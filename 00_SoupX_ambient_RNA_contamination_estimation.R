####################################################################################################
# Script: 00_SoupX_ambient_RNA_contamination_estimation.R
# Description: This script performs ambient RNA contamination estimation
#              for one example (WT_vehicle_01) of 18 samples.
# Author: [Wa Du]
# Affiliation: [Department of Surgery, Division of Vascular Diseases and Surgery]
# Date: 2025-10-22
####################################################################################################
# Load required packages
library(Seurat)
library(SoupX)

## Load the raw and filtered matrix ##############################################################

raw_counts <- Read10X("~path/sample_raw_feature_bc_matrix")
filtered_counts <- Read10X("~path/sample_filtered_feature_bc_matrix")

# Ensure genes match
common_genes <- intersect(rownames(raw_counts), rownames(filtered_counts))
raw_counts <- raw_counts[common_genes, ]
filtered_counts <- filtered_counts[common_genes, ]

# Make sure the gene order is identical
filtered_counts <- filtered_counts[rownames(raw_counts), ]

# Create SoupChannel object
sc <- SoupChannel(tod = raw_counts, toc = filtered_counts)

# Create and preprocess Seurat object
seurat_obj <- CreateSeuratObject(counts = filtered_counts)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)

# Optional: Check elbow plot
ElbowPlot(seurat_obj, ndims = 50)

# Clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.3)

# Assign clusters to SoupX
sc <- setClusters(sc, seurat_obj$seurat_clusters)

# Estimate contamination (gentle correction range)
sc <- autoEstCont(sc, contaminationRange = c(0.03, 0.15), forceAccept = TRUE)

# Inspect rho distribution
summary(sc$metaData$rho)
hist(sc$metaData$rho, breaks = 50)

# --- Get corrected counts and create final Seurat object ---
corrected_counts <- adjustCounts(sc)
WT_vehicle_01 <- CreateSeuratObject(counts = corrected_counts)
WT_vehicle_01

# Save processed Seurat object
saveRDS(WT_vehicle_01, file = "SoupX_WT_vehicle_01.rds")

##############################################################################################################################
