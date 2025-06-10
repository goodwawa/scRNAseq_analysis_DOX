####################################################################################################
# Script: 01_QC_and_filtering.R
# Description: This script performs Quality control, normalization, and doublet removal using DoubletFinder
#              for one example (WT_vehicle_01) of 18 samples.
# Author: [Wa Du]
# Affiliation: [Department of Surgery, Division of Vascular Diseases and Surgery]
# Date: 2025-06-10
####################################################################################################

# Load required packages
library(Seurat)
library(DoubletFinder)

# ----------------------------
# Load filtered 10X matrix
# ----------------------------
WT_vehicle_01 <- Read10X(data.dir = "~path/sample_filtered_feature_bc_matrix")

# Create Seurat object
WT_vehicle_01 <- CreateSeuratObject(counts = WT_vehicle_01)

# Add mitochondrial content as metadata
WT_vehicle_01[["percent.mt"]] <- PercentageFeatureSet(WT_vehicle_01, pattern = "^mt-")

# Filter out low-quality cells
WT_vehicle_01 <- subset(WT_vehicle_01, subset = percent.mt < 25 & nFeature_RNA > 200)

# ----------------------------
# Normalize and pre-process
# ----------------------------
WT_vehicle_01 <- NormalizeData(WT_vehicle_01)
WT_vehicle_01 <- FindVariableFeatures(WT_vehicle_01)
WT_vehicle_01 <- ScaleData(WT_vehicle_01)
WT_vehicle_01 <- RunPCA(WT_vehicle_01, npcs = 20)

# PCA inspection (interactive steps for user)
ElbowPlot(WT_vehicle_01, ndims = 20)
DimHeatmap(WT_vehicle_01, dims = 20, cells = 500, balanced = TRUE)

# ----------------------------
# DoubletFinder workflow
# ----------------------------

# Find optimal pK
sweep.res <- paramSweep(WT_vehicle_01, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Extract the optimal pK
best_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
write.csv(best_pK, file = "WT_vehicle_01_pK.csv")

# Estimate expected number of doublets (1% of cells)
nExp <- round(ncol(WT_vehicle_01) * 0.01)

# Run DoubletFinder
WT_vehicle_01 <- doubletFinder(
  seu = WT_vehicle_01,
  PCs = 1:20,
  pN = 0.25,
  pK = best_pK,
  nExp = nExp
)

# Check doublet classification
print(table(WT_vehicle_01@meta.data$DF.classifications_0.25_0.24_107))

# Remove predicted doublets
WT_vehicle_01 <- subset(WT_vehicle_01, subset = DF.classifications_0.25_0.24_107 == "Singlet")

# Save cleaned object
saveRDS(WT_vehicle_01, file = "cleaned_WT_vehicle_01.rds")
