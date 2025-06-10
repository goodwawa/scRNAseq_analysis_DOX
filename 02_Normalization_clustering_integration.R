####################################################################################################
# Scripts: 02_Normalization_clustering_integration.R
# Description: This script performs normalization, variable feature selection, scaling, PCA,
#              clustering, and RPCA-based integration across 18 scRNA-seq samples using Seurat.
# Author: [Wa Du]
# Affiliation: [Department of Surgery, Division of Vascular Diseases and Surgery]
# Date: 2025-06-10
####################################################################################################


# Load necessary libraries
library(Seurat)
library(tidyverse)
library(future)

# Set memory limit for integration
options(future.globals.maxSize = 100 * 1024^3)  # 100 GiB

#--------------------------#
# 1. Load Pre-processed Data
#--------------------------#
# Each sample is a cleaned Seurat object (e.g., after filtering out low-quality cells)
samples <- list(
  WT_vehicle_01 = readRDS("cleaned_WT_vehicle_01.rds"),
  WT_vehicle_02 = readRDS("cleaned_WT_vehicle_02.rds"),
  WT_vehicle_03 = readRDS("cleaned_WT_vehicle_03.rds"),
  WT_vehicle_04 = readRDS("cleaned_WT_vehicle_04.rds"),
  Tg_vehicle_01 = readRDS("cleaned_Tg_vehicle_01.rds"),
  Tg_vehicle_02 = readRDS("cleaned_Tg_vehicle_02.rds"),
  Tg_vehicle_03 = readRDS("cleaned_Tg_vehicle_03.rds"),
  Tg_vehicle_04 = readRDS("cleaned_Tg_vehicle_04.rds"),
  Tg_vehicle_05 = readRDS("cleaned_Tg_vehicle_05.rds"),
  WT_DOX_01     = readRDS("cleaned_WT_DOX_01.rds"),
  WT_DOX_02     = readRDS("cleaned_WT_DOX_02.rds"),
  WT_DOX_03     = readRDS("cleaned_WT_DOX_03.rds"),
  WT_DOX_04     = readRDS("cleaned_WT_DOX_04.rds"),
  Tg_DOX_01     = readRDS("cleaned_Tg_DOX_01.rds"),
  Tg_DOX_02     = readRDS("cleaned_Tg_DOX_02.rds"),
  Tg_DOX_03     = readRDS("cleaned_Tg_DOX_03.rds"),
  Tg_DOX_04     = readRDS("cleaned_Tg_DOX_04.rds"),
  Tg_DOX_05     = readRDS("cleaned_Tg_DOX_05.rds")
)

# Assign condition metadata
for (name in names(samples)) {
  if (grepl("WT_vehicle", name)) {
    samples[[name]]$condition <- "WT_vehicle"
  } else if (grepl("Tg_vehicle", name)) {
    samples[[name]]$condition <- "Tg_vehicle"
  } else if (grepl("WT_DOX", name)) {
    samples[[name]]$condition <- "WT_DOX"
  } else if (grepl("Tg_DOX", name)) {
    samples[[name]]$condition <- "Tg_DOX"
  }
}

#----------------------------------#
# 2. Remove Mitochondrial Genes
#----------------------------------#
for (i in seq_along(samples)) {
  mito_genes <- grep("^mt-", rownames(samples[[i]]), value = TRUE)
  samples[[i]] <- samples[[i]][!rownames(samples[[i]]) %in% mito_genes, ]
}

#------------------------------#
# 3. Merge Samples
#------------------------------#
DOX_combo <- merge(samples$WT_vehicle_01, 
                   y = samples[names(samples) != "WT_vehicle_01"], 
                   add.cell.ids = names(samples)[names(samples) != "WT_vehicle_01"])

# Confirm conditions
table(DOX_combo$condition)

#------------------------------#
# 4. Pre-integration Processing
#------------------------------#
DOX_combo <- NormalizeData(DOX_combo)
DOX_combo <- FindVariableFeatures(DOX_combo)
DOX_combo <- ScaleData(DOX_combo)
DOX_combo <- RunPCA(DOX_combo)

# Pre-integration clustering (optional visualization)
DOX_combo <- FindNeighbors(DOX_combo, dims = 1:30)
DOX_combo <- FindClusters(DOX_combo, resolution = 0.3, cluster.name = "unintegrated_clusters")
DOX_combo <- RunUMAP(DOX_combo, dims = 1:30, reduction.name = "umap.unintegrated")

# Save UMAP of unintegrated data
tiff("WT_vs_Tg_DOX_unintegrated.tiff", units = "in", width = 8, height = 6, res = 300)
DimPlot(DOX_combo, reduction = "umap.unintegrated", group.by = "condition")
dev.off()

#------------------------------#
# 5. Integration using RPCA
#------------------------------#
DefaultAssay(DOX_combo) <- "RNA"

DOX_combo <- IntegrateLayers(
  object = DOX_combo,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  k.anchor = 5,
  k.weight = 50,
  sd.weight = 1,
  reference = which.max(table(DOX_combo$condition))
)

#------------------------------#
# 6. Post-integration Clustering
#------------------------------#
DOX_combo <- FindNeighbors(DOX_combo, reduction = "integrated.rpca", dims = 1:20)
DOX_combo <- FindClusters(DOX_combo, resolution = 0.28, cluster.name = "rpca_clusters")
DOX_combo <- RunUMAP(DOX_combo, reduction = "integrated.rpca", dims = 1:20, reduction.name = "umap.rpca")

# Final check
table(DOX_combo$condition)

# Save final integrated object
saveRDS(DOX_combo, "DOX_combo_integrated.rds")

# Save integrated UMAP
tiff("WT_vs_Tg_DOX_integrated.tiff", units = "in", width = 8, height = 6, res = 300)
DimPlot(DOX_combo, reduction = "umap.rpca", group.by = "condition")
dev.off()
