####################################################################################################
# Script: 03_differential_expression.R
# Description: Perform cluster-wise differential gene expression analysis on endothelial cell (EC) subsets 
#              comparing WT_DOX vs Tg_DOX conditions using the MAST statistical test while controlling 
#              for batch effects.
# Author: Wa Du
# Affiliation: [Department of Surgery, Division of Vascular Diseases and Surgery]
# Date: 2025-06-10
####################################################################################################

# Load required packages
library(Seurat)
library(MAST)
library(dbplyr)

# Set default assay to RNA for differential expression testing
DefaultAssay(DOX_combo) <- "RNA"

# Define clusters and corresponding conditions for differential expression analysis
clusters <- c("0", "2", "5", "6")
conditions <- c("WT_DOX", "Tg_DOX")

# Loop through each cluster and perform differential expression analysis using MAST,
# controlling for batch effects. Results are saved as CSV files.
for (cluster_id in clusters) {
  
  ident_1 <- paste0(cluster_id, "_", conditions[1])
  ident_2 <- paste0(cluster_id, "_", conditions[2])
  
  de_results <- FindMarkers(
    object = DOX_combo,
    ident.1 = ident_1,
    ident.2 = ident_2,
    test.use = "MAST",
    latent.vars = "batch"
  )
  
  # Write DE results to CSV file with informative filename
  output_file <- paste0("cluster", cluster_id, "_EC_DOX_combo_DE.csv")
  write.csv(de_results, file = output_file)
}
