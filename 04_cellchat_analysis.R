####################################################################################################
# Scripts: 04_cellchat_analysis.R
# Description: This script performs Cellchat analysis.
# Author: [Wa Du]
# Affiliation: [Department of Surgery, Division of Vascular Diseases and Surgery]
# Date: 2025-06-10
####################################################################################################

# Load necessary libraries
library(Seurat)
library(CellChat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(future)
library(magrittr)
##NOTE: Ensure you have properly annotated metadata (e.g., condition, cluster identities) in your Seurat object

# Set default assay
DefaultAssay(DOX_combo_filtered) <- "RNA"

# Examine metadata and split object by condition
table(Idents(DOX_combo_filtered), DOX_combo_filtered$condition)
obj.list <- SplitObject(DOX_combo_filtered, split.by = "condition")

# Set database
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
CellChatDB.use <- CellChatDB

# Function to perform CellChat analysis
run_cellchat_analysis <- function(seurat_obj, name_prefix) {
  cellchat <- createCellChat(seurat_obj, assay = "RNA")
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  future::plan("multicore", workers = 4)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  
  tiff(paste0(name_prefix, "_number_of_interactions.tiff"), units = "in", width = 12, height = 6, res = 300)
  par(mfrow = c(1, 2), xpd = TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE,
                   title.name = paste0(name_prefix, "_number of interactions"))
  dev.off()
  
  tiff(paste0(name_prefix, "_interaction_weight-strength.tiff"), units = "in", width = 12, height = 6, res = 300)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE,
                   title.name = paste0(name_prefix, "_interaction weights/strength"))
  dev.off()
  
  saveRDS(cellchat, file = paste0(name_prefix, "_cellchat.rds"))
  write.csv(cellchat@netP$pathways, file = paste0(name_prefix, "_cell_cell_communication.csv"))
  
  return(cellchat)
}

# Run CellChat analysis
cellchat1 <- run_cellchat_analysis(obj.list$WT_vehicle, "WT_vehicle")
cellchat2 <- run_cellchat_analysis(obj.list$Tg_vehicle, "Tg_vehicle")
cellchat3 <- run_cellchat_analysis(obj.list$WT_DOX, "WT_DOX")
cellchat4 <- run_cellchat_analysis(obj.list$Tg_DOX, "Tg_DOX")

# Chord plots for EC-CM interactions
plot_chord <- function(cellchat, prefix) {
  tiff(paste0(prefix, "_EC_as_source.tiff"), units = "in", width = 12, height = 12, res = 300)
  netVisual_chord_gene(cellchat, sources.use = c(1, 2, 4, 9, 13, 17), targets.use = c(6, 8, 10),
                       legend.pos.x = 4, lab.cex = 0.6)
  dev.off()
  
  tiff(paste0(prefix, "_EC_as_target.tiff"), units = "in", width = 12, height = 12, res = 300)
  netVisual_chord_gene(cellchat, sources.use = c(6, 8, 10), targets.use = c(1, 2, 4, 9, 13, 17),
                       legend.pos.x = 4, lab.cex = 0.6)
  dev.off()
}

plot_chord(cellchat1, "WT_vehicle_EC_CM")
plot_chord(cellchat2, "Tg_vehicle_EC_CM")
plot_chord(cellchat3, "WT_DOX_EC_CM")
plot_chord(cellchat4, "Tg_DOX_EC_CM")

# Compare interactions across groups
object.list <- list(WT_vehicle = cellchat1, Tg_vehicle = cellchat2, WT_DOX = cellchat3, Tg_DOX = cellchat4)
combo_cellchat <- mergeCellChat(object.list)

tiff("WT_vs_Tg_DOX_interactions.tiff", units = "in", width = 8, height = 4, res = 300)
compareInteractions(combo_cellchat, show.legend = FALSE, group = c(1, 2, 3, 4)) +
  compareInteractions(combo_cellchat, show.legend = FALSE, group = c(1, 2, 3, 4), measure = "weight")
dev.off()

# Bubble plots for comparisons
bubble_plot <- function(object.list, filename, comp, sources, targets) {
  combo_cellchat <- mergeCellChat(object.list)
  tiff(filename, units = "in", width = 6, height = 7, res = 300)
  print(netVisual_bubble(combo_cellchat, sources.use = sources, targets.use = targets,
                         comparison = comp, max.dataset = 2, angle.x = 45))
  dev.off()
}

bubble_plot(list(WT_vehicle = cellchat1, Tg_vehicle = cellchat2), "upregulated_interactions_WT_vs_Tg_vehicle.tiff", c(1, 2),
            c("Capillary EC#1", "Capillary EC#2", "Arterial EC", "Venous EC"), c("CM#1", "CM#2", "CM#3"))

bubble_plot(list(WT_vehicle = cellchat1, WT_DOX = cellchat3), "upregulated_interaction_WT_vehicle_DOX.tiff", c(1, 2),
            c("Capillary EC#1", "Capillary EC#2", "Arterial EC", "Venous EC"), c("CM#1", "CM#2", "CM#3"))

bubble_plot(list(WT_DOX = cellchat3, Tg_DOX = cellchat4), "upregulated_interaction_WT_vs_Tg_DOX.tiff", c(1, 2),
            c("Capillary EC#1", "Capillary EC#2", "Arterial EC", "Venous EC"), c("CM#1", "CM#2", "CM#3"))