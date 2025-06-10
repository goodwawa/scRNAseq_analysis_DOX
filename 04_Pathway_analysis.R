####################################################################################################
# Title: Pathway Enrichment Analysis of Endothelial Cell Subset (DOX vs Vehicle Comparison)
# Description: This script performs pathway enrichment analysis using hallmark gene sets and individual 
#              individual genesets for integrated scRNA-seq data.
# Author: Wa Du
# Affiliation: Department of Surgery, Division of Vascular Diseases and Surgery
# Date: 2025-06-10
# Required Packages: SCPA, Seurat, msigdbr, dplyr, ggplot2, ggrepel, future, patchwork
####################################################################################################

# Load required packages
library(Seurat)
library(SCPA)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(future)
library(patchwork)


# ============================= Setup and Data Loading =============================================
# Increase global memory limit for future parallelization
options(future.globals.maxSize = 100 * 1024^3)  # 100 GiB

# Load the integrated, SCTransformed Seurat object
DOX_combo <- readRDS("DOX_combo_integrated_SCTransformed.rds")

# Subset endothelial cell clusters of interest (based on known cluster IDs)
EC_ALL <- subset(DOX_combo, 
                 subset = rpca_clusters %in% c(0, 2, 5, 6) & 
                   condition %in% c("WT_vehicle", "Tg_vehicle", "WT_DOX", "Tg_DOX"))
saveRDS(EC_ALL, file = "EC_All_subset.rds")

# Reload to confirm integrity
EC_ALL <- readRDS("EC_All_subset.rds")
DefaultAssay(EC_ALL) <- "RNA"

# Extract condition-specific subsets
EC_ALL_WT_DOX <- seurat_extract(EC_ALL, meta1 = "condition", value_meta1 = "WT_DOX", assay = "RNA")
EC_ALL_Tg_DOX <- seurat_extract(EC_ALL, meta1 = "condition", value_meta1 = "Tg_DOX", assay = "RNA")

# ============================= Hallmark Pathway Enrichment =========================================
# Keyword-based filtering of MSigDB gene sets
library(msigdbr)
library(dplyr)
library(ggplot2)
library(ggrepel)

pathways <- c("hallmark")
mkr_sets <- msigdbr(species = "Mus musculus") %>%
  filter(grepl(paste(pathways, collapse = "|"), gs_name, ignore.case = TRUE)) %>%
  format_pathways()

# Confirm pathway names loaded
head(mkr_sets)

# Perform comparative pathway enrichment
EC_ALL_hallmark <- compare_pathways(samples = list(EC_ALL_WT_DOX, EC_ALL_Tg_DOX),
                                    pathways = mkr_sets,
                                    downsample = 2000)

# Color-code significant pathways
EC_ALL_hallmark <- EC_ALL_hallmark %>% 
  mutate(color = case_when(
    FC > 1.0 & adjPval < 0.05 ~ '#6dbf88',             # upregulated
    FC < -1.0 & adjPval < 0.05 ~ 'mediumseagreen',     # downregulated
    abs(FC) <= 1.0 & adjPval < 0.05 ~ '#84b0f0',        # moderate fold-change
    adjPval >= 0.05 ~ 'black'                          # not significant
  ))

# Save results
saveRDS(EC_ALL_hallmark, file = "EC_ALL_hallmark.rds")
write.csv(EC_ALL_hallmark, file = "WT_vs_Tg_DOX_EC_ALL_hallmark.csv")

# ============================= Visualization of Key Pathways =======================================
# Focus on selected hallmark pathways of interest
highlight_paths <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                     "HALLMARK_HYPOXIA",
                     "HALLMARK_P53_PATHWAY",
                     "HALLMARK_APOPTOSIS",
                     "HALLMARK_APICAL_JUNCTION")

aa_path <- EC_ALL_hallmark %>%
  filter(grepl(paste(highlight_paths, collapse = "|"), Pathway, ignore.case = TRUE))

aa_path_label <- aa_path %>% filter(Pathway %in% highlight_paths)
label_left <- aa_path_label %>% filter(FC < -1.0)
label_right <- aa_path_label %>% filter(FC > 1.0)

# Generate volcano-style dot plot
tiff("WT_vs_Tg_DOX_EC_ALL_hallmark.tiff", units = "in", width = 5, height = 5, res = 300, bg = "transparent")
ggplot(EC_ALL_hallmark, aes(x = FC, y = qval)) +
  geom_vline(xintercept = c(-1.0, 1.0), linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_point(shape = 21, size = 2.6, fill = EC_ALL_hallmark$color, stroke = 0.6) +
  geom_point(data = aa_path, shape = 21, size = 2.8, fill = "mediumseagreen", color = "black", stroke = 0.6) +
  geom_text_repel(data = label_left, aes(label = Pathway),
                  size = 2, color = "black", box.padding = 0.4, point.padding = 0.3,
                  segment.color = "mediumseagreen", segment.size = 0.5, max.iter = 10000,
                  nudge_x = -4, nudge_y = 0.3, force_pull = 0.1, force = 2) +
  geom_text_repel(data = label_right, aes(label = Pathway),
                  size = 2, color = "black", box.padding = 0.4, point.padding = 0.3,
                  segment.color = "mediumseagreen", segment.size = 0.5, max.iter = 10000,
                  nudge_x = 4, nudge_y = 0.3, force_pull = 0.1, force = 2) +
  xlim(-13, 13) +
  ylim(0, 6) +
  xlab("Enrichment (Fold Change)") +
  ylab("Q-value") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        aspect.ratio = 1)
dev.off()
####################################################################################################
# ============================= Individual selected Pathway Enrichment =========================================
##single geneset for pathway analysis
pathways<-c("HALLMARK_APICAL_JUNCTION")
mkr_sets<-msigdbr("Mus musculus") %>%  filter(grepl(paste(pathways, collapse = "|"), gs_name, ignore.case = T)) %>%
  format_pathways()
head(pathways)
EC_ALL_HALLMARK_APICAL_JUNCTION<-compare_pathways(samples = list(EC_ALL_WT_DOX,EC_ALL_Tg_DOX),pathways = mkr_sets,downsample = 2000)
EC_ALL_HALLMARK_APICAL_JUNCTION<-EC_ALL_HALLMARK_APICAL_JUNCTION %>% 
  mutate(color=case_when(FC > 1.0 & adjPval < 0.05 ~ '#6dbf88',
                         FC < 1.0 & FC > -1.0 & adjPval < 0.05 ~ '#84b0f0',
                         FC < -1.0 & adjPval < 0.05 ~ 'mediumseagreen',
                         FC < 1.0 & FC > -1.0 & adjPval > 0.05 ~ 'black'))

# Save results
saveRDS(EC_ALL_HALLMARK_APICAL_JUNCTION,file = "EC_ALL_HALLMARK_APICAL_JUNCTION.rds")
write.csv(EC_ALL_HALLMARK_APICAL_JUNCTION,file = "WT_vs_Tg_DOX_EC_ALL_HALLMARK_APICAL_JUNCTION.csv")

# Generate volcano-style dot plot
tiff("WT_vs_Tg_DOX_EC_ALL_HALLMARK_APICAL_JUNCTION.tiff",units = "in",width =5,height = 5,res = 300,bg="transparent")
aa_path <- EC_ALL_HALLMARK_APICAL_JUNCTION %>% 
  filter(grepl(pattern = "HALLMARK_APICAL_JUNCTION", ignore.case = T, x = Pathway))
# Select specific pathways to label
label_pathways <- c("HALLMARK_APICAL_JUNCTION")  # <-- Adjust names as needed
aa_path_label <- aa_path %>% filter(Pathway %in% label_pathways)
# Split into left and right
label_left <- aa_path_label %>% filter(FC < -1)
label_right <- aa_path_label %>% filter(FC > 1)
ggplot(EC_ALL_HALLMARK_APICAL_JUNCTION, aes(FC, qval)) +
  geom_vline(xintercept = c(-1.0, 1.0), linetype = "dashed", col = 'black', lwd = 0.6) +
  geom_point(cex = 2.6, shape = 21, fill = EC_ALL_HALLMARK_APICAL_JUNCTION$color, stroke = 0.6) +
  geom_point(data = aa_path, shape = 21, cex = 2.8, fill = "mediumseagreen", color = "black", stroke = 0.6) +
  geom_text_repel(data = label_left,
                  aes(label = Pathway),
                  size = 2,
                  color = "black",
                  box.padding = 0.4,
                  point.padding = 0.3,
                  segment.color = "mediumseagreen",
                  max.overlaps = Inf,
                  min.segment.length = 0,
                  segment.size = 0.5,
                  force_pull = 0.1,
                  force = 2,
                  max.iter = 10000,
                  nudge_x = -6,
                  nudge_y = 0.3) +
  geom_text_repel(data = label_right,
                  aes(label = Pathway),
                  size = 2,
                  color = "black",
                  box.padding = 0.4,
                  point.padding = 0.3,
                  segment.color = "mediumseagreen",
                  max.overlaps = Inf,
                  min.segment.length = 0,
                  segment.size = 0.5,
                  force_pull = 0.1,
                  force = 2,
                  max.iter = 10000,
                  nudge_x = +6,
                  nudge_y = 0.3) +
  xlim(-14,14) +
  ylim(0, 4) +
  xlab("Enrichment") +
  ylab("Qval") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        aspect.ratio = 1)
dev.off()
##################################################################################################
