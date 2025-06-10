# ============================================================================
# Title: Automatic Cell Type Annotation with SingleR (Twigger Dataset)
# Author: elopezhuerta
# Description: This script performs reference-based cell annotation using 
#              SingleR and the HumanPrimaryCellAtlasData reference.
# ============================================================================

# --- Load Required Libraries ---
library(Seurat)         # version 4.3.0
library(celldex)        # version 1.0.0
library(SingleR)        # version 1.4.0
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)

# --- Input File: Post-QC Seurat object ---
integrated_seurat <- readRDS(file = "results/seurat_objects/TwiggerDataset/integrated_Twigger.rds")

# --- Load Annotation Reference (HPCA) ---
hpca.reference <- celldex::HumanPrimaryCellAtlasData()

# --- Extract normalized expression matrix ---
norm_expr <- GetAssayData(integrated_seurat, assay = "RNA", slot = "data")

# --- Filter out uninformative genes ---
bad_patterns <- "RP[0-9]|^MT|[0-9]{4}|^RPL|^RPS|[0-9]orf[0-9]|^ATP"
valid_genes <- rownames(norm_expr)[!grepl(bad_patterns, rownames(norm_expr))]
norm_expr <- norm_expr[valid_genes, ]

# --- Run SingleR ---
singleR_results <- SingleR(
  test = norm_expr,
  ref = hpca.reference,
  method = "single",
  assay.type.test = 1,
  labels = hpca.reference$label.main
)

# --- Remove low-confidence annotations (≤ 50 cells) ---
freq_table <- as.data.frame(table(singleR_results$pruned.labels))
low_conf <- freq_table %>%
  filter(Freq <= 50) %>%
  pull(Var1)

refined_labels <- ifelse(singleR_results$pruned.labels %in% low_conf | is.na(singleR_results$pruned.labels),
                         "unknown",
                         singleR_results$pruned.labels)

singleR_results$labels <- refined_labels

# --- Plot annotation heatmap ---
plotScoreHeatmap(singleR_results)

# --- Assign SingleR labels to Seurat metadata ---
label_df <- data.frame(Cell = rownames(singleR_results), Label = refined_labels)
cell_labels <- setNames(label_df$Label, label_df$Cell)

# Prefix to indicate automated annotation
integrated_seurat$Auto_Annotation <- paste0("SnR.", cell_labels)

# --- Remove 'unknown' labeled cells (optional step) ---
integrated_seurat <- subset(integrated_seurat, subset = Auto_Annotation != "SnR.unknown")

# --- Save integrated_Twigger Seurat object ---
saveRDS(integrated_seurat, file = "results/seurat_objects/TwiggerDataset/SnR.Annotated_Twigger.rds")

# ============================================================================
# Section: UMAP Visualization
# ============================================================================

# --- UMAP by Sample ---
umap_sample <- DimPlot(integrated_seurat, reduction = "umap", 
                       group.by = "orig.ident") +
  labs(title = "Merged Samples") +
  theme(legend.position = "bottom")

# --- UMAP by Clusters ---
umap_cluster <- DimPlot(integrated_seurat, reduction = "umap", 
                        group.by = "seurat_clusters", label = TRUE) +
  labs(title = "Clustering") +
  NoLegend()

# --- UMAP by Cell Types ---
umap_annotation <- DimPlot(integrated_seurat, reduction = "umap", 
                           group.by = "Auto_Annotation", label = TRUE, 
                           label.size = 4, repel = TRUE) +
  
  labs(title = "Automatic annotation") +
  NoLegend()

# --- Combine and Display ---
umap_sample + umap_cluster + umap_annotation

