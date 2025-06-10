# =============================================================================
# Title: Annotation Refinement (Twigger Dataset)
# Author: elopezhuerta
# Description: This script refines automatic cell-type annotations in 
#              a Seurat object using known marker genes and differential 
#              gene expression (DEG) analysis.
# =============================================================================

library(Seurat) # version 4.3.0
library(ggplot2)
library(dplyr)
library(plyr)

# Load previously annotated Seurat object (from SingleR)
annotated_seurat <- readRDS(file = "results/seurat_objects/WuDataset/SnR.Annotated_Wu.rds")

# Create a new metadata column to store refined annotations
annotated_seurat$Annotation <- annotated_seurat$Auto_Annotation

# Define a color palette for use in FeaturePlots
reds_palette <- RColorBrewer::brewer.pal(9, "YlOrRd")

# List of canonical marker genes to assist in manual annotation refinement
all_markers <- c(
  # T cells
  "CD3E", "CD3D", "TRAC",
  # CD8+ T cells
  "CD8A", "CD8B", "GZMA",
  # CD4+ T cells
  "IL7R", "CD4", "CCR7",
  # T regulatory cells
  "CTLA4", "FOXP3",
  # NK cells
  "NKG7", "KLRD1",
  # B cells
  "MS4A1", "CD19",
  # Plasma B cells
  "IGHG4", "IGHG1", "CD38", "SSR4", "IGHA1",
  # Myeloid lineage
  "LYZ",
  # Macrophages / Monocytes
  "CD68", "CD163", "CD14", "FCGR3A",
  # cDC1
  "CPVL", "CLEC9A", "LGALS2",
  # cDC2
  "CD1C", "FCER1A", "CLEC10A",
  # Mature DCs
  "CCR7", "LAMP3", "BIRC3",
  # Plasmacytoid DCs
  "CXCR3", "LILRA4", "JCHAIN",
  # Neutrophils
  "FCGR3B", "CXCR2",
  # Mast cells
  "TPSAB1", "CPA3",
  # Fibroblasts
  "COL1A1", "DCN",
  # Pericytes
  "RGS5", "ACTA2",
  # Epithelial cells
  "EPCAM", "KRT18",
  # Endothelial cells
  "PECAM1", "VWF",
  # Platelets
  "TUBB1", "PPBP"
)

# Calculate DEG to refine annotation manually
# Identify marker genes for each cluster (positive markers only)
cluster_markers <- FindAllMarkers(
  object = annotated_seurat,
  assay = "RNA",
  slot = "data",
  logfc.threshold = 1,
  min.pct = 0.3,
  only.pos = TRUE
)

# Select top 15 genes per cluster with log2FC > 2
top15_markers <- cluster_markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 2) %>%
  slice_head(n = 15) %>%
  ungroup()


# Calculate 75th pecentile upper quantile) of marker expression per cluster
# Split Seurat object by cluster identity
cluster_list <- SplitObject(annotated_seurat, split.by = "seurat_clusters")

# Initialize list to store per-cluster quantiles
quantile_list <- list()

for (ii in seq_along(cluster_list)) {
  # Get normalized expression matrix
  norm_matrix <- GetAssayData(cluster_list[[ii]], assay = "RNA", slot = "data")
  
  # Keep only marker genes that exist in current matrix
  current_markers <- all_markers[all_markers %in% rownames(norm_matrix)]
  
  # Compute the 75th percentile of expression per marker
  quantile_list[[ii]] <- apply(norm_matrix[current_markers, , drop = FALSE], 1, function(x) quantile(x, probs = 0.75))
}

# Assign cluster names to list elements
names(quantile_list) <- names(cluster_list)

# Combine into a single data frame (clusters in rows, genes in columns)
quantile_df <- plyr::ldply(quantile_list, rbind)

# ============================================================================
# Section 1: Manual Refinement of Lymphoid Cell Annotations
# ============================================================================

# --- Define canonical lymphoid markers ---
lymphoid_markers <- c(
  # T cells
  "CD3E", "CD3D", "TRAC",
  # CD8+ T cells
  "CD8A", "CD8B", "GZMA",
  # CD4+ T cells
  "IL7R", "CD4", "CCR7",
  # T regs
  "CTLA4", "FOXP3",
  # NK cells
  "NKG7", "KLRD1",
  # B cells
  "MS4A1", "CD19",
  # Plasma B cells
  "IGHG4", "IGHG1", "CD38", "SSR4", "IGHA1"
)

# --- Identify potential lymphoid clusters based on marker expression ---
likely_lymphoid <- unique(as.character(top15_markers$cluster[top15_markers$gene %in% lymphoid_markers]))
print(likely_lymphoid)
# --- Quick UMAP overview ---
p_annot <- DimPlot(annotated_seurat, group.by = "Auto_Annotation")
p_clust <- DimPlot(annotated_seurat, label = TRUE) + NoLegend()
p_annot + p_clust

# --- Create violin plots for each lymphoid marker ordered by expression ---
present_lymph_markers <- lymphoid_markers[lymphoid_markers %in% colnames(quantile_df)]

vln_plots_lymph <- lapply(present_lymph_markers, function(marker) {
  cluster_order <- quantile_df[, 1][order(quantile_df[[marker]], decreasing = TRUE)]
  VlnPlot(annotated_seurat, features = marker, pt.size = 0) +
    scale_x_discrete(limits = cluster_order) +
    labs(title = marker) +
    theme(
      axis.text.x = element_text(size = rel(0.5)),
      axis.text.y = element_text(size = rel(0.5)),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(size = rel(1), family = "serif")
    ) +
    NoLegend()
})
gridExtra::grid.arrange(grobs = vln_plots_lymph, ncol = 5)

# --- Extract lymphoid population for refinement ---
lymphoid_cells <- subset(
  annotated_seurat,
  subset = Annotation %in% c("SnR.T_cells", "SnR.B_cell", "SnR.NK_cell") &
    seurat_clusters %in% likely_lymphoid
)

# --- Check plasma B cell marker expression with low expression ---
FeaturePlot(lymphoid_cells, features = c("IGHG4", "IGHG1", "IGHA1"), slot = "data", cols = reds_palette)

# ============================================================================
# Section 1.1: CONCLUSION — Marker-Based Interpretation of Lymphoid Clusters
# ----------------------------------------------------------------------------
#  All lymphoid populations are located in cluster 1.
#  B cells are correctly annotated by initial automatic labels.
#  No plasma B cell population detected based on expression of IGHG/IGHA.
#  T cell cluster includes a mix of subtypes (CD4+, CD8+, Tregs, NK),
#  and requires subdivision based on canonical markers.
# ============================================================================

# --- Apply refined annotations ---

# T CD8
annotated_seurat$Annotation <- ifelse(
    annotated_seurat$seurat_clusters %in% c(2,25) & 
    annotated_seurat$Annotation == "SnR.T_cells",
  "T CD8",
  annotated_seurat$Annotation
)

# T CD4
cd4_cells <- WhichCells(lymphoid_cells, expression = IL7R >= 1)
annotated_seurat$Annotation <- ifelse(
  colnames(annotated_seurat) %in% cd4_cells & 
    annotated_seurat$seurat_clusters %in% c(3,13) & 
    annotated_seurat$Annotation == "SnR.T_cells",
  "T CD4",
  annotated_seurat$Annotation
)

# Regulatory T cells
treg_cells <- WhichCells(lymphoid_cells, expression = CTLA4 >= 1|FOXP3 >= 1)
annotated_seurat$Annotation <- ifelse(
  colnames(annotated_seurat) %in% treg_cells&
  annotated_seurat$seurat_clusters == 13, 
  "T regs",
  annotated_seurat$Annotation
)

# NK cells
nk_cells <- WhichCells(lymphoid_cells, expression = NKG7 >= 1)
annotated_seurat$Annotation <- ifelse(
  colnames(annotated_seurat) %in% nk_cells &
    annotated_seurat$Annotation == "SnR.NK_cell" &
    annotated_seurat$seurat_clusters == 2,
  "NK cell",
  annotated_seurat$Annotation
)

# Reinforce correct B cell label for cluster 
annotated_seurat$Annotation <- ifelse(
  annotated_seurat$seurat_clusters == 12 & annotated_seurat$Annotation == "SnR.B_cell",
  "B cell",
  annotated_seurat$Annotation
)

# ============================================================================
# Section 2: Manual Refinement of Myeloid Cell Annotations
# ============================================================================

# --- Define canonical myeloid markers ---
myeloid_markers <- c(
  "LYZ",                         # General myeloid
  "CD68", "CD163", "CD14", "FCGR3A",  # Macrophages / Monocytes (TAMs)
  "CPVL", "CLEC9A", "LGALS2",         # cDC1
  "CD1C", "FCER1A", "CLEC10A",        # cDC2
  "CCR7", "LAMP3", "BIRC3",           # Mature DC
  "CXCR3", "LILRA4", "JCHAIN",        # pDCs
  "FCGR3B", "CXCR2",                  # Neutrophils
  "TPSAB1", "CPA3"                    # Mast cells
)

# --- Identify candidate myeloid clusters ---
likely_myeloid_clusters <- unique(as.character(top15_markers$cluster[top15_markers$gene %in% myeloid_markers]))
print(likely_myeloid_clusters)

# --- Violin plots for marker expression ---
present_myeloid_markers <- myeloid_markers[myeloid_markers %in% colnames(quantile_df)]

vln_plots_myeloid <- lapply(present_myeloid_markers, function(marker) {
  cluster_order <- quantile_df[, 1][order(quantile_df[[marker]], decreasing = TRUE)]
  VlnPlot(annotated_seurat, features = marker, pt.size = 0) +
    scale_x_discrete(limits = cluster_order) +
    labs(title = marker) +
    theme(
      axis.text.x = element_text(size = rel(0.5)),
      axis.text.y = element_text(size = rel(0.5)),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(face = "bold", size = rel(1))
    ) +
    NoLegend()
})

gridExtra::grid.arrange(grobs = vln_plots_myeloid, ncol = 5)

# --- Subset myeloid cells for refinement ---
myeloid_subset <- subset(
  annotated_seurat,
  subset = Annotation %in% c("SnR.Macrophage", "SnR.Monocyte", "SnR.DC", "SnR.Neutrophils")&
    seurat_clusters %in% c(likely_myeloid_clusters,23)
)

# --- Visual confirmation of key markers ---
FeaturePlot(myeloid_subset, features = c(
  "CD68", "CD163", "CD14", "FCGR3A",      # TAMs
  "CPVL", "CLEC9A", "LGALS2",             # cDC1
  "CD1C", "FCER1A", "CLEC10A",            # cDC2
  "CCR7", "LAMP3", "BIRC3"                # Mature DC
), slot = "data", cols = reds_palette, ncol = 5)

# ============================================================================
# Section 2.1: CONCLUSION — Marker-Based Interpretation of Myeloid Clusters
# ----------------------------------------------------------------------------
#  Tissue macrophages are concentrated in clusters 4, 12, 17, partially in 10
#  CD14+ monocytes appear in cluster 15, possibly mixed with epithelial cells (doublets)
#  Conventional dendritic cells (cDCs) localized in cluster 12 and partially in 4
#  No clear neutrophil or mast cell populations detected
# ============================================================================
# --- Assign "Tissue macrophages" label to TAM-like clusters ---
annotated_seurat$Annotation <- ifelse(
  annotated_seurat$Annotation %in% c("SnR.Macrophage", "SnR.Monocyte", "SnR.DC") &
    annotated_seurat$seurat_clusters %in% c(7,11,17,18,23),
  "TAMs",
  annotated_seurat$Annotation
)

# --- Refine DC subtypes ---
#Defining pDC
annotated_seurat$Annotation<-ifelse(annotated_seurat$seurat_clusters=="26",
                             "pDC",
                             annotated_seurat$Annotation)

# cDC1
cDC1_cells <- WhichCells(myeloid_subset, expression = CLEC9A >= 1)
annotated_seurat$Annotation <- ifelse(
  colnames(annotated_seurat) %in% cDC1_cells & 
    colnames(annotated_seurat)%in%colnames(myeloid_subset),
  "cDC1",
  annotated_seurat$Annotation
)

# cDC2
cDC2_cells <- WhichCells(myeloid_subset, expression = CD1C >= 2 | FCER1A >= 2 | CLEC10A >= 2.5)
annotated_seurat$Annotation <- ifelse(
  colnames(annotated_seurat) %in% cDC2_cells & 
    colnames(annotated_seurat)%in%colnames(myeloid_subset),
  "cDC2",
  annotated_seurat$Annotation
)

# Mature DC
mature_dc_cells <- WhichCells(myeloid_subset, expression = CCR7>= 2.8 | BIRC3>= 2.8)
annotated_seurat$Annotation <- ifelse(
  colnames(annotated_seurat) %in% mature_dc_cells & 
    colnames(annotated_seurat)%in%colnames(myeloid_subset),
  "Mature DC",
  annotated_seurat$Annotation
)

# ============================================================================
# Section 3: Refinement of Stromal Cell (Fibroblast, Epithelial, Endothelial)
# ============================================================================

# --- Define markers for stromal cell types ---
stromal_markers <- c(
  "COL1A1", "DCN",               # Fibroblasts
  "RGS5", "ACTA2",               # Pericytes
  "EPCAM", "KRT18",              # Epithelial
  "PECAM1", "VWF",               # Endothelial
  "TUBB1", "PPBP"                # Platelets (not included here, but tracked)
)

# --- Identify candidate stromal clusters ---
stromal_clusters <- unique(as.character(top15_markers$cluster[top15_markers$gene %in% stromal_markers]))
print(stromal_clusters)

# --- Violin plots for expression of each stromal marker ---
present_stromal_markers <- stromal_markers[stromal_markers %in% colnames(quantile_df)]

violins_stromal <- lapply(present_stromal_markers, function(marker) {
  cluster_order <- quantile_df[, 1][order(quantile_df[[marker]], decreasing = TRUE)]
  VlnPlot(annotated_seurat, features = marker, pt.size = 0) +
    scale_x_discrete(limits = cluster_order) +
    labs(title = marker) +
    theme(
      axis.text.x = element_text(size = rel(0.8)),
      axis.text.y = element_text(size = rel(0.5)),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(family = "serif")
    ) +
    NoLegend()
})
gridExtra::grid.arrange(grobs = violins_stromal, ncol = 4)

# --- Confirm expression in subset of stromal-like annotated cells ---
stromal_subset <- subset(
  annotated_seurat,
  subset = Annotation %in% c(
    "SnR.Endothelial_cells", "SnR.Epithelial_cells",
    "SnR.Fibroblasts", "SnR.Smooth_muscle_cells",
    "SnR.Tissue_stem_cells", "SnR.Chondrocytes"
  )
)

# ============================================================================
# Section 3.1: CONCLUSION — Interpretation of Stromal Markers
# ----------------------------------------------------------------------------
#  Endothelial and epithelial cells are correctly annotated
#  Cluster 9 is fully epithelial
#  Pericytes found primarily in cluster 7 (and partially 17)
#  Smooth muscle, tissue stem, and chondrocytes should be merged into fibroblasts
# ============================================================================

# --- Refine annotations based on marker expression and clustering ---

# Fibroblasts 
annotated_seurat$Annotation <- ifelse(
  annotated_seurat$seurat_clusters %in% c(5,10),
  "Fibroblasts",
  annotated_seurat$Annotation
)

# Pericytes
annotated_seurat$Annotation <- ifelse(
  annotated_seurat$seurat_clusters %in% c(14,20)|
    annotated_seurat$Auto_Annotation %in% "SnR.Tissue_stem_cells",
  "Pericytes",
  annotated_seurat$Annotation
)

# Epithelial
annotated_seurat$Annotation <- ifelse(
  annotated_seurat$seurat_clusters %in% c(0,1,4,8,9,16,21,22,27) &
  annotated_seurat$Annotation == "SnR.Epithelial_cells",
  "Epithelial cells",
  annotated_seurat$Annotation
)

# Endothelial
annotated_seurat$Annotation <- ifelse(
  annotated_seurat$seurat_clusters %in% c(6,15,19,24) &
  annotated_seurat$Annotation == "SnR.Endothelial_cells",
  "Endothelial cells",
  annotated_seurat$Annotation
)

# ============================================================================
# Section 4: Identification of Cycling Cells (if any)
# ============================================================================

p_cycle <- FeaturePlot(annotated_seurat, features = c("MKI67", "CDK1"), cols = reds_palette)
p_clusters <- DimPlot(annotated_seurat, group.by = "seurat_clusters", label = TRUE) + NoLegend()
p_clusters | p_cycle

# ============================================================================
# Section 4.1: CONCLUSION — Cell Cycle Status
# ============================================================================
annotated_seurat$Annotation<-ifelse(annotated_seurat$seurat_clusters%in%c(21,25,23)&
                                      annotated_seurat$Annotation%in%c("TAMs","T CD8"),
                             paste("Cycling",annotated_seurat$Annotation, sep = " "),
                             annotated_seurat$Annotation)

# ============================================================================
# Section 5: Reassignment of Unresolved or Rare Cell Annotations
# ============================================================================

# --- Identify rare or unresolved annotations ---
rare_annotations <- unique(annotated_seurat$Annotation[grepl("^SnR\\.", annotated_seurat$Annotation)])
rare_cell_ids <- colnames(annotated_seurat)[annotated_seurat$Annotation %in% rare_annotations]

# --- Visualize distribution of rare/unresolved cells ---
DimPlot(annotated_seurat, cells.highlight = rare_cell_ids) +
  labs(title = "Cells with Rare or Unresolved Annotations")

# --- Reassign those cells using majority rule by cluster ---

# Step 1: Replace unresolved annotations with their seurat cluster label (temporarily)
resolved_annotation <- ifelse(
  annotated_seurat$Annotation %in% rare_annotations | is.na(annotated_seurat$Annotation),
  as.character(annotated_seurat$seurat_clusters),
  as.character(annotated_seurat$Annotation)
)

# Step 2: Build majority label map per cluster
cluster_levels <- levels(annotated_seurat$seurat_clusters)
majority_label_map <- setNames(
  object = sapply(cluster_levels, function(cluster) {
    cluster_annots <- annotated_seurat$Annotation[annotated_seurat$seurat_clusters == cluster]
    names(sort(table(cluster_annots), decreasing = TRUE))[1]
  }),
  nm = cluster_levels
)

# Step 3: Reassign cluster-based labels to unresolved cells
annotated_seurat$Annotation <- dplyr::recode(resolved_annotation, !!!majority_label_map)

# Step 4: Clean up metadata
Idents(annotated_seurat) <- annotated_seurat$Annotation
annotated_seurat$Auto_Annotation <- NULL  # Remove redundant column

# ============================================================================
# Section 6: Final Save
# ============================================================================
saveRDS(annotated_seurat, file = "results/seurat_objects/TwiggerDataset/Final_Annotated.rds")

# ============================================================================
# Section 7: Final Visualization
# ============================================================================
library(ggplot2)
library(cowplot)
# Create output directory for figures
fig_dir <- "results/figures/annotation_plots"

# --- UMAP of clustering ---
umap_cluster <- DimPlot(
  annotated_seurat, reduction = "umap", label = TRUE, group.by = "seurat_clusters"
) +
  labs(title = "Clustering") +
  theme(axis.title = element_text(size = rel(0.6)),
        axis.text = element_text(size = rel(0.6))) +
  NoLegend()

# --- UMAP of final annotation ---
Disc.Pal <- DiscretePalette(length(unique(annotated_seurat$Annotation)), shuffle = TRUE)

umap_annotation <- DimPlot(
  annotated_seurat, reduction = "umap", label = FALSE,
  group.by = "Annotation", cols = Disc.Pal
) +
  labs(title = "Annotation") +
  theme(axis.title = element_text(size = rel(0.6)),
        axis.text = element_text(size = rel(0.6)))

# --- Save side-by-side UMAPs ---
combined_umap <- umap_cluster + umap_annotation
combined_umap <- combined_umap+labs(title = "Wu dataset: Breast cancer")

ggsave(filename = file.path(fig_dir, "Wu_UMAP_Annotation.pdf"),
       plot = combined_umap, width = 10, height = 5)

# --- DotPlot of marker expression across final annotations ---

markers_to_display <- c(
  "CD3D", "CD3E", "CD8A", "CD4", "IL7R", "FOXP3", "CTLA4",   # T cell subtypes
  "NKG7", "KLRD1", "MS4A1", "CD19",                          # NK / B cells
  "CD68", "CD14", "FCGR3A",                                  # TAMs
  "CLEC9A", "CD1C", "CLEC10A", "LAMP3", "BIRC3",             # DCs
  "COL1A1", "DCN", "RGS5", "ACTA2",                          # Fibroblast / Pericytes
  "VWF", "PECAM1", "EPCAM", "KRT19"                          # Endo / Epithelial
)

annotation_order <- c(
  "T CD8", "T CD4", "T regs", "NK cell", "B cell",
  "Tissue macrophages", "cDC1", "cDC2", "Mature DC",
  "Fibroblasts", "Pericytes", "Endothelial cells", "Epithelial cells"
)

final_dotplot <- DotPlot(
  annotated_seurat, features = markers_to_display,
  assay = "RNA", cols = "Spectral"
) +
  scale_y_discrete(limits = rev(annotation_order)) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = rel(0.8)),
    axis.text.y = element_text(size = rel(0.8))
  )+labs(title = "Wu dataset: Breast cancer")

ggsave(filename = file.path(fig_dir, "Wu_DotPlot_Annotation.pdf"),
       plot = final_dotplot, width = 10, height = 6)
