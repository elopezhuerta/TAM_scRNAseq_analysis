# ============================================================================
# Title: Preprocessing and Integration of Wu scRNA-seq Dataset (GSE176078)
# Author: elopezhuerta
# Description: This script loads Matrix Market files from GSE176078,
#              creates Seurat object, performs batch correction and
#              integration. No QC was needed as the authors already provided
#              pre-filtered data. Subsampling was performed to reduce RAM usage.
# ============================================================================

# --- Load Required Libraries ---
library(dplyr)
library(Seurat)     # version 4.3.0
library(readr)
library(Matrix)
library(magrittr)
library(scuttle)
library(ggplot2)
library(patchwork)

# ============================================================================
# Section 1: Load Wu Dataset
# ============================================================================

# --- Define paths ---
output_dir <- "data/Wu_downloads_GSE176078"
data_path <- file.path(output_dir, "Wu_etal_2021_BRCA_scRNASeq")

# --- List files and match types ---
file_list <- list.files(data_path)
full_paths <- file.path(data_path, file_list)

# Specify path for each type of file
barcode_path <- full_paths[grep("barcodes",full_paths)]
feature_path <- full_paths[grep("genes",full_paths)]
matrix_path <- full_paths[grep("sparse",full_paths)]
metadata_path <- full_paths[grep("metadata",full_paths)]

# Read gene (feature) names — column 2 is the gene symbol
gene_info <- readr::read_tsv(feature_path, col_names = FALSE)
gene_names <- gene_info[[1]]

# Read barcode (cell ID) names
barcode_info <- readr::read_tsv(barcode_path, col_names = FALSE)
barcode_names <- barcode_info[[1]]

# Read the expression matrix (Matrix Market format)
expr_matrix <- Matrix::readMM(matrix_path) %>%
  magrittr::set_rownames(gene_names) %>%
  magrittr::set_colnames(barcode_names)

# Importing metadata provided by authors
Wu.metadata <- read.csv(metadata_path, header = T)
Wu.metadata <- tibble::column_to_rownames(Wu.metadata, "X")

# Create Seurat object (filtering low-quality cells)
SeuratObject <- CreateSeuratObject(counts = expr_matrix,
                                   meta.data = Wu.metadata,
                                   min.features = 200, min.cells = 3)

# --- Save Raw Merged Seurat Object ---
saveRDS(SeuratObject, file = "results/seurat_objects/WuDataset/Raw_WuDataset.rds")

# ============================================================================
# Section 2: Quality Control (QC)
# ============================================================================

# NOTE: Mitochondrial % and other QC filters were applied by the authors.
# We visualize QC metrics for verification only.

# --- Violin and scatter plots ---
VlnPlot(SeuratObject, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0)

plot_pre_qc <- FeatureScatter(
  SeuratObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# ============================================================================
# Section 3: Subsampling to Reduce RAM Usage
# ============================================================================

# LEVERAGING REPORTED ANNOTATION BY AUTHORS:
# The authors provided the celltype_major annotation.
# To reduce memory usage during integration, we subsample abundant cell types:
# - Targeted subsampling: Keep 16,000 cells total across T-cells, B-cells, Plasmablasts, PVL.
# - Remaining cell types (rare populations) are kept at full size.

# --- Subset cells of abundant lineages ---
cells_to_downsample <- subset(
  SeuratObject,
  subset = celltype_major%in%c("B-cells", "T-cells","Plasmablasts","PVL"))

# For reproducibility
set.seed(123)
# Subtypes have a total of 4700 cells. Eliminate 1600 of those cells
cells_to_remove <- sample(colnames(cells_to_downsample),ncol(cells_to_downsample)-16000)

# Subset final object
wu_subset <- SeuratObject[,!colnames(SeuratObject)%in%cells_to_remove]

#Free space
rm(cells_to_downsample,SeuratObject)

# ============================================================================
# Section 4: Batch Correction and Integration
# ============================================================================

# Split object by batch (based on 'lote' metadata)
batch_list <- SplitObject(wu_subset, split.by = "orig.ident")

# Normalize and find variable features for each batch
for (ii in seq_along(batch_list)) {
  batch_list[[ii]] <- NormalizeData(batch_list[[ii]], verbose = FALSE)
  batch_list[[ii]] <- FindVariableFeatures(
    batch_list[[ii]], selection.method = "vst",
    nfeatures = 2000, verbose = FALSE
  )
}

# Identify anchors for integration
anchors <- FindIntegrationAnchors(object.list = batch_list, dims = 1:30)

# Integrate data across batches (returns new "integrated" assay)
integrated_seurat <- IntegrateData(anchorset = anchors, dims = 1:30)

# Switch to integrated assay for downstream analysis
DefaultAssay(integrated_seurat) <- "integrated"

# Scale, PCA and UMAP (do not run tSNE yet)
integrated_seurat <- ScaleData(integrated_seurat, verbose = FALSE)
integrated_seurat <- RunPCA(integrated_seurat, npcs = 30, verbose = FALSE)
integrated_seurat <- RunUMAP(integrated_seurat, reduction = "pca", dims = 1:30)

# Clustering
integrated_seurat <- FindNeighbors(integrated_seurat)
integrated_seurat <- FindClusters(integrated_seurat, resolution = 1)
# Save
saveRDS(integrated_seurat, file = "results/seurat_objects/WuDataset/integrated_WuDataset.rds")