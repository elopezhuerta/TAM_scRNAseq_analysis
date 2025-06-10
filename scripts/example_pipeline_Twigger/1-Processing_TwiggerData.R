# ============================================================================
# Title: Preprocessing and Integration of Twigger scRNA-seq Datasets
# Author: elopezhuerta
# Description: This script loads Matrix Market files from selected Twigger samples,
#              creates Seurat objects, performs quality control (QC), batch
#              correction, and data integration.
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

# --- Define Input Directory and Sample IDs ---
data_path <- "data/Twigger_Downloads"
sample_ids <- c("NMC1", "NMC2", "NMC3", "AA4", "AB4", "AC4")

# --- List and Build Full File Paths ---
file_list <- list.files(data_path)
full_paths <- file.path(data_path, file_list)

# Initialize list to store Seurat objects
seurat_list <- list()

# Loop through each sample
for (sample in sample_ids) {
  
  # Get matching files for the current sample
  sample_files <- full_paths[grep(sample, full_paths)]
  
  # Sort to ensure order is always barcodes, features, matrix
  sample_files <- sample_files[order(sample_files)]
  
  barcode_path <- sample_files[1]
  feature_path <- sample_files[2]
  matrix_path  <- sample_files[3]
  
  # Read gene (feature) names — column 2 is the gene symbol
  gene_info <- readr::read_tsv(feature_path, col_names = FALSE)
  gene_names <- gene_info[[2]]
  
  # Read barcode (cell ID) names
  barcode_info <- readr::read_tsv(barcode_path, col_names = FALSE)
  barcode_names <- barcode_info[[1]]
  
  # Read the expression matrix (Matrix Market format)
  expr_matrix <- Matrix::readMM(matrix_path) %>%
    magrittr::set_rownames(gene_names) %>%
    magrittr::set_colnames(barcode_names)
  
  # Create Seurat object (filtering low-quality cells)
  seurat_list[[sample]] <- CreateSeuratObject(counts = expr_matrix,
                                              project = sample,
                                              min.features = 200, min.cells = 3)
}


# --- Merge All Seurat Objects ---
merged_seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  # Add unique sample identifiers to each cell using add.cell.ids
  add.cell.ids = sample_ids
)
merged_seurat$lote <- ifelse(grepl("NMC", merged_seurat$orig.ident), "Lote1", "Lote2")

# --- Save Raw Merged Seurat Object ---
saveRDS(merged_seurat, file = "results/seurat_objects/TwiggerDataset/Raw_Twigger.rds")

# ============================================================================
# Quality Control (QC)
# ============================================================================

# Calculate mitochondrial gene percentage using gene name pattern
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")

# Identify damaged cells with high mitochondrial content (outliers above threshold)
outlier_mt <- isOutlier(
  merged_seurat$percent.mt, 
  type = "higher",
  batch = merged_seurat$orig.ident)

# Identify cells with unusually high feature counts (genes detected)
outlier_features <- isOutlier(
  merged_seurat$nFeature_RNA, 
  type = "higher", 
  batch = merged_seurat$orig.ident,
  nmads = 6
)

# Identify cells with unusually high library sizes
outlier_counts <- isOutlier(
  merged_seurat$nCount_RNA, 
  type = "higher", 
  batch = merged_seurat$orig.ident,
  nmads = 6
)

# Filter out the identified outliers
qc_passed_seurat <- merged_seurat[, !(outlier_mt | outlier_features | outlier_counts)]

# ============================================================================
# Visualize QC metrics: pre- and post-QC
# ============================================================================

# Violin plots before QC
VlnPlot(merged_seurat, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0)

# Scatter plot before QC
plot_pre_qc <- FeatureScatter(
  merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"
) + labs(title = "Before QC")

# Violin plots after QC
VlnPlot(qc_passed_seurat, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0)

# Scatter plot after QC
plot_post_qc <- FeatureScatter(
  qc_passed_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"
) + labs(title = "After QC")

# Display pre- and post-QC scatter plots side-by-side
plot_pre_qc + plot_post_qc

# --------------------------------------------
# Batch correction and data integration
# --------------------------------------------
# Enrich myeloid cells before performing clustering to save computation resources
# Get cells with high myeloid markers 
Myeloid_cells <- colnames(subset(x = qc_passed_seurat,
                                 subset = PTPRC > 2|FCER1G>2|CD14>2|CD68>2|C1QB>2|`HLA-DRB1`>2|`HLA-DRA`>2|`HLA-DQB1`>2))

# Subset a myeloid-enriched seurat object
Mielod_enriched_seurat <- qc_passed_seurat[,Myeloid_cells]

# Split object by batch (based on 'lote' metadata)
batch_list <- SplitObject(Mielod_enriched_seurat, split.by = "lote")

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
integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:10)
integrated_seurat <- FindClusters(integrated_seurat, resolution = 0.5)

saveRDS(integrated_seurat, file = "results/seurat_objects/TwiggerDataset/integrated_Twigger.rds")