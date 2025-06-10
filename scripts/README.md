# Scripts Overview

This folder contains all scripts used in the `TAM_scRNAseq_analysis` project.  
Each script is organized by **dataset pipeline** or **analysis module**.

---

## Contents

### WuDataset_pipeline (Main dataset of the study — GSE176078)

- `0-Download_WuDataset_GSE176078.R`:  
  Downloads scRNA-seq dataset GSE176078 from GEO. The dataset is saved in `data/Wu_downloads_GSE176078/`.

- `1-Processing_WuDataset.R`:  
  Loads and processes the downloaded matrix files into Seurat objects. Performs batch correction and integration. Saves processed object in `results/seurat_objects/WuDataset/`.

- `2-AutomaticAnnotation_WuData.R`:  
  Performs automatic cell type annotation using **SingleR** and a reference dataset.

- `3-Manual Annotation Refinement_WuData.R`:  
  Refines annotations based on marker expression and cluster identity.  
  Saves the final annotated object and generates UMAP and DotPlot visualizations.

---

### example_pipeline_Twigger (Example additional dataset — Twigger, E-MTAB studies)

- `0-Download_Public_TwiggerData.R`:  
  Downloads scRNA-seq datasets from publicly available ArrayExpress studies (E-MTAB-9841, -10855, -10885).  
  The Twigger dataset (non-tumoral mammary tissue) is used as an example.  
  Data is saved in `data/Twigger_Downloads/`.

- `1-Processing_TwiggerData.R`:  
  Loads and processes the downloaded matrix files into Seurat objects.  
  Performs basic filtering, normalization, and batch correction.  
  Saves processed object in `results/seurat_objects/TwiggerDataset/`.

- `2-AutomaticAnnotation_TwiggerData.R`:  
  Performs automatic cell type annotation using **SingleR**.

- `3-Manual Annotation Refinement_Twigger.R`:  
  Refines annotations based on marker expression and cluster identity.  
  Saves the final annotated object and generates UMAP and DotPlot visualizations.

---

### in_process/ (Modular scripts for additional analyses)

- `1-QC Clustering Annotation Cleaning`: Scripts for clustering and annotation cleaning across various datasets.
- `2-TAM Modules`: Scripts to compute TAM-related gene modules and pathway enrichment.
- `3-Infiltration scores`: Script for computing immune infiltration scores.
- `4-Across tissues`: Scripts for cross-tissue integration and trajectory analysis.
- `5-Tissue vs Mono`: Scripts for comparative analysis of TAMs vs monocytes.
- `6-TCGA analysis`: Scripts for integration with TCGA datasets (BayesPrism draft, BP-TCGA).

---

## Notes

- All scripts are designed to work with **relative paths** to maintain reproducibility.
- You can run scripts independently if you already have the necessary input data in place.
- The `WuDataset_pipeline/` corresponds to the **main dataset analyzed in the study**.  
  The `example_pipeline_Twigger/` demonstrates the same pipeline applied to a complementary dataset.

---
