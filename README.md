# TAM_scRNAseq_analysis

This repository contains code and data analysis pipelines for processing and analyzing tumor-associated macrophage (TAM) populations and the tumor microenvironment using single-cell RNA-seq (scRNA-seq) datasets.

The main dataset analyzed in this study is the **Wu BRCA scRNA-seq dataset** (GSE176078).  
Additionally, an example pipeline is provided for a normal mammary tissue dataset (Twigger dataset, E-MTAB studies) used in complementary analyses.

---

## Project Structure

- `scripts/`: Analysis scripts organized by theme and processing component.
  - `in_process/`: Modular scripts used for cross-dataset analysis and development.
  - `WuDataset_pipeline/`: Full reproducible pipeline for the **Wu BRCA scRNA-seq dataset** (main dataset of the study).
  - `example_pipeline_Twigger/`: Example pipeline for an additional dataset used in the study (**Twigger normal mammary tissue**, E-MTAB series).
- `data/`: Input data and downloaded datasets.
  - `Twigger_Downloads/`: Downloads from ArrayExpress (E-MTAB studies).
  - `Wu_downloads_GSE176078/`: Downloads from GEO (GSE176078).
- `results/`: Processed results, figures, and tables.
  - `seurat_objects/`: Processed Seurat objects, organized by dataset (`WuDataset/`, `TwiggerDataset/`).
  - `figures/`: Figures generated during manual annotation and refinement.
  - `tables/`: Supplementary tables.

---

## Main Pipeline — Wu BRCA scRNA-seq Dataset (GSE176078)

The pipeline implemented in `scripts/WuDataset_pipeline/` was used to process the main dataset of this study:

1. Downloading dataset from GEO.
2. Pre-processing and data integration.
3. Automatic annotation using **SingleR**.
4. Manual annotation refinement based on marker expression and cluster identity.
5. Final visualization and saving of results.

---

## Example Pipeline — Twigger Normal Mammary Tissue Dataset (E-MTAB Series)

An example pipeline is provided in `scripts/example_pipeline_Twigger/` for the processing of a normal mammary tissue dataset, used as complementary analysis in the study.

This pipeline follows the same overall structure:

1. Downloading datasets.
2. Pre-processing and data integration.
3. Automatic annotation using **SingleR**.
4. Manual annotation refinement.
5. Final visualization and saving of results.

---

## Notes

- The `scripts/in_process/` folder contains additional modular scripts used in advanced analyses (e.g. TAM modules, infiltration scores, cross-tissue comparisons, etc.).
- The `figures/` and `tables/` folders contain figures and tables used in the final manuscript.
- All pipelines are organized for reproducibility and were used to generate key results reported in the study.

---

**Maintainer:** [elopezhuerta]

