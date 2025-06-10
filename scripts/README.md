# Scripts Overview

This folder contains all scripts used in the TAM_scRNAseq_analysis project. Each script is numbered in the order of execution and modularized by task.

## Contents

- `0-Download_Public_TwiggerData.R`:  
  Downloads scRNA-seq datasets from publicly available ArrayExpress studies (e.g., E-MTAB-9841, -10855, -10885). The datasets are from non tumoral mammary tissue and the Twigger dataset was used as example. Dataset will be saved in `data/Twigger_Downloads/`.

- `1-Process_TwiggerData.Rmd`:  
  Loads and processes the downloaded matrix files into Seurat objects. Performs basic filtering, normalization, and annotation, and stores results in `results/seurat_objects/`.

## Notes

- All scripts are designed to work with relative paths to maintain reproducibility.
- You can run scripts independently if you already have the necessary input data in place.
