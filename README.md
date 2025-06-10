# TAM_scRNAseq_analysis

This repository contains code and data analysis pipelines for processing and analyzing tumor-associated macrophage (TAM) populations and the tumor microenvironment using single-cell RNA-seq datasets.

## Project Structure

- `scripts/`: Analysis scripts organized by theme and pipeline component.
- `data/`: Input data and downloaded datasets.
- `results/`: Processed results, figures, and tables.
- `example_pipeline_Twigger/`: Example full pipeline applied to the Twigger dataset (E-MTAB studies).

## Example Pipelines

Several datasets were processed during this work, therefore, we chose two examples to fully document the pipeline for: tumor-associated macrophages (`example_pipeline_Wu`) and normal mammary tissue (`example_pipeline_Twigger/`) datasets.
1. Downloading datasets.
2. Pre-processing and integration.
3. Automatic annotation using SingleR.
4. Manual annotation refinement.
5. Final visualization and saving of results.
