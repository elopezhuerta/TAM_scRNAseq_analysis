# Results

This folder contains all processed outputs and visualizations generated during the `TAM_scRNAseq_analysis` project.

## Structure

- `seurat_objects/`: Processed Seurat objects (`.rds`), organized by dataset:
  - `WuDataset/`: Refined annotated Seurat object for the Wu dataset (main dataset of the study).
  - `TwiggerDataset/`: Refined annotated Seurat object for the Twigger example dataset.

- `figures/`: Publication-quality figures.
  - `annotation_plots/`: UMAPs and DotPlots of final refined annotations for each dataset (Wu and Twigger).
  - `manuscript/`: Scripts used to generate main and supplementary figures for the manuscript.

- `tables/`: Supplementary tables.

## Notes

- These files are **not tracked** in the repository via `.gitignore` to reduce repo size.
- You can regenerate all contents using the pipelines provided in `scripts/WuDataset_pipeline/` and `scripts/example_pipeline_Twigger/`, and the figure scripts provided in `figures/manuscript/`.
