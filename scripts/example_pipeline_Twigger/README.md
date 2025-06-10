# Example Pipeline - Twigger Dataset

This folder contains a reproducible pipeline to process and annotate the Twigger dataset.

## Steps

1. `0-Download_Public_TwiggerData.R`: Downloads raw data (.tsv and .mtx) from ArrayExpress database.
2. `1-Processing_TwiggerData.R`: Creates Seurat objects, performs QC, batch correction, and integration.
3. `2-AutomaticAnnotation_TwiggerData.R`: Performs automatic annotation using SingleR and HPCA reference.
4. `3-Manual_Annotation_Refinement_Twigger.R`: Refines annotations for lymphoid, myeloid, stromal, and other populations, saving the final Seurat object and generating visualization plots.

## Outputs

- `results/seurat_objects/TwiggerDataset/Final_Annotated.rds`: Final annotated Seurat object.
- `results/figures/final_annotation/`: UMAPs and DotPlots of final annotations.
