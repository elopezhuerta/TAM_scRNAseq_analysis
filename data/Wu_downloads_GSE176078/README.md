# 0-Download_WuDataset_GSE176078.R — README

## Purpose

This script downloads the supplementary file:

```
GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz
```

from the NCBI GEO repository and extracts its contents for further processing.

The **Wu BRCA scRNA-seq dataset (GSE176078)** is the **main dataset analyzed in this study**.  
The file contains preprocessed single-cell expression matrices and metadata, provided by the authors.

## Location of the downloaded data

The downloaded `.tar.gz` file will be stored in:

```
data/Wu_downloads_GSE176078/
```

The script also automatically extracts the `.tar.gz` archive into the same folder.  
After extraction, the expression matrices and metadata will be ready for loading in downstream scripts.

## How to run

Simply run:

```r
source("scripts/WuDataset_pipeline/0-Download_WuDataset_GSE176078.R")
```

## Notes

- The script uses `download.file()` with an increased timeout to handle large files.
- It leverages the predictable FTP structure of NCBI GEO:

```
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/
```

- After download, `untar()` is used to extract the contents.
- The next step in the pipeline is to run:

```
1-Processing_WuDataset.R
```

---

**Author:** elopezhuerta
