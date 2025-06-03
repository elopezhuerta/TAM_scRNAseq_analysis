# Twigger_Downloads

This folder stores scRNA-seq datasets downloaded from the following public ArrayExpress studies:

- [E-MTAB-9841](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-9841)
- [E-MTAB-10855](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10855)
- [E-MTAB-10885](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10885)

## Instructions

To populate this folder, run the script:

```r
scripts/0-Download_PublicData.R
```

The script will scrape the study directories and download only the selected sample files (e.g., NMC1–3, AA4–AC4).

## Notes

- This folder is excluded from version control via `.gitignore`.
- Files in this folder may be large and are not tracked in the repository.
