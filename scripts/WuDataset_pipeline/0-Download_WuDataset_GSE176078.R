# ============================================================================
# Title: Download Wu dataset GSE176078 BRCA scRNA-seq.
# Author: elopezhuerta
# Description: Downloads the compressed file GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz
#              from NCBI GEO FTP and stores it in data/GSE176078/.
# ============================================================================

# --- Load required library ---
# For download.file, base R is sufficient

# NCBI GEO always stores supplementary files in this predictable FTP structure.
# --- Set parameters ---
geo_accession <- "GSE176078"
ftp_base_url <- paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/", geo_accession, "/suppl/")
file_name <- "GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"

# --- Set output directory ---
output_dir <- "data/Wu_downloads_GSE176078"

# --- Build full URL ---
full_url <- paste0(ftp_base_url, file_name)
dest_file <- file.path(output_dir, file_name)

# --- Download the file ---
message("Downloading: ", file_name)
options(timeout = 1000)  # Increase timeout for large files

tryCatch({
  download.file(full_url, destfile = dest_file, mode = "wb")
  message("Download complete: ", dest_file)
}, error = function(e) {
  message("Failed to download: ", conditionMessage(e))
})

# --- Optional: print path ---
cat("File saved to: ", dest_file, "\n")

# Extract files from downloaded file from GEO
untar(dest_file,exdir = output_dir)
