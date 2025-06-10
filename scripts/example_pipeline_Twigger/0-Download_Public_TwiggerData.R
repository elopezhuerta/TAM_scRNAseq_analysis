# ============================================================================
# Title: Download Twigger Datasets from ArrayExpress
# Author: elopezhuerta
# Description: This script downloads selected samples (e.g., NMC1–3, AA4–AC4)
#              from three publicly available scRNA-seq datasets hosted on
#              ArrayExpress (E-MTAB-9841, -10855, -10885). The script uses
#              rvest to parse file listings and download selected files only.
# ============================================================================

# --- Load Required Libraries ---
library(rvest)   # version 1.0.4

# --- Set Output Directory ---
download_dir <- "data/Twigger_Downloads"

# --- Define Index URLs for ArrayExpress Studies ---
index_url <- list(
  "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/841/E-MTAB-9841/Files/",
  "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/855/E-MTAB-10855/Files/",
  "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/885/E-MTAB-10885/Files/"
)

# --- Define Sample IDs of Interest ---
sample_ids <- c("NMC1", "NMC2", "NMC3", "AA4", "AB4", "AC4")

# --- Loop Through Each Study and Download Matched Files ---
for (ii in seq_along(index_url)) {
  message("Checking: ", index_url[[ii]])
  
  # Parse HTML and extract matching file names
  page <- read_html(index_url[[ii]])
  files <- page %>%
    html_elements("a") %>%
    html_text() %>%
    grep(paste(sample_ids, collapse = "|"), ., value = TRUE)
  
  # Download each matched file
  for (file_name in files) {
    full_url <- paste0(index_url[[ii]], file_name)
    dest <- file.path(download_dir, file_name)
    message("Downloading: ", file_name)
    download.file(full_url, destfile = dest, mode = "wb")
  }
}
