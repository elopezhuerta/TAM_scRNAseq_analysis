library(rvest) # version 1.0.4
# Create local directory to store datasets
download_dir <- "data/Twigger_Downloads"

# URL of the index pages to download
index_url <- list()
# E-MTAB-9841
index_url[[1]] <- "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/841/E-MTAB-9841/Files/"
# E-MTAB-10855
index_url[[2]] <- "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/855/E-MTAB-10855/Files/"
# E-MTAB-10885
index_url[[3]] <- "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/885/E-MTAB-10885/Files/"

# Select only files of interest
sample_ids <- c("NMC1", "NMC2", "NMC3", "AA4", "AB4", "AC4")

# Create list to store files
for (ii in seq_along(index_url)) {
  # Read the page and extract all file links
  page <- read_html(index_url[[ii]])
  files <- page %>%
    html_elements("a") %>%
    html_text() %>%
    # Only download data files
    grep(paste(sample_ids, collapse = "|"), ., value = TRUE)
  # Download all files within same url
  for (file_name in files) {
    full_url <- paste0(index_url[[ii]], file_name)
    dest <- file.path(download_dir, file_name)
    message("Downloading: ", file_name)
    download.file(full_url, destfile = dest, mode = "wb")
  }
}
