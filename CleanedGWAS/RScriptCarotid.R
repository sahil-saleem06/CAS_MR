# Load required libraries
library(vroom)
library(data.table)
library(dplyr)
library(glue)

# Set base path
base_path <- "/project/voltron/Resources/GWAS-Summary-Statistics/MVP/Carotid/JULY2025"

# List subdirectories (e.g., Definition1/, Definition2/, etc.)
subdirs <- list.dirs(base_path, recursive = FALSE, full.names = TRUE)

# Outer loop: iterate over each subdirectory
for (subdir in subdirs) {
  message(glue("Entering directory: {subdir}"))
  
  # List GWAS files within this subdir
  gwas_files <- list.files(subdir, full.names = TRUE)
  
  # Inner loop: process each file
  for (gwas_file in gwas_files) {
    message(glue("Processing: {gwas_file}"))
    
    # Step 1: Load raw GWAS summary statistics
    file_raw <- tryCatch(
      fread(gwas_file),
      error = function(e) {
        message(glue("Failed to read: {gwas_file}"))
        return(NULL)
      }
    )
    if (is.null(file_raw)) next
    
    # Step 2: Clean and sort (robust version)
    if (!all(c("CHR", "POS") %in% colnames(file_raw))) {
      message(glue("Missing CHR or POS in: {gwas_file}, skipping."))
      next
    }
    
    # Ensure CHR is numeric
    file_raw[, CHR := gsub("^chr", "", CHR)]
    file_raw[, CHR := as.numeric(CHR)]
    
    # Drop columns if they exist
    cols_to_remove <- intersect(c("SNP", "BUILD"), colnames(file_raw))
    
    file_tidy <- file_raw %>%
      arrange(CHR, POS) %>%
      dplyr::select(-all_of(cols_to_remove)) %>%
      na.omit()
    
    # Step 2.5: Save cleaned pre-liftover data
    out_path <- glue(
      "/project/damrauer_scratch/Users/saleemsa/Cleaned GWAS/cleaned_{basename(gwas_file)}"
    )
    vroom::vroom_write(file_tidy, out_path, delim = "\t")
    message(glue("Saved cleaned file: {out_path}"))
    
    # Liftover steps are skipped as per your setup
  }
}
