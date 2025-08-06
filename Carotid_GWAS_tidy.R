library(dplyr)
library(vroom)

# Set the folder path
folder_path <- "/project/voltron/Resources/GWAS-Summary-Statistics/MVP/Carotid/JULY2025"

# List all files (add full.names = TRUE if you want full paths)
file_list <- list.files(folder_path, full.names = TRUE)

# Loop through each file
for (file in file_list) {
  definitions <- list.files(file, full.names = TRUE)
  
  for(gwas in definitions) {
    
    file_raw <- fread(gwas)
    
    # Fix CHR column
    file_raw[, CHR := gsub("^chr", "", CHR)]
    file_raw[, CHR := as.numeric(CHR)]
    
    # Sort efficiently
    file_sorted <- setorder(pad_raw, CHR, POS)
    
    
    # Drop any unnecessary columns (example: if 'BUILD' or 'SNP' exists, remove it)
    file_tidy <- file_sorted %>%
      select(-any_of(c("BUILD", "SNP"))) %>%  # safe way to remove if they exist
      na.omit()                               # drop any rows with NAs
    
    file_ready <- file_tidy %>%
      dplyr::mutate(N = (N_case+N_ctrl)) %>%
      select(CHROM = CHR, POS, MARKER = MarkerID, EFFECT_ALLELE = Allele1, OTHER_ALLELE = Allele2, EFF_ALL_FREQ = AF_Allele2, BETA, SE, PVAL = p.value, INFO = imputationInfo, N, N_CASE = N_case, N_CONTROL = N_ctrl)
    
    out_path <- paste0(
      "/project/damrauer_scratch/Users/saleemsa/CarotidStenosis/GWASInspectorInput/cleaned_",
      basename(gwas)
    )
    
    vroom_write(file_ready, out_path, delim = "\t")
    
  }
  
}
