library(R.utils)
library(ieugwasr)
library(data.table)
library(janitor)
library(tidyr)
library(dplyr)
library(tidyverse)

exposure_files <- list(
  #doneBMI = "/project/damrauer_scratch/Users/saleemsa/ExposureSumStats/BMI/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz",
  #CKD = "/project/damrauer_scratch/Users/saleemsa/ExposureSumStats/CKD/CKD_overall_ALL_JW_20180223_nstud30.dbgap.txt.gz"
  #done DBP = "/project/damrauer_scratch/Users/saleemsa/ExposureSumStats/DBP/GCST90310295.tsv.gz",
  #doneDB2 =  "/project/damrauer_scratch/Users/saleemsa/ExposureSumStats/Diabetes(Type2)/Suzuki.Nature2024.T2DGGI.EUR.sumstats.zip"
  #HDL = "/project/damrauer_scratch/Users/saleemsa/ExposureSumStats/HDL/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_HDL_INV_ALL_with_N_1.gz"
  LDL = "/project/damrauer_scratch/Users/saleemsa/ExposureSumStats/LDL/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_LDL_INV_ALL_with_N_1.gz"
  #doneLPA = "/project/damrauer_scratch/summary_stats/LPA/30790_raw.gwas.imputed_v3.both_sexes.varorder.tsv.bgz"
  #NONHDL = "/project/damrauer_scratch/Users/saleemsa/ExposureSumStats/Non-HDL/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_nonHDL_INV_ALL_with_N_1.gz"
  #doneSBP = "/project/damrauer_scratch/Users/saleemsa/ExposureSumStats/SBP/QC_SBP-meta-analysis_ICBP2024.tsv",
  #TC = "/project/damrauer_scratch/Users/saleemsa/ExposureSumStats/TC/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_TC_INV_ALL_with_N_1.gz"
  )

#Function to change to standard column names
standardize_columns <- function(df, exposure_name) {
  original_names <- names(df)
  cat("Original columns for", exposure_name, ":", paste(original_names, collapse = ", "), "\n")
  
  # Create a mapping based on common column name patterns
  standardized_df <- df
  
  # Common patterns for SNP identifier
  snp_patterns <- c("SNP", "rsid", "RSID", "rs", "MarkerName", "ID", "rs_id", "SNP_ID")
  snp_col <- original_names[original_names %in% snp_patterns][1]
  if (is.na(snp_col)) {
    # Try partial matching
    snp_col <- original_names[grepl("^(SNP|rsid|rs|ID)", original_names, ignore.case = TRUE)][1]
  }
  
  
  
  # Common patterns for p-value
  pval_patterns <- c("P", "p", "pval", "p_value", "P_VALUE", "PVAL", "p.value", "P.VALUE", "P-value", "Pval", "METAL_Pvalue")
  pval_col <- original_names[original_names %in% pval_patterns][1]
  if (is.na(pval_col)) {
    pval_col <- original_names[grepl("^(P|p)(_|\\.|$)", original_names, ignore.case = TRUE)][1]
  }
  
  # Common patterns for beta/effect size
  beta_patterns <- c("BETA", "beta", "Beta", "effect", "EFFECT", "b", "B", "logOR", "log_or", "odds_ratio", "METAL_Effect")
  beta_col <- original_names[original_names %in% beta_patterns][1]
  if (is.na(beta_col)) {
    beta_col <- original_names[grepl("^(BETA|beta|effect|b)", original_names, ignore.case = TRUE)][1]
  }
  
  if (!is.na(beta_col)) {
    if (beta_col == "odds_ratio") {
      cat("Converting odds_ratio to log scale for beta.\n")
      standardized_df$BETA <- log(as.numeric(standardized_df[[beta_col]]))
    } else {
      names(standardized_df)[names(standardized_df) == beta_col] <- "BETA"
    }
  }
  
  # Common patterns for standard error
  se_patterns <- c("SE", "se", "Se", "standard_error", "STANDARD_ERROR", "stderr", "SE_BETA", "StdErr", "METAL_StdErr")
  se_col <- original_names[original_names %in% se_patterns][1]
  if (is.na(se_col)) {
    se_col <- original_names[grepl("^(SE|se|standard)", original_names, ignore.case = TRUE)][1]
  }
  
  # Common patterns for effect allele
  ea_patterns <- c("effect_allele", "EFFECT_ALLELE", "A1", "a1", "Tested_Allele", "TESTED_ALLELE", "EA", "ea", "ALT", "Allele1", "allele1")
  ea_col <- original_names[original_names %in% ea_patterns][1]
  if (is.na(ea_col)) {
    ea_col <- original_names[grepl("^(A1|a1|effect|tested|EA|ea)", original_names, ignore.case = TRUE)][1]
  }
  
  # Common patterns for other allele
  oa_patterns <- c("other_allele", "OTHER_ALLELE", "A2", "a2", "Other_Allele", "OTHER_ALLELE", "OA", "oa", "ref_allele", "REF", "Allele2", "NonEffectAllele", "allele2")
  oa_col <- original_names[original_names %in% oa_patterns][1]
  if (is.na(oa_col)) {
    oa_col <- original_names[grepl("^(A2|a2|other|ref|OA|oa)", original_names, ignore.case = TRUE)][1]
  }
  
  # Common patterns for effect allele frequency
  eaf_patterns <- c("eaf", "EAF", "freq", "FREQ", "Freq_Tested_Allele", "FREQ_TESTED_ALLELE", "maf", "MAF", "effect_allele_frequency", "POOLED_ALT_AF")
  eaf_col <- original_names[original_names %in% eaf_patterns][1]
  if (is.na(eaf_col)) {
    eaf_col <- original_names[grepl("^(eaf|freq|maf)", original_names, ignore.case = TRUE)][1]
  }
  
  # Print what we found
  cat("Mapping for", exposure_name, ":\n")
  cat("  SNP:", snp_col, "\n")
  cat("  P-value:", pval_col, "\n")
  cat("  Beta:", beta_col, "\n")
  cat("  SE:", se_col, "\n")
  cat("  Effect allele:", ea_col, "\n")
  cat("  Other allele:", oa_col, "\n")
  cat("  EAF:", eaf_col, "\n")
  
  # Apply renaming
  if (!is.na(snp_col)) names(standardized_df)[names(standardized_df) == snp_col] <- "SNP"
  if (!is.na(pval_col)) names(standardized_df)[names(standardized_df) == pval_col] <- "P"
  if (!is.na(beta_col)) names(standardized_df)[names(standardized_df) == beta_col] <- "BETA"
  if (!is.na(se_col)) names(standardized_df)[names(standardized_df) == se_col] <- "SE"
  if (!is.na(ea_col)) names(standardized_df)[names(standardized_df) == ea_col] <- "Tested_Allele"
  if (!is.na(oa_col)) names(standardized_df)[names(standardized_df) == oa_col] <- "Other_Allele"
  if (!is.na(eaf_col)) names(standardized_df)[names(standardized_df) == eaf_col] <- "Freq_Tested_Allele"
  
  # Check if we have the required columns
  required_cols <- c("SNP", "P", "BETA", "SE", "Tested_Allele", "Other_Allele")
  missing_cols <- required_cols[!required_cols %in% names(standardized_df)]
  
  if (length(missing_cols) > 0) {
    cat("WARNING: Missing required columns for", exposure_name, ":", paste(missing_cols, collapse = ", "), "\n")
    cat("Available columns:", paste(names(standardized_df), collapse = ", "), "\n")
  }
  
  return(standardized_df)
}


for (exposure_name in names(exposure_files)) {
  file_path <- exposure_files[[exposure_name]]
  cat("Processing:", exposure_name, "\n")

  # 1. Load and filter exposure data
  sumStats <- fread(file_path)

  sumStats <- standardize_columns(sumStats, exposure_name)
  
  #For Odds Ratio to BETA
  #sumStats$BETA <- log(sumStats$OR)
  #sumStats$SE <- (log(sumStats$OR_95U) - log(sumStats$OR_95L)) / (2 * 1.96)
  #sumStats$Tested_Allele <- sumStats$RISK_ALLELE
  #sumStats$Freq_Tested_Allele <- NA
  
  #Create rsID
  #install.packages("/project/damrauer_scratch/Users/saleemsa/levinmisc/levinmisc-main", repos = NULL, type = "source")
  library(levinmisc)
  
  library(data.table)
  
  #sumStats[, c("chr", "pos") := tstrsplit(variant, ":", fixed = TRUE)[1:2]]
  #sumStats[, pos := as.integer(pos)]
  
  sumStats$P <- as.numeric(sumStats$P)
  sumStats <- sumStats[!is.na(sumStats$P) & sumStats$P < 5e-8]
  
  #sumStats <- annotate_rsids(sumStats, chrom_col = chr, pos_col = pos, dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh37::SNPlocs.Hsapiens.dbSNP144.GRCh37)
  
  
  sumStats$ID <- exposure_name
  sumStats$SNP_clean <- sub(":.*", "", sumStats$SNP)
  
  
  # 2. LD clump
  clumped <- ld_clump(
    dplyr::tibble(rsid = sumStats$SNP, pval = sumStats$P, id = sumStats$ID),
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = "/project/voltron/Resources/1000-Genomes-Phase3/MRC-IEU/EUR"
  )
  
  exposure_dat <- sumStats[sumStats$SNP %in% clumped$rsid, ]
  
  # 3. Format exposure and outcome
  exposure_dat <- exposure_dat %>%
    rename(
      #longSNP = SNP,
      #SNP = SNP_clean,
      beta.exposure = BETA,
      se.exposure = SE,
      effect_allele.exposure = Tested_Allele,
      other_allele.exposure = Other_Allele,
      eaf.exposure = Freq_Tested_Allele,
      pval.exposure = P
    ) %>%
    mutate(exposure = exposure_name, id.exposure = exposure_name)
  
  write.csv(exposure_dat, paste0("/project/damrauer_scratch/Users/saleemsa/ExposureSumStats/exposure_instrument_", exposure_name, ".csv"), row.names = FALSE)

}

