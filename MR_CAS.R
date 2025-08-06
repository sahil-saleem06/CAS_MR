# Load packages
require('tidyverse')
install.packages("TwoSampleMR", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
install.packages('R.utils')
install.packages("devtools")
devtools::install_github("rondolab/MR-PRESSO")

library(R.utils)
library(TwoSampleMR)
library(ieugwasr)
library(data.table)
library(MRPRESSO)
library(janitor)
library(ggplot2)
library(tidyr)
library(dplyr)

# Load outcome data
outcome_file <- "/project/damrauer_scratch/Users/saleemsa/CarotidStenosis/outcome_instrument_CAS"
sumStatsCAS <- fread(outcome_file)

# Define multiple exposure files
exposure_instrument_files <- list(
  #BMI = "/project/damrauer_scratch/Users/saleemsa/ExposureSumStats/exposure_instrument_BMI.csv",
  #CKD = "/project/damrauer_scratch/Users/saleemsa/ExposureSumStats/exposure_instrument_CKD.csv",
  DBP = "/project/damrauer_scratch/Users/saleemsa/ExposureSumStats/exposure_instrument_DBP.csv",
  T2D =  "/project/damrauer_scratch/Users/saleemsa/ExposureSumStats/exposure_instrument_DB2.csv",
  HDL = "/project/damrauer_scratch/Users/saleemsa/ExposureSumStats/exposure_instrument_HDL.csv",
  LDL = "/project/damrauer_scratch/Users/saleemsa/ExposureSumStats/exposure_instrument_LDL.csv",
  #LPA = "/project/damrauer_scratch/Users/saleemsa/ExposureSumStats/exposure_instrument_LPA.csv",
  NONHDL = "/project/damrauer_scratch/Users/saleemsa/ExposureSumStats/exposure_instrument_NONHDL.csv",
  SBP = "/project/damrauer_scratch/Users/saleemsa/ExposureSumStats/exposure_instrument_SBP.csv",
  TC = "/project/damrauer_scratch/Users/saleemsa/ExposureSumStats/exposure_instrument_TC.csv"
)

# Final result list
final_results_list <- list()

for (exposure_name in names(exposure_instrument_files)) {
  file_path <- exposure_instrument_files[[exposure_name]]
  cat("Processing:", exposure_name, "\n")
  
  # 1. Load and filter exposure data
  exposure_dat <- fread(file_path)
  
  #2. Format outcome data to exposure
  outcome_dat <- sumStatsCAS %>%
    filter(rsid %in% exposure_dat$SNP) %>%
    rename(
      SNP = rsid,
      beta.outcome = Effect,
      se.outcome = StdErr,
      effect_allele.outcome = EFFECT_ALLELE,
      other_allele.outcome = OTHER_ALLELE,
      pval.outcome = PVALUE,
      eaf.outcome = EFF_ALL_FREQ
    ) %>%
    mutate(id.outcome = "CAS", outcome = "CAS")
  
  # 4. Harmonize
  dat <- harmonise_data(exposure_dat, outcome_dat)
  
  # 5. MR
  res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
  res_single <- mr_singlesnp(dat)
  
  # 6 . PRESSO
  if (nrow(dat) >= 5 && length(unique(dat$SNP)) >= 5) {
    
    mr_presso_result <- tryCatch({
      mr_presso(
        BetaOutcome = "beta.outcome",
        BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome",
        SdExposure = "se.exposure",
        OUTLIERtest = TRUE,
        DISTORTIONtest = TRUE,
        data = dat,
        NbDistribution = 2000,
        SignifThreshold = 0.05
      )
    }, error = function(e) {
      warning(paste("MR-PRESSO failed for", exposure_name, ":", e$message))
      return(NULL)
    })
    
    if (!is.null(mr_presso_result)) {
      # Main MR PRESSO results
      main_df <- tryCatch({
        mr_presso_result[["Main MR results"]] %>%
          dplyr::rename(
            method = `MR Analysis`,
            b = `Causal Estimate`,
            se = Sd,
            pval = `P-value`
          ) %>%
          mutate(step = ifelse(method == "Raw", "Before PRESSO", "After PRESSO")) %>%
          select(step, method, b, se, pval)
      }, error = function(e) {
        warning(paste("Failed to process main PRESSO table for", exposure_name))
        NA
      })
      
      # Helper function
      parse_pval <- function(x) {
        x <- as.character(x)
        x <- ifelse(grepl("^<", x), sub("^<", "", x), x)
        suppressWarnings(as.numeric(x))
      }
      
      # Test results block
      tests_df <- tryCatch({
        global_test <- as.data.frame(mr_presso_result$`MR-PRESSO results`$`Global Test`)
        outlier_test <- as.data.frame(mr_presso_result$`MR-PRESSO results`$`Outlier Test`)
        distortion_test <- as.data.frame(mr_presso_result$`MR-PRESSO results`$`Distortion Test`)
        
        if (nrow(global_test) > 0) global_test$test <- "Global Test"
        if (nrow(outlier_test) > 0) outlier_test$test <- "Outlier Test"
        if (nrow(distortion_test) > 0) distortion_test$test <- "Distortion Test"
        
        # Clean p-values
        global_test$Pvalue <- parse_pval(global_test$Pvalue)
        outlier_test$Pvalue <- parse_pval(outlier_test$Pvalue)
        distortion_test$Pvalue <- parse_pval(distortion_test$Pvalue)
        
        bind_rows(global_test, outlier_test, distortion_test) %>%
          select(test, everything())
      }, error = function(e) {
        warning(paste("Failed to process PRESSO test tables for", exposure_name))
        NA
      })
      
    } else {
      main_df <- NA
      tests_df <- NA
    }
  } else {
    warning(paste("Skipping MR-PRESSO for", exposure_name, "due to insufficient SNPs"))
    main_df <- NA
    tests_df <- NA
  }
  
  # Save results
  final_results_list[[exposure_name]] <- list(
    mr_results = res,
    mr_single_snp = res_single,
    presso_main = main_df,
    presso_tests = tests_df
  )
}

#Summary Table
summary_table <- list()

for (exp_name in names(final_results_list)) {
  mr_res <- final_results_list[[exp_name]]$mr_results
  presso_main <- final_results_list[[exp_name]]$presso_main
  presso_tests <- final_results_list[[exp_name]]$presso_tests
  
  # Get IVW row
  ivw_row <- mr_res %>% filter(method == "Inverse variance weighted")
  
  if (nrow(ivw_row) == 1) {
    b <- ivw_row$b
    se <- ivw_row$se
    pval <- ivw_row$pval
    nsnp <- ivw_row$nsnp
    lo_ci <- b - 1.96 * se
    up_ci <- b + 1.96 * se
    or <- exp(b)
    or_lci95 <- exp(lo_ci)
    or_uci95 <- exp(up_ci)
    ci_range <- paste0(round(or_lci95, 2), "-", round(or_uci95, 2))
    
    # PRESSO
    global <- presso_tests %>% filter(test == "Global Test")
    distort <- presso_tests %>% filter(test == "Distortion Test")
    
    pval_presso_global <- ifelse(nrow(global) > 0, global$Pvalue[1], NA)
    pval_presso_distort <- ifelse(nrow(distort) > 0, distort$Pvalue[1], NA)
    
    # Add to summary table
    summary_table[[exp_name]] <- data.frame(
      outcome = "CAS",
      exposure = exp_name,
      method = "Inverse variance weighted",
      nsnp = nsnp,
      b = round(b, 3),
      se.mr = round(se, 4),
      pval.mr = signif(pval, 3),
      lo_ci = round(lo_ci, 3),
      up_ci = round(up_ci, 3),
      or = round(or, 3),
      or_lci95 = round(or_lci95, 2),
      or_uci95 = round(or_uci95, 2),
      ci_range = ci_range,
      pval.presso.global = signif(pval_presso_global, 3),
      pval.presso.distortion = signif(pval_presso_distort, 3),
      stringsAsFactors = FALSE
    )
  }
}

# Combine all exposures into a single table
final_summary_df <- bind_rows(summary_table)
final_summary_df$p.corrected <- signif(p.adjust(final_summary_df$pval.mr, method = "fdr"), 3)

# Save or view
#write.csv(final_summary_df, "mr_cas_summary.csv", row.names = FALSE)
print(final_summary_df)

#Forest Plot
if (nrow(final_summary_df) > 0) {
  library(ggplot2)
  
  forest_plot <- final_summary_df %>%
    mutate(
      exposure = factor(exposure, levels = exposure[order(or)]),
      significant = p.corrected < 0.05
    ) %>%
    ggplot(aes(x = or, y = exposure, color = significant)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = or_lci95, xmax = or_uci95), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "blue")) +
    scale_x_log10() +
    labs(
      title = "Mendelian Randomization Results: Risk Factors for Carotid Stenosis",
      x = "Odds Ratio (95% CI)",
      y = "Exposure",
      color = "FDR < 0.05"
    ) +
    theme_minimal()
  
  #ggsave("CAS_forest_plot.png", forest_plot, width = 10, height = 6)
  print(forest_plot)
}



