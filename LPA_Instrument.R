library(R.utils)
library(ieugwasr)
library(data.table)
library(janitor)
library(tidyr)
library(dplyr)
library(tidyverse)


file_path <- "/project/damrauer_scratch/summary_stats/LPA/30790_raw.gwas.imputed_v3.both_sexes.varorder.tsv.bgz"
exposure_name <- "LPA"

sumStats <- fread(file_path)

sumStats <- sumStats[
  !is.na(beta) & !is.na(se) & low_confidence_variant == FALSE
]

sumStats[, c("chr", "pos", "ref_allele", "alt_allele") := tstrsplit(variant, ":", fixed = TRUE)]

# Assuming `minor_allele` is always the alt_allele and also the effect allele
sumStats[, Tested_Allele := minor_allele]
sumStats[, Other_Allele := ifelse(ref_allele == minor_allele, alt_allele, ref_allele)]

setnames(sumStats, old = c("variant", "pval", "beta", "se", "minor_AF"), 
         new = c("SNP", "P", "BETA", "SE", "Freq_Tested_Allele"))

keep_cols <- c("SNP", "chr", "pos", "Tested_Allele", "Other_Allele", 
               "Freq_Tested_Allele", "BETA", "SE", "P")
sumStats <- sumStats[, ..keep_cols]

#create rsids
sumStats <- annotate_rsids(sumStats, chrom_col = chr, pos_col = pos, dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh37::SNPlocs.Hsapiens.dbSNP144.GRCh37)

n_before <- nrow(sumStats)

sumStats <- sumStats[
  !is.na(rsid)
]

n_after <- nrow(sumStats)

n_removed <- n_before - n_after
cat("Removed", n_removed, "rows with NA rsid values\n") #Removed 1150955 rows

#pvalue filter 
n_before <- nrow(sumStats)

sumStats$P <- as.numeric(sumStats$P)
sumStats <- sumStats[!is.na(sumStats$P) & sumStats$P < 5e-8]

n_after <- nrow(sumStats)

n_removed <- n_before - n_after
cat("Removed", n_removed, "rows with pvalue filter\n") #Removed 1150955 row

#LD clump
clumped <- ld_clump(
  dplyr::tibble(rsid = sumStats$rsid, pval = sumStats$P, id = sumStats$ID),
  plink_bin = genetics.binaRies::get_plink_binary(),
  bfile = "/project/voltron/Resources/1000-Genomes-Phase3/MRC-IEU/EUR"
)

exposure_dat <- sumStats[sumStats$rsid %in% clumped$rsid, ]

#Format exposure and outcome
exposure_dat <- exposure_dat %>%
  rename(
    #longSNP = SNP,
    SNP = rsid,
    beta.exposure = BETA,
    se.exposure = SE,
    effect_allele.exposure = Tested_Allele,
    other_allele.exposure = Other_Allele,
    eaf.exposure = Freq_Tested_Allele,
    pval.exposure = P
  ) %>%
  mutate(exposure = exposure_name, id.exposure = exposure_name)

write.csv(exposure_dat, paste0("/project/damrauer_scratch/Users/saleemsa/ExposureSumStats/exposure_instrument_", exposure_name, ".csv"), row.names = FALSE)





