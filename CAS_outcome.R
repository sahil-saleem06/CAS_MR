install.packages("/project/damrauer_scratch/Users/saleemsa/levinmisc/levinmisc-main", repos = NULL, type = "source")

library(levinmisc)
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(glue)
library(tibble)

library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

chainfile <- ("/project/voltron/Resources/liftOver/hg38ToHg19.over.chain")

carotidMeta <- fread("/project/damrauer_scratch/Users/nunez997/METAANALYSIS1.TBL")

#sumStats <- carotidMeta[!is.na(carotidMeta[["P-value"]]) & carotidMeta[["P-value"]] < 5e-8, ] %>% as.data.frame()


step1 <- as.data.frame(carotidMeta)

# 1. Rename chromosome and add start/end columns
step1 <- step1 %>%
  dplyr::rename(chrom = Chromosome) %>%
  dplyr::mutate(start = Position, end = start)
step1$chrom <- paste0("chr", step1$chrom)

print(head(step1))  # Check your dataframe here

# 2. Make GRanges object
step2 <- makeGRangesFromDataFrame(step1,
                                  keep.extra.columns = TRUE,
                                  seqnames.field = "chrom",
                                  start.field = "start",
                                  end.field = "end",
                                  starts.in.df.are.0based = FALSE)

print(step2)  # Check GRanges object summary

# 3. Import chain file
chain <- import.chain(chainfile)
print(chain)  # Check chain import

# 4. Perform liftOver
step4 <- liftOver(step2, chain)
print(step4)  # This is a list; check length and content

# 5. Unlist liftOver results (flatten list)
step5 <- unlist(step4)
print(step5)

# 6. Convert to tibble and add chr/pos columns
step6 <- as_tibble(step5) %>%
  mutate(CHR = as.character(seqnames), POSITION = start)

print(head(step6))

# 7. Select and rename columns
step7 <- step6 %>%
  select(CHR = seqnames,
         POSITION,
         EFFECT_ALLELE = Allele1,
         OTHER_ALLELE = Allele2,
         EFF_ALL_FREQ = Freq1,
         Effect,
         StdErr,
         PVALUE = `P.value`)
print(head(step7))

  
step7$CHR <- as.integer(gsub("^chr", "", step7$CHR))
  
carotid37_ID <- step7 %>%
  filter(!is.na(CHR) & !is.na(POSITION))

carotid37_ID <- levinmisc::annotate_rsids(carotid37_ID, chrom_col = CHR, pos_col = POSITION, dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh37::SNPlocs.Hsapiens.dbSNP144.GRCh37)

carotid37_ID <- carotid37_ID %>%
  filter(!is.na(rsid))
#Lost 11828404 SNPs

write.csv(carotid37_ID, "carotid37", row.names = FALSE)
