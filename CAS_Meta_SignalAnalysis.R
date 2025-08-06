require('tidyverse')
library(data.table)
library(httr)
library(jsonlite)
library(dplyr)
library(biomaRt)
library(ggplot2)

# Set input meta-analysis file
meta_file <- "/project/damrauer_scratch/Users/nunez997/METAANALYSIS1.TBL"
#Changed for CAS

# Gene list
gene_list <- c("NSD2", "KMT2A", "KDM5B", "KDM6B", "SETDB2", "SETD4", 
               "HDAC11", "HDAC3", "PRMT1", "KDM5C", "STAT3", "DLL4", 
               "NOTCH1", "NOTCH2", "JAG1", "JAG2", "PTGS2", "IL17A", "TGFBR2")

# Get coordinates of each gene
gene_coords <- fread("/project/damrauer_scratch/Users/saleemsa/PAD_Signals/gene_coords.csv")

# Get column names from meta-analysis file
header_line <- system(paste("head -n 1", meta_file), intern = TRUE)
colnames <- strsplit(header_line, "\t")[[1]]
print("Available columns:")
print(colnames)

# Function to extract gene region from meta-analysis file
extract_gene_region_meta <- function(gene, chr, start, end, meta_file_path, colnames) {
  # Remove "chr" prefix for matching with meta file (which uses just chromosome numbers)
  chr_num <- gsub("chr", "", chr)
  
  df <- fread(cmd = paste0("awk -F '\t' '$1==\"", chr_num, 
                           "\" && $2>=", start, " && $2<=", end, "' ", meta_file_path),
              col.names = colnames)
  
  if (nrow(df) > 0) {
    # Rename columns to match your plotting function
    df <- df %>%
      rename(
        CHR = Chromosome,
        POS = Position,
        SNP = MarkerName,
        p.value = `P-value`
      ) %>%
      mutate(
        CHR = paste0("chr", CHR),  # Add back "chr" prefix
        gene = gene
      )
  }
  return(df)
}

# Function to do Manhattan plot of single search (adapted for meta-analysis)
plot_gene_signal <- function(df, gene_label) {
  ggplot(df, aes(x = POS, y = -log10(p.value))) +
    geom_point(alpha = 0.6, color = "darkgreen") +
    geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed", size = 0.8) +
    labs(
      title = paste("Meta-Analysis Signal near", gene_label),
      x = "Genomic Position",
      y = "-log10(p-value)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_line(color = "gray90")
    )
}

# Create output directory(Changed for CAS)
base_output_dir <- "/project/damrauer_scratch/Users/saleemsa/CarotidStenosis/GeneSignals/"
dir.create(base_output_dir, showWarnings = FALSE, recursive = TRUE)

# Loop through genes and create plots
all_plots <- list()

for (i in 1:nrow(gene_coords)) {
  gene <- gene_coords$hgnc_symbol[i]
  chr <- gene_coords$chromosome_name[i]
  start <- gene_coords$region_start[i]
  end <- gene_coords$region_end[i]
  
  # Extract gene region from meta-analysis
  df <- extract_gene_region_meta(gene, chr, start, end, meta_file, colnames)
  
  # Create and save plot
  if (nrow(df) > 0) {
    gene_label <- paste(gene, "Meta", sep = "_")
    p <- plot_gene_signal(df, gene_label)
    all_plots[[gene]] <- p
    
    # Save individual plot
    ggsave(
      filename = file.path(base_output_dir, paste0(gene_label, "_plot.png")),
      plot = p,
      width = 8, height = 5, dpi = 300
    )
    
    message(paste("Created plot for", gene, "with", nrow(df), "SNPs"))
  } else {
    message(paste("No data for", gene, "in meta-analysis"))
  }
}

# Print summary
message(paste("Generated", length(all_plots), "plots in", base_output_dir))


#After looking at graph
#Only one signal in STAT3 gene

stat3_coords <- gene_coords[gene_coords$hgnc_symbol == "STAT3", ]
stat3_data <- extract_gene_region_meta("STAT3", 
                                       stat3_coords$chromosome_name, 
                                       stat3_coords$region_start, 
                                       stat3_coords$region_end, 
                                       meta_file, 
                                       colnames)

peak_variant <- stat3_data %>%
  filter(`p.value` == min(`p.value`, na.rm = TRUE))

print("Peak variant(s) in STAT3 region:")
print(peak_variant)

# Get variants within 1 order of magnitude of the peak
peak_pval <- min(stat3_data$`p.value`, na.rm = TRUE)
near_peak <- stat3_data %>%
  filter(`p.value` <= peak_pval * 10) %>%
  arrange(`p.value`)

print("Variants near the peak signal:")
print(near_peak)