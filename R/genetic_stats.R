library(gdsfmt)
library(GenomicRanges)
library(rtracklayer)

# --- 1. Configuration ---
# List of all your GDS files
gds_files <- c(
  "data/ACT197_SVs.gds",
  "data/ACT197_biallelic.gds",
  "data/ACT197_biallelic_PRUNED.gds",
  "data/ACT197_biallelic_codingonly.gds",
  "data/LDP324_nofiltered.gds",
  "data/Merged_Analysis_RealCoords.gds" 
)

# Reference Annotation for Gene Coverage
gff_fn <- "data/Lcu.1GRN.genes_description.sorted.gff3.gz"

# --- 2. Load Gene Models ---
cat("Loading GFF3 annotation...\n")
gff_data <- import.gff3(gff_fn)
gene_ranges <- gff_data[gff_data$type == "gene"]

# --- 3. Define the Summary Function ---
get_gds_stats <- function(file_path, gene_gr) {
  
  if(!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  cat(paste("Processing:", basename(file_path), "...\n"))
  
  # Open GDS agnostically (works for SeqArray and SNPRelate)
  f <- openfn.gds(file_path)
  
  # Find nodes (handles both naming conventions)
  node_id  <- index.gdsn(f, "snp.id", silent=TRUE)
  if (is.null(node_id)) node_id <- index.gdsn(f, "variant.id", silent=TRUE)
  
  node_chr <- index.gdsn(f, "snp.chromosome", silent=TRUE)
  if (is.null(node_chr)) node_chr <- index.gdsn(f, "chromosome", silent=TRUE)
  
  node_pos <- index.gdsn(f, "snp.position", silent=TRUE)
  if (is.null(node_pos)) node_pos <- index.gdsn(f, "position", silent=TRUE)
  
  # Extract data
  ids   <- read.gdsn(node_id)
  chrom <- read.gdsn(node_chr)
  pos   <- read.gdsn(node_pos)
  
  closefn.gds(f)
  
  # Total Variants
  total_vars <- length(ids)
  
  # Filter for the 7 main chromosomes 
  # (Matches "Lcu.1GRN.Chr1" to 7, or just "1" to "7")
  chrom_str <- as.character(chrom)
  chr_mask  <- grepl("Chr[1-7]$|^[1-7]$", chrom_str)
  
  vars_on_chr <- sum(chr_mask)
  avg_vars_chr <- vars_on_chr / 7
  
  # Subset to main chromosomes for spacing/density/coverage
  main_chroms <- chrom_str[chr_mask]
  main_pos    <- pos[chr_mask]
  
  # Calculate physical span of the data to estimate density
  # Span = sum of (max_position - min_position) per chromosome
  chr_spans <- tapply(main_pos, main_chroms, function(x) max(x, na.rm=TRUE) - min(x, na.rm=TRUE))
  total_span_bp <- sum(as.numeric(chr_spans), na.rm=TRUE)
  
  # Spacing and Density
  # Avoid division by zero if a file is empty
  if(vars_on_chr > 0) {
    avg_spacing_kb <- (total_span_bp / vars_on_chr) / 1000
    density_mb <- vars_on_chr / (total_span_bp / 1000000)
  } else {
    avg_spacing_kb <- NA; density_mb <- NA
  }
  
  # Gene Coverage
  # Ensure chromosome names match the GFF3 format exactly ("Lcu.1GRN.ChrX")
  mapped_chrs <- ifelse(grepl("^[1-7]$", main_chroms), 
                        paste0("Lcu.1GRN.Chr", main_chroms), 
                        main_chroms)
  
  variant_gr <- GRanges(seqnames = mapped_chrs,
                        ranges = IRanges(start = main_pos, width = 1))
  
  # Find variants inside genes
  overlaps <- findOverlaps(variant_gr, gene_gr)
  unique_genes_hit <- length(unique(subjectHits(overlaps)))
  
  gene_cov_pct <- (unique_genes_hit / length(gene_gr)) * 100
  
  # Return row
  return(data.frame(
    Dataset = basename(file_path),
    Total_Variants = total_vars,
    Variants_Chr1_7 = vars_on_chr,
    Avg_Vars_per_Chr = round(avg_vars_chr, 0),
    Avg_Spacing_kb = round(avg_spacing_kb, 2),
    Marker_Density_SNP_Mb = round(density_mb, 2),
    Gene_Coverage_Pct = round(gene_cov_pct, 2)
  ))
}

# --- 4. Run Analysis & Compile Table ---
results_list <- lapply(gds_files, get_gds_stats, gene_gr = gene_ranges)

# Combine into a single dataframe
final_table <- do.call(rbind, results_list)

# View the table
print(final_table)

# Save to CSV for easy inclusion in your manuscript
write.csv(final_table, "Results/Variant_Summary_Statistics.csv", row.names = FALSE)
cat("Done! Results saved to Variant_Summary_Statistics.csv\n")

#HQ Figures####
# =============================================================================
# ── PUBLICATION FIGURE: WINDOWED GENETIC STATS (CHR ~ METRIC) ────────────────
# =============================================================================
cat("Extracting windowed statistics for ACT and LDP panels...\n")

# Load required libraries
library(gdsfmt)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)    # Required for stat_compare_means
library(patchwork) # Required to stitch the columns together

# --- 1. Configuration ---
# Use the individual panel GDS files so we count SNPs actually present in each
files <- list(
  LDP = "data/LDP324_nofiltered.gds",
  ACT = "data/ACT197_biallelic.gds"
)

window_size <- 1000000 # 1 Mb bins
results <- list()

# --- 2. Calculate Stats per 1 Mb Window ---
for(panel in names(files)) {
  
  if(!file.exists(files[[panel]])) stop(paste("File not found:", files[[panel]]))
  cat(paste("  Processing", panel, "...\n"))
  
  f <- openfn.gds(files[[panel]])
  
  # Safely locate chromosome and position nodes
  chr_node <- index.gdsn(f, "snp.chromosome", silent=TRUE)
  if(is.null(chr_node)) chr_node <- index.gdsn(f, "chromosome")
  
  pos_node <- index.gdsn(f, "snp.position", silent=TRUE)
  if(is.null(pos_node)) pos_node <- index.gdsn(f, "position")
  
  chr <- read.gdsn(chr_node)
  pos <- read.gdsn(pos_node)
  closefn.gds(f)
  
  df <- data.frame(chr = as.character(chr), pos = pos, stringsAsFactors = FALSE)
  
  # Filter to main 7 chromosomes and calculate bins
  df_binned <- df %>%
    filter(grepl("Chr[1-7]$|^[1-7]$", chr)) %>%
    mutate(
      chr_num = as.numeric(gsub(".*Chr", "", chr)),
      chr_label = factor(paste0("Chr", chr_num), levels = paste0("Chr", 1:7)),
      bin_mb = floor(pos / window_size)
    ) %>%
    group_by(chr_label, bin_mb) %>%
    summarise(
      SNP_Density = n(),                               # SNPs per Mb
      Avg_Spacing_kb = (window_size / n()) / 1000,     # Spacing in Kb
      .groups = "drop"
    ) %>%
    mutate(Panel = panel)
  
  results[[panel]] <- df_binned
}

# --- 3. Prepare Data for Plotting ---
cat("Preparing data for faceted plotting...\n")

plot_df <- bind_rows(results)

# Ensure Panel is ordered with LDP on the left and ACT on the right
plot_df$Panel <- factor(plot_df$Panel, levels = c("LDP", "ACT"))

# Define standard high-contrast journal colors
panel_colors <- c("LDP" = "#e41a1c", "ACT" = "#377eb8")

# Define the pairwise comparison for the p-values
my_comparisons <- list(c("LDP", "ACT"))

# --- 4. Build Plot A: SNP Density Column ---
p_density <- ggplot(plot_df, aes(x = Panel, y = SNP_Density, fill = Panel)) +
  geom_violin(alpha = 0.4, color = NA, trim = TRUE) +
  geom_boxplot(width = 0.2, color = "black", outlier.size = 0.5, outlier.alpha = 0.3, alpha = 0.8) +
  
  # Facet by Chromosome (Rows)
  facet_grid(chr_label ~ .) +
  
  # Add Significance Brackets
  stat_compare_means(comparisons = my_comparisons, 
                     method = "wilcox.test", 
                     label = "p.format", 
                     size = 3.5, 
                     vjust = 0.5) +
  
  scale_fill_manual(values = panel_colors) +
  labs(title = "SNP Density", x = NULL, y = "SNPs per 1 Mb Window") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.background = element_rect(fill = "grey93", color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(face = "bold", color = "black"),
    panel.grid.minor = element_blank()
  )

# --- 5. Build Plot B: Average Spacing Column ---
# We use log10 scale here because spacing can have extreme outliers in centromeres
p_spacing <- ggplot(plot_df, aes(x = Panel, y = Avg_Spacing_kb, fill = Panel)) +
  geom_violin(alpha = 0.4, color = NA, trim = TRUE) +
  geom_boxplot(width = 0.2, color = "black", outlier.size = 0.5, outlier.alpha = 0.3, alpha = 0.8) +
  
  # Facet by Chromosome (Rows)
  facet_grid(chr_label ~ .) +
  
  # Add Significance Brackets
  stat_compare_means(comparisons = my_comparisons, 
                     method = "wilcox.test", 
                     label = "p.format", 
                     size = 3.5, 
                     vjust = 0.5) +
  
  scale_fill_manual(values = panel_colors) +
  scale_y_log10() + # Log scale is crucial for physical spacing visualization
  labs(title = "Average Spacing", x = NULL, y = "Log10 Spacing (Kb)") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.background = element_rect(fill = "grey93", color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(face = "bold", color = "black"),
    panel.grid.minor = element_blank()
  )

# --- 6. Assemble and Save Final Figure ---
cat("Stitching columns together with patchwork...\n")

# Combine them side-by-side
final_plot <- p_density | p_spacing

if(!dir.exists("Results")) dir.create("Results")

# Save as a tall, full-page A4 figure
ggsave("Results/Genetic_Stats_Faceted_Publication.png", 
       plot = final_plot, 
       width = 17, 
       height = 24, # 24 cm height gives all 7 chromosome rows plenty of vertical room
       units = "cm", 
       dpi = 600, 
       bg = "white")

cat("Done! High-resolution faceted figure saved to Results/Genetic_Stats_Faceted_Publication.png\n")