# Load necessary libraries
library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(dplyr)
library(tidyr)

cat("Starting Haplotype Block Analysis...\n")

#erosion vs introgression####
# --- 1. Configuration & Data Loading ---
gds_file <- "data/Merged_Analysis_RealCoords.gds"
erosion_file <- "Results/DeltaHe_Erosion_SignificantRegions.csv"
introgression_file <- "Results/DeltaHe_Introgression_SignificantRegions.csv"

# Check if files exist
if(!file.exists(gds_file) | !file.exists(erosion_file) | !file.exists(introgression_file)){
  stop("Input files not found. Check your paths.")
}

erosion_df <- read.csv(erosion_file, stringsAsFactors = FALSE)
introgression_df <- read.csv(introgression_file, stringsAsFactors = FALSE)

# Open GDS
f <- snpgdsOpen(gds_file)

# Define populations
samp_ids <- read.gdsn(index.gdsn(f, "sample.id"))
act_samples <- samp_ids[grepl("_ACT", samp_ids)]
act_samples <- act_samples[1:187] #Remove LR-86 lines!!
ldp_samples <- setdiff(samp_ids, act_samples) 

# --- 2. Define the Haplotype Block Function ---
# This function finds the maximum physical distance (bp) between any two SNPs 
# in the region that maintain strong Linkage Disequilibrium (r^2 > 0.8)
calc_max_ld_block <- function(gds_obj, target_samples, chrom, start_pos, end_pos, r2_threshold = 0.8) {
  
  # Get all SNPs in the GDS
  snp_chr <- read.gdsn(index.gdsn(gds_obj, "snp.chromosome"))
  snp_pos <- read.gdsn(index.gdsn(gds_obj, "snp.position"))
  snp_id  <- read.gdsn(index.gdsn(gds_obj, "snp.id"))
  
  # Filter for the specific window
  window_idx <- which(snp_chr == chrom & snp_pos >= start_pos & snp_pos <= end_pos)
  
  # Need at least 2 SNPs to form a block
  if(length(window_idx) < 2) return(0) 
  
  win_snp_ids <- snp_id[window_idx]
  win_snp_pos <- snp_pos[window_idx]
  
  # Calculate pairwise LD (r^2) matrix
  # suppressWarnings is used to keep console clean from monomorphic SNP warnings
  ld_res <- suppressWarnings(
    snpgdsLDMat(gds_obj, sample.id = target_samples, snp.id = win_snp_ids, method = "r", slide = -1, verbose = FALSE)
  )
  
  ld_mat <- ld_res$LD^2 # Square r to get r^2
  
  # Find the maximum physical distance between any two SNPs where r^2 > threshold
  max_block_size <- 0
  block_start <- NA
  block_end <- NA
  
  # Matrix is symmetric, we only need to check the upper triangle
  for (i in 1:(nrow(ld_mat)-1)) {
    for (j in (i+1):ncol(ld_mat)) {
      if (!is.na(ld_mat[i,j]) && ld_mat[i,j] >= r2_threshold) {
        dist <- win_snp_pos[j] - win_snp_pos[i]
        if (dist > max_block_size) {
          max_block_size <- dist
          block_start <- win_snp_pos[i]
          block_end <- win_snp_pos[j]
        }
      }
    }
  }
  
  return(list(size = max_block_size, start = block_start, end = block_end))
}

# --- 3. Process Regions ---
process_regions <- function(regions_df, region_type) {
  cat(sprintf("Processing %d %s regions...\n", nrow(regions_df), region_type))
  
  results <- list()
  
  for(i in 1:nrow(regions_df)) {
    r <- regions_df[i, ]
    
    # Calculate for ACTIVATE
    act_block <- calc_max_ld_block(f, act_samples, r$chr, r$start, r$stop)
    # Calculate for LDP
    ldp_block <- calc_max_ld_block(f, ldp_samples, r$chr, r$start, r$stop)
    
    results[[i]] <- data.frame(
      Region_ID = paste(r$chr, r$start, r$stop, sep="_"),
      Chromosome = r$chr,
      Window_Start = r$start,
      Window_Stop = r$stop,
      Region_Type = region_type,
      ACT_Block_Size_bp = act_block$size,
      ACT_Block_Start = act_block$start,
      ACT_Block_End = act_block$end,
      LDP_Block_Size_bp = ldp_block$size,
      LDP_Block_Start = ldp_block$start,
      LDP_Block_End = ldp_block$end
    )
  }
  return(do.call(rbind, results))
}

df_erosion <- process_regions(erosion_df, "Erosion (Top 1%)")
df_introgression <- process_regions(introgression_df, "Introgression (Bottom 1%)")

# Combine datasets
master_blocks <- rbind(df_erosion, df_introgression)

snpgdsClose(f)

# --- 4. Export CSV and BED Files ---
cat("Exporting data...\n")

if(!dir.exists("Results/Haplotypes")) dir.create("Results/Haplotypes", recursive = TRUE)

# Master CSV
write.csv(master_blocks, "Results/Haplotypes/Haplotype_Block_Lengths_erosion vs intro.csv", row.names = FALSE)

# Export BED file (Chromosome, Start, End, Name) for ACTIVATE blocks
bed_act <- master_blocks %>%
  filter(!is.na(ACT_Block_Start)) %>%
  mutate(Name = paste(Region_Type, "ACT_Block", sep="_")) %>%
  select(Chromosome, ACT_Block_Start, ACT_Block_End, Name)
write.table(bed_act, "Results/Haplotypes/ACTIVATE_Blocks.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# --- 5. Statistical Comparison ---
cat("\n--- Statistical Summary (Wilcoxon Rank-Sum Tests) ---\n")

# Restructure data for easier plotting and stats
plot_df <- master_blocks %>%
  select(Region_ID, Region_Type, ACT = ACT_Block_Size_bp, LDP = LDP_Block_Size_bp) %>%
  pivot_longer(cols = c("ACT", "LDP"), names_to = "Panel", values_to = "Block_Size_bp") %>%
  mutate(Block_Size_Kb = Block_Size_bp / 1000)

# 5a. Compare Panels within Eroded Regions
eroded_data <- plot_df %>% filter(Region_Type == "Erosion (Top 1%)")
test_eroded <- wilcox.test(Block_Size_Kb ~ Panel, data = eroded_data)
cat(sprintf("Erosion Regions: ACT median = %.1f Kb, LDP median = %.1f Kb (p-value = %.2e)\n", 
            median(eroded_data$Block_Size_Kb[eroded_data$Panel=="ACT"], na.rm=T),
            median(eroded_data$Block_Size_Kb[eroded_data$Panel=="LDP"], na.rm=T),
            test_eroded$p.value))

# 5b. Compare Panels within Introgressed Regions
intro_data <- plot_df %>% filter(Region_Type == "Introgression (Bottom 1%)")
test_intro <- wilcox.test(Block_Size_Kb ~ Panel, data = intro_data)
cat(sprintf("Introgression Regions: ACT median = %.1f Kb, LDP median = %.1f Kb (p-value = %.2e)\n", 
            median(intro_data$Block_Size_Kb[intro_data$Panel=="ACT"], na.rm=T),
            median(intro_data$Block_Size_Kb[intro_data$Panel=="LDP"], na.rm=T),
            test_intro$p.value))

# --- 6. Publication-Ready Plotting ---
cat("Generating plots...\n")

# Custom labels for the facets
facet_labels <- c(
  "Erosion (Top 1%)" = "Eroded Regions (Sweeps)",
  "Introgression (Bottom 1%)" = "Introgressed Regions (Novel Diversity)"
)

p_blocks <- ggplot(plot_df, aes(x = Panel, y = Block_Size_Kb, fill = Panel)) +
  # Violin plot for distribution shape
  geom_violin(alpha = 0.5, color = NA, trim = FALSE) +
  # Boxplot for median and quartiles
  geom_boxplot(width = 0.2, color = "black", alpha = 0.8, outlier.shape = NA) +
  # Jittered points to show actual data density
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.6, aes(color = Panel)) +
  facet_wrap(~ Region_Type, labeller = as_labeller(facet_labels)) +
  scale_fill_manual(values = c("ACT" = "darkred", "LDP" = "dodgerblue")) +
  scale_color_manual(values = c("ACT" = "darkred", "LDP" = "dodgerblue")) +
  labs(
    title = "Haplotype Block Expansion Driven by Selection",
    subtitle = "Maximum contiguous linkage (r² > 0.8) within the top 1% selective sweeps and introgressions",
    x = "Lentil Panel",
    y = "Haplotype Block Length (Kb)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "gray30"),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("Results/Haplotypes/Haplotype_Block_Comparison_erosion vs intro.png", plot = p_blocks, width = 9, height = 6, dpi = 600)

cat("Done! Check the 'Results/Haplotypes' directory for your figures, CSVs, and BED files.\n")

#breeding vs adaptation####
# --- 1. Configuration & Data Loading ---
gds_file <- "data/Merged_Analysis_RealCoords.gds"
breeding_file <- "Results/XPCLR_Scenario2_Breeding_SignificantRegions.csv"
adaptation_file <- "Results/XPCLR_Scenario1_Adaptation_SignificantRegions.csv"

# Check if files exist
if(!file.exists(gds_file) | !file.exists(breeding_file) | !file.exists(adaptation_file)){
  stop("Input files not found. Check your paths.")
}

breeding_df <- read.csv(breeding_file, stringsAsFactors = FALSE)
adaptation_df <- read.csv(adaptation_file, stringsAsFactors = FALSE)

# Open GDS to get the exact sample IDs present in the genotype matrix
f <- snpgdsOpen(gds_file)
gds_samples <- read.gdsn(index.gdsn(f, "sample.id"))

# Define the exact sample lists based on your XP-CLR M&M
# 1. ACT_187 (Remove the 10 RILs - assuming RILs have a specific naming pattern like "LR-86")
act_187_samples <- readLines("data/ACT187_samples.txt")

# 2. LDP_227 Non-Temperate (Adaptation Baseline)
# Assuming your metadata has a column 'Panel' and 'K_Group'
ldp_227_samples <- readLines("data/LDP227_Non_Temperate.txt")

# 3. LDP_97 Temperate (Breeding Baseline)
ldp_97_samples <- readLines("data/LDP97_Temperate.txt")

cat(sprintf("Sample sizes: ACT_187 (n=%d), LDP_227 (n=%d), LDP_97 (n=%d)\n", 
            length(act_187_samples), length(ldp_227_samples), length(ldp_97_samples)))

# --- 2. Define the Haplotype Block Function ---
# This function finds the maximum physical distance (bp) between any two SNPs 
# in the region that maintain strong Linkage Disequilibrium (r^2 > 0.8)
calc_max_ld_block <- function(gds_obj, target_samples, chrom, start_pos, end_pos, r2_threshold = 0.8) {
  
  # Get all SNPs in the GDS
  snp_chr <- read.gdsn(index.gdsn(gds_obj, "snp.chromosome"))
  snp_pos <- read.gdsn(index.gdsn(gds_obj, "snp.position"))
  snp_id  <- read.gdsn(index.gdsn(gds_obj, "snp.id"))
  
  # Filter for the specific window
  window_idx <- which(snp_chr == chrom & snp_pos >= start_pos & snp_pos <= end_pos)
  
  # Need at least 2 SNPs to form a block
  if(length(window_idx) < 2) return(0) 
  
  win_snp_ids <- snp_id[window_idx]
  win_snp_pos <- snp_pos[window_idx]
  
  # Calculate pairwise LD (r^2) matrix
  # suppressWarnings is used to keep console clean from monomorphic SNP warnings
  ld_res <- suppressWarnings(
    snpgdsLDMat(gds_obj, sample.id = target_samples, snp.id = win_snp_ids, method = "r", slide = -1, verbose = FALSE)
  )
  
  ld_mat <- ld_res$LD^2 # Square r to get r^2
  
  # Find the maximum physical distance between any two SNPs where r^2 > threshold
  max_block_size <- 0
  block_start <- NA
  block_end <- NA
  
  # Matrix is symmetric, we only need to check the upper triangle
  for (i in 1:(nrow(ld_mat)-1)) {
    for (j in (i+1):ncol(ld_mat)) {
      if (!is.na(ld_mat[i,j]) && ld_mat[i,j] >= r2_threshold) {
        dist <- win_snp_pos[j] - win_snp_pos[i]
        if (dist > max_block_size) {
          max_block_size <- dist
          block_start <- win_snp_pos[i]
          block_end <- win_snp_pos[j]
        }
      }
    }
  }
  
  return(list(size = max_block_size, start = block_start, end = block_end))
}

# --- 3. Process Regions with Specific Pairings ---
process_paired_regions <- function(regions_df, region_type, elite_samples, baseline_samples) {
  cat(sprintf("Processing %d %s regions...\n", nrow(regions_df), region_type))
  
  results <- list()
  
  for(i in 1:nrow(regions_df)) {
    r <- regions_df[i, ]
    
    # Calculate for Elite (ACT_187)
    act_block <- calc_max_ld_block(f, elite_samples, r$chr, r$start, r$stop)
    # Calculate for Specific Baseline (LDP_227 or LDP_97)
    ldp_block <- calc_max_ld_block(f, baseline_samples, r$chr, r$start, r$stop)
    
    results[[i]] <- data.frame(
      Region_ID = paste(r$chr, r$start, r$stop, sep="_"),
      Chromosome = r$chr,
      Window_Start = r$start,
      Window_Stop = r$stop,
      Region_Type = region_type,
      ACT_Block_Size_bp = act_block$size,
      LDP_Block_Size_bp = ldp_block$size
    )
  }
  return(do.call(rbind, results))
}

# Run Adaptation sweeps using ACT_187 vs LDP_227
df_adaptation_blocks <- process_paired_regions(
  regions_df = read.csv("Results/XPCLR_Scenario1_Adaptation_SignificantRegions.csv"), 
  region_type = "Adaptation Sweeps", 
  elite_samples = act_187_samples, 
  baseline_samples = ldp_227_samples
)

# Run Breeding sweeps using ACT_187 vs LDP_97
df_breeding_blocks <- process_paired_regions(
  regions_df = read.csv("Results/XPCLR_Scenario2_Breeding_SignificantRegions.csv"), 
  region_type = "Modern Breeding Sweeps", 
  elite_samples = act_187_samples, 
  baseline_samples = ldp_97_samples
)

# Combine datasets
master_blocks <- rbind(df_breeding_blocks, df_adaptation_blocks)

snpgdsClose(f)

# --- 4. Export CSV and BED Files ---
cat("Exporting data...\n")

if(!dir.exists("Results/Haplotypes")) dir.create("Results/Haplotypes", recursive = TRUE)

# Master CSV
write.csv(master_blocks, "Results/Haplotypes/Haplotype_Block_Lengths_Compared_breeding vs adapt.csv", row.names = FALSE)

# Export BED file (Chromosome, Start, End, Name) for ACTIVATE blocks
bed_act <- master_blocks %>%
  filter(!is.na(Window_Start)) %>%
  mutate(Name = paste(Region_Type, "ACT_Block", sep="_")) %>%
  select(Chromosome, Window_Start, Window_Stop, Name)
write.table(bed_act, "Results/Haplotypes/ACTIVATE_Blocks.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# --- 5. Statistical Comparison ---
cat("\n--- Statistical Summary (Wilcoxon Rank-Sum Tests) ---\n")

# Restructure data for easier plotting and stats
plot_df <- master_blocks %>%
  select(Region_ID, Region_Type, ACT = ACT_Block_Size_bp, LDP = LDP_Block_Size_bp) %>%
  pivot_longer(cols = c("ACT", "LDP"), names_to = "Panel", values_to = "Block_Size_bp") %>%
  mutate(Block_Size_Kb = Block_Size_bp / 1000)

# 5a. Compare Panels within Eroded Regions
breeding_data <- plot_df %>% filter(Region_Type == "Modern Breeding Sweeps")
test_breeding<- wilcox.test(Block_Size_Kb ~ Panel, data = breeding_data)
cat(sprintf("breeding Regions: ACT median = %.1f Kb, LDP median = %.1f Kb (p-value = %.2e)\n", 
            median(breeding_data$Block_Size_Kb[breeding_data$Panel=="ACT"], na.rm=T),
            median(breeding_data$Block_Size_Kb[breeding_data$Panel=="LDP"], na.rm=T),
            test_breeding$p.value))

# 5b. Compare Panels within Introgressed Regions
adaptation_data <- plot_df %>% filter(Region_Type == "Adaptation Sweeps")
test_adaptation <- wilcox.test(Block_Size_Kb ~ Panel, data = adaptation_data)
cat(sprintf("adaptation Regions: ACT median = %.1f Kb, LDP median = %.1f Kb (p-value = %.2e)\n", 
            median(adaptation_data$Block_Size_Kb[adaptation_data$Panel=="ACT"], na.rm=T),
            median(adaptation_data$Block_Size_Kb[adaptation_data$Panel=="LDP"], na.rm=T),
            test_intro$p.value))

# --- 6. Publication-Ready Plotting ---
cat("Generating plots...\n")

# Custom labels for the facets
facet_labels <- c(
  "Modern Breeding Sweeps" = "Improved Regions (breeding)",
  "Adaptation Sweeps" = "Adapted Regions (adaptation)"
)

p_blocks <- ggplot(plot_df, aes(x = Panel, y = Block_Size_Kb, fill = Panel)) +
  # Violin plot for distribution shape
  geom_violin(alpha = 0.5, color = NA, trim = FALSE) +
  # Boxplot for median and quartiles
  geom_boxplot(width = 0.2, color = "black", alpha = 0.8, outlier.shape = NA) +
  # Jittered points to show actual data density
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.6, aes(color = Panel)) +
  facet_wrap(~ Region_Type, labeller = as_labeller(facet_labels)) +
  scale_fill_manual(values = c("ACT" = "darkred", "LDP" = "dodgerblue")) +
  scale_color_manual(values = c("ACT" = "darkred", "LDP" = "dodgerblue")) +
  labs(
    title = "Haplotype Block Expansion Driven by Selection",
    subtitle = "Maximum contiguous linkage (r² > 0.8) within the top 1% selective sweeps and adaptations",
    x = "Lentil Panel",
    y = "Haplotype Block Length (Kb)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "gray30"),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("Results/Haplotypes/Haplotype_Block_Comparison_breeding vs adaptation.png", plot = p_blocks, width = 9, height = 6, dpi = 600)

cat("Done! Check the 'Results/Haplotypes' directory for your figures, CSVs, and BED files.\n")
