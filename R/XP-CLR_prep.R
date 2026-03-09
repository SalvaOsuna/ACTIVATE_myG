#Prepare the Files in R for XP-CLR
library(gdsfmt)
library(SeqArray)

# Define file paths
gds_file <- "data/Merged_Analysis_RealCoords.gds"
temp_seq_gds <- "data/Temp_Merged_SeqArray.gds"
out_vcf <- "data/Merged_Analysis.vcf.gz"

# 1. Extract Sample IDs
f <- openfn.gds(gds_file)
sample_ids <- read.gdsn(index.gdsn(f, "sample.id"))
closefn.gds(f)

# Identify the two populations (Adjust the grep patterns if your IDs differ)
act_samples <- sample_ids[grepl("_ACT$", sample_ids)]
ldp_samples <- sample_ids[grepl("_AGL$", sample_ids)]

# Save to text files (one ID per line, no headers)
writeLines(act_samples, "data/ACT_samples.txt")
writeLines(ldp_samples, "data/LDP_samples.txt")

cat(paste("Saved", length(act_samples), "ACT samples and", length(ldp_samples), "LDP samples.\n"))

# --- 2. Convert and Export to VCF ---
cat("\nConverting SNPRelate GDS to SeqArray GDS...\n")
# Convert the SNP GDS file to a temporary SeqArray GDS file
seqSNP2GDS(gds_file, temp_seq_gds, verbose = FALSE)

# Open the newly converted file
genofile <- seqOpen(temp_seq_gds)

cat("Exporting merged dataset to VCF...\n")
# Export to a compressed VCF
seqGDS2VCF(genofile, out_vcf, use_Rsamtools = TRUE, verbose = TRUE)

# Close the file
seqClose(genofile)

# Clean up the temporary file
unlink(temp_seq_gds)
cat("VCF Export Complete! Temporary files cleaned up.\n")


#Formulas in R####
# Load necessary libraries
library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(dplyr)
library(data.table)

# --- 1. Configuration ---
gds_file <- "data/Merged_Analysis_RealCoords.gds"
window_size <- 50000  # 50 kb windows
step_size <- 10000    # 10 kb slide
thin_bp <- 2000       # Keep max 1 SNP per 2kb to correct for LD inflation (Eq. 7)

cat("Opening GDS file and extracting allele frequencies...\n")
f <- snpgdsOpen(gds_file)

# Identify samples
samp_ids <- read.gdsn(index.gdsn(f, "sample.id"))
act_samples <- samp_ids[grepl("_ACT", samp_ids)]
ldp_samples <- setdiff(samp_ids, act_samples) # LDP is everything else

# --- 2. Calculate Allele Frequencies natively (Fast) ---
cat("Calculating allele frequencies for ACT and LDP populations...\n")
freq_act <- snpgdsSNPRateFreq(f, sample.id=act_samples, with.id=TRUE)
freq_ldp <- snpgdsSNPRateFreq(f, sample.id=ldp_samples, with.id=TRUE)

df <- data.frame(
  snp.id = read.gdsn(index.gdsn(f, "snp.id")),
  chr = read.gdsn(index.gdsn(f, "snp.chromosome")),
  pos = read.gdsn(index.gdsn(f, "snp.position")),
  p1 = freq_act$AlleleFreq,  # Object Population (ACT)
  p2 = freq_ldp$AlleleFreq,  # Reference Population (LDP)
  stringsAsFactors = FALSE
)
snpgdsClose(f)

# --- 3. Statistical Prep (Chen et al. 2010) ---
cat("Applying Gaussian variance model (Eq. 1)...\n")
# Filter out SNPs that are fixed in the reference population to avoid division by zero
df <- df %>% filter(p2 > 0.02 & p2 < 0.98 & !is.na(p1) & !is.na(p2))

# Calculate expected variance parameter under neutrality
df$var_expected <- df$p2 * (1 - df$p2)

# Estimate global drift parameter (w) across the genome
w_global <- mean((df$p1 - df$p2)^2 / df$var_expected, na.rm=TRUE)

# Standardize standard deviation (sigma) for the Gaussian model
df$sigma <- sqrt(w_global * df$var_expected)
df$sigma[df$sigma < 1e-4] <- 1e-4 # Prevent infinite likelihoods

# Direction of sweep: if p1 > p2, allele is sweeping toward 1. Else toward 0.
df$target <- ifelse(df$p1 > df$p2, 1, 0)

# Restrict to main chromosomes
main_chrs <- unique(df$chr)[grepl("Chr[1-7]$|^[1-7]$", unique(df$chr))]

# --- 4. Sliding Window XP-CLR Optimization ---
cat("Running XP-CLR sliding window likelihood optimization...\n")
all_results <- list()

for(c in main_chrs) {
  cat(paste("  Processing", c, "...\n"))
  df_c <- df %>% filter(chr == c)
  
  starts <- seq(1, max(df_c$pos), by=step_size)
  chr_results <- list()
  
  for(s in starts) {
    e <- s + window_size
    win_snps <- df_c %>% filter(pos >= s & pos <= e)
    
    # Require at least 5 SNPs in a window to compute a valid likelihood
    if(nrow(win_snps) < 5) next
    
    # Eq 7: Down-weight correlated SNPs by spatial thinning
    win_snps$bin <- floor(win_snps$pos / thin_bp)
    win_snps <- win_snps %>% distinct(bin, .keep_all=TRUE)
    
    # Null Model (Eq. 1): No selection shift
    ll_null <- sum(dnorm(win_snps$p1, mean = win_snps$p2, sd = win_snps$sigma, log = TRUE))
    
    # Alternative Model (Eq. 4): Shift frequencies toward fixation target
    # We optimize the local selection shift coefficient 'c' (0 = no shift, 1 = total sweep)
    alt_ll_func <- function(c_shift) {
      mu_alt <- win_snps$p2 + c_shift * (win_snps$target - win_snps$p2)
      sum(dnorm(win_snps$p1, mean = mu_alt, sd = win_snps$sigma, log = TRUE))
    }
    
    # Optimize bounded between 0 and 1
    opt <- optimize(alt_ll_func, interval = c(0, 1), maximum = TRUE)
    ll_alt <- opt$objective
    
    # Likelihood Ratio Test Statistic (Eq. 6)
    score <- 2 * (ll_alt - ll_null)
    
    chr_results[[length(chr_results)+1]] <- data.frame(
      chr = c,
      start = s,
      stop = e,
      mid_pos_mb = (s + window_size/2) / 1e6,
      n_snps = nrow(win_snps),
      xpclr = max(0, score) # Score cannot be theoretically negative
    )
  }
  all_results[[c]] <- do.call(rbind, chr_results)
}

clean_df <- do.call(rbind, all_results)

# --- 5. Normalize to Z-scores ---
cat("Calculating Z-scores...\n")
global_mean <- mean(clean_df$xpclr, na.rm=TRUE)
global_sd <- sd(clean_df$xpclr, na.rm=TRUE)

clean_df <- clean_df %>%
  mutate(
    z_score = (xpclr - global_mean) / global_sd,
    chr_num = as.numeric(gsub(".*Chr", "", chr))
  )

threshold <- quantile(clean_df$z_score, 0.99, na.rm=TRUE)

# --- 6. Plotting the Manhattan-like Plot ---
cat("Generating Manhattan Plot...\n")
clean_df$color_group <- ifelse(clean_df$chr_num %% 2 == 0, "Even", "Odd")

p <- ggplot(clean_df, aes(x = mid_pos_mb, y = z_score)) +
  geom_point(aes(color = color_group), alpha = 0.7, size = 0.9) +
  geom_point(data = filter(clean_df, z_score > threshold), color = "firebrick", size = 1.3) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 0.8) +
  facet_grid(. ~ chr_num, scales = "free_x", space = "free_x", switch = "x") +
  scale_color_manual(values = c("Even" = "steelblue4", "Odd" = "skyblue3")) +
  labs(
    title = "XP-CLR Test of Selection (ACT vs LDP)",
    subtitle = paste("Genome-wide Z-score threshold (Top 1%):", round(threshold, 2)),
    x = "Physical Distance (Mb)",
    y = "XP-CLR (Z-score)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text.x = element_text(face = "bold", size = 10)
  )

if(!dir.exists("Results")) dir.create("Results")
ggsave("Results/ACTIVATE_Native_XPCLR.png", plot = p, width = 12, height = 6, dpi = 600)

significant_regions <- clean_df %>% filter(z_score > threshold) %>% arrange(desc(z_score))
write.csv(significant_regions, "Results/Significant_XPCLR_Regions.csv", row.names = FALSE)

cat("Done! Plot and data saved in Results/.\n")

#scenarios####
# Load necessary libraries
library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(dplyr)

# --- 1. Define the XP-CLR Wrapper Function ---
run_custom_xpclr <- function(gds_file, ref_file, obj_file, out_prefix, 
                             window_size = 50000, step_size = 10000, thin_bp = 2000) {
  
  cat(sprintf("\n======================================================\n"))
  cat(sprintf("Running XP-CLR: %s\n", out_prefix))
  cat(sprintf("======================================================\n"))
  
  # Read sample lists
  ref_samples <- readLines(ref_file)
  obj_samples <- readLines(obj_file)
  
  cat(paste("Loaded", length(ref_samples), "Reference samples and", length(obj_samples), "Objective samples.\n"))
  
  # Open GDS
  f <- snpgdsOpen(gds_file)
  
  # --- 2. Calculate Allele Frequencies ---
  cat("Calculating allele frequencies...\n")
  freq_obj <- snpgdsSNPRateFreq(f, sample.id=obj_samples, with.id=TRUE)
  freq_ref <- snpgdsSNPRateFreq(f, sample.id=ref_samples, with.id=TRUE)
  
  df <- data.frame(
    snp.id = read.gdsn(index.gdsn(f, "snp.id")),
    chr = read.gdsn(index.gdsn(f, "snp.chromosome")),
    pos = read.gdsn(index.gdsn(f, "snp.position")),
    p1 = freq_obj$AlleleFreq,  # Objective Population
    p2 = freq_ref$AlleleFreq,  # Reference Population
    stringsAsFactors = FALSE
  )
  snpgdsClose(f)
  
  # --- 3. Statistical Preparation (Chen et al. 2010) ---
  cat("Applying Gaussian variance model...\n")
  # Filter out SNPs that are fixed in the reference population to avoid division by zero
  df <- df %>% filter(p2 > 0.02 & p2 < 0.98 & !is.na(p1) & !is.na(p2))
  
  # Calculate expected variance parameter under neutrality
  df$var_expected <- df$p2 * (1 - df$p2)
  
  # Estimate global drift parameter (w)
  w_global <- mean((df$p1 - df$p2)^2 / df$var_expected, na.rm=TRUE)
  
  # Standardize standard deviation (sigma) for the Gaussian model
  df$sigma <- sqrt(w_global * df$var_expected)
  df$sigma[df$sigma < 1e-4] <- 1e-4 # Prevent infinite likelihoods
  
  # Direction of sweep: if p1 > p2, allele is sweeping toward 1. Else toward 0.
  df$target <- ifelse(df$p1 > df$p2, 1, 0)
  
  # Restrict to main chromosomes
  main_chrs <- unique(df$chr)[grepl("Chr[1-7]$|^[1-7]$", unique(df$chr))]
  
  # --- 4. Sliding Window Likelihood Optimization ---
  cat("Running sliding window likelihood optimization...\n")
  all_results <- list()
  
  for(c in main_chrs) {
    df_c <- df %>% filter(chr == c)
    starts <- seq(1, max(df_c$pos), by=step_size)
    chr_results <- list()
    
    for(s in starts) {
      e <- s + window_size
      win_snps <- df_c %>% filter(pos >= s & pos <= e)
      
      # Require at least 5 SNPs
      if(nrow(win_snps) < 5) next
      
      # Down-weight correlated SNPs by spatial thinning
      win_snps$bin <- floor(win_snps$pos / thin_bp)
      win_snps <- win_snps %>% distinct(bin, .keep_all=TRUE)
      
      # Null Model: No selection shift
      ll_null <- sum(dnorm(win_snps$p1, mean = win_snps$p2, sd = win_snps$sigma, log = TRUE))
      
      # Alternative Model: Shift frequencies toward fixation target
      alt_ll_func <- function(c_shift) {
        mu_alt <- win_snps$p2 + c_shift * (win_snps$target - win_snps$p2)
        sum(dnorm(win_snps$p1, mean = mu_alt, sd = win_snps$sigma, log = TRUE))
      }
      
      # Optimize bounded between 0 and 1
      opt <- optimize(alt_ll_func, interval = c(0, 1), maximum = TRUE)
      
      score <- 2 * (opt$objective - ll_null)
      
      chr_results[[length(chr_results)+1]] <- data.frame(
        chr = c,
        start = s,
        stop = e,
        mid_pos_mb = (s + window_size/2) / 1e6,
        n_snps = nrow(win_snps),
        xpclr = max(0, score) 
      )
    }
    all_results[[c]] <- do.call(rbind, chr_results)
  }
  
  clean_df <- do.call(rbind, all_results)
  
  # --- 5. Normalize to Z-scores ---
  cat("Calculating Z-scores...\n")
  global_mean <- mean(clean_df$xpclr, na.rm=TRUE)
  global_sd <- sd(clean_df$xpclr, na.rm=TRUE)
  
  clean_df <- clean_df %>%
    mutate(
      z_score = (xpclr - global_mean) / global_sd,
      chr_num = as.numeric(gsub(".*Chr", "", chr))
    )
  
  threshold <- quantile(clean_df$z_score, 0.99, na.rm=TRUE)
  
  # --- 6. Plot and Save ---
  cat("Exporting data and generating plot...\n")
  clean_df$color_group <- ifelse(clean_df$chr_num %% 2 == 0, "Even", "Odd")
  
  plot_title <- gsub("_", " ", out_prefix)
  
  p <- ggplot(clean_df, aes(x = mid_pos_mb, y = z_score)) +
    geom_point(aes(color = color_group), alpha = 0.7, size = 0.9) +
    geom_point(data = filter(clean_df, z_score > threshold), color = "firebrick", size = 1.3) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 0.8) +
    facet_grid(. ~ chr_num, scales = "free_x", space = "free_x", switch = "x") +
    scale_color_manual(values = c("Even" = "steelblue4", "Odd" = "skyblue3")) +
    labs(
      title = paste("XP-CLR Test:", plot_title),
      subtitle = paste("Genome-wide Z-score threshold (Top 1%):", round(threshold, 2)),
      x = "Physical Distance (Mb)",
      y = "XP-CLR (Z-score)"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text.x = element_text(face = "bold", size = 10)
    )
  
  if(!dir.exists("Results")) dir.create("Results")
  
  # Save Plot
  ggsave(sprintf("Results/XPCLR_%s_Manhattan.png", out_prefix), plot = p, width = 12, height = 6, dpi = 600)
  
  # Save Significant Regions
  sig_regions <- clean_df %>% filter(z_score > threshold) %>% arrange(desc(z_score))
  write.csv(sig_regions, sprintf("Results/XPCLR_%s_SignificantRegions.csv", out_prefix), row.names = FALSE)
  
  # Save Raw Scores (optional, but good for archiving)
  write.csv(clean_df, sprintf("Results/XPCLR_%s_AllScores.csv", out_prefix), row.names = FALSE)
  
  cat(sprintf("Success! Results saved in the 'Results' folder with prefix '%s'.\n", out_prefix))
}

# --- EXECUTE SCENARIOS ---
# Define the path to your merged GDS file
gds_dataset <- "data/Merged_Analysis_RealCoords.gds"

# SCENARIO 1: The Adaptation Scan (Global Diversity -> Canadian Adaptation)
run_custom_xpclr(
  gds_file   = gds_dataset,
  ref_file   = "data/LDP227_Non_Temperate.txt",
  obj_file   = "data/ACT187_samples.txt",
  out_prefix = "Scenario1_Adaptation"
)

# SCENARIO 2: The Breeding/Improvement Scan (Temperate Baseline -> Elite Canadian Breeding)
run_custom_xpclr(
  gds_file   = gds_dataset,
  ref_file   = "data/LDP97_Temperate.txt",
  obj_file   = "data/ACT187_samples.txt",
  out_prefix = "Scenario2_Breeding"
)
