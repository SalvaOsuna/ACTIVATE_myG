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
                             window_size = 500000, step_size = 50000, thin_bp = 2000) {
  
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

#plotting from python####
library(ggplot2)
library(dplyr)
library(data.table)

# --- 1. Load and Merge Raw Python Outputs ---
cat("Loading raw XP-CLR files...\n")

# Define the prefix for the scenario you want to plot
# Change this to "Py_Scenario2_Breeding" or "Py_Scenario1_Adaptation" when you want to plot the second run
target_scenario <- "Py_Scenario1_Adaptation"

# Find all 7 chromosome files generated by Python (ignoring the lack of extension)
file_pattern <- paste0("^", target_scenario, "_Chr[1-7]$")
xpclr_files <- list.files("Results", pattern = file_pattern, full.names = TRUE)

if(length(xpclr_files) == 0) {
  stop("No files found! Check your Results folder path.")
}

# Read and combine all files into one dataframe. 
# fread() automatically handles tab-delimited text, even without a .txt extension.
results_list <- lapply(xpclr_files, fread)
xpclr_df <- rbindlist(results_list)

# Rename columns to standard lowercase 
colnames(xpclr_df) <- tolower(colnames(xpclr_df))

# --- 2. Clean Data and Calculate Z-Scores ---
cat("Calculating Z-scores...\n")

# The Python tool outputs the score in a column usually named 'xpclr'
# We filter out windows that failed calculation (NaN or Inf)
clean_df <- xpclr_df %>%
  filter(!is.na(xpclr), is.finite(xpclr)) %>%
  mutate(
    # Ensure xp-clr scores are treated as numeric
    xpclr = as.numeric(xpclr),
    # Calculate midpoint of the window in Megabases for plotting
    mid_pos_mb = ((start + stop) / 2) / 1000000,
    # Extract chromosome number
    chr_num = as.numeric(gsub(".*Chr", "", chrom))
  )

# Calculate genome-wide Z-scores
global_mean <- mean(clean_df$xpclr, na.rm = TRUE)
global_sd <- sd(clean_df$xpclr, na.rm = TRUE)

clean_df <- clean_df %>%
  mutate(z_score = (xpclr - global_mean) / global_sd)

# Define Significance Threshold (Top 1%)
threshold <- quantile(clean_df$z_score, 0.99, na.rm = TRUE)
cat(sprintf("Top 1%% Significance Threshold (Z-score): %.3f\n", threshold))

# --- 3. Plot the Manhattan Plot ---
cat("Generating Manhattan plot...\n")
clean_df$color_group <- ifelse(clean_df$chr_num %% 2 == 0, "Even", "Odd")

# Format title based on scenario
plot_title <- gsub("_", " ", target_scenario)

p <- ggplot(clean_df, aes(x = mid_pos_mb, y = z_score)) +
  geom_point(aes(color = color_group), alpha = 0.7, size = 0.9) +
  geom_point(data = filter(clean_df, z_score > threshold), color = "firebrick", size = 1.3) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 0.8) +
  facet_grid(. ~ chr_num, scales = "free_x", space = "free_x", switch = "x") +
  scale_color_manual(values = c("Even" = "steelblue4", "Odd" = "skyblue3")) +
  labs(
    title = paste("Python XP-CLR:", plot_title),
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

# --- 4. Export Results ---
# Save the plot
ggsave(sprintf("Results/%s_Manhattan.png", target_scenario), plot = p, width = 12, height = 6, dpi = 600)

# Save a CSV of ONLY the significant regions for gene annotation
significant_regions <- clean_df %>% filter(z_score > threshold) %>% arrange(desc(z_score))
write.csv(significant_regions, sprintf("Results/%s_SignificantRegions.csv", target_scenario), row.names = FALSE)

cat("Done! Plot and significant regions saved.\n")

#a different plot option####
library(ggplot2)
library(dplyr)
library(data.table)
library(patchwork)
library(ggrepel)

cat("Processing XP-CLR data for lollipop plots...\n")

# --- 1. Define a function to prepare the data ---
prep_lollipop_data <- function(scenario_prefix, threshold_quantile = 0.99) {
  # Load all 7 chromosome files
  files <- list.files("Results", pattern = paste0("^", scenario_prefix, "_Chr[1-7]$"), full.names = TRUE)
  if(length(files) == 0) stop(paste("No files found for", scenario_prefix))
  
  df <- rbindlist(lapply(files, fread))
  colnames(df) <- tolower(colnames(df))
  
  # Clean and format
  clean_df <- df %>%
    filter(!is.na(xpclr), is.finite(xpclr)) %>%
    mutate(
      xpclr = as.numeric(xpclr),
      chr_num = as.numeric(gsub(".*Chr", "", chrom)),
      bp_mid = (start + stop) / 2
    )
  
  # Calculate cumulative positions for a continuous X-axis
  chr_lens <- clean_df %>% group_by(chr_num) %>% summarize(chr_max = max(bp_mid))
  chr_lens$offset <- cumsum(as.numeric(chr_lens$chr_max)) - chr_lens$chr_max
  
  clean_df <- clean_df %>% 
    left_join(chr_lens, by = "chr_num") %>%
    mutate(bp_cum = bp_mid + offset)
  
  # Calculate X-axis breaks (center of each chromosome)
  axis_df <- clean_df %>% group_by(chr_num) %>% summarize(center = mean(bp_cum))
  
  # Calculate threshold based on raw XP-CLR scores
  thresh_val <- quantile(clean_df$xpclr, threshold_quantile, na.rm = TRUE)
  
  return(list(data = clean_df, axis = axis_df, threshold = thresh_val))
}

# --- 2. Load Data for Both Scenarios ---
data_adapt <- prep_lollipop_data("Py_Scenario1_Adaptation")
data_breed <- prep_lollipop_data("Py_Scenario2_Breeding") 
# --- 3. Define the Color Palette ---
# Approximating the colors from your provided image
chr_colors <- c(
  "1" = "#1b9e77", "2" = "#d95f02", "3" = "#7570b3", 
  "4" = "#e7298a", "5" = "#66a61e", "6" = "#e6ab02", "7" = "#a6761d"
)

# --- 4. Define the Plotting Function ---
make_lollipop_plot <- function(prep_data, title_label, show_x_labels = FALSE) {
  df <- prep_data$data
  ax <- prep_data$axis
  thr <- prep_data$threshold
  
  p <- ggplot(df, aes(x = bp_cum, y = xpclr, color = as.character(chr_num))) +
    # The Lollipop sticks
    geom_segment(aes(xend = bp_cum, yend = 0), alpha = 0.8, linewidth = 0.4) +
    # The Lollipop heads
    geom_point(size = 1.2, alpha = 0.9) +
    # The Threshold line
    geom_hline(yintercept = thr, color = "firebrick", linetype = "dashed", linewidth = 0.8) +
    # Add the threshold value in red on the y-axis
    annotate("text", x = min(df$bp_cum) - (max(df$bp_cum)*0.02), y = thr, 
             label = round(thr, 2), color = "firebrick", hjust = 1, vjust = -0.5, size = 3.5, fontface = "bold") +
    scale_color_manual(values = chr_colors) +
    scale_x_continuous(breaks = ax$center, labels = paste("LG", ax$chr_num, sep=""), expand = c(0.02, 0.02)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    # Clip = "off" allows the red threshold text to sit outside the plot area
    coord_cartesian(clip = "off") + 
    labs(x = NULL, y = "XP-CLR", title = title_label) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 16),
      axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 15)),
      axis.text.y = element_text(size = 11, color = "black"),
      axis.line.y = element_line(color = "black"),
      # Toggle X-axis labels based on whether it's the top or bottom panel
      axis.text.x = if(show_x_labels) element_text(size = 11, face = "bold", color = "black") else element_blank(),
      axis.ticks.x = if(show_x_labels) element_line(color = "black") else element_blank(),
      axis.line.x = if(show_x_labels) element_line(color = "black") else element_blank(),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 40) # Add left margin for the red text
    )
  
  return(p)
}

# --- 5. Generate and Combine Plots ---
cat("Building final stacked plot...\n")

# Panel A: Adaptation (Hide X-axis to keep it clean)
plot_a <- make_lollipop_plot(data_adapt, "a", show_x_labels = FALSE)

# Panel B: Breeding (Show X-axis)
plot_b <- make_lollipop_plot(data_breed, "b", show_x_labels = TRUE)

# Use patchwork to stack them cleanly
final_plot <- plot_a / plot_b

# Save the high-resolution image
ggsave("Results/XPCLR_Stacked_Lollipop.png", plot = final_plot, width = 12, height = 8, dpi = 600)

cat("Done! Check 'XPCLR_Stacked_Lollipop.png' in your Results folder.\n")

#lollipop plot with R results####
library(ggplot2)
library(dplyr)
library(data.table)
library(patchwork)

cat("Loading native R XP-CLR data for lollipop plots...\n")

# --- 1. Define a function to prepare the R CSV data ---
prep_r_lollipop_data <- function(file_path, threshold_quantile = 0.99) {
  
  if(!file.exists(file_path)) stop(paste("File not found:", file_path))
  
  # Load the single CSV file
  df <- fread(file_path)
  
  # Normalize column names to lowercase to make them easy to target
  colnames(df) <- tolower(colnames(df))
  
  # The native R output might name the position 'pos', 'position', or 'bp'
  # We dynamically grab the score and position columns
  clean_df <- df %>%
    dplyr::rename(
      bp_mid = mid_pos_mb,
      xpclr_val = xpclr,
      chr_name = chr
    ) %>%
    filter(!is.na(xpclr_val), is.finite(xpclr_val)) %>%
    mutate(
      xpclr_val = as.numeric(xpclr_val),
      # Safely extract the chromosome number (handles "Lcu.1GRN.Chr1" or just "1")
      chr_num = as.numeric(gsub(".*Chr", "", chr_name))
    )
  
  # Calculate cumulative positions for a continuous X-axis
  chr_lens <- clean_df %>% group_by(chr_num) %>% summarize(chr_max = max(bp_mid, na.rm = TRUE))
  chr_lens$offset <- cumsum(as.numeric(chr_lens$chr_max)) - chr_lens$chr_max
  
  clean_df <- clean_df %>% 
    left_join(chr_lens, by = "chr_num") %>%
    mutate(bp_cum = bp_mid + offset)
  
  # Calculate X-axis breaks (center of each chromosome)
  axis_df <- clean_df %>% group_by(chr_num) %>% summarize(center = mean(bp_cum, na.rm = TRUE))
  
  # Calculate threshold based on raw XP-CLR scores (Top 1%)
  thresh_val <- quantile(clean_df$xpclr_val, threshold_quantile, na.rm = TRUE)
  
  return(list(data = clean_df, axis = axis_df, threshold = thresh_val))
}

# --- 2. Load Data for Both Scenarios ---
# Make sure these point to the folder where your CSVs are located (e.g., "Results/" or "data/")
file_adapt <- "Results/XPCLR_Scenario1_Adaptation_AllScores.csv"
file_breed <- "Results/XPCLR_Scenario2_Breeding_AllScores.csv"

data_adapt <- prep_r_lollipop_data(file_adapt)
data_breed <- prep_r_lollipop_data(file_breed)

# --- 3. Define the Color Palette ---
chr_colors <- c(
  "1" = "#1b9e77", "2" = "#d95f02", "3" = "#7570b3", 
  "4" = "#e7298a", "5" = "#66a61e", "6" = "#e6ab02", "7" = "#a6761d"
)

# --- 4. Define the Plotting Function ---
make_lollipop_plot <- function(prep_data, title_label, show_x_labels = FALSE) {
  df <- prep_data$data
  ax <- prep_data$axis
  thr <- prep_data$threshold
  
  p <- ggplot(df, aes(x = bp_cum, y = xpclr_val, color = as.character(chr_num))) +
    geom_segment(aes(xend = bp_cum, yend = 0), alpha = 0.8, linewidth = 0.4) +
    geom_point(size = 1.2, alpha = 0.9) +
    geom_hline(yintercept = thr, color = "firebrick", linetype = "dashed", linewidth = 0.8) +
    annotate("text", x = min(df$bp_cum) - (max(df$bp_cum)*0.02), y = thr, 
             label = round(thr, 2), color = "firebrick", hjust = 1, vjust = -0.5, size = 3.5, fontface = "bold") +
    scale_color_manual(values = chr_colors) +
    scale_x_continuous(breaks = ax$center, labels = paste("LG", ax$chr_num, sep=""), expand = c(0.02, 0.02)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    coord_cartesian(clip = "off") + 
    labs(x = NULL, y = "XP-CLR", title = title_label) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 16),
      axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 15)),
      axis.text.y = element_text(size = 11, color = "black"),
      axis.line.y = element_line(color = "black"),
      axis.text.x = if(show_x_labels) element_text(size = 11, face = "bold", color = "black") else element_blank(),
      axis.ticks.x = if(show_x_labels) element_line(color = "black") else element_blank(),
      axis.line.x = if(show_x_labels) element_line(color = "black") else element_blank(),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 40)
    )
  
  return(p)
}

# --- 5. Generate and Combine Plots ---
cat("Building final stacked R plot...\n")

# Panel A: Adaptation (Native R)
plot_a <- make_lollipop_plot(data_adapt, "a", show_x_labels = FALSE)

# Panel B: Breeding (Native R)
plot_b <- make_lollipop_plot(data_breed, "b", show_x_labels = TRUE)

# Stack with patchwork
final_plot <- plot_a / plot_b

# Save the high-resolution image
ggsave("Results/R_XPCLR_Stacked_Lollipop.png", plot = final_plot, width = 12, height = 8, dpi = 600)

cat("Done! Check 'R_XPCLR_Stacked_Lollipop.png' in your working directory.\n")