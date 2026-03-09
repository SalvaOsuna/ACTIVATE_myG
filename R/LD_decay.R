# Load necessary libraries
library(gdsfmt)
library(SeqArray)
library(SNPRelate)
library(ggplot2)
library(dplyr)

#all chromosomes####
# --- 1. Configuration ---
gds_files <- c(
  "data/ACT197_SVs.gds",
  "data/ACT197_biallelic.gds",
  "data/ACT197_biallelic_PRUNED.gds",
  "data/ACT197_biallelic_codingonly.gds"
)

dataset_labels <- c(
  "SVs",
  "All Biallelic SNPs ",
  "LD-Pruned SNPs",
  "Coding SNPs"
)

max_distance_bp <- 10000000  # Max distance to plot (10 Mb)
max_snps_per_chr <- 10000    # Sample size per chromosome to prevent memory crashes

# --- 2. Define LD Calculation Function ---
calculate_ld_decay <- function(file_path, label) {
  
  if(!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  cat(paste("\nProcessing:", label, "...\n"))
  
  # Smart File Format Detection
  f_check <- openfn.gds(file_path)
  file_format <- get.attr.gdsn(f_check$root)$FileFormat
  closefn.gds(f_check)
  
  needs_cleanup <- FALSE
  target_file <- file_path
  
  # Convert if it is a SeqArray file
  if (!is.null(file_format) && file_format == "SEQ_ARRAY") {
    cat("  Converting SeqArray to SNP GDS format...\n")
    target_file <- paste0("temp_LD_", basename(file_path))
    seqGDS2SNP(file_path, target_file, verbose = FALSE)
    needs_cleanup <- TRUE
  }
  
  # Open with SNPRelate
  genofile <- snpgdsOpen(target_file)
  
  # Filter for informative SNPs (MAF > 5%, Missing < 20%)
  cat("  Filtering for MAF > 0.05...\n")
  snpset <- snpgdsSelectSNP(genofile, maf=0.05, missing.rate=0.2, autosome.only=FALSE, verbose=FALSE)
  
  # Get coordinate metadata
  snp_id <- read.gdsn(index.gdsn(genofile, "snp.id"))
  chr <- as.character(read.gdsn(index.gdsn(genofile, "snp.chromosome")))
  pos <- read.gdsn(index.gdsn(genofile, "snp.position"))
  
  df_meta <- data.frame(snp_id = snp_id, chr = chr, pos = pos, stringsAsFactors=FALSE)
  df_meta <- df_meta[df_meta$snp_id %in% unlist(snpset), ]
  
  # Restrict to main 7 chromosomes (ignores unanchored scaffolds)
  main_chrs <- unique(df_meta$chr)[grepl("Chr[1-7]$|^[1-7]$", unique(df_meta$chr))]
  
  ld_results <- list()
  
  # Calculate LD iteratively per chromosome to save RAM
  cat("  Calculating pairwise r^2...\n")
  for(c in main_chrs) {
    chr_snps <- df_meta[df_meta$chr == c, ]
    
    # Randomly subsample if there are too many SNPs (crucial for the 17.3M dataset)
    if(nrow(chr_snps) > max_snps_per_chr) {
      chr_snps <- chr_snps[sort(sample(1:nrow(chr_snps), max_snps_per_chr)), ]
    }
    
    # Calculate LD (r) matrix
    ld_mat <- snpgdsLDMat(genofile, snp.id=chr_snps$snp_id, method="r", slide=-1, verbose=FALSE)$LD
    
    # Keep only the upper triangle of the matrix to avoid duplicates
    ld_mat[lower.tri(ld_mat, diag=TRUE)] <- NA
    
    # Calculate physical distances between all sampled pairs
    pos_mat <- abs(outer(chr_snps$pos, chr_snps$pos, "-"))
    
    # Flatten matrices into vectors
    r2_vec <- as.vector(ld_mat^2)
    dist_vec <- as.vector(pos_mat)
    
    # Filter out NAs and restrict to our max distance limit (e.g., 2 Mb)
    valid <- !is.na(r2_vec) & dist_vec > 0 & dist_vec <= max_distance_bp
    
    ld_results[[c]] <- data.frame(
      Distance = dist_vec[valid],
      R2 = r2_vec[valid]
    )
  }
  
  snpgdsClose(genofile)
  if(needs_cleanup) unlink(target_file)
  
  # Combine chromosome results
  chr_ld_df <- do.call(rbind, ld_results)
  chr_ld_df$Dataset <- label
  
  # Subsample the final datapoints for plot rendering speed (10,000 points is plenty for a smooth curve) or change it if increase max_distance_bp
  if(nrow(chr_ld_df) > 20000) {
    chr_ld_df <- chr_ld_df[sample(1:nrow(chr_ld_df), 20000), ]
  }
  
  return(chr_ld_df)
}

# --- 3. Run the Function Across All Datasets ---
all_ld_data <- mapply(calculate_ld_decay, gds_files, dataset_labels, SIMPLIFY = FALSE)

# Combine everything into one giant dataframe
plot_data <- do.call(rbind, all_ld_data)

# --- 4. Plot the LD Decay Curves ---
cat("\nGenerating LD Decay plot...\n")

# Use GAM (Generalized Additive Model) to smooth the millions of points into clear trend lines
p <- ggplot(plot_data, aes(x = Distance / 1000000, y = R2, color = Dataset)) +
  # Optional: include raw points with high transparency. 
  # Uncomment the next line if you want to see the "cloud" of points behind the lines
  # geom_point(alpha = 0.05, size = 0.5) +
  
  # Fit the decay curves
  geom_smooth(method = "gam", formula = y ~ log(x), se = FALSE, linewidth = 1.5) +
  
  # Threshold line (commonly r^2 = 0.2 or half-decay is used to define "LD block size")
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "darkgrey") +
  
  labs(
    title = "LD Decay Comparison in ACTIVATE Lentil Panel",
    x = "Physical Distance (Megabases)",
    y = expression("Linkage Disequilibrium (" * r^2 * ")"),
    color = "Dataset"
  ) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face="bold")
  )

print(p)

# Save the plot
ggsave("Results/ACTIVATE_LD_Decay_Comparison.png", plot = p, width = 9, height = 6, dpi = 600)
cat("Done! Plot saved to ACTIVATE_LD_Decay_Comparison.png\n")

#by chromosome####
# Load necessary libraries
library(gdsfmt)
library(SeqArray)
library(SNPRelate)
library(ggplot2)
library(dplyr)

# --- 1. Configuration ---
gds_files <- c(
  "data/ACT197_SVs.gds",
  "data/ACT197_biallelic.gds",
  "data/ACT197_biallelic_PRUNED.gds",
  "data/ACT197_biallelic_codingonly.gds"
)

dataset_labels <- c(
  "SVs",
  "All Biallelic SNPs",
  "LD-Pruned SNPs",
  "Coding SNPs"
)

max_distance_bp <- 10000000  # Max distance to plot (10 Mb)
max_snps_per_chr <- 10000    # Sample size per chromosome to prevent memory crashes
points_to_plot_per_chr <- 5000 # Subsample size per chromosome for rendering speed

# Ensure output directory exists
if(!dir.exists("Results")) dir.create("Results")

# --- 2. Define LD Calculation Function ---
calculate_ld_decay <- function(file_path, label) {
  
  if(!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  cat(paste("\nProcessing:", label, "...\n"))
  
  # Smart File Format Detection
  f_check <- openfn.gds(file_path)
  file_format <- get.attr.gdsn(f_check$root)$FileFormat
  closefn.gds(f_check)
  
  needs_cleanup <- FALSE
  target_file <- file_path
  
  # Convert if it is a SeqArray file
  if (!is.null(file_format) && file_format == "SEQ_ARRAY") {
    cat("  Converting SeqArray to SNP GDS format...\n")
    target_file <- paste0("temp_LD_", basename(file_path))
    seqGDS2SNP(file_path, target_file, verbose = FALSE)
    needs_cleanup <- TRUE
  }
  
  # Open with SNPRelate
  genofile <- snpgdsOpen(target_file)
  
  # Filter for informative SNPs (MAF > 5%, Missing < 20%)
  cat("  Filtering for MAF > 0.05...\n")
  snpset <- snpgdsSelectSNP(genofile, maf=0.05, missing.rate=0.2, autosome.only=FALSE, verbose=FALSE)
  
  # Get coordinate metadata
  snp_id <- read.gdsn(index.gdsn(genofile, "snp.id"))
  chr <- as.character(read.gdsn(index.gdsn(genofile, "snp.chromosome")))
  pos <- read.gdsn(index.gdsn(genofile, "snp.position"))
  
  df_meta <- data.frame(snp_id = snp_id, chr = chr, pos = pos, stringsAsFactors=FALSE)
  df_meta <- df_meta[df_meta$snp_id %in% unlist(snpset), ]
  
  # Restrict to main 7 chromosomes (ignores unanchored scaffolds)
  main_chrs <- unique(df_meta$chr)[grepl("Chr[1-7]$|^[1-7]$", unique(df_meta$chr))]
  
  ld_results <- list()
  
  # Calculate LD iteratively per chromosome to save RAM
  cat("  Calculating pairwise r^2...\n")
  for(c in main_chrs) {
    chr_snps <- df_meta[df_meta$chr == c, ]
    
    # Randomly subsample if there are too many SNPs
    if(nrow(chr_snps) > max_snps_per_chr) {
      chr_snps <- chr_snps[sort(sample(1:nrow(chr_snps), max_snps_per_chr)), ]
    }
    
    # Calculate LD (r) matrix
    ld_mat <- snpgdsLDMat(genofile, snp.id=chr_snps$snp_id, method="r", slide=-1, verbose=FALSE)$LD
    
    # Keep only the upper triangle of the matrix to avoid duplicates
    ld_mat[lower.tri(ld_mat, diag=TRUE)] <- NA
    
    # Calculate physical distances between all sampled pairs
    pos_mat <- abs(outer(chr_snps$pos, chr_snps$pos, "-"))
    
    # Flatten matrices into vectors
    r2_vec <- as.vector(ld_mat^2)
    dist_vec <- as.vector(pos_mat)
    
    # Filter out NAs and restrict to our max distance limit
    valid <- !is.na(r2_vec) & dist_vec > 0 & dist_vec <= max_distance_bp
    
    # ADDED: Include the Chromosome column so we can group by it later
    ld_results[[c]] <- data.frame(
      Chromosome = c,
      Distance = dist_vec[valid],
      R2 = r2_vec[valid],
      stringsAsFactors = FALSE
    )
  }
  
  snpgdsClose(genofile)
  if(needs_cleanup) unlink(target_file)
  
  # Combine chromosome results
  chr_ld_df <- do.call(rbind, ld_results)
  chr_ld_df$Dataset <- label
  
  # ADDED: Subsample the datapoints PER CHROMOSOME for plot rendering speed
  # Base R split-apply-combine to safely subsample without dplyr version issues
  split_df <- split(chr_ld_df, chr_ld_df$Chromosome)
  subsampled_list <- lapply(split_df, function(df) {
    if(nrow(df) > points_to_plot_per_chr) {
      df[sample(1:nrow(df), points_to_plot_per_chr), ]
    } else {
      df
    }
  })
  
  final_chr_ld_df <- do.call(rbind, subsampled_list)
  rownames(final_chr_ld_df) <- NULL
  
  return(final_chr_ld_df)
}

# --- 3. Run the Function Across All Datasets ---
all_ld_data <- mapply(calculate_ld_decay, gds_files, dataset_labels, SIMPLIFY = FALSE)

# Combine everything into one giant dataframe
plot_data <- do.call(rbind, all_ld_data)

# --- 4. Plot the LD Decay Curves per Chromosome ---
cat("\nGenerating per-chromosome LD Decay plots...\n")

# Get unique chromosomes from the dataset
unique_chrs <- unique(plot_data$Chromosome)

for (chr_name in unique_chrs) {
  cat(paste("  Plotting", chr_name, "...\n"))
  
  # Subset data for this specific chromosome
  chr_plot_data <- plot_data[plot_data$Chromosome == chr_name, ]
  
  # Create the plot
  p <- ggplot(chr_plot_data, aes(x = Distance / 1000000, y = R2, color = Dataset)) +
    # Fit the decay curves
    geom_smooth(method = "gam", formula = y ~ log(x), se = FALSE, linewidth = 1.5) +
    # Threshold line (r^2 = 0.2)
    geom_hline(yintercept = 0.2, linetype = "dashed", color = "darkgrey") +
    labs(
      title = paste("LD Decay Comparison -", chr_name),
      x = "Physical Distance (Megabases)",
      y = expression("Linkage Disequilibrium (" * r^2 * ")"),
      color = "Dataset"
    ) +
    scale_color_brewer(palette = "Set1") +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_text(face="bold")
    )
  
  # Clean up the chromosome name for the filename (e.g., handles "Lcu.1GRN.Chr1")
  clean_chr_name <- gsub("[^A-Za-z0-9_]", "_", chr_name)
  out_file <- paste0("Results/ACTIVATE_LD_Decay_", clean_chr_name, ".png")
  
  # Save the plot
  ggsave(out_file, plot = p, width = 9, height = 6, dpi = 600)
}

cat("Done! All chromosome plots saved in the Results folder.\n")