# Load necessary libraries
library(SeqArray)
library(SNPRelate)
library(ggplot2)
library(gridExtra)
library(gdsfmt)

# --- 1. Define Input Files and Titles ---
gds_files <- c(
  "data/ACT187_SVs.gds",
  "data/ACT187_biallelic.gds",
  "data/ACT187_biallelic_PRUNED.gds",
  "data/ACT187_biallelic_codingonly.gds"
)

plot_titles <- c(
  "Structural Variants (SVs)",
  "All Biallelic SNPs",
  "LD-Pruned SNPs",
  "Coding-Region SNPs"
)

# --- 2. Define Smart PCA & Plotting Function ---
generate_pca_plot <- function(file_path, title) {
  
  if(!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  cat(paste("\nProcessing:", title, "...\n"))
  
  # Step A: Smart File Format Detection using gdsfmt
    f_check <- openfn.gds(file_path)
  # Extract the FileFormat attribute from the root of the GDS file
  file_format <- get.attr.gdsn(f_check$root)$FileFormat
  closefn.gds(f_check)
  
  needs_cleanup <- FALSE
  target_file <- file_path
  
  # Step B: Convert only if necessary
  if (!is.null(file_format) && file_format == "SEQ_ARRAY") {
    cat("  Detected SeqArray format. Converting to SNP GDS format...\n")
    target_file <- paste0("temp_SNPRelate_", basename(file_path))
    seqGDS2SNP(file_path, target_file, verbose = FALSE)
    needs_cleanup <- TRUE
  } else {
    cat("  Detected SNP GDS format. No conversion needed...\n")
  }
  
  # Step C: Run PCA
  cat("  Running PCA...\n")
  # Now SNPRelate can safely open it regardless of its original format
  genofile <- snpgdsOpen(target_file)
  
  # Run PCA (autosome.only = FALSE is critical for non-human genomes)
  pca <- snpgdsPCA(genofile, autosome.only = FALSE, verbose = FALSE)
  
  # Calculate the percentage of variance explained by each PC
  pc_percent <- pca$varprop * 100
  
  # Extract PC1 and PC2 into a dataframe
  pca_df <- data.frame(
    Sample = pca$sample.id,
    PC1 = pca$eigenvect[, 1],
    PC2 = pca$eigenvect[, 2]
  )
  
  # Step D: Create the ggplot
  p <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
    geom_point(alpha = 0.7, color = "dodgerblue4", size = 2) +
    labs(
      title = title,
      x = sprintf("PC1 (%.1f%%)", pc_percent[1]),
      y = sprintf("PC2 (%.1f%%)", pc_percent[2])
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  
  # Close the GDS file
  snpgdsClose(genofile)
  
  # Step E: Clean up
  # Delete the temporary file only if we created one
  if (needs_cleanup) {
    unlink(target_file)
  }
  
  return(p)
}

# --- 3. Execute and Plot ---
# Apply the function to all files
pca_plots <- mapply(generate_pca_plot, gds_files, plot_titles, SIMPLIFY = FALSE)

# Remove any NULLs in case a file was missing
pca_plots <- Filter(Negate(is.null), pca_plots)

# Arrange the plots in a 2x2 grid
cat("\nGenerating combined plot...\n")
grid.arrange(grobs = pca_plots, ncol = 2)

# Save the combined plot to a high-resolution PNG
ggsave("Results/ACTIVATE_PCA_Comparisons_NOLR86.png", 
       arrangeGrob(grobs = pca_plots, ncol = 2), 
       width = 12, height = 10, dpi = 600)
