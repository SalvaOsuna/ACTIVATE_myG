# Load necessary libraries
library(SNPRelate)
library(ggplot2)
library(gridExtra)
library(SeqArray)


# --- 1. Define Input Files and Titles ---
gds_files <- c(
  "data/ACT197_SVs.gds",
  "data/ACT197_biallelic.gds",
  "data/ACT197_biallelic_PRUNED.gds",
  "data/ACT197_biallelic_codingonly.gds"
)

# Titles for the plots
plot_titles <- c(
  "Structural Variants (SVs)",
  "All Biallelic SNPs",
  "LD-Pruned SNPs",
  "Coding-Region SNPs"
)

# --- 2. Define PCA & Plotting Function ---
generate_pca_plot <- function(file_path, title) {
  
  if(!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  cat(paste("Running PCA for:", title, "...\n"))
  
  # Open the GDS file
  genofile <- snpgdsOpen(file_path)
  
  # Run PCA
  # autosome.only = FALSE is critical for non-human genomes like lentil
  pca <- snpgdsPCA(genofile, autosome.only = FALSE, verbose = FALSE)
  
  # Calculate the percentage of variance explained by each PC
  pc_percent <- pca$varprop * 100
  
  # Extract PC1 and PC2 into a dataframe
  pca_df <- data.frame(
    Sample = pca$sample.id,
    PC1 = pca$eigenvect[, 1],
    PC2 = pca$eigenvect[, 2]
  )
  
  # Create the ggplot
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
  
  return(p)
}

# --- 3. Execute and Plot ---
# Apply the function to all files
pca_plots <- mapply(generate_pca_plot, gds_files, plot_titles, SIMPLIFY = FALSE)

# Remove any NULLs in case a file was missing
pca_plots <- Filter(Negate(is.null), pca_plots)

# Arrange the plots in a 2x2 grid
cat("Generating combined plot...\n")
grid.arrange(grobs = pca_plots, ncol = 2)

# Optional: Save the combined plot to a high-resolution PDF or PNG
ggsave("ACTIVATE_PCA_Comparisons.png", 
       arrangeGrob(grobs = pca_plots, ncol = 2), 
       width = 12, height = 10, dpi = 300)