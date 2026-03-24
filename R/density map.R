#SNP density map
# Load necessary libraries
library(gdsfmt)
library(ggplot2)
library(dplyr)

#single panel####

cat("Extracting SNP positions from GDS file...\n")

# --- 1. Extract Data ---
gds_file <- "data/Merged_Analysis_RealCoords.gds"
f <- snpgdsOpen(gds_file)

snp_df <- data.frame(
  chr = read.gdsn(index.gdsn(f, "snp.chromosome")),
  pos = read.gdsn(index.gdsn(f, "snp.position")),
  stringsAsFactors = FALSE
)
snpgdsClose(f)

# --- 2. Clean and Bin the Data ---
cat("Binning genome into 1 Mb windows...\n")

window_size <- 1000000 # 1 Mb bins

snp_density <- snp_df %>%
  # Keep only the main 7 lentil linkage groups
  filter(grepl("Chr[1-7]$|^[1-7]$", chr)) %>%
  mutate(
    chr_num = as.numeric(gsub(".*Chr", "", chr)),
    # Calculate the starting point of the bin for each SNP
    bin_start = floor(pos / window_size) * window_size,
    # Convert to Megabases for a cleaner X-axis
    bin_mb = bin_start / 1e6 
  ) %>%
  # Count the number of SNPs in each 1 Mb bin
  group_by(chr_num, bin_mb) %>%
  summarise(snp_count = n(), .groups = "drop")

# --- 3. Plotting the Density Map ---
cat("Generating SNP Density Plot...\n")

# Reverse the chromosome order so LG1 is at the top of the plot
snp_density$chr_factor <- factor(paste("Chr", snp_density$chr_num, sep=""), 
                                 levels = paste("Chr", 7:1, sep=""))

p_density <- ggplot(snp_density, aes(x = bin_mb, y = chr_factor, fill = snp_count)) +
  # Use tiles to create the chromosomal heatmap blocks
  geom_tile(width = 1, height = 0.6, color = NA) +
  
  # A high-contrast color palette (Plasma goes from dark purple to bright yellow)
  # You can change "plasma" to "viridis" or "magma" if you prefer
  scale_fill_viridis_c(option = "C", name = "SNPs per 1 Mb\n") +
  
  labs(
    #title = "Genome-Wide Marker Distribution",
    #subtitle = "SNP density across 1 Mb sliding windows for the merged dataset",
    x = "Physical Distance (Mb)",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "gray20", margin = margin(b = 15)),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
    legend.title = element_text(face = "bold", size = 11),
    legend.position = "right",
    # Remove the gridlines for a cleaner, floating chromosome look
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

# Save the plot
if(!dir.exists("Results")) dir.create("Results")
ggsave("Results/Merged_SNP_Density_Map.png", plot = p_density, width = 9, height = 5, dpi = 400)

cat("Done! SNP density map saved to the Results folder.\n")


#multi panel####
# Load necessary libraries
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales) # For comma formatting in legends

cat("Starting Multi-Panel SNP Density Extraction...\n")

# --- 1. Configuration ---
gds_files <- c(
  "data/LDP324_nofiltered.gds",
  "data/ACT197_biallelic.gds",
  "data/Merged_Analysis_RealCoords.gds"
)

dataset_labels <- c(
  "LDP324_SNPs",
  "ACT197_SNPs",
  "LDP324+ACT197"
)

window_size <- 1000000 # 1 Mb bins

# --- 2. Define Data Extraction Function ---
get_snp_density <- function(file_path, label) {
  
  if(!file.exists(file_path)) stop(paste("File not found:", file_path))
  cat(paste("  Processing:", label, "...\n"))
  
  # Smart File Format Detection (crucial for the ACT197 dataset)
  f_check <- openfn.gds(file_path)
  file_format <- get.attr.gdsn(f_check$root)$FileFormat
  closefn.gds(f_check)
  
  needs_cleanup <- FALSE
  target_file <- file_path
  
  if (!is.null(file_format) && file_format == "SEQ_ARRAY") {
    cat("    Converting SeqArray to SNP GDS format...\n")
    target_file <- paste0("temp_Density_", basename(file_path))
    seqGDS2SNP(file_path, target_file, verbose = FALSE)
    needs_cleanup <- TRUE
  }
  
  # Extract Data
  f <- snpgdsOpen(target_file)
  df <- data.frame(
    chr = read.gdsn(index.gdsn(f, "snp.chromosome")),
    pos = read.gdsn(index.gdsn(f, "snp.position")),
    stringsAsFactors = FALSE
  )
  snpgdsClose(f)
  if(needs_cleanup) unlink(target_file)
  
  # Clean and Bin
  cat("    Binning genome into 1 Mb windows...\n")
  density_df <- df %>%
    filter(grepl("Chr[1-7]$|^[1-7]$", chr)) %>%
    mutate(
      chr_num = as.numeric(gsub(".*Chr", "", chr)),
      bin_start = floor(pos / window_size) * window_size,
      bin_mb = bin_start / 1e6 
    ) %>%
    group_by(chr_num, bin_mb) %>%
    summarise(snp_count = n(), .groups = "drop") %>%
    mutate(Dataset = label)
  
  return(density_df)
}

# --- 3. Run Extraction Across All Datasets ---
all_density_data <- mapply(get_snp_density, gds_files, dataset_labels, SIMPLIFY = FALSE)

# ADD THIS LINE: Explicitly name the list elements so the loop can find them
names(all_density_data) <- dataset_labels

# --- 4. Plotting the Density Maps ---
cat("\nGenerating publication-ready faceted density plots...\n")

plots <- list()
plot_index <- 1

for (ds in dataset_labels) {
  
  df_sub <- all_density_data[[ds]]
  
  # Standardize chromosome factor (Chr1 to Chr7 from left to right)
  df_sub$chr_factor <- factor(paste("Chr", df_sub$chr_num, sep=""), 
                              levels = paste("Chr", 1:7, sep=""))
  
  # Base Plotting: X = Chromosome, Y = Distance
  p <- ggplot(df_sub, aes(x = chr_factor, y = bin_mb, fill = snp_count)) +
    geom_tile(width = 0.7, height = 1, color = NA) +
    
    # Use comma formatting for the legend
    scale_fill_viridis_c(option = "C", name = "SNPs/Mb", labels = comma) +
    
    # Reverse Y axis so 0 Mb is at the top (standard biological formatting)
    scale_y_reverse(expand = c(0, 0)) +
    
    # Apply standard size 10 font universally
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = NULL, #element_text(size = 10, face = "bold", color = "black", angle = 45, hjust = 1),
      legend.title = element_text(face = "bold", size = 11),
      
      # Move legend to the bottom so the 3 plots aren't horizontally crushed
      legend.position = "bottom",
      legend.key.width = unit(0.9, "cm"),
      legend.key.height = unit(0.2, "cm"),
      
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    ) +
    
    # Format the legend to sit nicely below the plot
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5,angle = 45))
  
  # Handle Y-axis for side-by-side layout (only panel 'a' gets the Y-axis labels)
  if (plot_index == 1) {
    p <- p + labs(title = NULL, x = NULL, y = "Physical Distance (Mb)") +
      theme(
        axis.title.y = element_text(size = 10, face = "bold", margin = margin(r = 10)),
        axis.text.y = element_text(size = 10, color = "black")
      )
  } else {
    p <- p + labs(title = NULL, x = NULL, y = NULL) +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }
  
  plots[[ds]] <- p
  plot_index <- plot_index + 1
}

# Arrange plots side-by-side (1 row, 3 columns) and add a), b), c)
final_plot <- wrap_plots(plots, ncol = 3) + 
  plot_annotation(tag_levels = 'a', tag_suffix = ')') & 
  theme(plot.tag = element_text(face = 'bold', size = 10))

print(final_plot)

# Save explicitly mapped to an A4 layout (17 cm wide)
if(!dir.exists("Results")) dir.create("Results")
ggsave("Results/Merged_SNP_Density_Map_Publication_Vertical.png", 
       plot = final_plot, 
       width = 18, 
       height = 10,  # Adjusted height to look proportional with vertical chromosomes
       units = "cm", 
       dpi = 600)

cat("Done! Clean publication plot saved to Results/Merged_SNP_Density_Map_Publication_Vertical.png\n")