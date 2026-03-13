#SNP density map
# Load necessary libraries
library(gdsfmt)
library(ggplot2)
library(dplyr)

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