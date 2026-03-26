# Script 3: Analyze Haplotype Blocks in Eroded and Introgressed Regions
# Reads the PLINK blocks, maps them to regions, and computes summary metrics.

if (!requireNamespace("GenomicRanges", quietly = TRUE))
    BiocManager::install("GenomicRanges", update = FALSE)
if (!requireNamespace("dplyr", quietly = TRUE))
    install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE))
    install.packages("tidyr")

library(GenomicRanges)
library(dplyr)
library(tidyr)

# Define paths
base_dir <- "c:/Users/Salva/OneDrive - University of Saskatchewan/UsasK/github/ACTIVATE_myG"
results_dir <- file.path(base_dir, "Results")
hap_dir <- file.path(results_dir, "Haplotype_anti")

act_blocks_file <- file.path(hap_dir, "activate_blocks.blocks.det")
ldp_blocks_file <- file.path(hap_dir, "ldp_blocks.blocks.det")

eroded_file <- file.path(results_dir, "DeltaHe_Erosion_SignificantRegions.csv")
introgressed_file <- file.path(results_dir, "Significant_Introgression_DeltaHe.csv")
breeding_file <- file.path(results_dir, "XPCLR_Scenario2_Breeding_SignificantRegions.csv")
adaptation_file <- file.path(results_dir, "XPCLR_Scenario1_Adaptation_SignificantRegions.csv")

# 1. Load data
cat("Loading region and block data...\n")
act_blocks <- read.table(act_blocks_file, header = TRUE, stringsAsFactors = FALSE)
ldp_blocks <- read.table(ldp_blocks_file, header = TRUE, stringsAsFactors = FALSE)

eroded_regions <- read.csv(eroded_file, stringsAsFactors = FALSE)
introg_regions <- read.csv(introgressed_file, stringsAsFactors = FALSE)
breeding_regions <- read.csv(breeding_file, stringsAsFactors = FALSE)
adapt_regions <- read.csv(adaptation_file, stringsAsFactors = FALSE)

# Clean up chromosome names if necessary (e.g. from PLINK)
act_blocks$CHR <- as.character(act_blocks$CHR)
ldp_blocks$CHR <- as.character(ldp_blocks$CHR)

# Make GRanges objects for regions
eroded_gr <- GRanges(
  seqnames = eroded_regions$chr,
  ranges = IRanges(start = eroded_regions$start, end = eroded_regions$stop),
  region_type = "Eroded"
)

introg_gr <- GRanges(
  seqnames = introg_regions$chr,
  ranges = IRanges(start = introg_regions$start, end = introg_regions$stop),
  region_type = "Introgressed"
)

breed_gr <- GRanges(
  seqnames = breeding_regions$chr,
  ranges = IRanges(start = breeding_regions$start, end = breeding_regions$stop),
  region_type = "Breeding"
)

adapt_gr <- GRanges(
  seqnames = adapt_regions$chr,
  ranges = IRanges(start = adapt_regions$start, end = adapt_regions$stop),
  region_type = "Adaptation"
)

regions_gr <- c(eroded_gr, introg_gr, breed_gr, adapt_gr)

# 2. Function to map blocks and label them
map_blocks <- function(blocks_df, pop_name) {
  # Create GRanges for blocks
  blocks_gr <- GRanges(
    seqnames = blocks_df$CHR,
    ranges = IRanges(start = blocks_df$BP1, end = blocks_df$BP2),
    block_id = paste0(blocks_df$CHR, "_", blocks_df$BP1, "_", blocks_df$BP2),
    kb = blocks_df$KB,
    nsnps = blocks_df$NSNPS,
    population = pop_name,
    snps_list = blocks_df$SNPS
  )
  
  # Find overlaps (blocks that fall within or overlap the target regions)
  # We use type="any" meaning any block overlapping the region is tagged
  overlaps <- findOverlaps(blocks_gr, regions_gr)
  
  # Extract the indices
  q_hits <- queryHits(overlaps)
  s_hits <- subjectHits(overlaps)
  
  # Assign region types
  blocks_df$Region_Type <- "Genome-wide"
  blocks_df$Region_Type[q_hits] <- regions_gr$region_type[s_hits]
  
  blocks_df$Population <- pop_name
  return(blocks_df)
}

cat("Mapping blocks to functional regions...\n")
act_blocks_mapped <- map_blocks(act_blocks, "ACTIVATE")
ldp_blocks_mapped <- map_blocks(ldp_blocks, "LDP")

# Combine datasets
all_blocks <- bind_rows(act_blocks_mapped, ldp_blocks_mapped)

# 3. Compute Summary Statistics
cat("Computing block summaries...\n")
summary_stats <- all_blocks %>%
  group_by(Population, Region_Type) %>%
  summarise(
    Total_Blocks = n(),
    Avg_Size_kb = mean(KB, na.rm=TRUE),
    Median_Size_kb = median(KB, na.rm=TRUE),
    Max_Size_kb = max(KB, na.rm=TRUE),
    Avg_SNPs = mean(NSNPS, na.rm=TRUE),
    Median_SNPs = median(NSNPS, na.rm=TRUE),
    .groups = "drop"
  )

print(summary_stats)

# Save the mapped block sets
out_csv <- file.path(hap_dir, "All_Mapped_Haplotype_Blocks.csv")
write.csv(all_blocks, out_csv, row.names = FALSE)

out_summary <- file.path(hap_dir, "Haplotype_Blocks_Summary.csv")
write.csv(summary_stats, out_summary, row.names = FALSE)

cat(sprintf("Saved detailed mapping to %s\n", out_csv))
cat(sprintf("Saved summary statistics to %s\n", out_summary))
