#required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("SeqArray", "SNPRelate", "vcfR", "ggplot2"))

#VCF to GDS####
library(SeqArray)
library(SNPRelate)

# Define your file names
 #setwd("data/")
vcf_fn <- "data/Lcu.1GRN-ACTIVATE-lightly filtered.vcf.gz"
gds_fn <- "data/Lcu_ACTIVATE_biallelic.gds"

# 1. Convert VCF to GDS
# parallel = 4 uses 4 cores to speed this up. Adjust based on your PC.

seqVCF2GDS(vcf_fn, 
           gds_fn, 
           method = "biallelic.only")

# 2. Open the connection to the new GDS file
genofile <- seqOpen(gds_fn)

# 3. Quick Summary of what's inside
# This will tell you exactly how many samples and variants R can see
seqSummary(genofile)

# The file is now "open" and ready for analysis. 
# Remember to close it when done: seqClose(genofile)

#quick PCA####
# Run PCA using SNPRelate (very fast on GDS files)
# We use 'autosome.only=FALSE' since plant chromosomes might not be numbered standardly
pca <- snpgdsPCA(genofile, autosome.only = FALSE, num.thread = 4)

# Calculate the % variance explained by PC1 and PC2
pc.percent <- pca$varprop * 100
head(round(pc.percent, 2))

# Make a dataframe for plotting
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

# Plot with ggplot2
library(ggplot2)
ggplot(tab, aes(x=EV1, y=EV2)) +
  geom_point(color="blue", alpha=0.6) +
  labs(title="PCA of 197 Lentil Breeding Lines",
       x=paste0("PC1 (", round(pc.percent[1], 2), "%)"),
       y=paste0("PC2 (", round(pc.percent[2], 2), "%)")) +
  theme_minimal()

#heterozygosity####
# Calculate heterozygosity for each individual
het_res <- snpgdsIndivBeta(genofile, autosome.only = FALSE, num.thread = 4, inbreeding = T)

# Create a data frame
het_df <- data.frame(sample.id = het_res$sample.id, 
                     inbreeding_coef = diag(het_res$beta),  # Extract the diagonal (This is the Inbreeding Coef for each individual)
                     heterozygosity = 1 - het_res$beta) # Approximation of Het

# Look at the distribution
summary(het_df$inbreeding_coef)

# Histogram of Inbreeding
hist(het_df$inbreeding_coef, 
     main="Distribution of Inbreeding Coefficients", 
     xlab="Inbreeding Coefficient (Beta)", col="lightblue")

#Load LDP####
# Convert the Diversity Panel file too
# Define file names
setwd("C:/Users/tlv329/OneDrive - University of Saskatchewan/UsasK/ACTIVATE_project/data analysis/VariantSet/lentil/LDP/")
ldp_vcf <- "LDP-Lcu.1GRN.vcf.gz"
ldp_gds <- "LDP_Diversity_Panel.gds"

# Convert VCF to GDS (SNPRelate Format)
snpgdsVCF2GDS(vcf.fn = ldp_vcf, 
              out.fn = ldp_gds, 
              method = "biallelic.only", # Usually recommended to ignore complex multi-allelic sites
              verbose = TRUE)

# Verify by opening it immediately
genofile <- snpgdsOpen(ldp_gds)
snpgdsSummary(genofile)
snpgdsClose(genofile)

# I now have two .gds files:
# 1. Lcu_Breeding_Lines.gds
# 2. LDP_Diversity_Panel.gds

#Filter GDS using GFF3####
 #this is to keep just the variants within coding regions in the ACT197 panel
library(SNPRelate)
library(GenomicRanges)
library(rtracklayer)
library(gdsfmt) # Needed for index.gdsn

# Files
gds_fn <- "ACT197_biallelic.gds"
gff_fn <- "Lcu.1GRN.genes_description.sorted.gff3.gz"
out_fn <- "ACT197_codingonly.gds" # New filename

# 1. Load GFF3 Annotation
cat("Reading GFF3 annotation...\n")
gff_data <- import.gff3(gff_fn)
# Filter for genes (adjust 'type' if your GFF uses 'CDS' or 'mRNA')
gene_ranges <- gff_data[gff_data$type == "gene"]


# 2. Open GDS in SNPRelate Mode (As requested)
cat("Opening GDS file...\n")
genofile <- snpgdsOpen(gds_fn, allow.duplicate=TRUE)


# 3. Extract Chromosome, Position, and IDs
# Note: We need to handle potential node name differences if the file was made with SeqArray.
# This block attempts to find the correct nodes automatically.

# Try standard SNPRelate names first
if (!is.null(index.gdsn(genofile, "snp.id", silent=TRUE))) {
  all_ids <- read.gdsn(index.gdsn(genofile, "snp.id"))
  all_chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
  all_pos <- read.gdsn(index.gdsn(genofile, "snp.position"))
} else {
  # Fallback for SeqArray files opened in SNPRelate
  cat("Note: Using SeqArray node names (variant.id)...\n")
  all_ids <- read.gdsn(index.gdsn(genofile, "variant.id"))
  all_chr <- read.gdsn(index.gdsn(genofile, "chromosome"))
  all_pos <- read.gdsn(index.gdsn(genofile, "position"))
}

# 4. Find Overlaps
# Create GRanges for the SNPs
cat("Mapping SNPs to Genomic Ranges...\n")
snp_gr <- GRanges(seqnames = all_chr, 
                  ranges = IRanges(start = all_pos, width = 1))

cat("Finding overlaps with genes...\n")
overlaps <- findOverlaps(snp_gr, gene_ranges)

# Get the INDICES (1, 2, 5...) of the SNPs that overlap
valid_indices <- unique(queryHits(overlaps))

# Convert Indices to actual IDs
# We need the actual ID numbers to pass to snpgdsCreateGenoSet
ids_to_keep <- all_ids[valid_indices]

cat(paste("Total SNPs:", length(all_ids), "\n"))
cat(paste("Coding SNPs:", length(ids_to_keep), "\n"))


# 5. Close the File
# We MUST close it because snpgdsCreateGenoSet needs to read the file from disk fresh.
snpgdsClose(genofile)


# 6. Create the New Subsetted GDS File
if(length(ids_to_keep) > 0) {
  cat("Creating new GDS file...\n")
  snpgdsCreateGenoSet(src.fn = gds_fn, 
                      dest.fn = out_fn, 
                      snp.id = ids_to_keep, 
                      verbose = TRUE)
  cat(paste("Success! Saved as:", out_fn))
} else {
  cat("Error: No SNPs found in coding regions. Check your GFF3 chromosome names match the VCF.")
}

#GDS to VCF####
library(SeqArray)
library(SNPRelate)

# 1. Define filenames
input_gds  <- "ACT197_codingonly.gds"
temp_gds   <- "ACT197_temp.gds"   # Intermediate SeqArray format
output_vcf <- "ACT197_codingonly.vcf"

# 2. Convert SNP GDS (SNPRelate format) to Seq GDS (SeqArray format)
# This step is necessary because seqGDS2VCF only accepts SeqArray files
seqSNP2GDS(input_gds, temp_gds)

# 3. Open the new SeqArray GDS file
genofile <- seqOpen(temp_gds)

# 4. Export to VCF
# 'info.var=NULL' and 'fmt.var=NULL' exports all available fields
# Use 'gzip=TRUE' if you want a .vcf.gz file directly (saves space)
seqGDS2VCF(genofile, output_vcf, info.var=NULL, fmt.var=NULL, verbose=TRUE)

# 5. Close the file and clean up
seqClose(genofile)

# Optional: Remove the intermediate file if you don't need it
file.remove(temp_gds)

print("Conversion complete!")

#LD decay####
  # 1. SETUP & LIBRARIES
library(SeqArray)
library(SNPRelate)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(data.table)

# --- CONFIGURATION (Adjust these paths) ---

# GDS Output Names (Intermediate binary files)
gds_breeding      <- "C:/Users/tlv329/OneDrive - University of Saskatchewan/UsasK/ACTIVATE_project/data analysis/VariantSet/lentil/genetic_diversity/ACT197_biallelic.gds"
gds_ldp           <- "C:/Users/tlv329/OneDrive - University of Saskatchewan/UsasK/ACTIVATE_project/data analysis/VariantSet/lentil/genetic_diversity/LDP324_nofiltered.gds"
gds_act_filtered  <- "C:/Users/tlv329/OneDrive - University of Saskatchewan/UsasK/ACTIVATE_project/data analysis/VariantSet/lentil/genetic_diversity/ACT197_codingonly.gds"
gds_act_LDpruned  <- "C:/Users/tlv329/OneDrive - University of Saskatchewan/UsasK/ACTIVATE_project/data analysis/VariantSet/lentil/genetic_diversity/ACT197_biallelic_PRUNED.gds"

# LD Calculation Parameters
# These are conservative to save time. Increase 'slide.max.bp' if you expect very long LD.
slide_max_bp      <- 5000000  # Max distance to calculate LD (500kb). Adjust based on Lentil genome size.

  # 3. LD CALCULATION FUNCTION
    #change seqOpen() or snpgdsOpen()

calc_ld_decay <- function(gds_file, snp_subset = NULL, dataset_name) {
  cat(paste0("\n--- Processing: ", dataset_name, " ---\n"))
  
  # 1. Open File
  genofile <- snpgdsOpen(gds_file, allow.duplicate = TRUE)
  on.exit(snpgdsClose(genofile), add = TRUE) # Safety latch to always close file
  
  # 2. Define Window Size (in SNPs)
  # We use 5000 SNPs to ensure we cover 1Mb even in dense regions.
  slide_n <- 4000
  
  cat(paste0("Calculating LD with SNP window = ", slide_n, "...\n"))
  
  # 3. Calculate LD Matrix
  ld_res <- snpgdsLDMat(genofile, 
                        sample.id = NULL, 
                        snp.id = snp_subset, 
                        method = "corr", 
                        slide = slide_n,
                        verbose = FALSE)
  
  # 4. Extract Distance vs R2
  # --- FIX: DETECT GDS FORMAT (SeqArray vs Standard) ---
  
  # Check if "variant.id" exists (SeqArray format)
  node_variant <- index.gdsn(genofile, "variant.id", silent=TRUE)
  
  if (!is.null(node_variant)) {
    # It is a SeqArray file (Created by seqVCF2GDS)
    all_snp_ids <- read.gdsn(node_variant)
    all_positions <- read.gdsn(index.gdsn(genofile, "position"))
  } else {
    # It is a standard SNPRelate file (Created by snpgdsVCF2GDS)
    all_snp_ids <- read.gdsn(index.gdsn(genofile, "snp.id"))
    all_positions <- read.gdsn(index.gdsn(genofile, "snp.position"))
  }
  
  used_snp_ids <- ld_res$snp.id
  
  # Map subset IDs to positions
  match_idx <- match(used_snp_ids, all_snp_ids)
  snp_pos <- all_positions[match_idx]
  
  r2_mat <- ld_res$LD^2
  n_snps <- length(snp_pos)
  
  cat("Extracting and filtering by physical distance...\n")
  
  out_list <- list()
  
  # Loop over Lags (rows of the LD matrix)
  for(k in 1:nrow(r2_mat)) {
    valid_idx <- 1:(n_snps - k)
    
    # Calculate Physical Distances
    dists <- snp_pos[valid_idx + k] - snp_pos[valid_idx]
    r2_vals <- r2_mat[k, valid_idx]
    
    # FILTER: Keep only pairs < slide_max_bp bp
    keep <- which(dists <= slide_max_bp & !is.na(r2_vals))
    
    if(length(keep) > 0) {
      # Downsample to save RAM (max 500 points per lag)
      if(length(keep) > 500) keep <- sample(keep, 500)
      
      chunk <- data.table(dist = dists[keep], r2 = r2_vals[keep])
      out_list[[k]] <- chunk
    }
  }
  
  result_dt <- rbindlist(out_list)
  result_dt[, group := dataset_name]
  
  return(result_dt)
}

  # 4. RUN ANALYSES

# A. LDP (Diversity Panel) - Full (because it is exome already)
# Note: LDP is already "coding", so we don't filter it further usually.
df_LDP324 <- calc_ld_decay(gds_ldp, dataset_name = "LDP324")
write.csv(df_LDP324, "LD_decay_LDP324.csv")
 df_LDP324 <- data.table(read.csv("LD_decay_LDP324.csv"))[,-1]

# B. Breeding Lines - Filtered (Coding Only)
df_ACT197_coding <- calc_ld_decay(gds_act_filtered, dataset_name = "ACT197_coding")
write.csv(df_ACT197_coding, "df_ACT197_coding.csv")
 #df_ACT197_coding <- data.table(read.csv("LD_decay_ACT197_LDpruned.csv", header = T))[-1]

# C. Breeding Lines - Full set (>17M)
df_ACT197 <- calc_ld_decay(gds_act, dataset_name = "ACT197")
write.csv(df_ACT197, "LD_decay_ACT197_LDpruned.csv")
#df_ACT197 <- data.table(read.csv("LD_decay_ACT197_coding.csv", header = T))[-1]

  # 5. COMBINE AND SMOOTH

cat("Combining data...\n")
full_data <- rbind(df_LDP324, df_ACT197_coding)

# Bin the data for plotting (Average r2 every 1kb) to make the line smooth
full_data[, dist_bin := cut(dist, breaks=seq(0, slide_max_bp, by=1000), labels=FALSE)]
full_data[, dist_real := dist_bin * 1000]

summary_data <- full_data[, .(mean_r2 = mean(r2, na.rm=TRUE)), by=.(dist_real, group)]

  # 6. PLOT

cat("Plotting...\n")
pdf("Three_Curve_LD_Decay.pdf", width=8, height=6)

p <- ggplot(summary_data, aes(x=dist_real/1000, y=mean_r2, color=group, linetype=group)) +
  geom_smooth(method = "loess", formula = y~log(x), se = T, size=1.2) + # Loess smoothing line
  scale_color_manual(values=c("ACT197_LDpruned"="grey50", 
                              "ACT197"="blue", 
                              "LDP324"="red")) +
  scale_linetype_manual(values=c("ACT197_LDpruned"="dashed", 
                                 "ACT197"="solid", 
                                 "LDP324"="solid")) +
  labs(title="LD Decay Comparison: Breeding Lines vs. Diversity Panel",
       subtitle="Comparison of LD decay rates across genomic contexts",
       x="Distance (kb)",
       y=expression(Linkage~Disequilibrium~(r^2))) +
  theme_classic() +
  theme(legend.position = c(0.9, 0.9))

print(p)
dev.off()


cat("Done! Check 'Three_Curve_LD_Decay.pdf'\n")

#average markers####
library(SNPRelate)

# 1. Open the pruned file
f_pruned <- snpgdsOpen(gds_ldp, allow.duplicate = T) # Update if your filename differs

# 2. Load the raw chromosome and position data
all_chr_ids <- read.gdsn(index.gdsn(f_pruned, "snp.chromosome"))
all_positions <- read.gdsn(index.gdsn(f_pruned, "snp.position"))

# 3. Define the Valid Chromosomes
valid_chromosomes <- c("Lcu.1GRN.Chr1", "Lcu.1GRN.Chr2", "Lcu.1GRN.Chr3", 
                       "Lcu.1GRN.Chr4", "Lcu.1GRN.Chr5", "Lcu.1GRN.Chr6", 
                       "Lcu.1GRN.Chr7")

# 4. Filter the data to keep only these 7
# We find the indices where the chromosome name matches our list
valid_indices <- which(all_chr_ids %in% valid_chromosomes)

clean_chr <- all_chr_ids[valid_indices]
clean_pos <- all_positions[valid_indices]

cat(paste("Total SNPs in file:", length(all_chr_ids), "\n"))
cat(paste("SNPs on Main Chromosomes:", length(clean_chr), "\n"))


# --- METRIC A: Average Markers per Chromosome ---
# Now calculate table only on the clean data
chr_counts <- table(clean_chr)

cat("\n--- Counts per Main Chromosome ---\n")
print(chr_counts)

avg_per_chr <- mean(chr_counts)
cat(paste0("\nAverage SNPs per Chromosome: ", round(avg_per_chr, 0), "\n"))


# --- METRIC B: Average Spacing (Distance between SNPs) ---

# Function to get spacing for one chromosome
get_spacing <- function(indices) {
  # Extract positions for this specific group of indices
  p <- clean_pos[indices]
  
  # Need at least 2 points to measure distance
  if(length(p) < 2) return(NA)
  
  # Calculate distance between neighbors (pos2-pos1, pos3-pos2...)
  return(mean(diff(p))) 
}

# Apply to the clean data
# tapply works by splitting the indices (1:N) by the chromosome name factor
# We pass the indices to the function so it can look up the positions
spacings <- tapply(seq_along(clean_chr), clean_chr, get_spacing)

# Calculate the global average from the 7 chromosomes
global_avg_spacing_bp <- mean(spacings, na.rm=TRUE)

cat(paste0("Average Spacing: 1 SNP every ", round(global_avg_spacing_bp/1000, 2), " kb\n"))

# --- METRIC C: Density (SNPs per Mb) ---
density_per_mb <- 1000000 / global_avg_spacing_bp
cat(paste0("Marker Density: ~", round(density_per_mb, 2), " SNPs/Mb\n"))

# Close the file
snpgdsClose(f_pruned)

#average markers in common regions####
library(SNPRelate)

# --- CONFIGURATION ---
# 1. Load the Common Keys first
# We need to know which markers are "Common". We saved these in the merged file.
f_merge <- snpgdsOpen("Merged_Analysis.gds")
common_keys <- read.gdsn(index.gdsn(f_merge, "snp.rs.id"))
snpgdsClose(f_merge)

# 2. Open the file you want to analyze (LDP is fine since positions are the same)
f_target <- snpgdsOpen(gds_ldp, allow.duplicate = TRUE) 

# 3. Load the raw chromosome and position data
all_chr_ids <- read.gdsn(index.gdsn(f_target, "snp.chromosome"))
all_positions <- read.gdsn(index.gdsn(f_target, "snp.position"))

# 4. Define the Valid Chromosomes
valid_chromosomes <- c("Lcu.1GRN.Chr1", "Lcu.1GRN.Chr2", "Lcu.1GRN.Chr3", 
                       "Lcu.1GRN.Chr4", "Lcu.1GRN.Chr5", "Lcu.1GRN.Chr6", 
                       "Lcu.1GRN.Chr7")

# --- FILTERING STEP (ADAPTED) ---

# Create keys for the current file
current_keys <- paste(all_chr_ids, all_positions, sep=":")

# Condition A: Is it on a main chromosome?
is_valid_chr <- all_chr_ids %in% valid_chromosomes

# Condition B: Is it a common marker?
is_common <- current_keys %in% common_keys

# Combine Filters
valid_indices <- which(is_valid_chr & is_common)

# Apply Filter
clean_chr <- all_chr_ids[valid_indices]
clean_pos <- all_positions[valid_indices]

cat(paste("Total SNPs in file:      ", length(all_chr_ids), "\n"))
cat(paste("Common SNPs on Main Chrs:", length(clean_chr), "\n"))


# --- METRIC A: Average Markers per Chromosome ---
chr_counts <- table(clean_chr)

cat("\n--- Counts per Main Chromosome (Common Set) ---\n")
print(chr_counts)

avg_per_chr <- mean(chr_counts)
cat(paste0("\nAverage SNPs per Chromosome: ", round(avg_per_chr, 0), "\n"))


# --- METRIC B: Average Spacing (Distance between SNPs) ---
get_spacing <- function(indices) {
  # Extract positions for this specific group of indices
  p <- clean_pos[indices]
  
  # Need at least 2 points to measure distance
  if(length(p) < 2) return(NA)
  
  # Calculate distance between neighbors
  # Note: This assumes positions are sorted (standard in GDS)
  return(mean(diff(p))) 
}

# Apply to the clean data
spacings <- tapply(seq_along(clean_chr), clean_chr, get_spacing)

# Calculate the global average
global_avg_spacing_bp <- mean(spacings, na.rm=TRUE)

cat(paste0("Average Spacing: 1 SNP every ", round(global_avg_spacing_bp/1000, 2), " kb\n"))

# --- METRIC C: Density (SNPs per Mb) ---
density_per_mb <- 1000000 / global_avg_spacing_bp
cat(paste0("Marker Density: ~", round(density_per_mb, 2), " SNPs/Mb\n"))

# Close the file
snpgdsClose(f_target)


#common SNPs LDP324 vs ACT197####
 
library(SNPRelate)
library(gdsfmt)

# --- CONFIGURATION ---
# Update these filenames to match your files
file_LDP      <- "LDP324_nofiltered.gds"          # Your 324 lines file
file_Breeding <- "ACT197_biallelic.gds"           # Your 197 lines file (Full 17M)

# --- HELPER FUNCTION: Get Keys (Chr:Pos) ---
get_keys <- function(gds_fn, name_tag) {
  cat(paste0("Opening ", name_tag, "...\n"))
  
  # Open File
  f <- snpgdsOpen(gds_fn, allow.duplicate=TRUE)
  on.exit(snpgdsClose(f)) # Close when done
  
  # Smart Node Detection (SeqArray vs SNPRelate naming)
  # Check if "snp.position" exists. If not, try "position".
  if (!is.null(index.gdsn(f, "snp.position", silent=TRUE))) {
    # SNPRelate Format
    chr <- read.gdsn(index.gdsn(f, "snp.chromosome"))
    pos <- read.gdsn(index.gdsn(f, "snp.position"))
  } else {
    # SeqArray Format
    chr <- read.gdsn(index.gdsn(f, "chromosome"))
    pos <- read.gdsn(index.gdsn(f, "position"))
  }
  
  # Create unique keys "Chr:Position"
  cat(paste0("Creating identifier keys for ", length(pos), " variants...\n"))
  keys <- paste(chr, pos, sep=":")
  return(keys)
}

# --- MAIN EXECUTION ---

# 1. Get LDP Keys
keys_LDP <- get_keys(file_LDP, "LDP (Diversity)")

# 2. Get Breeding Keys
keys_Breed <- get_keys(file_Breeding, "Breeding (HiFi)")

# 3. Find Intersection
cat("\nCalculating intersection (this might take a minute)...\n")
common_markers <- intersect(keys_LDP, keys_Breed)
n_common <- length(common_markers)

# 4. Report Results
cat("\n==============================================\n")
cat(paste("LDP Total Markers:      ", length(keys_LDP), "\n"))
cat(paste("Breeding Total Markers: ", length(keys_Breed), "\n"))
cat("----------------------------------------------\n")
cat(paste("COMMON MARKERS:         ", n_common, "\n"))
cat("==============================================\n")

# Calculate % overlap relative to LDP (the smaller set)
overlap_pct <- (n_common / length(keys_LDP)) * 100
cat(paste0("Coverage: You have data for ", round(overlap_pct, 2), "% of the LDP markers in your Breeding Set.\n"))


#Common list####
library(SNPRelate)

# 1. Define filenames
fn_LDP   <- "LDP324_nofiltered.gds"
fn_ACT <- "ACT197_biallelic.gds"
fn_Merge <- "Merged_Panel_LDP&ACT.gds"

# 2. Get the Intersection Keys again (re-running logic from before)
# (Assuming you ran the previous code and have 'keys_LDP', 'keys_Breed')
common_keys <- intersect(keys_LDP, keys_Breed)

# 3. Get the SNP IDs for these keys in BOTH files
# We need to know: "Key 'Chr1:1005' is ID #5 in LDP but ID #500 in Breeding"

# Helper to map Keys -> IDs
get_ids_for_keys <- function(gds_fn, target_keys, name_tag) {
  f <- snpgdsOpen(gds_fn, allow.duplicate=TRUE)
  on.exit(snpgdsClose(f))
  
  # Detect Format
  if (!is.null(index.gdsn(f, "snp.position", silent=TRUE))) {
    chr <- read.gdsn(index.gdsn(f, "snp.chromosome"))
    pos <- read.gdsn(index.gdsn(f, "snp.position"))
    ids <- read.gdsn(index.gdsn(f, "snp.id"))
  } else {
    chr <- read.gdsn(index.gdsn(f, "chromosome"))
    pos <- read.gdsn(index.gdsn(f, "position"))
    ids <- read.gdsn(index.gdsn(f, "variant.id"))
  }
  
  current_keys <- paste(chr, pos, sep=":")
  
  # Match
  # mapping is a vector where mapping[i] is the index in 'current_keys' 
  # that matches target_keys[i]
  mapping <- match(target_keys, current_keys)
  
  # Return the IDs
  return(ids[mapping])
}

cat("Mapping common keys to IDs...\n")
ids_in_LDP   <- get_ids_for_keys(fn_LDP, common_keys, "LDP")
ids_in_ACT <- get_ids_for_keys(fn_ACT, common_keys, "ACT")

# Check for NAs (Should be none if intersect worked correctly)
if(any(is.na(ids_in_LDP)) | any(is.na(ids_in_ACT))) {
  stop("Error: key mismatch during ID retrieval.")
}

#The reason Merged_Analysis.gds contains dummy positions is that when we created 
#it (specifically for PCA), we used 1:ncol(geno_Combined) for the positions to simply 
#get the matrix into a file format snpgdsPCA could accept without errors.

#merge files####
# 4. Merge
# Note: snpgdsCombineGeno usually requires the SNP IDs to match exactly or be aligned.
# A safer way for potentially mismatched IDs is to use "snpgdsCombineGeno" 
# but we must ensure we only pass the overlap.

# Actually, the easiest way to PCA grouped by panel is NOT to physically merge the huge files 
# (which can be messy with different ID schemes), but to run PCA on the smaller dataset (LDP) 
# and PROJECT the breeding lines onto it.

# OR: We extract the genotype matrix for the common SNPs from both and rbind() them.
# Given you have 175k markers, this fits in RAM (175k * 500 samples is small).

cat("Extracting Genotype Matrices for Common SNPs...\n")

# A. Extract LDP Matrix
f_ldp <- snpgdsOpen(fn_LDP)
geno_LDP <- snpgdsGetGeno(f_ldp, snp.id=ids_in_LDP, verbose=FALSE) 
# Result: [Samples x SNPs]
samples_LDP <- read.gdsn(index.gdsn(f_ldp, "sample.id"))
snpgdsClose(f_ldp)

# B. Extract Breeding Matrix
f_ACT <- snpgdsOpen(fn_ACT)
geno_ACT <- snpgdsGetGeno(f_ACT, snp.id=ids_in_ACT, verbose=FALSE)
samples_ACT <- read.gdsn(index.gdsn(f_ACT, "sample.id")) # or "sample.id" node
snpgdsClose(f_ACT)

# C. Combine Matrices
# We must ensure columns (SNPs) are in same order. 
# Since we used the same 'common_keys' vector to fetch IDs, they are aligned!
geno_Combined <- rbind(geno_LDP, geno_ACT)
all_samples <- c(samples_LDP, samples_ACT)

# Create Group Labels
groups <- c(rep("LDP Panel", length(samples_LDP)), 
            rep("ACT Panel", length(samples_ACT)))

cat("Combined Matrix Dimensions:", dim(geno_Combined), "\n")

#PCA merged data####
# 5. Save to a temporary GDS for analysis
# We create a new GDS with the combined data
snpgdsCreateGeno("Merged_Analysis.gds", 
                 genmat = geno_Combined, 
                 sample.id = all_samples,
                 snpfirstdim=F,
                 snp.id = 1:ncol(geno_Combined), # Dummy IDs 1..N
                 snp.rs.id = common_keys,        # Use keys as names
                 snp.chromosome = rep(1, ncol(geno_Combined)), # Dummy chr
                 snp.position = 1:ncol(geno_Combined))         # Dummy pos

# 6. Run PCA
f_merge <- snpgdsOpen("Merged_Analysis.gds")
pca <- snpgdsPCA(f_merge, autosome.only=FALSE)
snpgdsClose(f_merge)

# 7. Plot
pca_df <- data.frame(sample.id = pca$sample.id,
                     EV1 = pca$eigenvect[,1],
                     EV2 = pca$eigenvect[,2],
                     Group = groups)

pc.percent <- pca$varprop * 100
pc1_lab <- paste0("PC1 (", round(pc.percent[1], 2), "%)")
pc2_lab <- paste0("PC2 (", round(pc.percent[2], 2), "%)")

library(ggplot2)
p <- ggplot(pca_df, aes(x=EV1, y=EV2, color=Group)) +
  geom_point(alpha=0.6) +
  labs(title="PCA: Breeding Lines vs Diversity Panel",
       subtitle=paste0("Based on ", n_common, " common SNPs"),
       x=pc1_lab,
       y=pc2_lab) +
  theme_minimal()

ggsave("PCA_PC1-PC2.png", p, width = 6, height = 4, bg = "White")

#make it interactive!
# 1. Install plotly if you haven't already
if (!require("plotly")) install.packages("plotly")

library(ggplot2)
library(plotly)

# 2. Prepare your Data Frame
# Ensure your 'pca_df' from the previous step exists and has these columns:
# sample.id, EV1, EV2, Group

# If you need to recreate the dataframe, here is the structure:
# pca_df <- data.frame(sample.id = pca$sample.id,
#                      EV1 = pca$eigenvect[,1],
#                      EV2 = pca$eigenvect[,2],
#                      Group = groups)

# 3. Create the Static ggplot first (slightly modified for plotly)
# Note the 'text' aesthetic - this tells plotly what to show when you hover!

p <- ggplot(pca_df, aes(x=EV1, y=EV2, color=Group, text=sample.id)) +
  geom_point(alpha=0.7, size=2) +
  labs(title="Interactive PCA: Breeding vs Diversity",
       x=paste0("PC1 (", round(pc.percent[1], 2), "%)"),
       y=paste0("PC2 (", round(pc.percent[2], 2), "%)")) +
  theme_minimal() +
  scale_color_manual(values=c("ACT Panel"="blue", "LDP Panel"="red"))

# 4. Convert to Interactive Plot
# tooltip="text" tells it to display the 'text' aesthetic we mapped to sample.id
interactive_plot <- ggplotly(p, tooltip="text")

# 5. Show it
interactive_plot

# 6. (Optional) Save as a standalone HTML file
# You can open this file in any web browser later without opening R
htmlwidgets::saveWidget(interactive_plot, "Interactive_Lentil_PCA.html")

#more dimmensions? I got you
# 1. Install/Load GGally if needed
if (!require("GGally")) install.packages("GGally")
library(GGally)
library(SNPRelate)

# 2. Run PCA (if not already in memory)
# We assume 'Merged_Analysis.gds' exists from previous steps
f_merge <- snpgdsOpen("Merged_Analysis.gds")

# Calculate more components (we need at least 5)
pca <- snpgdsPCA(f_merge, autosome.only=FALSE, eigen.cnt = 10) 
snpgdsClose(f_merge)

# 3. Create the Extended Data Frame (PC1 to PC5)
# Note: Ensure the 'groups' vector (Breeding vs Diversity) matches sample order!
# If 'groups' isn't in memory, recreate it based on sample IDs or previous logic.

pca_df_multi <- data.frame(
  sample.id = pca$sample.id,
  PC1 = pca$eigenvect[,1],
  PC2 = pca$eigenvect[,2],
  PC3 = pca$eigenvect[,3],
  PC4 = pca$eigenvect[,4],
  PC5 = pca$eigenvect[,5],
  Group = as.factor(groups) # Ensure this exists
)

# 4. Create the Labels (Variance Explained)
pc.percent <- pca$varprop * 100
pc_labels <- paste0("PC", 1:5, "\n(", round(pc.percent[1:5], 1), "%)")

# Assign these labels to the column names for the plot
colnames(pca_df_multi)[2:6] <- pc_labels

# 5. Generate the Pairs Plot
# We use 'ggpairs' to create the matrix.
# - upper="null": Hides the upper triangle (redundant) to save space/ink
# - diag="density": Shows distribution curves on the diagonal
# - lower: Shows the scatter plots

p <- ggpairs(pca_df_multi, 
             columns = 2:6,        # Columns containing PC data
             aes(color = Group, alpha = 0.5), # Color by Group
             
             # Customize the "Lower" triangle (Scatterplots)
             lower = list(continuous = wrap("points", size=1.5)),
             
             # Customize the "Diagonal" (Density plots)
             diag = list(continuous = wrap("densityDiag", alpha=0.3)),
             
             # Hide the "Upper" triangle (Correlations are 0 by definition in PCA)
             upper = list(continuous = "blank")
) +
  scale_color_manual(values = c("ACT Panel"="blue", "LDP Panel"="red")) +
  scale_fill_manual(values = c("ACT Panel"="blue", "LDP Panel"="red")) +
  labs(title = "PCA Matrix: PC1 to PC5",
       subtitle = paste("Based on", length(pca$snp.id), "common SNPs")) +
  theme_bw()

# 6. Show Plot
print(p)

# 7. Save (It needs to be large to see details)
ggsave("PCA_Pairs_Plot_PC1-PC5.png", p, width = 10, height = 7.5, bg = "White")

#PCA adding origin in the LDP####
library(ggplot2)
library(dplyr)

# 1. Prepare Base Data
# Assuming 'pca_df' already exists from your previous PCA run. 
# If not, recreate it from the 'pca' object:
# pca_df <- data.frame(sample.id = pca$sample.id, EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2])

# 2. Create Group 1: Panel (Based on Suffix)
# We use grepl to search for the specific endings
pca_df$Panel <- case_when(
  grepl("_ACT$", pca_df$sample.id) ~ "ACTIVATE (Breeding)",
  grepl("_AGL$", pca_df$sample.id) ~ "LDP (Diversity)",
  TRUE ~ "Unknown" # Safety catch
)

# 3. Create Group 2: Climate (Merge with Metadata)
# Load the origin file
ldp_meta <- read.csv("LDP324_origin.csv", stringsAsFactors = FALSE)

# Merge: pca_df (x) + ldp_meta (y)
# all.x = TRUE ensures we keep the ACT lines even though they aren't in the CSV
pca_df <- merge(pca_df, ldp_meta, by.x = "sample.id", by.y = "Name", all.x = TRUE)

# Rename 'Classification' to 'Climate' for clarity
# ACT lines will have NA in this column initially
pca_df$Climate <- pca_df$Classification

# 4. Define Shape Categories
# The user requested: Breeding + Other + NA -> Circle
# Mediterranean -> Square
# Temperate -> Triangle
# South Asian -> X

# We create a new column specifically for controlling the SHAPE aesthetic
pca_df$Shape_Group <- pca_df$Climate

# Fill NAs (Breeding lines) with a label that maps to Circle
pca_df$Shape_Group[is.na(pca_df$Shape_Group)] <- "Breeding/Other"
# Map "Other" to the same label so they look the same
pca_df$Shape_Group[pca_df$Shape_Group == "Other"] <- "Breeding/Other"

# 5. Define the Plot
# Check variance for labels (if pca object is available)
if(exists("pca")) {
  pc_percent <- pca$varprop * 100
  xlab_txt <- paste0("PC1 (", round(pc_percent[1], 2), "%)")
  ylab_txt <- paste0("PC2 (", round(pc_percent[2], 2), "%)")
} else {
  xlab_txt <- "PC1"; ylab_txt <- "PC2"
}

ggplot(pca_df, aes(x = EV1, y = EV2, color = Panel, shape = Shape_Group)) +
  
  geom_point(size = 3, alpha = 0.7) +
  
  # --- COLOR: Panel ---
  scale_color_manual(values = c("ACTIVATE (Breeding)" = "blue", 
                                "LDP (Diversity)" = "firebrick")) +
  
  # --- SHAPE: Climate ---
  # R Shape Codes: 16=Circle, 15=Square, 17=Triangle, 4=X (Cross)
  scale_shape_manual(values = c("Breeding/Other" = 16,  # Circle
                                "Mediterranean"  = 15,  # Square
                                "Temperate"      = 17,  # Triangle
                                "South Asian"    = 4)) + # X
  
  labs(title = "PCA: Breeding vs Diversity Panel by Climate",
       subtitle = "Shapes indicate agro-climatic adaptation (LDP only)",
       x = xlab_txt, 
       y = ylab_txt,
       shape = "Climate Origin",
       color = "Panel") +
  
  theme_bw() +
  theme(legend.position = "right")

#Make it interactive!
library(ggplot2)
library(plotly)

# 1. Ensure Data is Ready (Same logic as previous step)
# Make sure pca_df has columns: EV1, EV2, Panel, Shape_Group, sample.id

# 2. Define the ggplot object (assign to 'p')
# Note: We add the 'text' aesthetic. This is what Plotly uses for the hover box.
p <- ggplot(pca_df, aes(x = EV1, y = EV2, 
                        color = Panel, 
                        shape = Shape_Group,
                        text = paste("Line:", sample.id, 
                                     "<br>Panel:", Panel,
                                     "<br>Climate:", Shape_Group,
                                     "<br>Country:", Origin,
                                     "<br>GenGroup:", GenGroup))) +
  
  geom_point(size = 3, alpha = 0.7) +
  
  # --- COLORS ---
  scale_color_manual(values = c("ACTIVATE (Breeding)" = "blue", 
                                "LDP (Diversity)" = "firebrick")) +
  
  # --- SHAPES ---
  # Note: Plotly does its best to match these R shapes to WebGL shapes
  scale_shape_manual(values = c("Breeding/Other" = 16,  # Circle
                                "Mediterranean"  = 15,  # Square
                                "Temperate"      = 17,  # Triangle
                                "South Asian"    = 4)) + # X
  
  labs(title = "Interactive PCA: Breeding vs Diversity by Climate",
       x = "PC1", 
       y = "PC2",
       shape = "Climate",
       color = "Panel") +
  
  theme_bw()

# 3. Convert to Interactive Plotly Object
# tooltip="text" forces it to display ONLY what we wrote in the paste() command above
interactive_pca <- ggplotly(p, tooltip = "text")

# 4. Display it in RStudio Viewer
interactive_pca

# 5. Save as HTML (Optional)
# This creates a file you can open in Chrome/Edge to share with colleagues
htmlwidgets::saveWidget(interactive_pca, "Interactive_PCA_Climate.html")

#PCA by GenGroup####
library(ggplot2)
library(plotly)
library(dplyr)
library(RColorBrewer) # Good for distinct palettes

# 1. Prepare Base Data (Recreate if needed)
# Ensure pca_df has sample.id, EV1, EV2
# if(!exists("pca_df")) pca_df <- data.frame(sample.id = pca$sample.id, EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2])

# 2. Define Panel (Shape)
pca_df$Panel <- case_when(
  grepl("_ACT$", pca_df$sample.id) ~ "ACTIVATE (Breeding)",
  grepl("_AGL$", pca_df$sample.id) ~ "LDP (Diversity)",
  TRUE ~ "Unknown"
)

# 3. Merge GenGroup Metadata
# Assuming GenGroup is in your origin file. Update filename if it's different.
ldp_meta <- read.csv("LDP324_origin.csv", stringsAsFactors = FALSE)

# Check if 'GenGroup' column exists, if not, adjust column name
# Merge
pca_df_merged <- merge(pca_df, ldp_meta, by.x="sample.id", by.y="Name", all.x=TRUE)

# 4. Handle Colors (GenGroup)
# The LDP lines will have GenGroup (e.g., "Cluster1", "Cluster2").
# The ACT lines will be NA. Let's label them clearly.
pca_df_merged$GenGroup[is.na(pca_df_merged$GenGroup)] <- "Breeding Program"

# Ensure GenGroup is a factor (optional, helps with legend order)
# pca_df_merged$GenGroup <- as.factor(pca_df_merged$GenGroup)

# 5. Define Plot
# We use a custom color palette because 8 clusters + 1 Breeding group = 9 colors.
# Paired or Set1 are good qualitative palettes.

p <- ggplot(pca_df_merged, aes(x = EV1, y = EV2, 
                               color = GenGroup, 
                               shape = Panel,
                               text = paste("Line:", sample.id, 
                                            "<br>Panel:", Panel,
                                            "<br>Country:", Origin,
                                            "<br>GenGroup:", GenGroup))) +
  
  geom_point(size = 3, alpha = 0.8) +
  
  # --- SHAPES ---
  # LDP = X (4), ACT = Circle (16)
  scale_shape_manual(values = c("LDP (Diversity)" = 16, 
                                "ACTIVATE (Breeding)" = 4)) +
  
  # --- COLORS ---
  # Assign "Breeding Program" to Black/Blue so it pops out against the 8 clusters
  # We let ggplot pick the rest, or define a manual scale if you have specific cluster colors.
  scale_color_manual(values = c("Breeding Program" = "black", 
                                # Define colors for your 8 clusters (1-8) if known, or use defaults
                                "1"="#E41A1C", "2"="#377EB8", "3"="#4DAF4A", "4"="#984EA3",
                                "5"="#FF7F00", "6"="#FFFF33", "7"="#A65628", "8"="#F781BF",
                                # Fallback if names are "Cluster 1" etc
                                "Cluster 1"="#E41A1C", "Cluster 2"="#377EB8")) + 
  # Note: If your GenGroups are named differently (e.g., "K1", "K2"), remove the manual scale 
  # or update the names above. Removing 'scale_color_manual' lets R pick auto colors.
  
  labs(title = "PCA: Breeding Lines vs LDP by Genetic Clusters",
       subtitle = "Based on 136957 common SNPs",
       x = pc1_lab, 
       y = pc2_lab,
       color = "Genetic Group (K)",
       shape = "Panel") +
  
  theme_bw()

ggsave("PCA_PC1-PC2_GenGroup.png", p, width = 8, height = 5, bg = "White")

# 6. Interactive
interactive_pca_gen <- ggplotly(p, tooltip = "text")
interactive_pca_gen

# 7. Save
htmlwidgets::saveWidget(interactive_pca_gen, "Interactive_PCA_GenGroups.html")

#Calculate Fst (LDPvsACT)####
library(SNPRelate)

# 1. Open the Merged File
f_merge <- snpgdsOpen("Merged_Analysis.gds")

# 2. Define Populations (Re-creating the list from your previous steps)
# Ensure 'groups' vector exists from the PCA step! 
# groups <- c(rep("Diversity Panel", 324), rep("Breeding Panel", 197))

pop_list <- list(
  "Breeding" = read.gdsn(index.gdsn(f_merge, "sample.id"))[groups == "ACT Panel"],
  "Diversity" = read.gdsn(index.gdsn(f_merge, "sample.id"))[groups == "LDP Panel"]
)

# 3. Calculate Fst (FIXED)
cat("Calculating genome-wide Fst...\n")

fst_res <- snpgdsFst(f_merge, 
                     sample.id = unlist(pop_list), 
                     population = as.factor(groups),
                     method = "W&C84",       # <--- UPDATED HERE
                     autosome.only = FALSE,
                     verbose = TRUE)

# 4. Report Results
cat("\n==============================================\n")
cat(paste("Global Fst (Weighted):", round(fst_res$Fst, 4), "\n"))
cat(paste("Mean Fst (Unweighted):", round(fst_res$MeanFst, 4), "\n"))
cat("==============================================\n")

# 5. Histogram
hist(fst_res$FstSNP, 
     breaks=50, 
     col="cadetblue", 
     main="Distribution of per-SNP Fst", 
     xlab="Fst per Locus")

snpgdsClose(f_merge)

#genetic diversity stats####
library(SNPRelate)
library(dplyr)
library(tidyr)

library(SNPRelate)
library(dplyr)

calc_diversity_stats <- function(gds_filename, label) {
  
  cat(paste0("\n--- Analyzing: ", label, " ---\n"))
  f <- snpgdsOpen(gds_filename, allow.duplicate=TRUE)
  on.exit(snpgdsClose(f))
  
  # 1. Get SNP Info
  snp_ids_all <- read.gdsn(index.gdsn(f, "snp.id"))
  chroms_all  <- read.gdsn(index.gdsn(f, "snp.chromosome"))
  pos_all     <- read.gdsn(index.gdsn(f, "snp.position"))
  
  # --- FILTERING STEP ---
  # Define valid chromosomes
  valid_chrs <- c("Lcu.1GRN.Chr1", "Lcu.1GRN.Chr2", "Lcu.1GRN.Chr3", 
                  "Lcu.1GRN.Chr4", "Lcu.1GRN.Chr5", "Lcu.1GRN.Chr6", 
                  "Lcu.1GRN.Chr7")
  
  # Find indices of valid SNPs
  valid_idx <- which(chroms_all %in% valid_chrs)
  
  # Subset the metadata
  snp_ids <- snp_ids_all[valid_idx]
  chroms  <- chroms_all[valid_idx]
  pos     <- pos_all[valid_idx]
  n_snps  <- length(snp_ids)
  
  cat(paste0("Filtering: Kept ", n_snps, " variants on main chromosomes (Discarded ", length(snp_ids_all) - n_snps, " scaffolds).\n"))
  
  # 2. Calculate Allele Frequencies (For He & PIC)
  cat("Calculating Allele Frequencies...\n")
  # We pass 'snp.id' to calculate ONLY for the valid subset
  freq_res <- snpgdsSNPRateFreq(f, snp.id=snp_ids)
  p <- freq_res$AlleleFreq
  q <- 1 - p
  
  # Expected Het & PIC
  He <- 2 * p * q
  PIC <- He - (2 * p^2 * q^2)
  
  # 3. Calculate Observed Heterozygosity (Ho) - MANUALLY
  cat("Calculating Observed Heterozygosity (Ho)...\n")
  
  Ho <- numeric(n_snps)
  chunk_size <- 10000
  
  # Loop guarantees alignment with our subsetted 'snp_ids'
  for(i in seq(1, n_snps, by=chunk_size)) {
    
    end <- min(i + chunk_size - 1, n_snps)
    idx_range <- i:end
    
    # Get IDs for this chunk (from our VALID subset)
    chunk_ids <- snp_ids[idx_range]
    
    # Get Genotypes
    mat <- snpgdsGetGeno(f, snp.id=chunk_ids)
    
    # Calculate
    n_het <- colSums(mat == 1, na.rm=TRUE)
    n_valid <- colSums(!is.na(mat))
    
    Ho[idx_range] <- n_het / n_valid
    
    if(i %% 100000 == 1) cat(".")
  }
  cat("\n")
  
  # 4. Assemble Data Frame
  df <- data.frame(
    SNP_ID = snp_ids,
    Chr = chroms,
    Pos = pos,
    Ho = Ho,
    He = He,
    PIC = PIC,
    Dataset = label
  )
  
  # 5. Clean Up
  df <- df %>% mutate(
    He = ifelse(He < 1e-6, NA, He), 
    Fis = 1 - (Ho / He),
    Ratio_HoHe = Ho / He
  )
  
  return(df)
}

#run the diversity stats####
# --- CONFIGURATION ---
file_LDP      <- "LDP324_nofiltered.gds"
file_Breeding <- "ACT197_codingonly.gds"
file_Breeding_total <-"ACT197_biallelic.gds"
bin_size_bp   <- 1000000  # 1 Mb Window size (Adjust as needed)

# 1. Run the Calculator
stats_LDP   <- calc_diversity_stats(file_LDP, "LDP (Diversity)")
stats_Breed <- calc_diversity_stats(file_Breeding, "Breeding (Elite)")

# 2. ANALYSIS A: By Chromosome
# We group by Chr and take the median of all SNPs
avg_by_chr_LDP <- stats_LDP %>% 
  group_by(Chr) %>% 
  summarise(
    Mean_Ho = median(Ho, na.rm=TRUE),
    Mean_He = median(He, na.rm=TRUE),
    Mean_PIC = median(PIC, na.rm=TRUE),
    Mean_Fis = median(Fis, na.rm=TRUE)
  ) %>% mutate(Panel = "LDP")

avg_by_chr_Breed <- stats_Breed %>% 
  group_by(Chr) %>% 
  summarise(
    Mean_Ho = median(Ho, na.rm=TRUE),
    Mean_He = median(He, na.rm=TRUE),
    Mean_PIC = median(PIC, na.rm=TRUE),
    Mean_Fis = median(Fis, na.rm=TRUE)
  ) %>% mutate(Panel = "Breeding")

# Combine and View
final_chr_stats <- rbind(avg_by_chr_LDP, avg_by_chr_Breed)
print(final_chr_stats)


# 3. ANALYSIS B: By Genomic Bin (Sliding Window)
# This creates a "Manhattan-style" dataset for diversity

calc_bins <- function(df, bin_size) {
  df %>%
    mutate(Bin_Start = floor(Pos / bin_size) * bin_size) %>%
    group_by(Chr, Bin_Start, Dataset) %>%
    summarise(
      n_SNPs = n(),
      Bin_Ho = mean(Ho, na.rm=TRUE),
      Bin_He = mean(He, na.rm=TRUE),
      Bin_PIC = mean(PIC, na.rm=TRUE),
      Bin_Fis = mean(Fis, na.rm=TRUE),
      .groups = 'drop'
    )
}

bins_LDP   <- calc_bins(stats_LDP, bin_size_bp)
bins_Breed <- calc_bins(stats_Breed, bin_size_bp)

# Merge Bins for Direct Comparison (Erosion Scan)
# We want to compare the SAME bin across both datasets
comparison_bins <- merge(bins_LDP, bins_Breed, 
                         by=c("Chr", "Bin_Start"), 
                         suffixes=c(".LDP", ".Breed"))

# Calculate "Delta He" (Loss of Diversity)
# Positive value = Diversity lost in breeding
comparison_bins$Diversity_Loss <- comparison_bins$Bin_He.LDP - comparison_bins$Bin_He.Breed

# Find the Top 10 "Erosion Hotspots"
hotspots <- comparison_bins %>% 
  arrange(desc(Diversity_Loss)) %>% 
  head(10)

cat("\n--- Top 10 Regions of Genetic Erosion (High Diversity in LDP, Low in Breeding) ---\n")
print(hotspots[, c("Chr", "Bin_Start", "n_SNPs.LDP", "Bin_He.LDP", "Bin_He.Breed", "Diversity_Loss")])

#plotting####
 #data transformation
library(dplyr)
library(ggplot2)

# 1. Prepare the Data
# We work with 'comparison_bins' from the previous step.
# Let's clean the chromosome names for plotting (remove "Lcu.1GRN.Chr")
plot_data <- comparison_bins %>%
  mutate(Chr_Num = as.numeric(gsub("Lcu.1GRN.Chr", "", Chr))) %>%
  arrange(Chr_Num, Bin_Start)

# 2. Calculate Chromosome Offsets
# We determine the length of each chromosome to stack them side-by-side
data_cum <- plot_data %>%
  group_by(Chr_Num) %>%
  summarise(max_bp = max(Bin_Start)) %>%
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>%
  select(Chr_Num, bp_add)

# 3. Add Cumulative Position to the Main Data
plot_data <- plot_data %>%
  inner_join(data_cum, by = "Chr_Num") %>%
  mutate(bp_cum = Bin_Start + bp_add)

# 4. Calculate Axis Label Positions (Center of each Chr)
axis_set <- plot_data %>%
  group_by(Chr_Num) %>%
  summarize(center = (max(bp_cum) + min(bp_cum)) / 2)

#ploting genetic erosion
# Define a custom color palette for alternating chromosomes (Classic style)
# Grey and Blue allows red "peaks" to stand out
chr_colors <- c("grey70", "grey40")

ggplot(plot_data, aes(x = bp_cum, y = Diversity_Loss)) +
  
  # A. The Main Data Points
  # Color points by Chromosome (Alternating)
  geom_point(aes(color = as.factor(Chr_Num)), alpha = 0.8, size = 1.5) +
  scale_color_manual(values = rep(chr_colors, 7)) +
  
  # B. Highlight the "Zero" line
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  
  # C. Highlight "Erosion Peaks" (Top 5% of loss)
  # We add a red layer for the most eroded bins
  geom_point(data = subset(plot_data, Diversity_Loss > quantile(Diversity_Loss, 0.95)),
             aes(x = bp_cum, y = Diversity_Loss),
             color = "firebrick", size = 2) +
  
  # D. Custom X-Axis (Show Chr Numbers instead of huge BP numbers)
  scale_x_continuous(label = axis_set$Chr_Num, breaks = axis_set$center) +
  
  # E. Formatting
  labs(title = "Genetic Erosion Scan: Diversity Loss in Breeding Lines",
       subtitle = "Red peaks indicate loci where Breeding Diversity << LDP Diversity",
       x = "Chromosome",
       y = expression(Delta ~ H[e] ~ (LDP - Breeding))) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


#multifacet
# 1. Reshape Data to "Long" format for Faceting
# We want separate rows for "Bin_He.Breed", "Bin_He.LDP", etc.
library(tidyr)

long_data <- plot_data %>%
  select(Chr_Num, bp_cum, Bin_He.LDP, Bin_He.Breed, Bin_Ho.LDP, Bin_Ho.Breed, Bin_PIC.LDP, Bin_PIC.Breed) %>%
  pivot_longer(cols = starts_with("Bin"), 
               names_to = "Metric_Panel", 
               values_to = "Value") %>%
  separate(Metric_Panel, into = c("Stat", "Panel"), sep = "\\.")

# 2. Plot Everything Stacked
ggplot(long_data, aes(x = bp_cum, y = Value, color = Panel)) +
  
  # Use lines instead of points for a "smoother" look across bins
  geom_point(alpha = 0.8) +
  
  # Facet by Statistic (He, Ho, PIC)
  facet_grid(Stat ~ Chr_Num, scales = "free") +
  
  scale_x_continuous(label = axis_set$Chr_Num, breaks = axis_set$center) +
  scale_color_manual(values = c("LDP" = "forestgreen", "Breed" = "navyblue")) +
  
  labs(title = "Genome-Wide Diversity Metrics by Bin",
       x = "Chromosome", 
       y = "Value") +
  
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

#same but in common markers####
library(SNPRelate)

# Open the merged file where we saved the common set
f_merge <- snpgdsOpen("Merged_Analysis.gds")

# Retrieve the "keys" (Chr:Pos strings) we saved in snp.rs.id
common_keys <- read.gdsn(index.gdsn(f_merge, "snp.rs.id"))
snpgdsClose(f_merge)

cat(paste("Targeting", length(common_keys), "common markers.\n"))

calc_diversity_stats_filtered <- function(gds_filename, label, target_keys=NULL) {
  
  cat(paste0("\n--- Analyzing: ", label, " ---\n"))
  f <- snpgdsOpen(gds_filename, allow.duplicate=TRUE)
  on.exit(snpgdsClose(f))
  
  # 1. Read Metadata
  snp_ids_all <- read.gdsn(index.gdsn(f, "snp.id"))
  
  # Handle Format Differences (SeqArray vs SNPRelate)
  if (!is.null(index.gdsn(f, "snp.chromosome", silent=TRUE))) {
    chroms_all <- read.gdsn(index.gdsn(f, "snp.chromosome"))
    pos_all    <- read.gdsn(index.gdsn(f, "snp.position"))
  } else {
    chroms_all <- read.gdsn(index.gdsn(f, "chromosome"))
    pos_all    <- read.gdsn(index.gdsn(f, "position"))
  }
  
  # 2. CONSTRUCT FILTERS
  
  # Filter A: Chromosomes Only (Remove Scaffolds)
  valid_chrs <- c("Lcu.1GRN.Chr1", "Lcu.1GRN.Chr2", "Lcu.1GRN.Chr3", 
                  "Lcu.1GRN.Chr4", "Lcu.1GRN.Chr5", "Lcu.1GRN.Chr6", 
                  "Lcu.1GRN.Chr7")
  
  # Boolean Vector: TRUE if chromosome is valid
  is_valid_chr <- chroms_all %in% valid_chrs
  
  # Filter B: Common Markers (Intersection)
  if(!is.null(target_keys)) {
    current_keys <- paste(chroms_all, pos_all, sep=":")
    is_common_marker <- current_keys %in% target_keys
    
    # COMBINE: Must be Valid Chr AND Common Marker
    keep_mask <- is_valid_chr & is_common_marker
    cat(paste0("Filtering: Applying Chromosome + Intersection Filter...\n"))
    
  } else {
    # Only Chromosome Filter
    keep_mask <- is_valid_chr
    cat(paste0("Filtering: Applying Chromosome Filter only...\n"))
  }
  
  # 3. APPLY SUBSET
  # Get the indices of TRUE values
  keep_idx <- which(keep_mask)
  
  snp_ids <- snp_ids_all[keep_idx]
  chroms  <- chroms_all[keep_idx]
  pos     <- pos_all[keep_idx]
  
  n_snps <- length(snp_ids)
  cat(paste0("Final Set: ", n_snps, " markers retained.\n"))
  
  if(n_snps == 0) stop("No markers retained after filtering!")
  
  # 4. Calculate Allele Frequencies
  cat("Calculating Allele Frequencies...\n")
  freq_res <- snpgdsSNPRateFreq(f, snp.id=snp_ids)
  p <- freq_res$AlleleFreq
  q <- 1 - p
  
  # Stats
  He <- 2 * p * q
  PIC <- He - (2 * p^2 * q^2)
  
  # 5. Calculate Ho (Manual Loop)
  cat("Calculating Observed Heterozygosity (Ho)...\n")
  Ho <- numeric(n_snps)
  chunk_size <- 10000
  
  for(i in seq(1, n_snps, by=chunk_size)) {
    end <- min(i + chunk_size - 1, n_snps)
    idx_range <- i:end
    
    chunk_ids <- snp_ids[idx_range]
    mat <- snpgdsGetGeno(f, snp.id=chunk_ids)
    
    n_het <- colSums(mat == 1, na.rm=TRUE)
    n_valid <- colSums(!is.na(mat))
    
    Ho[idx_range] <- n_het / n_valid
    
    if(i %% 50000 == 1) cat(".")
  }
  cat("\n")
  
  # 6. Output
  df <- data.frame(
    SNP_ID = snp_ids,
    Chr = chroms,
    Pos = pos,
    Ho = Ho,
    He = He,
    PIC = PIC,
    Dataset = label
  )
  
  df <- df %>% mutate(
    He = ifelse(He < 1e-6, NA, He), 
    Fis = 1 - (Ho / He),
    Ratio_HoHe = Ho / He
  )
  
  return(df)
}

stats_LDP   <- calc_diversity_stats_filtered(file_LDP, "LDP (Diversity)", target_keys = common_keys)
stats_Breed <- calc_diversity_stats_filtered(file_Breeding_total, "Breeding (Elite)", target_keys = common_keys)

#Proceed with Binning/Plotting as before...

#Common markers gene coverage####
 #Create the "Merged" File with REAL Coordinates
library(SNPRelate)
library(gdsfmt)

# --- CONFIGURATION ---
fn_LDP    <- "LDP324_nofiltered.gds"
fn_Breed  <- "ACT197_biallelic.gds"
fn_Output <- "Merged_Analysis_RealCoords.gds" # New filename to avoid confusion

# --- HELPER: Get Keys (Chr:Pos) ---
get_keys <- function(gds_fn) {
  f <- snpgdsOpen(gds_fn, allow.duplicate=TRUE)
  on.exit(snpgdsClose(f))
  
  # Handle SeqArray vs SNPRelate nodes
  if (!is.null(index.gdsn(f, "snp.position", silent=TRUE))) {
    chr <- read.gdsn(index.gdsn(f, "snp.chromosome"))
    pos <- read.gdsn(index.gdsn(f, "snp.position"))
    ids <- read.gdsn(index.gdsn(f, "snp.id"))
  } else {
    chr <- read.gdsn(index.gdsn(f, "chromosome"))
    pos <- read.gdsn(index.gdsn(f, "position"))
    ids <- read.gdsn(index.gdsn(f, "variant.id"))
  }
  return(data.frame(ID=ids, Key=paste(chr, pos, sep=":"), Chr=chr, Pos=pos, stringsAsFactors=FALSE))
}

# 1. Get Map Info from Both Files
cat("Reading map info from LDP...\n")
map_LDP <- get_keys(fn_LDP)

cat("Reading map info from Breeding Panel...\n")
map_Breed <- get_keys(fn_Breed)

# 2. Find Common Markers (Based on Chr:Pos Key)
common_keys <- intersect(map_LDP$Key, map_Breed$Key)
n_common <- length(common_keys)
cat(paste("Found", n_common, "common markers.\n"))

# 3. Filter Maps to Common Set (to get IDs)
# We match against the common keys to keep the order identical
# We use the LDP map as the "Master Map" for coordinates since they are the same
map_LDP_subset <- map_LDP[match(common_keys, map_LDP$Key), ]
map_Breed_subset <- map_Breed[match(common_keys, map_Breed$Key), ]

# Verify alignment
if(!all(map_LDP_subset$Key == map_Breed_subset$Key)) stop("CRITICAL ERROR: Keys not aligned!")

# 4. Extract Genotypes
cat("Extracting Genotype Matrices...\n")

# A. LDP Matrix
f_ldp <- snpgdsOpen(fn_LDP, allow.duplicate=TRUE)
geno_LDP <- snpgdsGetGeno(f_ldp, snp.id=map_LDP_subset$ID, verbose=FALSE)
samples_LDP <- read.gdsn(index.gdsn(f_ldp, "sample.id"))
snpgdsClose(f_ldp)

# B. Breeding Matrix
f_breed <- snpgdsOpen(fn_Breed, allow.duplicate=TRUE)
geno_Breed <- snpgdsGetGeno(f_breed, snp.id=map_Breed_subset$ID, verbose=FALSE)
samples_Breed <- read.gdsn(index.gdsn(f_breed, "sample.id")) # or "sample.id" if standard
snpgdsClose(f_breed)

# 5. Combine Data
geno_Combined <- rbind(geno_LDP, geno_Breed)
all_samples <- c(samples_LDP, samples_Breed)

# 6. Create the GDS with REAL COORDINATES
cat("Creating new GDS with real physical positions...\n")

snpgdsCreateGeno(fn_Output, 
                 genmat = geno_Combined, 
                 sample.id = all_samples, 
                 snp.id = 1:n_common,          # Sequential IDs for the new file
                 snp.rs.id = common_keys,      # Store keys here
                 
                 # --- THIS IS THE FIX ---
                 snp.chromosome = map_LDP_subset$Chr,  # Real Chromosomes
                 snp.position = map_LDP_subset$Pos,    # Real Positions
                 
                 
                 snpfirstdim = FALSE)

cat(paste("Success! Saved as:", fn_Output, "\n"))

#Run the Gene Coverage Analysis
library(SNPRelate)
library(GenomicRanges)
library(rtracklayer)

# --- CONFIGURATION ---
gff_file <- "Lcu.1GRN.genes_description.sorted.gff3.gz"
gds_file <- "Merged_Analysis_RealCoords.gds" # Use the NEW file

# 1. Load Genes
cat("Reading GFF3 annotation...\n")
gff_data <- import.gff3(gff_file)
all_genes_gr <- gff_data[gff_data$type == "gene"]
total_genes_count <- length(all_genes_gr)

# 2. Load SNP Positions from GDS
cat("Reading SNP positions...\n")
f <- snpgdsOpen(gds_file)
snp_chr <- read.gdsn(index.gdsn(f, "snp.chromosome"))
snp_pos <- read.gdsn(index.gdsn(f, "snp.position"))
snpgdsClose(f)

# 3. Create GRanges for SNPs
# Important: Check if GFF and GDS use same chr naming (e.g. "Lcu.1GRN.Chr1")
# If one uses just "1" and other "Chr1", you get 0 overlaps.
# Let's inspect the first few to be safe:
cat("Sample GDS Chromosomes:", head(snp_chr), "\n")
cat("Sample GFF Chromosomes:", head(seqlevels(all_genes_gr)), "\n")

snp_gr <- GRanges(seqnames = snp_chr, 
                  ranges = IRanges(start = snp_pos, width = 1))

# 4. Find Overlaps
cat("Calculating overlaps...\n")
overlaps <- findOverlaps(query = all_genes_gr, subject = snp_gr)
covered_genes_count <- length(unique(queryHits(overlaps)))

# 5. Results
cat("\n==============================================\n")
cat(paste("Total Genes in Reference:       ", total_genes_count, "\n"))
cat(paste("Genes represented by Common SNPs:", covered_genes_count, "\n"))
cat(paste("Percentage of Total Genome:      ", round((covered_genes_count/total_genes_count)*100, 2), "%\n"))
cat("==============================================\n")

# Count how many times each gene index appears in the overlap list
snps_per_gene <- table(queryHits(overlaps))
cat("\n--- SNP Density per Gene ---\n")
cat(paste("Mean SNPs per covered gene:", round(mean(snps_per_gene), 1), "\n"))
cat(paste("Median SNPs per covered gene:", median(snps_per_gene), "\n"))

# Simple Histogram
hist(as.numeric(snps_per_gene), 
     breaks = 50, col = "orange", main = "Distribution of SNPs per Gene",
     xlab = "Number of SNPs", xlim = c(0, 50))

#Gene Coverage for Original Panels####
library(SNPRelate)
library(GenomicRanges)
library(rtracklayer)

# --- CONFIGURATION ---
gff_file     <- "Lcu.1GRN.genes_description.sorted.gff3.gz"
file_LDP     <- "LDP324_nofiltered.gds"
file_Breed   <- "ACT197_biallelic.gds" # Assuming this is your full/biallelic breeding file

# 1. Load Reference Genes (Once)
cat("Reading GFF3 annotation...\n")
gff_data <- import.gff3(gff_file)
all_genes_gr <- gff_data[gff_data$type == "gene"]
total_genes <- length(all_genes_gr)

cat(paste("Total Annotated Genes:", total_genes, "\n"))

# --- HELPER FUNCTION ---
calc_gene_coverage <- function(gds_path, label) {
  cat(paste0("\n--- Analyzing: ", label, " ---\n"))
  
  # Open File
  f <- snpgdsOpen(gds_path, allow.duplicate=TRUE)
  on.exit(snpgdsClose(f))
  
  # Smart Node Detection (SeqArray vs SNPRelate)
  if (!is.null(index.gdsn(f, "snp.position", silent=TRUE))) {
    # SNPRelate Format
    chr <- read.gdsn(index.gdsn(f, "snp.chromosome"))
    pos <- read.gdsn(index.gdsn(f, "snp.position"))
  } else {
    # SeqArray Format
    chr <- read.gdsn(index.gdsn(f, "chromosome"))
    pos <- read.gdsn(index.gdsn(f, "position"))
  }
  
  n_markers <- length(pos)
  cat(paste("Total Markers in file:", n_markers, "\n"))
  
  # Create GRanges
  snp_gr <- GRanges(seqnames = chr, 
                    ranges = IRanges(start = pos, width = 1))
  
  # Overlap
  cat("Calculating overlaps with genes...\n")
  hits <- findOverlaps(query = all_genes_gr, subject = snp_gr)
  
  # Count unique genes hit
  covered_genes <- length(unique(queryHits(hits)))
  pct <- (covered_genes / total_genes) * 100
  
  cat(paste("Genes Covered: ", covered_genes, "\n"))
  cat(paste("Coverage %:    ", round(pct, 2), "%\n"))
  
  # Optional: SNP Density in genes
  # Average SNPs per covered gene
  snps_per_gene <- table(queryHits(hits))
  cat(paste("Avg SNPs/Gene: ", round(mean(snps_per_gene), 1), "\n"))
  
  return(covered_genes)
}

# 2. Run Analysis
cov_LDP   <- calc_gene_coverage(file_LDP, "LDP (Exome Capture)")
cov_Breed <- calc_gene_coverage(file_Breed, "Breeding Panel (HiFi)")


#Filter ALS Genes####
library(SNPRelate)
library(GenomicRanges)
library(ggplot2)
library(gdsfmt)

# --- CONFIGURATION ---
input_gds  <- "ACT197_biallelic.gds"
output_gds <- "ACT197_ALS_Genes.gds"

# 1. Define the ALS Gene Regions manually
# I have cleaned the numbers (removed commas) from your list.
als_df <- data.frame(
  GeneID = c("Lcu.1GRN.1g006000", "Lcu.1GRN.1g009150", "Lcu.1GRN.1g018620",
             "Lcu.1GRN.2g003830", "Lcu.1GRN.2g031220", "Lcu.1GRN.3g036140",
             "Lcu.1GRN.3g053870", "Lcu.1GRN.4g029260", "Lcu.1GRN.5g068330",
             "Lcu.1GRN.6g009300", "Lcu.1GRN.7g062030"),
  Chr = c("Lcu.1GRN.Chr1", "Lcu.1GRN.Chr1", "Lcu.1GRN.Chr1",
          "Lcu.1GRN.Chr2", "Lcu.1GRN.Chr2", "Lcu.1GRN.Chr3",
          "Lcu.1GRN.Chr3", "Lcu.1GRN.Chr4", "Lcu.1GRN.Chr5",
          "Lcu.1GRN.Chr6", "Lcu.1GRN.Chr7"),
  Start = c(33155801, 61505633, 167559577, 
            7996094, 185989575, 286141769, 
            379219354, 245281691, 502872300, 
            73416912, 596103509),
  End = c(33156692, 61506547, 167560329, 
          8013970, 185990081, 286142257, 
          379221345, 245283665, 502872872, 
          73417475, 596104837)
)

# Convert genes to GRanges
genes_gr <- GRanges(seqnames = als_df$Chr,
                    ranges = IRanges(start = als_df$Start, end = als_df$End))

# 2. Open the Source GDS
f_in <- snpgdsOpen(input_gds)

# 3. Load SNP Coordinates
# Handle SeqArray vs SNPRelate node names if necessary
if (!is.null(index.gdsn(f_in, "snp.position", silent=TRUE))) {
  ids <- read.gdsn(index.gdsn(f_in, "snp.id"))
  chr <- read.gdsn(index.gdsn(f_in, "snp.chromosome"))
  pos <- read.gdsn(index.gdsn(f_in, "snp.position"))
} else {
  ids <- read.gdsn(index.gdsn(f_in, "variant.id"))
  chr <- read.gdsn(index.gdsn(f_in, "chromosome"))
  pos <- read.gdsn(index.gdsn(f_in, "position"))
}

# 4. Find SNPs within ALS Genes
snps_gr <- GRanges(seqnames = chr, ranges = IRanges(start = pos, width = 1))

cat("Finding overlaps between SNPs and ALS genes...\n")
overlaps <- findOverlaps(query = snps_gr, subject = genes_gr)
target_indices <- unique(queryHits(overlaps))
target_ids <- ids[target_indices]

cat(paste("Found", length(target_ids), "markers within the 11 ALS genes.\n"))

# 5. Create the Subset GDS
snpgdsClose(f_in) # Close before creating new file

if(length(target_ids) > 0) {
  snpgdsCreateGenoSet(src.fn = input_gds, 
                      dest.fn = output_gds, 
                      snp.id = target_ids, 
                      verbose = TRUE)
  cat(paste("Success! Saved:", output_gds, "\n"))
} else {
  stop("No SNPs found in these gene regions!")
}

# PART B: PCA Analysis on ALS Genes

# 1. Open the new ALS-specific GDS
f_als <- snpgdsOpen(output_gds)

# 2. Run PCA
# Note: With very few markers, we set eigen.cnt lower and allow monomorphic removal
pca_als <- snpgdsPCA(f_als, autosome.only=FALSE, remove.monosnp=TRUE)

# 3. Prepare Plotting Data
pc_percent <- pca_als$varprop * 100
pca_df <- data.frame(sample.id = pca_als$sample.id,
                     PC1 = pca_als$eigenvect[,1],
                     PC2 = pca_als$eigenvect[,2])

# 4. Plot
# If you have metadata (like 'Resistant' vs 'Susceptible'), merge it here!
ggplot(pca_df, aes(x=PC1, y=PC2)) +
  geom_point(color="purple", alpha=0.7, size=3) +
  labs(title = "Local PCA: Acetolactate Synthase (ALS) Genes",
       subtitle = paste("Structure based on", length(target_ids), "markers in 11 ALS loci"),
       x = paste0("PC1 (", round(pc_percent[1], 2), "%)"),
       y = paste0("PC2 (", round(pc_percent[2], 2), "%)")) +
  theme_bw()

snpgdsClose(f_als)

#grouping by ACT100####
library(ggplot2)
library(dplyr)

# --- CONFIGURATION ---
subset_file <- "ACT100_cot_color.txt"

# 1. Load the Subset List
# Assuming tab-delimited text. If CSV, use read.csv()
subset_info <- read.table(subset_file, header=TRUE, stringsAsFactors=FALSE)

# Extract the list of IDs
# Ensure the column name matches exactly what is in your file
target_ids <- subset_info$SampleID 

cat(paste("Loaded", length(target_ids), "IDs from the rotation trial list.\n"))

# 2. Update the PCA Data Frame
# We assume 'pca_df' already exists from your previous PCA run.
# If not, re-run the "pca_df <- data.frame(...)" block from the specific PCA you want to plot.

pca_df$Rotation_Group <- ifelse(pca_df$sample.id %in% target_ids, 
                                "Selected for Trial", 
                                "Not in Trial")

# Check overlap counts
print(table(pca_df$Rotation_Group))

# 3. Plot
# We use Grey for non-trial lines to make the Trial lines (Green) pop out.

# Retrieve variance percentages if available (or use generic labels)
if(exists("pca")) {
  pc_percent <- pca$varprop * 100
  xlab_txt <- paste0("PC1 (", round(pc_percent[1], 2), "%)")
  ylab_txt <- paste0("PC2 (", round(pc_percent[2], 2), "%)")
} else {
  xlab_txt <- "PC1"
  ylab_txt <- "PC2"
}

ggplot(pca_df, aes(x=PC1, y=PC2, color=Rotation_Group, alpha=Rotation_Group)) +
  
  # Plot points
  geom_point(size=3) +
  
  # Custom Colors: Grey for background, Green for your trial
  scale_color_manual(values = c("Not in Trial" = "purple", 
                                "Selected for Trial" = "forestgreen")) +
  
  # Custom Transparency: Make background lines fainter so trial lines stand out
  scale_alpha_manual(values = c("Not in Trial" = 0.4, 
                                "Selected for Trial" = 1.0)) +
  
  labs(title = "Genetic Representation of Rotation Trial",
       subtitle = "Comparison of selected subset (N=100) vs. full breeding panel",
       x = paste0("PC1 (", round(pc_percent[1], 2), "%)"),
       y = paste0("PC2 (", round(pc_percent[2], 2), "%)")) +
  
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank())

#XP-CLR approach####
library(SNPRelate)
library(ggplot2)
library(dplyr)

# --- CONFIGURATION ---
gds_file <- "Merged_Analysis.gds"
window_size_bp <- 50000   # 50kb window (Adjust as needed)
step_size_bp   <- 10000   # 10kb step
ld_threshold   <- 0.7     # Correlation cutoff (Chen et al. suggest high LD)

# 1. Open Data & Define Populations
f <- snpgdsOpen(gds_file)
sample_ids <- read.gdsn(index.gdsn(f, "sample.id"))

# Define Groups (Adjust patterns if needed)
# Object = Breeding (Selection); Reference = LDP (Diversity)
group_obj <- sample_ids[grepl("_ACT$", sample_ids)]
group_ref <- sample_ids[grepl("_AGL$", sample_ids)]

# 2. Calculate Per-SNP Fst (The "Differentiation" Signal)
cat("Calculating genome-wide Fst...\n")

# Load all SNP IDs first to be explicit
all_snps <- read.gdsn(index.gdsn(f, "snp.rs.id"))

fst_res <- snpgdsFst(f, 
                     sample.id = c(group_obj, group_ref),
                     # Explicitly pass the SNPs we want to analyze
                     population = factor(c(rep("Obj", length(group_obj)), 
                                           rep("Ref", length(group_ref)))),
                     method = "W&C84", # Weir & Cockerham 1984
                     autosome.only = FALSE,
                     verbose = FALSE)

cat("Parsing SNP IDs for map info...\n")
split_ids <- strsplit(all_snps, ":")
# Extract Raw Chromosome String (e.g., "Lcu.1GRN.Chr1")
raw_chr <- sapply(split_ids, `[`, 1)


# Extract Chromosome (1st part) and Position (2nd part)
# unlist to convert list to vector
chr_vec <- gsub("Lcu.1GRN.Chr", "", raw_chr)
pos_vec <- as.numeric(sapply(split_ids, `[`, 2))

# Extract into dataframe
snp_df <- data.frame(
  snp.id = all_snps,
  chr = chr_vec,
  pos = pos_vec,
  
  # IMPORTANT: Use $FstSNP for the per-site values
  fst = fst_res$FstSNP 
)

# Remove NAs (invariant sites or missing data)
snp_df <- snp_df[!is.na(snp_df$fst), ]

# Sanity Check: Fst should not be negative in this context (though W&C can technically be slightly neg)
# We often clamp negatives to 0 for selection scans
snp_df$fst[snp_df$fst < 0] <- 0

cat(paste("Valid Fst SNPs:", nrow(snp_df), "\n"))


# 3. Calculate "Chen et al." LD Weights
# The paper weights SNP i by 1 / (1 + sum of correlations > cutoff)
# We calculate this per chromosome to save memory.

cat("Calculating LD weights (Chen et al. 2010 method)...\n")

snp_df <- snp_df |>
  filter( chr %in% c("1","2","3","4","5","6","7"))

# Initialize weights column
snp_df$weight <- 1 

chroms <- unique(snp_df$chr)

for(c in chroms) {
  cat(paste("  Processing Chr:", c, "\n"))
  
  # Get SNPs for this chromosome
  chr_snps <- snp_df$snp.id[snp_df$chr == c]
  
  # Calculate LD Matrix (Correlation r)
  # sliding.window=500 keeps it fast (only looks at neighbors)
  ld_mat <- snpgdsLDMat(f,  slide=500, method="corr", verbose=FALSE)$LD
  
  # Thresholding: Convert to 1 if r > threshold, 0 otherwise
  # We subtract 1 because the diagonal (self-correlation) is always 1
  high_ld_counts <- colSums(abs(ld_mat) > ld_threshold, na.rm=TRUE) - 1
  
  # Calculate Weight: 1 / (1 + neighbors_in_LD)
  # Eq 7 in Chen et al. 2010
  weights <- 1 / (1 + high_ld_counts)
  
  # Store back in dataframe
  snp_df$weight[snp_df$snp.id %in% chr_snps] <- weights
}

# 4. Calculate Sliding Window Scores
cat("Calculating Sliding Window Scores...\n")

# Define Windows
# We'll create a grid of windows for the whole genome
make_windows <- function(df, size, step) {
  wins <- list()
  for(c in unique(df$chr)) {
    sub <- df[df$chr == c, ]
    starts <- seq(min(sub$pos), max(sub$pos), by=step)
    
    for(s in starts) {
      e <- s + size
      idx <- which(sub$pos >= s & sub$pos < e)
      
      # Only keep windows with sufficient data
      if(length(idx) > 5) { 
        
        # FIX 1: Use na.rm=TRUE in sums
        w_fst <- sum(sub$fst[idx] * sub$weight[idx], na.rm = TRUE)
        total_w <- sum(sub$weight[idx], na.rm = TRUE)
        
        # FIX 2: Avoid division by zero
        if(total_w > 0) {
          score <- w_fst / total_w
          wins[[length(wins)+1]] <- data.frame(chr=c, start=s, end=e, 
                                               n_snps=length(idx), score=score)
        }
      }
    }
  }
  return(do.call(rbind, wins))
}

# Run the windowing again
win_res <- make_windows(snp_df, window_size_bp, step_size_bp)
win_res <- win_res[is.finite(win_res$score), ]

# 5. Normalize Scores (Z-score)
# Chen et al. normalize by genome-wide mean/sd to make it comparable
global_mean <- mean(win_res$score, na.rm = TRUE)
global_sd   <- sd(win_res$score, na.rm = TRUE)

win_res$z_score <- (win_res$score - global_mean) / global_sd
head(win_res)

# 6. Plotting
# Define Top 1% Threshold
thresh <- quantile(win_res$z_score, 0.99, na.rm = T)

p <- ggplot(win_res, aes(x=start/1e6, y=z_score)) +
  geom_point(color="grey", alpha=0.5, size=0.8) +
  geom_point(data=subset(win_res, z_score > thresh), color="firebrick", size=1) +
  geom_hline(yintercept=thresh, linetype="dashed", color="blue") +
  facet_grid(.~chr, scales="free_x", space="free_x") +
  labs(title="XP-CLR: LD-Weighted Differentiation",
       subtitle=paste("Breeding vs Diversity | Window:", window_size_bp/1000, "kb"),
       x="Position (Mb)", 
       y="Weighted Differentiation (Z-Score)") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

print(p)
ggsave("XP-CLR LD-wheighted differenciation.png", height = 5, width = 10, bg= "white")
# Save Results
write.csv(win_res, "XPCLR_Emulation_Results.csv", row.names=FALSE)
cat("Done! Results saved to XPCLR_Emulation_Results.csv\n")
snpgdsClose(f)
