#required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("SeqArray", "SNPRelate", "vcfR", "ggplot2"))

#get SV####
library(SeqArray)
library(SNPRelate)
library(vcfR)

# --- 1. Convert VCF to GDS ---
vcf_fn <- "data/Lcu.1GRN-ACTIVATE-lightly filtered.vcf.gz"
gds_fn <- "data/ACTIVATE_SVs_only.gds"

# Parse the VCF into a SeqArray GDS file (this handles the INFO column parsing)
# (Only run this once, it takes a few minutes for a large VCF)
seqVCF2GDS(vcf_fn, gds_fn, verbose = TRUE)
genofile <- seqOpen(gds_fn)

# --- 2. ROBUST Identification of SVs ---
# A. Get total number of variants
all_var_ids <- seqGetData(genofile, "variant.id")
n_vars <- length(all_var_ids)

# B. Safely check for SVTYPE in the INFO node
info_nodes <- ls.gdsn(index.gdsn(genofile, "annotation/info", silent=TRUE))

if(!is.null(info_nodes) && "SVTYPE" %in% info_nodes) {
  sv_type <- seqGetData(genofile, "annotation/info/SVTYPE")
} else {
  sv_type <- rep(NA, n_vars) # Fill with NAs if missing
}

# C. Calculate sequence length differences
# FIX: Extract $ref and $alt separately
ref_alleles <- seqGetData(genofile, "$ref")
alt_alleles <- seqGetData(genofile, "$alt")

# Calculate length of REF
ref_len <- nchar(ref_alleles)

# Extract first ALT allele and calculate its length
first_alt <- sapply(strsplit(alt_alleles, ","), `[`, 1)
alt_len <- nchar(first_alt)

# Calculate the size of the indel
indel_length <- abs(ref_len - alt_len)

# D. Check for symbolic alleles (e.g., <DUP>, <INV>) 
is_symbolic <- grepl("<", first_alt)

# E. Combine logic: It's an SV if:
# 1. SVTYPE is present, OR
# 2. Length difference is >= 50bp, OR
# 3. It's a symbolic allele
sv_indices <- which(!is.na(sv_type) | indel_length >= 50 | is_symbolic)

sv_var_ids <- all_var_ids[sv_indices]

cat(paste("Identified", length(sv_var_ids), "Structural Variants out of", n_vars, "total variants.\n"))


# --- 3. Filter for Chromosomes 1 to 7 ---
sv_chroms <- seqGetData(genofile, "chromosome")[sv_indices]
target_chrs <- paste0("Lcu.1GRN.Chr", 1:7)

valid_sv_mask <- sv_chroms %in% target_chrs
final_sv_ids <- sv_var_ids[valid_sv_mask]

cat(paste("SVs mapped to Chr 1-7:", length(final_sv_ids), "\n"))

# Set the active filter
seqSetFilter(genofile, variant.id = final_sv_ids)

# --- 4. Create GM Object (Genotype Map) ---
raw_chr <- seqGetData(genofile, "chromosome")
pos <- seqGetData(genofile, "position")
clean_chr <- as.numeric(gsub("Lcu.1GRN.Chr", "", raw_chr))
sv_names <- paste0("SV_", clean_chr, "_", pos)

myGM <- data.frame(
  Name = sv_names,
  Chromosome = clean_chr,
  Position = pos,
  stringsAsFactors = FALSE
)


# --- 5. Create GD Object (Genotype Data) ---
# seqGetData returns alleles as a matrix of (Samples x Variants)
geno_matrix <- seqGetData(genofile, "$dosage_alt")

cat("Matrix dimensions (Samples x SVs):", dim(geno_matrix), "\n")

# Create the GD dataframe: First column MUST be taxa
# We use geno_matrix directly, no transposition needed!
myGD <- data.frame(taxa = seqGetData(genofile, "sample.id"), 
                   geno_matrix, 
                   stringsAsFactors = FALSE)

# Rename the marker columns in GD to match the Name column in GM
colnames(myGD)[-1] <- myGM$Name


# --- 6. Save ---
save(myGD, myGM, file = "data/ACTIVATE_GAPIT_SV_Input.RData")
cat("Success! GAPIT input files saved.\n")

# Or save as CSVs if you prefer inspecting them first
write.csv(myGM, "data/GAPIT_GM_SVs.csv", row.names = FALSE, quote = FALSE)
write.csv(myGD, "data/GAPIT_GD_SVs.csv", row.names = FALSE, quote = FALSE)

seqClose(genofile)

#impute####
# --- 1. Load Data ---
# Load the objects created in the previous step
load("data/ACTIVATE_GAPIT_SV_Input.RData")
# This loads 'myGD' and 'myGM' into your environment

# --- 2. Filter by Missing Rate ---
max_missing_rate <- 0.20 # 20% threshold

cat("Original GD dimensions:", dim(myGD), "\n")

# Calculate the proportion of NAs for each marker (ignoring column 1 'taxa')
missing_rates <- colMeans(is.na(myGD[, -1]))

# Identify markers that pass the filter
markers_to_keep <- names(missing_rates)[missing_rates <= max_missing_rate]

cat("Filtering out", sum(missing_rates > max_missing_rate), 
    "markers with >", max_missing_rate*100, "% missing data.\n")

# Subset GD to keep only 'taxa' and the valid markers
myGD_clean <- myGD[, c("taxa", markers_to_keep)]

# Subset GM so it perfectly matches the new GD
myGM_clean <- myGM[myGM$Name %in% markers_to_keep, ]


# --- 3. Impute Remaining NAs ---
cat("Imputing remaining missing values...\n")

# Choice A: Mean Imputation (Recommended for GAPIT/numeric formats)
# Replaces NA with the average dosage of that marker across all other lines.
impute_mean <- function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
}

# Choice B: Mode (Major Allele) Imputation
# Replaces NA with the most frequent genotype (0, 1, or 2) for that marker.
impute_mode <- function(x) {
  # Find the most common value
  ux <- unique(x[!is.na(x)])
  mode_val <- ux[which.max(tabulate(match(x, ux)))]
  x[is.na(x)] <- mode_val
  return(x)
}

# Apply the imputation across all marker columns
# (Using impute_mean here, swap to impute_mode if preferred)
myGD[, -1] <- lapply(myGD[, -1], impute_mean)

# Verify no NAs remain
remaining_nas <- sum(is.na(myGD[, -1]))
cat("Total NAs remaining:", remaining_nas, "\n")
cat("Cleaned GD dimensions:", dim(myGD), "\n")


# --- 4. Save Cleaned Objects ---
# Overwrite the old objects in memory with the clean ones
myGD_clean <- myGD
myGM_clean <- myGM

# Save to a new RData file ready for GAPIT
save(myGD_clean, myGM_clean, file = "data/ACTIVATE_GAPIT_SV_Cleaned.RData")
cat("Success! Cleaned data saved.\n")

#Export SVs to GDS and VCF####
# Make sure the genofile is still open and the seqSetFilter(genofile, variant.id = final_sv_ids) 
# from the previous script is still active!

# --- 1. Export to a new, clean GDS file ---
# This file will ONLY contain the SVs, making it much smaller and faster
out_gds <- "data/ACTIVATE_197_SVs.gds"
cat("Exporting pure SV dataset to GDS...\n")
seqExport(genofile, out_gds, optimize = TRUE, verbose = TRUE)


# --- 2. Export to a new compressed VCF file ---
# This writes the filtered SVs back into standard VCF text format
out_vcf <- "data/ACTIVATE_Final_SVs.vcf.gz"
cat("Exporting pure SV dataset to VCF...\n")

# Note: We keep all INFO and FORMAT tags intact by default
seqGDS2VCF(genofile, out_vcf, use_Rsamtools = TRUE, verbose = TRUE)

cat("Success! You now have standalone SV files.\n")

# You can now safely close the original working file
seqClose(genofile)

#Remove RILs####
# Load necessary libraries
library(gdsfmt)
library(SeqArray)
library(SNPRelate)

# --- 1. Define Inputs, Outputs, and Samples to Remove ---
gds_files <- c(
  "data/ACT197_SVs.gds",
  "data/ACT197_biallelic.gds",
  "data/ACT197_biallelic_PRUNED.gds",
  "data/ACT197_biallelic_codingonly.gds"
)

out_files <- c(
  "data/ACT187_SVs.gds",
  "data/ACT187_biallelic.gds",
  "data/ACT187_biallelic_PRUNED.gds",
  "data/ACT187_biallelic_codingonly.gds"
)

gen_to_remove <- c(
  "LR-86-10_ACT", "LR-86-11_ACT", "LR-86-12_ACT",
  "LR-86-18_ACT", "LR-86-2_ACT",  "LR-86-45_ACT",
  "LR-86-48_ACT", "LR-86-50_ACT", "LR-86-57_ACT",
  "LR-86-73_ACT"
)

# --- 2. Define the Smart Subsetting Function ---
filter_gds_samples <- function(in_file, out_file, remove_ids) {
  
  if(!file.exists(in_file)) {
    warning(paste("File not found:", in_file))
    return(NULL)
  }
  
  cat(paste("\nProcessing:", basename(in_file), "...\n"))
  
  # Check File Format
  f_check <- openfn.gds(in_file)
  file_format <- get.attr.gdsn(f_check$root)$FileFormat
  closefn.gds(f_check)
  
  if (!is.null(file_format) && file_format == "SEQ_ARRAY") {
    cat("  Format: SeqArray\n")
    
    # Open with SeqArray
    f <- seqOpen(in_file)
    
    # Get all samples and determine which to keep
    all_samples <- seqGetData(f, "sample.id")
    keep_samples <- setdiff(all_samples, remove_ids)
    
    cat(paste("  Original samples:", length(all_samples), "\n"))
    cat(paste("  Samples to keep:", length(keep_samples), "\n"))
    
    # Apply filter
    seqSetFilter(f, sample.id = keep_samples, verbose=FALSE) 
    
    # Export to new file
    cat(paste("  Exporting to:", basename(out_file), "\n"))
    seqExport(f, out_file, verbose=FALSE)
    
    seqClose(f)
    
  } else {
    cat("  Format: SNP GDS (SNPRelate)\n")
    
    # Open with SNPRelate
    f <- snpgdsOpen(in_file)
    
    # Get all samples and determine which to keep
    all_samples <- read.gdsn(index.gdsn(f, "sample.id"))
    keep_samples <- setdiff(all_samples, remove_ids)
    
    cat(paste("  Original samples:", length(all_samples), "\n"))
    cat(paste("  Samples to keep:", length(keep_samples), "\n"))
    
    # Create the new subsetted file
    cat(paste("  Exporting to:", basename(out_file), "\n"))
    snpgdsCreateGenoSet(in_file, out_file, sample.id = keep_samples, verbose=FALSE)
    
    snpgdsClose(f)
  }
  
  cat("  Done.\n")
}

# --- 3. Execute Subsetting ---
# Apply the function to all matching input/output pairs
mapply(filter_gds_samples, gds_files, out_files, MoreArgs = list(remove_ids = gen_to_remove))