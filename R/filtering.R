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

# --- 2. Identify SVs ---
sv_type <- seqGetData(genofile, "annotation/info/SVTYPE")
ref_alt <- seqGetData(genofile, "$ref,$alt")
indel_length <- abs(nchar(ref_alt$ref) - nchar(sapply(strsplit(ref_alt$alt, ","), `[`, 1)))

# Get indices and IDs of variants that are SVs
sv_indices <- which(!is.na(sv_type) | indel_length >= 50)
sv_var_ids <- seqGetData(genofile, "variant.id")[sv_indices]

# --- 3. Filter for Chromosomes 1 to 7 ---
# Get the chromosome assignments just for our SVs
sv_chroms <- seqGetData(genofile, "chromosome")[sv_indices]

# Define exactly what we want to keep
target_chrs <- paste0("Lcu.1GRN.Chr", 1:7)

# Keep only SVs that are on the target chromosomes
valid_sv_mask <- sv_chroms %in% target_chrs
final_sv_ids <- sv_var_ids[valid_sv_mask]

cat(paste("Total SVs found:", length(sv_var_ids), "\n"))
cat(paste("SVs mapped to Chr 1-7:", length(final_sv_ids), "\n"))

# Set the active filter to these specific variants
seqSetFilter(genofile, variant.id = final_sv_ids)

# --- 4. Create the GM Object (Genotype Map) ---
# Now, when we pull data, it ONLY pulls from Chr 1-7
raw_chr <- seqGetData(genofile, "chromosome")
pos <- seqGetData(genofile, "position")

# Because we already filtered out scaffolds, this clean-up will work perfectly
clean_chr <- as.numeric(gsub("Lcu.1GRN.Chr", "", raw_chr))

# Create unique IDs
sv_names <- paste0("SV_", clean_chr, "_", pos)

myGM <- data.frame(
  Name = sv_names,
  Chromosome = clean_chr,
  Position = pos,
  stringsAsFactors = FALSE
)

# Check to ensure no NAs were created
if(any(is.na(myGM$Chromosome))) {
  warning("Warning: Some chromosomes failed to parse as numeric!")
} else {
  cat("Chromosomes successfully cleaned (1 to 7).\n")
}

# --- 5. Create the GD Object (Genotype Data) ---
geno_matrix <- seqGetData(genofile, "$dosage_alt")
t_geno <- t(geno_matrix)

myGD <- data.frame(taxa = seqGetData(genofile, "sample.id"), t_geno, stringsAsFactors = FALSE)
colnames(myGD)[-1] <- myGM$Name

# --- 6. Save ---
save(myGD, myGM, file = "data/ACTIVATE_GAPIT_SV_Input.RData")

# Or save as CSVs if you prefer inspecting them first
write.csv(myGM, "data/GAPIT_GM_SVs.csv", row.names = FALSE, quote = FALSE)
write.csv(myGD, "data/GAPIT_GD_SVs.csv", row.names = FALSE, quote = FALSE)

seqClose(genofile)