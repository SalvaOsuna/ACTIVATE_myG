# Script 1: Subset GDS and Export to VCF for ACTIVATE and LDP populations
# This script loads the merged GDS file and subsets it by population, exporting to VCF

# Install missing packages if necessary
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
if (!requireNamespace("SNPRelate", quietly = TRUE))
  BiocManager::install("SNPRelate", update = FALSE)

library(SNPRelate)

# Close any previously opened GDS files (left open by previous script crashes)
if (!requireNamespace("gdsfmt", quietly = TRUE))
  BiocManager::install("gdsfmt", update = FALSE)
gdsfmt::showfile.gds(closeall=TRUE)

# Define paths
base_dir <- "c:/Users/Salva/OneDrive - University of Saskatchewan/UsasK/github/ACTIVATE_myG"
gds_file <- file.path(base_dir, "data/Merged_Analysis_RealCoords.gds")
act_samples_file <- file.path(base_dir, "data/ACT187_samples.txt")
ldp_samples_file <- file.path(base_dir, "data/LDP324_samples.txt")

# Convert the SNPGDS to a SeqArray GDS (SeqArray handles non-human chromosomes perfectly)
temp_seq_gds <- file.path(base_dir, "Results/tmp_merged_seqarray.gds")
if(!file.exists(temp_seq_gds)) {
  cat("Converting SNPGDS to SeqArray GDS format...\n")
  if (!requireNamespace("SeqArray", quietly = TRUE)) {
    BiocManager::install("SeqArray", update = FALSE)
  }
  library(SeqArray)
  seqSNP2GDS(gds_file, temp_seq_gds)
} else {
  library(SeqArray)
}

cat("Opening SeqArray GDS...\n")
seqfile <- seqOpen(temp_seq_gds)
seq_samples <- seqGetData(seqfile, "sample.id")

# Filter ACTIVATE samples strictly to those present in the GDS
act_samples_keep <- intersect(act_samples, seq_samples)
cat(sprintf("Exporting ACTIVATE subset: %d samples matched in GDS.\n", length(act_samples_keep)))
seqSetFilter(seqfile, sample.id = act_samples_keep)
seqGDS2VCF(seqfile, vcf.fn = "c:/Users/Salva/OneDrive - University of Saskatchewan/UsasK/github/ACTIVATE_myG/Results/activate_merged_snps.vcf")
system("gzip -f c:/Users/Salva/OneDrive - University of Saskatchewan/UsasK/github/ACTIVATE_myG/Results/activate_merged_snps.vcf")

# Filter LDP samples strictly to those present in the GDS
ldp_samples_keep <- intersect(ldp_samples, seq_samples)
cat(sprintf("Exporting LDP subset: %d samples matched in GDS.\n", length(ldp_samples_keep)))
# Reset filter first, then filter sequentially
seqResetFilter(seqfile)
seqSetFilter(seqfile, sample.id = ldp_samples_keep)
seqGDS2VCF(seqfile, vcf.fn = "c:/Users/Salva/OneDrive - University of Saskatchewan/UsasK/github/ACTIVATE_myG/Results/ldp_merged_snps.vcf")
system("gzip -f c:/Users/Salva/OneDrive - University of Saskatchewan/UsasK/github/ACTIVATE_myG/Results/ldp_merged_snps.vcf")

# Close the file and clean up
seqClose(seqfile)
# Uncomment to delete the temporary SeqArray file after execution
# file.remove(temp_seq_gds)
cat("Subsetting and export complete.\n")
