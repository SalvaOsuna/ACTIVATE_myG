library(gdsfmt)
library(GenomicRanges)
library(rtracklayer)

# --- 1. Configuration ---
# List of all your GDS files
gds_files <- c(
  "data/ACT197_SVs.gds",
  "data/ACT197_biallelic.gds",
  "data/ACT197_biallelic_PRUNED.gds",
  "data/ACT197_biallelic_codingonly.gds",
  "data/LDP324_nofiltered.gds",
  "data/Merged_Analysis_RealCoords.gds" 
)

# Reference Annotation for Gene Coverage
gff_fn <- "data/Lcu.1GRN.genes_description.sorted.gff3.gz"

# --- 2. Load Gene Models ---
cat("Loading GFF3 annotation...\n")
gff_data <- import.gff3(gff_fn)
gene_ranges <- gff_data[gff_data$type == "gene"]

# --- 3. Define the Summary Function ---
get_gds_stats <- function(file_path, gene_gr) {
  
  if(!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  cat(paste("Processing:", basename(file_path), "...\n"))
  
  # Open GDS agnostically (works for SeqArray and SNPRelate)
  f <- openfn.gds(file_path)
  
  # Find nodes (handles both naming conventions)
  node_id  <- index.gdsn(f, "snp.id", silent=TRUE)
  if (is.null(node_id)) node_id <- index.gdsn(f, "variant.id", silent=TRUE)
  
  node_chr <- index.gdsn(f, "snp.chromosome", silent=TRUE)
  if (is.null(node_chr)) node_chr <- index.gdsn(f, "chromosome", silent=TRUE)
  
  node_pos <- index.gdsn(f, "snp.position", silent=TRUE)
  if (is.null(node_pos)) node_pos <- index.gdsn(f, "position", silent=TRUE)
  
  # Extract data
  ids   <- read.gdsn(node_id)
  chrom <- read.gdsn(node_chr)
  pos   <- read.gdsn(node_pos)
  
  closefn.gds(f)
  
  # Total Variants
  total_vars <- length(ids)
  
  # Filter for the 7 main chromosomes 
  # (Matches "Lcu.1GRN.Chr1" to 7, or just "1" to "7")
  chrom_str <- as.character(chrom)
  chr_mask  <- grepl("Chr[1-7]$|^[1-7]$", chrom_str)
  
  vars_on_chr <- sum(chr_mask)
  avg_vars_chr <- vars_on_chr / 7
  
  # Subset to main chromosomes for spacing/density/coverage
  main_chroms <- chrom_str[chr_mask]
  main_pos    <- pos[chr_mask]
  
  # Calculate physical span of the data to estimate density
  # Span = sum of (max_position - min_position) per chromosome
  chr_spans <- tapply(main_pos, main_chroms, function(x) max(x, na.rm=TRUE) - min(x, na.rm=TRUE))
  total_span_bp <- sum(as.numeric(chr_spans), na.rm=TRUE)
  
  # Spacing and Density
  # Avoid division by zero if a file is empty
  if(vars_on_chr > 0) {
    avg_spacing_kb <- (total_span_bp / vars_on_chr) / 1000
    density_mb <- vars_on_chr / (total_span_bp / 1000000)
  } else {
    avg_spacing_kb <- NA; density_mb <- NA
  }
  
  # Gene Coverage
  # Ensure chromosome names match the GFF3 format exactly ("Lcu.1GRN.ChrX")
  mapped_chrs <- ifelse(grepl("^[1-7]$", main_chroms), 
                        paste0("Lcu.1GRN.Chr", main_chroms), 
                        main_chroms)
  
  variant_gr <- GRanges(seqnames = mapped_chrs,
                        ranges = IRanges(start = main_pos, width = 1))
  
  # Find variants inside genes
  overlaps <- findOverlaps(variant_gr, gene_gr)
  unique_genes_hit <- length(unique(subjectHits(overlaps)))
  
  gene_cov_pct <- (unique_genes_hit / length(gene_gr)) * 100
  
  # Return row
  return(data.frame(
    Dataset = basename(file_path),
    Total_Variants = total_vars,
    Variants_Chr1_7 = vars_on_chr,
    Avg_Vars_per_Chr = round(avg_vars_chr, 0),
    Avg_Spacing_kb = round(avg_spacing_kb, 2),
    Marker_Density_SNP_Mb = round(density_mb, 2),
    Gene_Coverage_Pct = round(gene_cov_pct, 2)
  ))
}

# --- 4. Run Analysis & Compile Table ---
results_list <- lapply(gds_files, get_gds_stats, gene_gr = gene_ranges)

# Combine into a single dataframe
final_table <- do.call(rbind, results_list)

# View the table
print(final_table)

# Save to CSV for easy inclusion in your manuscript
write.csv(final_table, "Results/Variant_Summary_Statistics.csv", row.names = FALSE)
cat("Done! Results saved to Variant_Summary_Statistics.csv\n")