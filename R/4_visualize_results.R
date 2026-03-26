# Script 4: Visualize Haplotype Block Analysis Results

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("SeqArray", quietly = TRUE)) BiocManager::install("SeqArray", update = FALSE)

library(ggplot2)
library(dplyr)
library(tidyr)
library(SeqArray)

base_dir <- "c:/Users/Salva/OneDrive - University of Saskatchewan/UsasK/github/ACTIVATE_myG"
hap_dir <- file.path(base_dir, "Results", "Haplotype_anti")

blocks <- read.csv(file.path(hap_dir, "All_Mapped_Haplotype_Blocks.csv"), stringsAsFactors = FALSE)

region_colors <- c(
  "Adaptation"   = "#1b9e77",
  "Breeding"     = "#d95f02",
  "Eroded"       = "darkorange",
  "Introgressed" = "dodgerblue",
  "Genome-wide"  = "grey80"
)

# Filter out very large blocks (e.g. > 99th percentile) strictly for visualization clarity in boxplots
limit_kb <- quantile(blocks$KB, 0.95, na.rm=TRUE)

cat("Generating Plot A: Block Size Comparisons...\n")
# Fig A: Block Size Comparisons (Boxplots)
pA <- ggplot(blocks, aes(x = Region_Type, y = KB, fill = Region_Type)) +
  geom_boxplot(outlier.shape = NA, alpha=0.7) +
  scale_y_continuous(limits = c(0, limit_kb)) +
  scale_fill_manual(values = region_colors) +
  facet_wrap(~ Population) +
  theme_minimal(base_size = 14) +
  labs(title = "Haplotype Block Size Comparison",
       subtitle = "ACTIVATE vs LDP across functional regions",
       x = "Region Type",
       y = "Block Size (kb)") +
  theme(plot.title = element_text(face="bold"),
        axis.text.x = element_text(angle=45, hjust=1))

ggsave(file.path(hap_dir, "FigA_BlockSize_Boxplots.png"), pA, width=10, height=5, bg="white")

cat("Generating Plot B: ECDF of Block Sizes...\n")
# Fig B: ECDF Plot of Block Sizes
# Compare ACTIVATE vs LDP within Introgressed, Eroded, Breeding, Adaptation
target_blocks <- blocks %>% filter(Region_Type %in% c("Introgressed", "Eroded", "Breeding", "Adaptation"))

pB <- ggplot(target_blocks, aes(x = KB, color = Region_Type)) +
  stat_ecdf(geom = "step", linewidth=1) +
  scale_x_log10() +
  facet_wrap(~ Population) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = region_colors) +
  labs(title = "Empirical Cumulative Distribution Function (ECDF)",
       subtitle = "Haplotype Block Sizes in Target Regions",
       x = "Block Size (kb) [log10 scale]",
       y = "Cumulative Probability") +
  theme(plot.title = element_text(face="bold"))

ggsave(file.path(hap_dir, "FigB_BlockSize_ECDF.png"), pB, width=10, height=5, bg="white")


cat("Generating Plot C: Proxy Haplotype Diversity (# of SNPs per block)...\n")
# Fig C: Complexity/Diversity proxy (Number of SNPs per block)
# Violin plot of NSNPS
pC <- ggplot(target_blocks, aes(x = Region_Type, y = NSNPS, fill = Region_Type)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width=0.2, outlier.shape=NA, color="black", alpha=0.5) +
  scale_y_continuous(limits = c(0, quantile(target_blocks$NSNPS, 0.95))) +
  scale_fill_manual(values = region_colors) +
  facet_wrap(~ Population) +
  theme_minimal(base_size = 14) +
  labs(title = "Block Complexity (SNPs per Block)",
       subtitle = "Proxy for haplotype diversity within functional regions",
       x = "Region Type",
       y = "Number of Variants (NSNPS)") +
  theme(plot.title = element_text(face="bold"),
        axis.text.x = element_text(angle=45, hjust=1))

ggsave(file.path(hap_dir, "FigC_Diversity_NSNPS.png"), pC, width=10, height=5, bg="white")


cat("Generating Plot D: Regional Zoom-in (Exact Haplotypes)...\n")
# Fig D: Exact Haplotype Frequencies for the largest block in Eroded & Introgressed
tryCatch({
  gds_path <- file.path(base_dir, "Results", "tmp_merged_seqarray.gds")
  if(!file.exists(gds_path)) {
     cat("SeqArray GDS not found. Re-converting from SNP GDS to do haplotype extraction...\n")
     orig_gds <- file.path(base_dir, "data", "Merged_Analysis_RealCoords.gds")
     seqSNP2GDS(orig_gds, gds_path, verbose=FALSE)
  }
  
  genofile <- seqOpen(gds_path)
  
  # Helper to get haplotype frequencies for a given block string and population
  get_hap_freqs <- function(snp_str, pop_samples) {
    if(is.na(snp_str) || nchar(snp_str) < 5) return(NULL)
    
    # Parse CHR:POS out of the string "chr:pos|chr:pos..."
    variant_ids <- strsplit(snp_str, "\\|")[[1]]
    chr_pos <- do.call(rbind, strsplit(variant_ids, ":"))
    chrs <- chr_pos[,1]
    pos  <- as.numeric(chr_pos[,2])
    
    seqSetFilter(genofile, sample.id=pop_samples, verbose=FALSE)
    
    # Find matching variants in GDS (using position and chromosome)
    gds_chr <- seqGetData(genofile, "chromosome")
    gds_pos <- seqGetData(genofile, "position")
    
    # We want exactly these positions on that chromosome
    sel_chr <- chrs[1] # all SNPs in a block share a chromosome
    match_idx <- which(gds_chr == sel_chr & (gds_pos %in% pos))
    
    if(length(match_idx) == 0) return(NULL)
    
    # Filter to specific SNPs
    seqSetFilter(genofile, variant.sel=match_idx, action="intersect", verbose=FALSE)
    
    # Extract genotypes
    # Format: variants x samples x ploidy
    geno <- seqGetData(genofile, "$dosage") # this is sample x variant matrix (0,1,2,NA)
    
    if(is.null(geno)) return(NULL)
    
    # Create haplotype string for each sample
    haps <- apply(geno, 1, function(row) {
      if(any(is.na(row))) return(NA)
      paste(row, collapse="")
    })
    
    haps <- haps[!is.na(haps)]
    if(length(haps) == 0) return(NULL)
    
    freq_df <- as.data.frame(table(haps))
    freq_df$Freq <- freq_df$Freq / sum(freq_df$Freq) * 100
    # Keep top 10 haplotypes to avoid clutter
    freq_df <- freq_df %>% arrange(desc(Freq)) %>% head(10)
    return(freq_df)
  }
  
  act_samples <- read.table(file.path(base_dir, "data", "ACT187_samples.txt"), stringsAsFactors=FALSE)$V1
  ldp_samples <- read.table(file.path(base_dir, "data", "LDP324_samples.txt"), stringsAsFactors=FALSE)$V1
  
  # Find largest eroded block in ACTIVATE
  top_eroded_blk <- target_blocks %>% 
                     filter(Region_Type == "Eroded", Population == "ACTIVATE") %>% 
                     arrange(desc(KB)) %>% head(1)
                     
  # Find largest introgressed block in ACTIVATE
  top_introg_blk <- target_blocks %>% 
                     filter(Region_Type == "Introgressed", Population == "ACTIVATE") %>% 
                     arrange(desc(KB)) %>% head(1)

  # Find largest breeding block in ACTIVATE
  top_breed_blk <- target_blocks %>% 
                     filter(Region_Type == "Breeding", Population == "ACTIVATE") %>% 
                     arrange(desc(KB)) %>% head(1)

  # Find largest adaptation block in ACTIVATE
  top_adapt_blk <- target_blocks %>% 
                     filter(Region_Type == "Adaptation", Population == "ACTIVATE") %>% 
                     arrange(desc(KB)) %>% head(1)
  
  plot_hap_freq <- function(block_row, title_label) {
     if(nrow(block_row) == 0) return()
     
     hap_act <- get_hap_freqs(block_row$SNPS, act_samples)
     hap_ldp <- get_hap_freqs(block_row$SNPS, ldp_samples)
     
     if(is.null(hap_act) && is.null(hap_ldp)) return()
     
     if(!is.null(hap_act)) hap_act$Population <- "ACTIVATE"
     if(!is.null(hap_ldp)) hap_ldp$Population <- "LDP"
     
     comb <- bind_rows(hap_act, hap_ldp)
     
     pH <- ggplot(comb, aes(x=reorder(haps, -Freq), y=Freq, fill=Population)) +
       geom_bar(stat="identity", position="dodge") +
       scale_fill_manual(values=c("ACTIVATE"="#E69F00", "LDP"="#56B4E9")) +
       theme_minimal() +
       theme(axis.text.x = element_text(angle=45, hjust=1)) +
       labs(title=title_label,
            subtitle=paste0(block_row$CHR, ": ", block_row$BP1, " - ", block_row$BP2, " (", block_row$KB, " kb)"),
            x="Haplotype Sequence (Top 10)",
            y="Frequency (%)")
     return(pH)
  }
  
  cat("Plotting top Eroded block...\n")
  pD1 <- plot_hap_freq(top_eroded_blk, "Haplotype Frequencies: Top Eroded Block")
  if(!is.null(pD1)) ggsave(file.path(hap_dir, "FigD1_Top_Eroded_Haps.png"), pD1, width=8, height=5, bg="white")
  
  cat("Plotting top Introgressed block...\n")
  pD2 <- plot_hap_freq(top_introg_blk, "Haplotype Frequencies: Top Introgressed Block")
  if(!is.null(pD2)) ggsave(file.path(hap_dir, "FigD2_Top_Introg_Haps.png"), pD2, width=8, height=5, bg="white")

  cat("Plotting top Breeding block...\n")
  pD3 <- plot_hap_freq(top_breed_blk, "Haplotype Frequencies: Top Breeding Block")
  if(!is.null(pD3)) ggsave(file.path(hap_dir, "FigD3_Top_Breeding_Haps.png"), pD3, width=8, height=5, bg="white")
  
  cat("Plotting top Adaptation block...\n")
  pD4 <- plot_hap_freq(top_adapt_blk, "Haplotype Frequencies: Top Adaptation Block")
  if(!is.null(pD4)) ggsave(file.path(hap_dir, "FigD4_Top_Adaptation_Haps.png"), pD4, width=8, height=5, bg="white")
  
  seqClose(genofile)
  
}, error = function(e) {
  cat("Warning: Could not complete exact haplotype extraction for Fig D.\n")
  cat(e$message, "\n")
})

cat("Visualization script completed successfully!\n")
