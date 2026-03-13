#genetic erosion####
# Load necessary libraries
library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(dplyr)
library(data.table)

# --- 1. Configuration ---
gds_file <- "data/Merged_Analysis_RealCoords.gds"
window_size <- 1000000  # 500 Kb windows (bins)
step_size <- 500000    # 100 Kb slide (Change to 500000 for non-overlapping rigid bins)

cat("Opening GDS file and defining populations...\n")
f <- snpgdsOpen(gds_file)

# Identify samples
samp_ids <- read.gdsn(index.gdsn(f, "sample.id"))
act_samples <- samp_ids[grepl("_ACT", samp_ids)]
act_samples <- act_samples[-c(188:197)]
ldp_samples <- setdiff(samp_ids, act_samples) 

# --- 2. Calculate SNP-level Diversity Metrics ---
cat("Extracting genotypes and calculating frequencies...\n")

# A function to calculate Ho, He, and PIC for a given population
calc_diversity_metrics <- function(gds_obj, target_samples) {
  
  # 1. Get Allele Frequencies (p)
  freq_info <- snpgdsSNPRateFreq(gds_obj, sample.id=target_samples, with.id=TRUE)
  p <- freq_info$AlleleFreq
  
  # Note: 2p(1-p) is symmetric, so it doesn't matter if p is the major or minor allele!
  He <- 2 * p * (1 - p)
  
  # PIC formula: 2p(1-p) - 2[2p(1-p)]^2  ==  He - 2*(He^2)
  PIC <- He - ((He^2) / 2)
  
  # 2. Get Observed Genotypes to calculate Ho
  # snpgdsGetGeno returns a matrix where 1 = heterozygous
  geno_matrix <- snpgdsGetGeno(gds_obj, sample.id=target_samples, verbose=FALSE)
  
  # Ho = (count of 1s) / (count of non-missing genotypes)
  Ho <- colSums(geno_matrix == 1, na.rm=TRUE) / colSums(!is.na(geno_matrix))
  
  # Return a dataframe of metrics
  data.frame(
    snp.id = freq_info$snp.id,
    p = p,
    Ho = Ho,
    He = He,
    PIC = PIC
  )
}

cat("Calculating metrics for ACTIVATE panel...\n")
metrics_act <- calc_diversity_metrics(f, act_samples)

cat("Calculating metrics for LDP diversity panel...\n")
metrics_ldp <- calc_diversity_metrics(f, ldp_samples)

# Combine with physical positions
df_snps <- data.frame(
  snp.id = read.gdsn(index.gdsn(f, "snp.id")),
  chr = read.gdsn(index.gdsn(f, "snp.chromosome")),
  pos = read.gdsn(index.gdsn(f, "snp.position")),
  stringsAsFactors = FALSE
) %>%
  left_join(metrics_act %>% rename_with(~paste0(., "_act"), -snp.id), by="snp.id") %>%
  left_join(metrics_ldp %>% rename_with(~paste0(., "_ldp"), -snp.id), by="snp.id")

snpgdsClose(f)

# --- 3. Sliding Window Genome Scan ---
cat("Running genome-wide 500 Kb sliding window scan...\n")

main_chrs <- unique(df_snps$chr)[grepl("Chr[1-7]$|^[1-7]$", unique(df_snps$chr))]
all_windows <- list()

for(c in main_chrs) {
  cat(paste("  Scanning", c, "...\n"))
  df_c <- df_snps %>% filter(chr == c)
  
  starts <- seq(1, max(df_c$pos, na.rm=TRUE), by=step_size)
  chr_results <- list()
  
  for(s in starts) {
    e <- s + window_size
    win_snps <- df_c %>% filter(pos >= s & pos <= e)
    
    # Require at least 5 SNPs in a window to calculate a reliable average
    if(nrow(win_snps) < 5) next
    
    chr_results[[length(chr_results)+1]] <- data.frame(
      chr = c,
      start = s,
      stop = e,
      mid_pos_mb = (s + window_size/2) / 1e6,
      n_snps = nrow(win_snps),
      # Average the metrics across the window
      avg_Ho_act = mean(win_snps$Ho_act, na.rm=TRUE),
      avg_He_act = mean(win_snps$He_act, na.rm=TRUE),
      avg_PIC_act = mean(win_snps$PIC_act, na.rm=TRUE),
      avg_Ho_ldp = mean(win_snps$Ho_ldp, na.rm=TRUE),
      avg_He_ldp = mean(win_snps$He_ldp, na.rm=TRUE),
      avg_PIC_ldp = mean(win_snps$PIC_ldp, na.rm=TRUE)
    )
  }
  all_windows[[c]] <- do.call(rbind, chr_results)
}

clean_df <- do.call(rbind, all_windows)

# --- 4. Calculate Delta He and Two-Tailed Thresholds ---
cat("Calculating Reduction of Diversity (Delta He) and thresholds...\n")
clean_df <- clean_df %>%
  mutate(
    chr_num = as.numeric(gsub(".*Chr", "", chr)),
    # Reduction of Diversity: He(LDP) - He(ACTIVATE)
    delta_He = avg_He_ldp - avg_He_act
  )

# Define top 1% (Erosion) and bottom 1% (Introgression) thresholds
threshold_erosion <- quantile(clean_df$delta_He, 0.99, na.rm=TRUE)
threshold_introgression <- quantile(clean_df$delta_He, 0.01, na.rm=TRUE)
# --- 5. Export Metrics ---
if(!dir.exists("Results")) dir.create("Results")

# Export Global Averages
global_metrics <- data.frame(
  Panel = c("ACTIVATE", "LDP"),
  Global_Ho = c(mean(df_snps$Ho_act, na.rm=TRUE), mean(df_snps$Ho_ldp, na.rm=TRUE)),
  Global_He = c(mean(df_snps$He_act, na.rm=TRUE), mean(df_snps$He_ldp, na.rm=TRUE)),
  Global_PIC = c(mean(df_snps$PIC_act, na.rm=TRUE), mean(df_snps$PIC_ldp, na.rm=TRUE))
)
write.csv(global_metrics, "Results/Global_Diversity_Metrics.csv", row.names = FALSE)

# Export the top 1% (Genetic Erosion / Sweeps)
significant_erosion <- clean_df %>% 
  filter(delta_He >= threshold_erosion) %>% 
  arrange(desc(delta_He))
write.csv(significant_erosion, "Results/Significant_Erosion_Sweeps_DeltaHe.csv", row.names = FALSE)

# Export the bottom 1% (Introgression / Diversification)
significant_introgression <- clean_df %>% 
  filter(delta_He <= threshold_introgression) %>% 
  arrange(delta_He) # Sorts from most negative upward
write.csv(significant_introgression, "Results/Significant_Introgression_DeltaHe.csv", row.names = FALSE)

# Group the SNP dataframe by chromosome and calculate the means
chr_metrics <- df_snps %>%
  # Filter out any unplaced scaffolds if necessary, keeping only main chromosomes
  filter(grepl("Chr[1-7]$|^[1-7]$", chr)) %>%
  group_by(chr) %>%
  summarize(
    ACTIVATE_Ho = mean(Ho_act, na.rm = TRUE),
    ACTIVATE_He = mean(He_act, na.rm = TRUE),
    ACTIVATE_PIC = mean(PIC_act, na.rm = TRUE),
    LDP_Ho = mean(Ho_ldp, na.rm = TRUE),
    LDP_He = mean(He_ldp, na.rm = TRUE),
    LDP_PIC = mean(PIC_ldp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Sort them nicely by chromosome number
  mutate(chr_num = as.numeric(gsub(".*Chr", "", chr))) %>%
  arrange(chr_num) %>%
  select(-chr_num) # Remove the sorting column for a cleaner table

# Save to a CSV for your Supplementary Tables
write.csv(chr_metrics, "Results/Chromosome_Diversity_Metrics.csv", row.names = FALSE)
cat("Chromosome metrics saved to 'Results/Chromosome_Diversity_Metrics.csv'!\n")

# --- 6. Plotting the Two-Tailed Manhattan Plot ---
cat("Generating Two-Tailed Delta He Manhattan Plot...\n")
clean_df$color_group <- ifelse(clean_df$chr_num %% 2 == 0, "Even", "Odd")

p_delta <- ggplot(clean_df, aes(x = mid_pos_mb, y = delta_He)) +
  # Background chromosomes
  geom_point(aes(color = color_group), alpha = 0.6, size = 1.2) +
  
  # Highlight Top 1% (Erosion) in Dark Orange
  geom_point(data = filter(clean_df, delta_He >= threshold_erosion), color = "darkorange", size = 1.6) +
  
  # Highlight Bottom 1% (Introgression) in Dodger Blue
  geom_point(data = filter(clean_df, delta_He <= threshold_introgression), color = "dodgerblue", size = 1.6) +
  
  # Add both threshold lines
  geom_hline(yintercept = threshold_erosion, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkred", linewidth = 0.5) +
  geom_hline(yintercept = threshold_introgression, linetype = "dashed", color = "black", linewidth = 0.8) +
  
  # Add text labels for the exact threshold values on the Y-axis
  #annotate("text", x = min(clean_df$mid_pos_mb), y = threshold_erosion, 
  #         label = round(threshold_erosion, 3), vjust = -0.5, hjust = 0, color = "darkorange", fontface = "bold") +
  #annotate("text", x = min(clean_df$mid_pos_mb), y = threshold_introgression, 
  #         label = round(threshold_introgression, 3), vjust = 1.5, hjust = 0, color = "dodgerblue", fontface = "bold") +
  
  facet_grid(. ~ chr_num, scales = "free_x", space = "free_x", switch = "x") +
  scale_color_manual(values = c("Even" = "gray50", "Odd" = "gray75")) +
  labs(
  #  title = "Genome-Wide Scan for Genetic Erosion and Introgression",
  #  subtitle = "Reduction of Diversity (\u0394He) highlighting extreme bottlenecks (top 1%) and introgressions (bottom 1%)",
    x = "Chromosome",
    y = expression(paste(Delta, "He (LDP - ACTIVATE)"))
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text.x = element_text(face = "bold", size = 10)
  )

ggsave("Results/Delta_He_TwoTailed_Plot_1Percent.png", plot = p_delta, width = 11, height = 5, dpi = 600)

cat("Done! 1% Two-tailed diversity scan complete. Check the 'Results' folder.\n")

#Extracting Proteins from Delta He Sweeps####
library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(dplyr)
library(stringr)

cat("Loading Lentil Annotation and Proteome...\n")

# --- 1. Load the Genome Data ---
# Update these paths to your actual lentil GFF3 and Proteome FASTA files
gff_file <- "data/Lcu.1GRN.genes_description.sorted.gff3.gz" 
proteome_file <- "data/Lcu.1GRN.pep.fasta"

annotation <- import(gff_file)
genes <- annotation[annotation$type == "gene"]
proteome <- readAAStringSet(proteome_file)

# --- 2. Define a Function to Extract and Chunk Proteins ---
extract_and_chunk <- function(sweep_file, scenario_name, chunk_size = 100) {
  cat(sprintf("\nProcessing %s...\n", scenario_name))
  
  if(!file.exists(sweep_file)) {
    warning(paste("File not found:", sweep_file))
    return(NULL)
  }
  
  # Load the Delta He sweep regions
  sweeps_df <- read.csv(sweep_file, stringsAsFactors = FALSE)
  
  if(nrow(sweeps_df) == 0) {
    cat("No significant regions found in this file. Skipping.\n")
    return(NULL)
  }
  
  # Create a GRanges object for the sweeps
  sweeps_gr <- GRanges(
    seqnames = sweeps_df$chr,
    ranges = IRanges(start = sweeps_df$start, end = sweeps_df$stop)
  )
  
  # Find overlapping genes
  overlaps <- findOverlaps(genes, sweeps_gr)
  candidate_genes <- genes[queryHits(overlaps)]
  
  # Extract Gene IDs
  gene_ids <- candidate_genes$ID
  gene_ids <- str_remove_all(gene_ids, "^gene:") # Clean the ID prefix
  
  cat(sprintf("Found %d candidate genes in the sweep regions.\n", length(unique(gene_ids))))
  
  # Filter the proteome: Match gene IDs and keep ONLY the .1 primary transcript
  # We construct a regex pattern to match "GeneID.1"
  search_patterns <- paste0(gene_ids, "\\.1\\b") 
  
  # Find matching sequences in the proteome
  matched_indices <- unlist(lapply(search_patterns, function(pat) grep(pat, names(proteome))))
  candidate_peps <- proteome[unique(matched_indices)]
  
  num_seqs <- length(candidate_peps)
  cat(sprintf("Extracted %d primary transcript (.1) protein sequences.\n", num_seqs))
  
  # Save the full candidate gene list for your records
  out_csv <- sprintf("Results/%s_Candidate_Genes.csv", scenario_name)
  write.csv(data.frame(Gene_ID = gene_ids), out_csv, row.names = FALSE)
  
  # Chunk and save for InterPro
  chunk_dir <- sprintf("Results/%s_InterPro_Chunks", scenario_name)
  if (!dir.exists(chunk_dir)) dir.create(chunk_dir, recursive = TRUE)
  
  if(num_seqs > 0) {
    chunk_indices <- ceiling(seq_len(num_seqs) / chunk_size)
    pep_chunks <- split(candidate_peps, chunk_indices)
    
    for (i in seq_along(pep_chunks)) {
      chunk_filename <- file.path(chunk_dir, sprintf("chunk_%02d.fasta", i))
      writeXStringSet(pep_chunks[[i]], chunk_filename)
    }
    cat(sprintf("Saved %d chunks to %s\n", length(pep_chunks), chunk_dir))
  }
}

# --- 3. Run the Extraction for Both Datasets ---

# Run for Genetic Erosion (Top 1% Sweeps)
extract_and_chunk(
  sweep_file = "Results/Significant_Erosion_Sweeps_DeltaHe.csv", 
  scenario_name = "DeltaHe_Erosion"
)

# Run for Introgression (Bottom 1% Diversification)
# (Make sure this filename matches exactly what the previous script output)
extract_and_chunk(
  sweep_file = "Results/Significant_Introgression_DeltaHe.csv", 
  scenario_name = "DeltaHe_Introgression"
)

cat("\nAll extractions complete! Ready for InterPro upload.\n")

#Candidate genes####
# --- 1. Load Data ---
cat("Loading GFF3 annotation and Significant regions...\n")

# Load the Lentil Annotation File
gff_file <- "data/Lcu.1GRN.genes_description.sorted.gff3.gz"
gff <- import(gff_file)

# Keep only "gene" or "mRNA" features to avoid counting exons multiple times
genes <- gff[gff$type %in% c("gene", "mRNA")]

# Load the Significant Regions from your target scenario
target_scenario <- "DeltaHe_Erosion" # Change to DeltaHe_Introgression or DeltaHe_Erosion for the other run
sig_regions <- read.csv(sprintf("Results/%s_SignificantRegions.csv", target_scenario))

# Convert to GenomicRanges
sig_ranges <- GRanges(
  seqnames = sig_regions$chr,
  ranges = IRanges(start = sig_regions$start, end = sig_regions$stop),
  z_score = sig_regions$z_score
)

# --- 2. Intersect Windows with Genes ---
cat("Finding candidate genes within selected regions...\n")

overlaps <- findOverlaps(sig_ranges, genes)
candidate_genes <- genes[subjectHits(overlaps)]
candidate_genes$z_score <- sig_ranges$z_score[queryHits(overlaps)]

# Use any_of() to safely select columns regardless of the specific GFF3 format
candidates_df <- as.data.frame(candidate_genes) %>%
  dplyr::select(any_of(c(
    "seqnames", "start", "end", "ID", "Name", 
    "Description", "Note", "Ontology_term", "Dbxref", "z_score"
  ))) %>% 
  distinct(ID, .keep_all = TRUE) 

candidates_df <- candidates_df %>%
  mutate(
    # Use sapply to look at every single row in the Description column
    Description = sapply(Description, function(x) {
      # If the list is empty (character(0)), return NA
      if (length(x) == 0) {
        return(NA_character_)
      } else {
        # If it has text, combine it and strip out the backslashes and quotes
        clean_text <- paste(x, collapse = "; ")
        clean_text <- gsub("\"", "", clean_text)
        return(clean_text)
      }
    })
  )

# Now it is safe to save!
write.csv(candidates_df, sprintf("Results/%s_Candidate_Genes.csv", target_scenario), row.names = FALSE)
cat(sprintf("Found %d candidate genes! Saved to Results.\n", nrow(candidates_df)))

#generate and plot GO####
library(dplyr)
library(tidyr)
library(ggplot2)
library(GO.db)
library(AnnotationDbi)

cat("Loading InterPro annotations and candidate genes...\n")

# 1. Load your candidate genes
target_scenario <- "DeltaHe_Erosion" # Change to DeltaHe_Introgression or DeltaHe_Erosion for the other run
candidates_df <- read.csv(sprintf("Results/%s_Candidate_Genes.csv", target_scenario))

# Clean the IDs in your candidates dataframe so they match the InterPro output exactly
candidates_df <- candidates_df %>%
  mutate(clean_id = gsub("^(gene:|mRNA:|transcript:)", "", ID))

cat("Finding and stitching InterPro TSV chunks...\n")

# --- 2. Load and Stitch the Headerless InterPro TSVs ---

# Define the folder where you saved all those downloaded TSV chunks
# (Update this path if you saved them somewhere else!)
tsv_folder <- "Results/DeltaHe_Erosion_InterPro_Chunks" 

# Search the folder for all files matching the chunk naming pattern
file_pattern <- paste0("^iprscan5-", target_scenario, ".*\\.tsv$")
tsv_files <- list.files(path = tsv_folder, pattern = file_pattern, full.names = TRUE)

if(length(tsv_files) == 0) {
  stop("No TSV files found! Check your folder path and target_scenario name.")
}

cat(sprintf("Found %d InterPro TSV files. Combining them now...\n", length(tsv_files)))

# Read all files into a list, then bind them into one large dataframe
interpro_list <- lapply(tsv_files, function(file) {
  read.delim(file, header = FALSE, stringsAsFactors = FALSE)
})
interpro_df <- bind_rows(interpro_list)

# Assign the proper standard InterPro column names
colnames(interpro_df) <- c("protein_acc", "md5", "length", "analysis", 
                           "sig_acc", "sig_desc", "start", "stop", "evalue", 
                           "status", "date", "ipr_acc", "ipr_desc", 
                           "go_terms", "pathways")

# --- 3. Extract and Clean the GO terms ---
cat("Extracting and cleaning GO terms...\n")

go_mapping <- interpro_df %>%
  dplyr::select(protein_acc, go_terms) %>%
  dplyr::rename(clean_id = protein_acc) %>%
  # Remove rows where GO terms are missing or just a hyphen
  filter(go_terms != "" & !is.na(go_terms) & go_terms != "-") %>%
  # Split multiple GO terms on the pipe character
  separate_rows(go_terms, sep = "\\|") %>%
  # Strip out the database tags in parentheses (e.g., "(InterPro)" or "(PANTHER)")
  mutate(go_terms = gsub("\\(.*\\)", "", go_terms)) %>%
  # Remove whitespace
  mutate(go_terms = trimws(go_terms)) %>%
  # Snip off the transcript suffix (e.g., turning .1 into just the gene ID)
  mutate(clean_id = sub("\\.[0-9]+$", "", clean_id)) %>%
  # Keep only unique combinations of gene + GO term
  distinct(clean_id, go_terms)

# 4. Merge with your candidate genes
final_go_df <- candidates_df %>%
  inner_join(go_mapping, by = "clean_id")

cat(sprintf("Successfully mapped GO terms to %d candidate genes!\n", length(unique(final_go_df$clean_id))))

# 5. Translate GO IDs to readable biological terms
cat("Translating GO IDs...\n")

mapped_terms <- suppressMessages(
  AnnotationDbi::select(GO.db, keys = unique(final_go_df$go_terms), columns=c("TERM","ONTOLOGY"), keytype = "GOID")
)

final_go_df <- final_go_df %>%
  left_join(mapped_terms, by = c("go_terms" = "GOID")) %>%
  # Keep only valid BP, MF, CC ontologies (filters out missing data)
  filter(ONTOLOGY %in% c("BP", "MF", "CC")) %>%
  mutate(
    term_name = ifelse(is.na(TERM), go_terms, TERM),
    # Expand the acronyms for beautiful facet labels
    ontology_full = case_when(
      ONTOLOGY == "BP" ~ "Biological Process",
      ONTOLOGY == "MF" ~ "Molecular Function",
      ONTOLOGY == "CC" ~ "Cellular Component"
    )
  )
write.csv(final_go_df, "Results/DeltaHe_Erosion_Final_GO.csv", row.names = FALSE)

# 6. Count and rank WITHIN each Ontology category
go_counts <- final_go_df %>%
  count(ontology_full, term_name, name = "gene_count") %>%
  group_by(ontology_full) %>%
  # Grab the top 10 terms per category to keep the plot balanced
  slice_max(order_by = gene_count, n = 10, with_ties = FALSE) %>%
  ungroup()

# 7. Generate the Faceted Lollipop Plot
cat("Generating Faceted GO Term plot...\n")

plot_title <- gsub("_", " ", target_scenario)

p_go_faceted <- ggplot(go_counts , aes(x = reorder(term_name, gene_count), y = gene_count, color = ontology_full)) +
  # The Lollipop sticks
  geom_segment(aes(xend = term_name, yend = 0), linewidth = 1) +
  # The Lollipop heads
  geom_point(size = 4, alpha = 0.9) +
  coord_flip() + 
  # Facet by Ontology, letting the Y-axis (which is flipped to X) scale independently
  # space = "free_y" ensures categories with fewer than 10 terms don't look awkwardly stretched
  facet_grid(ontology_full ~ ., scales = "free_y", space = "free_y") +
  # Assign distinct colors to the three categories
  scale_color_manual(values = c("Biological Process" = "#1b9e77", 
                                "Molecular Function" = "#d95f02", 
                                "Cellular Component" = "#7570b3")) +
  labs(
    title = paste("Functional Classification of Selected Genes:", plot_title),
    subtitle = "Top GO terms faceted by Ontology",
    x = "Gene Ontology (GO) Term",
    y = "Number of Selected Genes"
  ) +
  theme_bw() +
  theme(
    legend.position = "none", # Hide legend since the facet strips act as labels
    axis.text.y = element_text(size = 11, color = "black"),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text.y = element_text(size = 11, face = "bold", angle = 270), # Rotate text for readability
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold")
  )
# Save the plot
ggsave(sprintf("Results/%s_InterPro_GO_Faceted.png", target_scenario), plot = p_go_faceted, width = 12, height = 9, dpi = 600)
cat("Done! Beautiful GO Term plot saved.\n")

#Venn diagram####
library(dplyr)
library(ggplot2)
library(ggVennDiagram)

cat("Preparing Biological Process lists for Venn Diagram...\n")

# --- 1. Load your final GO dataframes ---
# (Assuming you saved your final_go_df for each scenario, or you have them loaded in your environment)
# For example, you might have saved them like this in your previous scripts:
# write.csv(final_go_df, "Results/Py_Scenario1_Adaptation_Final_GO.csv", row.names = FALSE)

df_eros <- read.csv("Results/DeltaHe_Erosion_Final_GO.csv", stringsAsFactors = FALSE)
df_intro <- read.csv("Results/DeltaHe_Introgression_Final_GO.csv", stringsAsFactors = FALSE)

# --- 2. Filter for Biological Processes and extract unique terms ---
# We use unique() because we just want to know IF a pathway was targeted, 
# not how many times it was targeted.
bp_eros <- df_eros %>%
  filter(ONTOLOGY == "BP") %>%
  pull(term_name) %>%
  unique()

bp_intro <- df_intro %>%
  filter(ONTOLOGY == "BP") %>%
  pull(term_name) %>%
  unique()

# --- 3. Create the List Object for the Venn Diagram ---
venn_list <- list(
  "Erosion" = bp_eros,
  "Introgression" = bp_intro
)

# --- 4. Plot the Venn Diagram ---
cat("Generating plot...\n")

p_venn <- ggVennDiagram(venn_list, 
                        label_alpha = 0, # Removes the ugly white box behind the numbers
                        category.names = c("Erosion", "Introgression"),
                        set_color = "black", # Outline of the circles
                        set_size = 4) +      # Size of the title text
  # Use a clean color gradient (e.g., from light grey to a deep blue/green)
  scale_fill_gradient(low = "#F4FAFE", high = "#2B83BA") +
  labs(
    title = "Overlap of Selected Biological Processes",
    subtitle = "Shared vs. Unique pathways driven by evolutionary pressure",
    fill = "Number of\nGO Terms"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    legend.position = "right"
  )

# Save the plot
ggsave("Results/Erosion_introg_Venn_Diagram.png", plot = p_venn, width = 8, height = 8, dpi = 600, bg = "white")
cat("Venn diagram saved to Results folder!\n\n")

# --- 5. Extract the specific terms for your manuscript ---
cat("--- Extracted Pathways for Manuscript Results ---\n")

# Shared by both
shared_bp <- intersect(bp_eros, bp_intro)
cat(sprintf("Number of Shared BP Terms: %d\n", length(shared_bp)))

# Unique to Adaptation
unique_adapt <- setdiff(bp_eros, bp_intro)
cat(sprintf("Number of Unique Adaptation BP Terms: %d\n", length(unique_adapt)))

# Unique to Breeding
unique_breed <- setdiff(bp_intro, bp_eros)
cat(sprintf("Number of Unique Breeding BP Terms: %d\n", length(unique_breed)))

# Save these lists so you can read them easily
write.csv(data.frame(Term = shared_bp), "Results/Venn_Shared_Erosion_introg.csv", row.names = FALSE)
write.csv(data.frame(Term = unique_adapt), "Results/Venn_Unique_Erosion_BP.csv", row.names = FALSE)
write.csv(data.frame(Term = unique_breed), "Results/Venn_Unique_Introgression_BP.csv", row.names = FALSE)

cat("CSV files containing the exact overlapping and unique terms have been saved.\n")

#lolipop plot with just unique terms####
library(dplyr)
library(ggplot2)

cat("Categorizing GO Terms across scenarios...\n")

# --- 1. Load your final GO dataframes ---
# (Make sure these point to the dataframes you generated in the GO plotting step)
df_eros <- read.csv("Results/DeltaHe_Erosion_Final_GO.csv", stringsAsFactors = FALSE)
df_intro <- read.csv("Results/DeltaHe_Introgression_Final_GO.csv", stringsAsFactors = FALSE)

# --- 2. Identify Unique and Shared Terms ---
terms_eros <- unique(df_eros$term_name)
terms_intro <- unique(df_intro$term_name)

shared_terms <- intersect(terms_eros, terms_intro)
unique_adapt <- setdiff(terms_eros, terms_intro)
unique_breed <- setdiff(terms_intro, terms_eros)

# --- 3. Create the Master Comparison CSV ---
# Build a simple dataframe tagging each term
presence_df <- bind_rows(
  data.frame(term_name = unique_adapt, Presence = "Erosion"),
  data.frame(term_name = unique_breed, Presence = "Introgression"),
  data.frame(term_name = shared_terms, Presence = "Both Scenarios")
)

# Grab the Ontology categories to make the CSV more informative
ontology_map <- bind_rows(
  df_eros %>% dplyr::select(term_name, ONTOLOGY, ontology_full),
  df_intro %>% dplyr::select(term_name, ONTOLOGY, ontology_full)
) %>% distinct(term_name, .keep_all = TRUE)

master_csv <- presence_df %>%
  left_join(ontology_map, by = "term_name") %>%
  arrange(Presence, ONTOLOGY, term_name)

write.csv(master_csv, "Results/GO_Terms_Intro_Eros_Comparison_Master.csv", row.names = FALSE)
cat("Master comparison CSV saved to 'Results/GO_Terms_Intro_Eros_Comparison_Master.csv'!\n\n")

# --- 4. Define a Function to Plot Unique Lollipops ---
plot_unique_lollipop <- function(df, target_terms, title_label, filename_out) {
  
  # Filter the dataframe to ONLY include terms from our target list
  unique_df <- df %>%
    filter(term_name %in% target_terms) %>%
    filter(ONTOLOGY %in% c("BP", "MF", "CC"))
  
  # Count and rank the top 10 UNIQUE terms per ontology
  go_counts <- unique_df %>%
    count(ontology_full, term_name, name = "gene_count") %>%
    group_by(ontology_full) %>%
    slice_max(order_by = gene_count, n = 10, with_ties = FALSE) %>%
    ungroup()
  
  # Generate the plot
  p <- ggplot(go_counts, aes(x = reorder(term_name, gene_count), y = gene_count, color = ontology_full)) +
    geom_segment(aes(xend = term_name, yend = 0), linewidth = 1) +
    geom_point(size = 4, alpha = 0.9) +
    coord_flip() + 
    facet_grid(ontology_full ~ ., scales = "free_y", space = "free_y") +
    scale_color_manual(values = c("Biological Process" = "#1b9e77", 
                                  "Molecular Function" = "#d95f02", 
                                  "Cellular Component" = "#7570b3")) +
    labs(
      title = paste("Unique Functional Drivers:", title_label),
      subtitle = "Excluding GO terms shared across both evolutionary phases",
      x = "Gene Ontology (GO) Term",
      y = "Number of Selected Genes"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 11, color = "black"),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text.y = element_text(size = 11, face = "bold", angle = 270),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold")
    )
  
  ggsave(filename_out, plot = p, width = 12, height = 9, dpi = 600)
}

# --- 5. Generate the Unique Plots ---
cat("Generating unique-only lollipop plots...\n")

plot_unique_lollipop(df_eros, unique_adapt, "Erosion", "Results/Erosion_Unique_GO_Faceted.png")
plot_unique_lollipop(df_intro, unique_breed, "Introgression", "Results/Introgression_Unique_GO_Faceted.png")

cat("Done! Both unique plots saved successfully.\n")

