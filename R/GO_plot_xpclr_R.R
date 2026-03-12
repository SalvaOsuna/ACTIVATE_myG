#Gene Extraction and GO Plotting
BiocManager::install(c("GenomicRanges", "rtracklayer", "GO.db", "AnnotationDbi"))
BiocManager::install("Biostrings")
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(tidyr)
BiocManager::install("GO.db")
library(GO.db)
library(AnnotationDbi)
library(Biostrings)

# --- 1. Load Data ---
cat("Loading GFF3 annotation and Significant XP-CLR regions...\n")

# Load the Lentil Annotation File
gff_file <- "data/Lcu.1GRN.genes_description.sorted.gff3.gz"
gff <- import(gff_file)

# Keep only "gene" or "mRNA" features to avoid counting exons multiple times
genes <- gff[gff$type %in% c("gene", "mRNA")]

# Load the Significant Regions from your target scenario
target_scenario <- "XPCLR_Scenario2_Breeding" # Change to XPCLR_Scenario2_Breeding or XPCLR_Scenario1_Adaptation for the other run
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
  arrange(desc(z_score)) %>%
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

#generate GO####
library(Biostrings)
library(dplyr)
library(stringr)

cat("Loading protein sequences and candidate gene list...\n")

# 1. Define the scenario
target_scenario <- "XPCLR_Scenario2_Breeding" # Change to XPCLR_Scenario2_Breeding or XPCLR_Scenario1_Adaptation for the other run

# 2. Load your candidate genes
candidates_df <- read.csv(sprintf("Results/%s_Candidate_Genes.csv", target_scenario))

# Clean up the IDs (sometimes GFF3 adds prefixes like "gene:" or "mRNA:" that aren't in the FASTA)
query_ids <- str_remove_all(candidates_df$ID, "^(gene:|mRNA:|transcript:)")

# 3. Load the massive lentil protein FASTA file
fasta_file <- "data/Lcu.1GRN.pep.fasta"
pep_seqs <- readAAStringSet(fasta_file)

cat(sprintf("Loaded %d total protein sequences. Searching for our candidates...\n", length(pep_seqs)))

# 4. Find the exact matching sequences
# We create a logical vector to flag any FASTA header that contains one of our query IDs
keep_seqs <- rep(FALSE, length(pep_seqs))

for(q_id in query_ids) {
  # fixed(q_id) ensures exact string matching without regex confusion
  keep_seqs <- keep_seqs | str_detect(names(pep_seqs), fixed(q_id))
}

# Subset the FASTA file
candidate_peps <- pep_seqs[keep_seqs]

cat("Creating InterPro upload chunks...\n")

# 5. Define chunk size and create a specific folder
chunk_size <- 100
chunk_dir <- sprintf("Results/%s_InterPro_Chunks", target_scenario)

# Create the folder if it doesn't already exist
if (!dir.exists(chunk_dir)) {
  dir.create(chunk_dir, recursive = TRUE)
}

# 6. Split the protein sequences and save them
num_seqs <- length(candidate_peps)
# Create a vector that assigns each sequence to a chunk number (1, 1, ..., 2, 2, ...)
chunk_indices <- ceiling(seq_len(num_seqs) / chunk_size)

# Split the AAStringSet into a list of smaller AAStringSets
pep_chunks <- split(candidate_peps, chunk_indices)

cat(sprintf("Splitting %d sequences into %d chunks...\n", num_seqs, length(pep_chunks)))

# Loop through the list and save each chunk
for (i in seq_along(pep_chunks)) {
  # Format the filename with leading zeros (e.g., chunk_01.fasta, chunk_02.fasta)
  chunk_filename <- file.path(chunk_dir, sprintf("chunk_%02d.fasta", i))
  writeXStringSet(pep_chunks[[i]], chunk_filename)
}

cat(sprintf("Success! All chunks saved to the folder: %s\n", chunk_dir))

#ploting GO####
library(dplyr)
library(tidyr)
library(ggplot2)
library(GO.db)
library(AnnotationDbi)

cat("Loading InterPro annotations and candidate genes...\n")

# 1. Load your candidate genes
target_scenario <- "XPCLR_Scenario1_Adaptation" # Change to XPCLR_Scenario2_Breeding or XPCLR_Scenario1_Adaptation for the other run
candidates_df <- read.csv(sprintf("Results/%s_Candidate_Genes.csv", target_scenario))

# Clean the IDs in your candidates dataframe so they match the InterPro output exactly
candidates_df <- candidates_df %>%
  mutate(clean_id = gsub("^(gene:|mRNA:|transcript:)", "", ID))

cat("Finding and stitching InterPro TSV chunks...\n")

# --- 2. Load and Stitch the Headerless InterPro TSVs ---

# Define the folder where you saved all those downloaded TSV chunks
# (Update this path if you saved them somewhere else!)
tsv_folder <- "Results/XPCLR_Scenario1_Adaptation_InterPro_Chunks" 

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

p_go_faceted <- ggplot(go_counts, aes(x = reorder(term_name, gene_count), y = gene_count, color = ontology_full)) +
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


