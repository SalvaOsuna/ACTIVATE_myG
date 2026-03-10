#Gene Extraction and GO Plotting
BiocManager::install(c("GenomicRanges", "rtracklayer", "GO.db", "AnnotationDbi"))

library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(tidyr)
library(GO.db)
library(AnnotationDbi)

# --- 1. Load Data ---
cat("Loading GFF3 annotation and Significant XP-CLR regions...\n")

# Load the Lentil Annotation File
gff_file <- "data/Lcu.1GRN.genes_description.sorted.gff3.gz"
gff <- import(gff_file)

# Keep only "gene" or "mRNA" features to avoid counting exons multiple times
genes <- gff[gff$type %in% c("gene", "mRNA")]

# Load the Significant Regions from your target scenario
target_scenario <- "Py_Scenario2_Breeding" # Change to Py_Scenario2_Breeding or Py_Scenario1_Adaptation for the other run
sig_regions <- read.csv(sprintf("Results/%s_SignificantRegions.csv", target_scenario))

# Convert to GenomicRanges
sig_ranges <- GRanges(
  seqnames = sig_regions$chrom,
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

write.csv(candidates_df, sprintf("Results/%s_Candidate_Genes.csv", target_scenario), row.names = FALSE)
cat(sprintf("Found %d candidate genes! Saved to Results.\n", nrow(candidates_df)))

# --- 3. Process and Translate Gene Ontology (GO) Terms ---
cat("Translating GO IDs to readable biological terms...\n")

# Separate comma-separated GO terms and extract the raw IDs
go_df <- candidates_df %>%
  filter(!is.na(Ontology_term) & Ontology_term != "") %>%
  separate_rows(Ontology_term, sep = ",") %>%
  mutate(go_id = trimws(Ontology_term)) %>%
  # Ensure we only try to translate properly formatted GO IDs (e.g., "GO:0009628")
  filter(grepl("^GO:\\d+", go_id))

# Look up the English terms in the GO.db dictionary
# suppressMessages keeps the console output clean during the lookup
mapped_terms <- suppressMessages(
  select(GO.db, keys = unique(go_df$go_id), columns = "TERM", keytype = "GOID")
)

# Join the English terms back to our genes
go_df <- go_df %>%
  left_join(mapped_terms, by = c("go_id" = "GOID")) %>%
  # Fallback: if a GO ID is obsolete/missing in the database, keep the raw ID
  mutate(term_name = ifelse(is.na(TERM), go_id, TERM))

# Count the frequency of each biological term
go_counts <- go_df %>%
  count(term_name, name = "gene_count") %>%
  arrange(desc(gene_count)) %>%
  # Filter out highly generic terms if desired (optional)
  # filter(!term_name %in% c("biological_process", "molecular_function")) %>% 
  head(20) # Keep top 20 for the plot

# --- 4. Plot the Top Translated GO Terms ---
cat("Generating GO Term plot...\n")

plot_title <- gsub("_", " ", target_scenario)

p_go <- ggplot(go_counts, aes(x = reorder(term_name, gene_count), y = gene_count)) +
  geom_segment(aes(x = reorder(term_name, gene_count), xend = term_name, y = 0, yend = gene_count), color = "steelblue") +
  geom_point(color = "firebrick", size = 4, alpha = 0.8) +
  coord_flip() + # Flip coordinates for readable text
  labs(
    title = paste("Top 20 Biological Functions Selected in", plot_title),
    subtitle = "Gene Ontology (GO) terms extracted from top 1% XP-CLR sweeps",
    x = "Biological Function",
    y = "Number of Selected Genes"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 11, color = "black"), # Slightly larger text for readability
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold")
  )

# Save the plot
ggsave(sprintf("Results/%s_GO_Terms.png", target_scenario), plot = p_go, width = 11, height = 7, dpi = 600)
cat("Done! Beautiful GO Term plot saved.\n")

#word minning option####
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

cat("Extracting keywords from Gene Descriptions...\n")

# 1. Clean up the Description column
desc_df <- candidates_df %>%
  # Filter out missing descriptions
  filter(!is.na(Description) & Description != "character(0)" & Description != "") %>%
  # Coerce to character and make lowercase
  mutate(clean_desc = tolower(as.character(Description))) %>%
  # Remove special characters and punctuation
  mutate(clean_desc = str_replace_all(clean_desc, "[^a-z\\s]", ""))

# 2. Split sentences into individual words
word_df <- desc_df %>%
  separate_rows(clean_desc, sep = "\\s+") %>%
  filter(clean_desc != "")

# 3. Define "stop words" to ignore so we only get meaningful biology
stop_words <- c("protein", "putative", "uncharacterized", "domain", 
                "containing", "family", "type", "like", "a", "of", "and", 
                "the", "in", "is", "with", "factor", "subunit", "class")

# 4. Count and rank the keywords
word_counts <- word_df %>%
  filter(!clean_desc %in% stop_words) %>%
  count(clean_desc, name = "word_count") %>%
  arrange(desc(word_count)) %>%
  head(20) # Keep top 20

# 5. Plot the Keyword frequencies
p_keywords <- ggplot(word_counts, aes(x = reorder(clean_desc, word_count), y = word_count)) +
  geom_segment(aes(x = reorder(clean_desc, word_count), xend = clean_desc, y = 0, yend = word_count), color = "seagreen") +
  geom_point(color = "forestgreen", size = 4, alpha = 0.8) +
  coord_flip() + 
  labs(
    title = paste("Top 20 Functional Keywords in", plot_title),
    subtitle = "Mined from gene descriptions in top 1% XP-CLR sweeps",
    x = "Biological Keyword",
    y = "Number of Occurrences"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 11, color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold")
  )

# Save the plot
ggsave(sprintf("Results/%s_Keywords.png", target_scenario), plot = p_keywords, width = 10, height = 7, dpi = 600)
cat("Done! Keyword plot saved to Results folder.\n")