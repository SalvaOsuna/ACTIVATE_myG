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
target_scenario <- "Py_Scenario1_Adaptation" # Change to Py_Scenario2_Breeding or Py_Scenario1_Adaptation for the other run
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


#generate GO####
library(Biostrings)
library(dplyr)
library(stringr)

cat("Loading protein sequences and candidate gene list...\n")

# 1. Define the scenario
target_scenario <- "Py_Scenario1_Adaptation" # Change to Py_Scenario2_Breeding or Py_Scenario1_Adaptation for the other run

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

# 5. Export the new, tiny FASTA file!
out_fasta <- sprintf("Results/%s_Candidate_Proteins.fasta", target_scenario)
writeXStringSet(candidate_peps, out_fasta)

cat(sprintf("Success! Extracted %d candidate protein sequences to %s\n", length(candidate_peps), out_fasta))

#ploting GO####
library(dplyr)
library(tidyr)
library(ggplot2)
library(GO.db)
library(AnnotationDbi)

cat("Loading InterPro annotations and candidate genes...\n")

# 1. Load your candidate genes
target_scenario <- "Py_Scenario2_Breeding" # Change to Py_Scenario2_Breeding or Py_Scenario1_Adaptation for the other run
candidates_df <- read.csv(sprintf("Results/%s_Candidate_Genes.csv", target_scenario), sep = ";")

# Clean the IDs in your candidates dataframe so they match the InterPro output exactly
candidates_df <- candidates_df %>%
  mutate(clean_id = gsub("^(gene:|mRNA:|transcript:)", "", ID))

# Reading TSV file using read_tsv()
library(readr)
# I did this intermediate step because in interpro I cannot scan more than 100 proteins
df <- read_tsv('Results/iprscan5-Breeding_1-100.tsv',col_names = F)
df2 <- read_tsv('Results/iprscan5_breeding_100-179.tsv',col_names = F)
interpro_file <- rbind(df, df2)
write_tsv(interpro_file, 'Results/Py_Scenario2_Breeding_interpro_results.tsv',col_names = F)

# 2. Load the headerless InterPro TSV
interpro_file <- "Results/Py_Scenario2_Breeding_interpro_results.tsv"

interpro_df <- read.delim(interpro_file, header = FALSE, stringsAsFactors = FALSE)

# Assign the proper standard InterPro column names
colnames(interpro_df) <- c("ID", "md5", "length", "analysis", 
                           "sig_acc", "sig_desc", "start", "stop", "evalue", 
                           "status", "date", "ipr_acc", "ipr_desc", 
                           "go_terms", "pathways")

# 3. Extract and Clean the GO terms
go_mapping <- interpro_df %>%
  dplyr::select(ID, go_terms) %>%
  # Remove rows where GO terms are missing or just a hyphen
  filter(go_terms != "" & !is.na(go_terms) & go_terms != "-") %>%
  # Split multiple GO terms on the pipe character
  separate_rows(go_terms, sep = "\\|") %>%
  # Strip out the database tags in parentheses (e.g., "(InterPro)" or "(PANTHER)")
  mutate(go_terms = gsub("\\(.*\\)", "", go_terms)) %>%
  # Remove the whitespace just in case
  mutate(go_terms = trimws(go_terms)) %>%
  # Keep only unique combinations of gene + GO term
  distinct(ID, go_terms)

go_mapping <- go_mapping %>%
  # Remove a literal dot (\\.) followed by digits ([0-9]+) at the end of the string ($)
  mutate(ID = sub("\\.[0-9]+$", "", ID)) %>%
  # Re-run distinct just in case multiple transcripts collapsed into the same gene ID
  distinct(ID, go_terms)

# 4. Merge with your candidate genes
final_go_df <- candidates_df %>%
  inner_join(go_mapping, by = "ID")

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


