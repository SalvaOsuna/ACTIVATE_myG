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
target_scenario <- "XPCLR_Scenario1_Adaptation" # Change to XPCLR_Scenario2_Breeding or XPCLR_Scenario1_Adaptation for the other run

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
target_scenario <- "XPCLR_Scenario2_Breeding" # Change to XPCLR_Scenario2_Breeding or XPCLR_Scenario1_Adaptation for the other run
candidates_df <- read.csv(sprintf("Results/%s_Candidate_Genes.csv", target_scenario))

# Clean the IDs in your candidates dataframe so they match the InterPro output exactly
candidates_df <- candidates_df %>%
  mutate(clean_id = gsub("^(gene:|mRNA:|transcript:)", "", ID))

cat("Finding and stitching InterPro TSV chunks...\n")

# --- 2. Load and Stitch the Headerless InterPro TSVs ---

# Define the folder where you saved all those downloaded TSV chunks
# (Update this path if you saved them somewhere else!)
tsv_folder <- "Results/XPCLR_Scenario2_Breeding_InterPro_Chunks" 

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
write.csv(final_go_df, "Results/XPCLR_Scenario2_Breeding_Final_GO.csv", row.names = FALSE)

# 6. Count and rank WITHIN each Ontology category
go_counts <- final_go_df %>%
  count(ontology_full, term_name, name = "gene_count") %>%
  group_by(ontology_full) %>%
  # Grab the top 10 terms per category to keep the plot balanced
  slice_max(order_by = gene_count, n = 20, with_ties = FALSE) %>%
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

df_adapt <- read.csv("Results/XPCLR_Scenario1_Adaptation_Final_GO.csv", stringsAsFactors = FALSE)
df_breed <- read.csv("Results/XPCLR_Scenario2_Breeding_Final_GO.csv", stringsAsFactors = FALSE)

# --- 2. Filter for Biological Processes and extract unique terms ---
# We use unique() because we just want to know IF a pathway was targeted, 
# not how many times it was targeted.
bp_adapt <- df_adapt %>%
  filter(ONTOLOGY == "BP") %>%
  pull(term_name) %>%
  unique()

bp_breed <- df_breed %>%
  filter(ONTOLOGY == "BP") %>%
  pull(term_name) %>%
  unique()

# --- 3. Create the List Object for the Venn Diagram ---
venn_list <- list(
  "Adaptation (Scenario 1)" = bp_adapt,
  "Breeding (Scenario 2)" = bp_breed
)

# --- 4. Plot the Venn Diagram ---
cat("Generating plot...\n")

p_venn <- ggVennDiagram(venn_list, 
                        label_alpha = 0, # Removes the ugly white box behind the numbers
                        category.names = c("Adaptation", "Breeding"),
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
ggsave("Results/BP_Venn_Diagram.png", plot = p_venn, width = 8, height = 8, dpi = 600, bg = "white")
cat("Venn diagram saved to Results folder!\n\n")

# --- 5. Extract the specific terms for your manuscript ---
cat("--- Extracted Pathways for Manuscript Results ---\n")

# Shared by both
shared_bp <- intersect(bp_adapt, bp_breed)
cat(sprintf("Number of Shared BP Terms: %d\n", length(shared_bp)))

# Unique to Adaptation
unique_adapt <- setdiff(bp_adapt, bp_breed)
cat(sprintf("Number of Unique Adaptation BP Terms: %d\n", length(unique_adapt)))

# Unique to Breeding
unique_breed <- setdiff(bp_breed, bp_adapt)
cat(sprintf("Number of Unique Breeding BP Terms: %d\n", length(unique_breed)))

# Save these lists so you can read them easily
write.csv(data.frame(Term = shared_bp), "Results/Venn_Shared_BP.csv", row.names = FALSE)
write.csv(data.frame(Term = unique_adapt), "Results/Venn_Unique_Adaptation_BP.csv", row.names = FALSE)
write.csv(data.frame(Term = unique_breed), "Results/Venn_Unique_Breeding_BP.csv", row.names = FALSE)

cat("CSV files containing the exact overlapping and unique terms have been saved.\n")

#lolipop plot with just unique terms####
library(dplyr)
library(ggplot2)

cat("Categorizing GO Terms across scenarios...\n")

# --- 1. Load your final GO dataframes ---
# (Make sure these point to the dataframes you generated in the GO plotting step)
df_adapt <- read.csv("Results/XPCLR_Scenario1_Adaptation_Final_GO.csv", stringsAsFactors = FALSE)
df_breed <- read.csv("Results/XPCLR_Scenario2_Breeding_Final_GO.csv", stringsAsFactors = FALSE)

# --- 2. Identify Unique and Shared Terms ---
terms_adapt <- unique(df_adapt$term_name)
terms_breed <- unique(df_breed$term_name)

shared_terms <- intersect(terms_adapt, terms_breed)
unique_adapt <- setdiff(terms_adapt, terms_breed)
unique_breed <- setdiff(terms_breed, terms_adapt)

# --- 3. Create the Master Comparison CSV ---
# Build a simple dataframe tagging each term
presence_df <- bind_rows(
  data.frame(term_name = unique_adapt, Presence = "Scenario 1 (Adaptation Only)"),
  data.frame(term_name = unique_breed, Presence = "Scenario 2 (Breeding Only)"),
  data.frame(term_name = shared_terms, Presence = "Both Scenarios")
)

# Grab the Ontology categories to make the CSV more informative
ontology_map <- bind_rows(
  df_adapt %>% dplyr::select(term_name, ONTOLOGY, ontology_full),
  df_breed %>% dplyr::select(term_name, ONTOLOGY, ontology_full)
) %>% distinct(term_name, .keep_all = TRUE)

master_csv <- presence_df %>%
  left_join(ontology_map, by = "term_name") %>%
  arrange(Presence, ONTOLOGY, term_name)

write.csv(master_csv, "Results/GO_Terms_Comparison_Master.csv", row.names = FALSE)
cat("Master comparison CSV saved to 'Results/GO_Terms_Comparison_Master.csv'!\n\n")

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

plot_unique_lollipop(df_adapt, unique_adapt, "Adaptation Phase", "Results/Py_Scenario1_Unique_GO_Faceted.png")
plot_unique_lollipop(df_breed, unique_breed, "Modern Breeding Phase", "Results/Py_Scenario2_Unique_GO_Faceted.png")

cat("Done! Both unique plots saved successfully.\n")
