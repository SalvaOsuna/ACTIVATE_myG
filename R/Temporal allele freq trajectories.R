# =============================================================================
# Temporal Allele Frequency Trajectories — Top XP-CLR Sweep Candidates
# =============================================================================
# Analysis #3 from peer review recommendations (Caballero et al.)
#
# WHAT THIS SCRIPT DOES:
#   1. Reads the ACTIVATE panel VCF + release-year metadata
#   2. Parses XP-CLR output for both sweep scenarios (adaptation / breeding)
#   3. Stratifies lines into decadal cohorts
#   4. Calculates alt allele frequency per candidate SNP per cohort
#   5. Fits a linear regression (frequency ~ mean cohort year) per SNP
#      with Benjamini-Hochberg FDR correction
#   6. Produces publication-ready ggplot2 figures + CSV outputs
#
# REQUIRED INPUT FILES:
#   A) VCF (.vcf or .vcf.gz)  — Merged two-population VCF (Merged_Analysis.vcf.gz)
#                                The script automatically subsets to ACTIVATE-only
#                                samples using the sample IDs in the metadata file,
#                                so LDP samples are excluded before any calculation.
#   B) Metadata CSV            — columns: sample_id, year_release (numeric)
#                                Must contain ACTIVATE sample IDs only.
#   C) XP-CLR CSVs             — one AllScores file per scenario:
#                                  XPCLR_Scenario1_Adaptation_AllScores.csv
#                                  XPCLR_Scenario2_Breeding_AllScores.csv
#                                columns: chr, start, stop, mid_pos_mb,
#                                n_snps, xpclr, z_score, chr_num, color_group
#                                AllScores (not ficantRegions) is used so the
#                                script can rank all windows and pick the true top N.
#
# OUTPUT FILES (written to CONFIG$output_dir):
#   allele_freq_by_cohort.csv       long-format frequency table per cohort
#   regression_results.csv          slope / p-value / FDR / R² per SNP
#   plot_trajectories_adaptation.pdf
#   plot_trajectories_breeding.pdf
#   plot_trajectories_combined.pdf  side-by-side faceted overview
#
# =============================================================================


# =============================================================================
# ── 1. CONFIGURATION — edit this block before running ────────────────────────
# =============================================================================

CONFIG <- list(
  
  # ── Input file paths ─────────────────────────────────────────────────────────
  # Merged two-population VCF (ACTIVATE + LDP combined).
  # The script will automatically restrict genotype extraction to ACTIVATE
  # samples only, using the sample IDs present in your metadata file.
  vcf_file      = "data/Merged_Analysis.vcf.gz",
  metadata_file = "data/activate_metadata.csv",
  
  # XP-CLR AllScores files — one per scenario.
  # Using AllScores (not ficantRegions) so the script ranks all windows
  # genome-wide and selects the true top N by XP-CLR score.
  xpclr_file_adaptation = "Results/XPCLR_Scenario1_Adaptation_AllScores.csv",
  xpclr_file_breeding   = "Results/XPCLR_Scenario2_Breeding_AllScores.csv",
  xpclr_file_combined   = NULL,   # <- set a path here if both scenarios are in one file
  
  # If using a combined file, what is the column name that identifies the scenario?
  # Values in that column must be exactly "adaptation" and "breeding".
  xpclr_scenario_col    = "scenario",
  
  # ── Metadata columns ─────────────────────────────────────────────────────────
  sample_id_col = "sample_id",    # column matching VCF sample names exactly
  year_col      = "year_release", # numeric year column (e.g. 2003)
  
  # ── Analysis parameters ──────────────────────────────────────────────────────
  top_n_windows    = 30,     # top N XP-CLR windows per scenario to use
  min_maf_global   = 0.05,   # discard SNPs with global MAF below this
  pvalue_threshold = 0.05,   # FDR-adjusted ficance threshold
  min_slope        = 0.003,  # minimum |Δfreq/year| — filters trivial trends
  
  # ── Output ──────────────────────────────────────────────────────────────────
  output_dir    = "Results/output_temporal_trajectories",
  figure_width  = 12,   # inches
  figure_height = 8,
  dpi           = 300
)


# =============================================================================
# ── 2. PACKAGES ───────────────────────────────────────────────────────────────
# =============================================================================

required_pkgs <- c(
  "SeqArray", "SeqVarTools",
  "dplyr", "tidyr", "purrr", "tibble",
  "ggplot2", "ggrepel", "scales",
  "RColorBrewer", "broom"
)

missing_pkgs <- required_pkgs[
  !sapply(required_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(missing_pkgs, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(SeqArray);  library(SeqVarTools)
  library(dplyr);     library(tidyr); library(purrr); library(tibble)
  library(ggplot2);   library(ggrepel); library(scales)
  library(RColorBrewer); library(broom)
})


# =============================================================================
# ── 3. HELPER FUNCTIONS ───────────────────────────────────────────────────────
# =============================================================================

# ── 3a. XP-CLR loader ─────────────────────────────────────────────────────────
# Reads your specific format:
#   chr | start | stop | mid_pos_mb | n_snps | xpclr | z_score | chr_num | color_group
# The 'xpclr' column is used as the score; 'start' and 'stop' define windows.

load_xpclr <- function(cfg) {
  
  read_one <- function(path, scenario_name) {
    df <- read.csv(path, stringsAsFactors = FALSE)
    colnames(df) <- tolower(trimws(colnames(df)))
    
    # Verify required columns are present
    required <- c("chr", "start", "stop", "xpclr")
    missing  <- setdiff(required, colnames(df))
    if (length(missing) > 0)
      stop("XP-CLR file '", basename(path), "' is missing columns: ",
           paste(missing, collapse = ", "),
           "\n  Found: ", paste(colnames(df), collapse = ", "))
    
    df %>%
      transmute(
        chr          = as.character(chr),
        window_start = as.integer(start),
        window_end   = as.integer(stop),
        xpclr_score  = as.numeric(xpclr),
        z_score      = if ("z_score" %in% colnames(df)) as.numeric(z_score) else NA_real_,
        n_snps_win   = if ("n_snps"  %in% colnames(df)) as.integer(n_snps)  else NA_integer_,
        scenario     = scenario_name
      ) %>%
      filter(!is.na(xpclr_score))
  }
  
  # Option A: two separate files
  if (!is.null(cfg$xpclr_file_adaptation) && !is.null(cfg$xpclr_file_breeding)) {
    adapt  <- read_one(cfg$xpclr_file_adaptation, "adaptation")
    breed  <- read_one(cfg$xpclr_file_breeding,   "breeding")
    return(bind_rows(adapt, breed))
  }
  
  # Option B: single combined file with scenario column
  if (!is.null(cfg$xpclr_file_combined)) {
    df <- read.csv(cfg$xpclr_file_combined, stringsAsFactors = FALSE)
    colnames(df) <- tolower(trimws(colnames(df)))
    scen_col <- tolower(cfg$xpclr_scenario_col)
    if (!scen_col %in% colnames(df))
      stop("Scenario column '", scen_col, "' not found in combined XP-CLR file.\n",
           "  Found: ", paste(colnames(df), collapse = ", "))
    return(df %>%
             transmute(
               chr          = as.character(chr),
               window_start = as.integer(start),
               window_end   = as.integer(stop),
               xpclr_score  = as.numeric(xpclr),
               z_score      = if ("z_score" %in% colnames(df)) as.numeric(z_score) else NA_real_,
               n_snps_win   = if ("n_snps"  %in% colnames(df)) as.integer(n_snps)  else NA_integer_,
               scenario     = .data[[scen_col]]
             ) %>%
             filter(!is.na(xpclr_score)))
  }
  
  stop("No XP-CLR file(s) specified. ",
       "Set either xpclr_file_adaptation + xpclr_file_breeding, ",
       "or xpclr_file_combined in CONFIG.")
}

# ── 3b. Cohort assignment ─────────────────────────────────────────────────────

assign_decadal_cohort <- function(years) {
  decade_start  <- floor(years / 5) * 5
  cohort_labels <- paste0(decade_start, "s")
  factor(cohort_labels, levels = unique(sort(cohort_labels)))
}

# ── 3c. Allele frequency helpers ──────────────────────────────────────────────

# Alternative allele frequency (directional — not folded at 0.5)
# This keeps trajectories monotone so rising/falling trends are visually clear.
calc_alt_freq <- function(d) {
  d <- d[!is.na(d)]
  if (!length(d)) return(NA_real_)
  sum(d) / (2 * length(d))
}

calc_maf <- function(d) {
  af <- calc_alt_freq(d)
  if (is.na(af)) return(NA_real_)
  min(af, 1 - af)
}

# ── 3d. Regression ────────────────────────────────────────────────────────────

run_regression <- function(df) {
  df <- df %>% filter(!is.na(mean_year), !is.na(alt_freq))
  if (nrow(df) < 3) return(NULL)     # need at least 3 cohorts to fit a line
  fit <- lm(alt_freq ~ mean_year, data = df)
  broom::tidy(fit) %>%
    filter(term == "mean_year") %>%
    rename(slope = estimate) %>%
    mutate(r_squared = summary(fit)$r.squared)
}


# =============================================================================
# ── 4. LOAD & VALIDATE DATA ───────────────────────────────────────────────────
# =============================================================================

if (!dir.exists(CONFIG$output_dir))
  dir.create(CONFIG$output_dir, recursive = TRUE)

# ── 4a. Metadata ──────────────────────────────────────────────────────────────
message("\n[1/6] Loading metadata...")

meta <- read.csv(CONFIG$metadata_file, stringsAsFactors = FALSE, sep = ";") %>%
  rename_with(tolower) %>%
  rename(sample_id    = all_of(tolower(CONFIG$sample_id_col)),
         year_release = all_of(tolower(CONFIG$year_col))) %>%
  mutate(year_release = as.integer(year_release)) %>%
  filter(!is.na(year_release)) %>%
  mutate(cohort = assign_decadal_cohort(year_release))

cohort_summary <- meta %>%
  group_by(cohort) %>%
  summarise(n         = n(),
            year_min  = min(year_release),
            year_max  = max(year_release),
            mean_year = mean(year_release),
            .groups   = "drop")

message("  Lines with valid release year: ", nrow(meta))
message("  Cohort breakdown:")
print(cohort_summary)

# Warn if any cohort is very small
small_cohorts <- cohort_summary %>% filter(n < 5)
if (nrow(small_cohorts) > 0)
  warning("Cohorts with < 5 lines — regression will be unreliable for these:\n  ",
          paste(paste0(small_cohorts$cohort, " (n=", small_cohorts$n, ")"),
                collapse = ", "))

# ── 4b. XP-CLR ────────────────────────────────────────────────────────────────
message("\n[2/6] Loading XP-CLR results...")

xpclr_all <- load_xpclr(CONFIG)

detected_scenarios <- unique(xpclr_all$scenario)
message("  Scenarios detected: ", paste(detected_scenarios, collapse = ", "))
message("  Total windows loaded: ", nrow(xpclr_all))

if (!all(c("adaptation", "breeding") %in% detected_scenarios))
  warning("Expected scenarios 'adaptation' and 'breeding'.\n",
          "  Got: ", paste(detected_scenarios, collapse = ", "),
          "\n  Check file paths or scenario column values in your CSVs.")

# Select top N windows per scenario by raw XP-CLR score
top_windows <- xpclr_all %>%
  group_by(scenario) %>%
  slice_max(xpclr_score, n = CONFIG$top_n_windows, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(scenario, chr, window_start)

message("\n  Top ", CONFIG$top_n_windows,
        " windows selected per scenario (ranked by XP-CLR score):")
print(top_windows %>%
        select(scenario, chr, window_start, window_end, xpclr_score, z_score, n_snps_win))


# =============================================================================
# ── 5. VCF → GDS CONVERSION + GENOTYPE EXTRACTION ────────────────────────────
# =============================================================================

message("\n[3/6] Preparing GDS file from merged VCF...")
message("  Source: Merged_Analysis.vcf.gz (contains both ACTIVATE and LDP samples)")
message("  Genotype extraction will be restricted to ACTIVATE samples only")
message("  (i.e. samples present in your metadata file).")

gds_path <- file.path(CONFIG$output_dir, "merged_analysis.gds")

if (!file.exists(gds_path)) {
  message("  Converting VCF to GDS (one-time step, may take several minutes")
  message("  for a merged two-population VCF)...")
  seqVCF2GDS(CONFIG$vcf_file, gds_path, verbose = FALSE)
  message("  Created: ", gds_path)
} else {
  message("  GDS already exists — skipping conversion.")
}

gds <- seqOpen(gds_path)
gds_samples <- seqGetData(gds, "sample.id")

# Keep only ACTIVATE samples — those present in the metadata.
# LDP samples are in the VCF but absent from the metadata, so they are
# automatically excluded here. No manual filtering step is needed.
common_samples <- intersect(gds_samples, meta$sample_id)

if (length(common_samples) == 0)
  stop(
    "No ACTIVATE sample IDs matched between the merged VCF and metadata.\n",
    "  VCF example IDs  : ", paste(head(gds_samples,      3), collapse = ", "), "\n",
    "  Meta example IDs : ", paste(head(meta$sample_id,   3), collapse = ", "), "\n",
    "  Check that sample_id_col = '", CONFIG$sample_id_col,
    "' matches the VCF sample names exactly."
  )

n_ldp_excluded <- length(gds_samples) - length(common_samples)
message("  Total samples in merged VCF : ", length(gds_samples))
message("  ACTIVATE samples kept       : ", length(common_samples))
message("  LDP samples excluded        : ", n_ldp_excluded,
        "  (absent from metadata — correctly excluded)")

# Restrict GDS and metadata to ACTIVATE samples only
seqSetFilter(gds, sample.id = common_samples)
meta <- meta %>% filter(sample_id %in% common_samples)

# ── Extract genotypes for every top window ────────────────────────────────────
message("\n[4/6] Extracting genotypes from top XP-CLR windows...")

# NOTE on chromosome names:
#   Your VCF likely uses "Lcu.1GRN.Chr1" … "Lcu.1GRN.Chr7".
#   Your XP-CLR file may use the same names OR shorter forms like "Chr1" or "1".
#   The script tries to match as-is first; if nothing is found it will warn you.
#   If names don't match, add a chr_name_map below (e.g. "1" -> "Lcu.1GRN.Chr1").
#   Example: chr_name_map <- setNames(paste0("Lcu.1GRN.Chr", 1:7), as.character(1:7))

# Retrieve chromosome names actually present in the GDS for reference
gds_chromosomes <- unique(seqGetData(gds, "chromosome"))
message("  Chromosomes in GDS: ", paste(sort(gds_chromosomes), collapse = ", "))
message("  Chromosomes in XP-CLR: ",
        paste(sort(unique(top_windows$chr)), collapse = ", "))

# Optional name mapping — edit if your names differ (set to NULL to skip)
# Format: named vector where names = XP-CLR names, values = GDS/VCF names
# Example for short -> full names:
#   chr_name_map <- setNames(paste0("Lcu.1GRN.Chr", 1:7), paste0("Chr", 1:7))
chr_name_map <- NULL

extract_window <- function(chr, window_start, window_end,
                           xpclr_score, z_score, scenario, ...) {
  # Apply chromosome name mapping if provided
  chr_gds <- if (!is.null(chr_name_map) && chr %in% names(chr_name_map))
    chr_name_map[[chr]] else chr
  
  # Step 1: reset filters, then restrict to ACTIVATE samples only
  seqResetFilter(gds)
  seqSetFilter(gds, sample.id = common_samples, verbose = FALSE)
  
  # Step 2: build a logical variant selector for this genomic window.
  # We read chromosome + position for all variants, then filter by range.
  # This avoids seqSetFilterPos() which is not available in all SeqArray versions.
  all_chr <- seqGetData(gds, "chromosome")
  all_pos <- seqGetData(gds, "position")
  var_sel  <- all_chr == as.character(chr_gds) &
    all_pos >= window_start &
    all_pos <= window_end
  
  if (!any(var_sel)) {
    message("  ⚠  No variants found in ", chr, ":", window_start, "-", window_end,
            "\n     If this is unexpected, check chromosome name format.",
            "\n     GDS uses: ", paste(head(gds_chromosomes, 3), collapse = ", "))
    return(NULL)
  }
  
  # Step 3: apply both sample and variant filters together
  seqSetFilter(gds, sample.id = common_samples,
               variant.sel = var_sel, verbose = FALSE)
  
  n_var <- sum(var_sel)
  
  # Dosage matrix: rows = samples, cols = variants (0 = hom-ref, 1 = het, 2 = hom-alt)
  geno    <- seqGetData(gds, "$dosage")
  pos     <- seqGetData(gds, "position")
  chrom   <- seqGetData(gds, "chromosome")
  snp_ids <- paste0(chrom, ":", pos)
  colnames(geno) <- snp_ids
  
  # Explicitly assign sample names as rownames.
  # seqGetData("$dosage") returns a numeric matrix whose rownames are
  # positional integers (1, 2, 3 ...), not sample ID strings.
  # We retrieve the currently active sample IDs from the GDS filter
  # and assign them so the subsequent left_join on sample_id works correctly.
  rownames(geno) <- seqGetData(gds, "sample.id")
  
  # Apply global MAF filter
  global_maf <- apply(geno, 2, calc_maf)
  geno <- geno[, global_maf >= CONFIG$min_maf_global, drop = FALSE]
  if (ncol(geno) == 0) return(NULL)
  
  as.data.frame(geno) %>%
    rownames_to_column("sample_id") %>%
    pivot_longer(-sample_id, names_to = "snp_id", values_to = "dosage") %>%
    mutate(
      chr          = chr,
      window_start = window_start,
      window_end   = window_end,
      xpclr_score  = xpclr_score,
      z_score      = z_score,
      scenario     = scenario
    )
}

all_geno <- pmap_dfr(top_windows, extract_window)
seqClose(gds)

# --- NEW: Handle Overlapping Windows ---
# Sort by scenario, sample, and SNP, prioritizing the highest XP-CLR score.
# Then use distinct() to drop the duplicates, keeping the SNP assigned 
# to the window with the strongest selection signal.
all_geno <- all_geno %>%
  arrange(scenario, sample_id, snp_id, desc(xpclr_score)) %>%
  distinct(scenario, sample_id, snp_id, .keep_all = TRUE)

n_snps <- n_distinct(all_geno$snp_id)
message("  Extracted ", n_snps, " UNIQUE candidate SNPs across ",
        nrow(top_windows), " windows (duplicates from overlapping regions removed).")

# Sanity check: confirm sample IDs in all_geno match the metadata
n_matched <- n_distinct(all_geno$sample_id[all_geno$sample_id %in% meta$sample_id])
n_unmatched <- n_distinct(all_geno$sample_id[!all_geno$sample_id %in% meta$sample_id])
message("  Sample ID check — matched to metadata : ", n_matched,
        "  |  unmatched : ", n_unmatched)
if (n_unmatched > 0) {
  unmatched_examples <- head(
    unique(all_geno$sample_id[!all_geno$sample_id %in% meta$sample_id]), 5)
  warning("Some sample IDs in the genotype matrix did not match metadata.\n",
          "  Examples: ", paste(unmatched_examples, collapse = ", "), "\n",
          "  These rows will be dropped in the left_join and excluded from",
          " frequency calculations.\n",
          "  Check that your metadata sample_id column matches VCF sample names exactly.")
}

if (n_snps == 0)
  stop(
    "No SNPs were extracted from any window.\n",
    "  Most likely cause: chromosome name mismatch between VCF and XP-CLR output.\n",
    "  VCF chromosomes  : ", paste(sort(gds_chromosomes), collapse = ", "), "\n",
    "  XP-CLR chromosomes: ", paste(sort(unique(top_windows$chr)), collapse=", "), "\n",
    "  Fix: set chr_name_map in section [4/6] of this script."
  )


# =============================================================================
# ── 6. ALLELE FREQUENCIES PER COHORT ─────────────────────────────────────────
# =============================================================================

message("\n[5/6] Calculating cohort allele frequencies and running regressions...")

freq_table <- all_geno %>%
  left_join(meta %>% select(sample_id, cohort, year_release),
            by = "sample_id") %>%
  filter(!is.na(cohort), !is.na(dosage)) %>%
  group_by(scenario, snp_id, chr, window_start, window_end,
           xpclr_score, z_score, cohort) %>%
  summarise(
    n_lines  = n(),
    alt_freq = calc_alt_freq(dosage),
    maf      = calc_maf(dosage),
    .groups  = "drop"
  ) %>%
  left_join(cohort_summary %>% select(cohort, mean_year, n),
            by = "cohort") %>%
  rename(cohort_size = n)

write.csv(freq_table,
          file.path(CONFIG$output_dir, "allele_freq_by_cohort.csv"),
          row.names = FALSE)
message("  Saved: allele_freq_by_cohort.csv  (", nrow(freq_table), " rows)")

# ── Linear regression: alt_freq ~ mean_year, one model per SNP ───────────────

reg_results <- freq_table %>%
  group_by(scenario, snp_id, chr, window_start, xpclr_score, z_score) %>%
  group_modify(~ {
    res <- run_regression(.x)
    if (is.null(res)) tibble() else res
  }) %>%
  ungroup() %>%
  mutate(
    # BH correction across all tests (both scenarios together)
    fdr_p       = p.adjust(p.value, method = "BH"),
    significant = fdr_p < CONFIG$pvalue_threshold &
      abs(slope) >= CONFIG$min_slope,
    direction   = case_when(
      slope >  CONFIG$min_slope  ~ "increasing",
      slope < -CONFIG$min_slope  ~ "decreasing",
      TRUE                       ~ "stable"
    )
  ) %>%
  arrange(scenario, fdr_p)

write.csv(reg_results,
          file.path(CONFIG$output_dir, "regression_results.csv"),
          row.names = FALSE)

n_sig <- sum(reg_results$significant, na.rm = TRUE)
message("  Saved: regression_results.csv")
message("  Significant trajectories: ", n_sig,
        "  (FDR < ", CONFIG$pvalue_threshold,
        ", |slope| >= ", CONFIG$min_slope, " /yr)")

if (n_sig > 0) {
  message("  Summary:")
  print(reg_results %>%
          filter(significant) %>%
          select(scenario, snp_id, chr, xpclr_score, z_score,
                 slope, p.value, fdr_p, r_squared, direction))
}


# =============================================================================
# ── 7. FIGURES ────────────────────────────────────────────────────────────────
# =============================================================================

message("\n[6/6] Generating figures...")

DIRECTION_PALETTE <- c(
  increasing = "#D6604D",   # warm red — allele increasing under selection
  decreasing = "#4393C3",   # blue     — allele decreasing
  stable     = "grey70"
)

# For each window, pick the single best SNP to label (lowest FDR, then |slope|)
label_ids <- reg_results %>%
  filter(significant) %>%
  group_by(scenario, window_start) %>%
  slice_min(fdr_p, n = 1, with_ties = FALSE) %>%
  slice_max(abs(slope), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  pull(snp_id)

# Merge regression info into frequency table
plot_df <- freq_table %>%
  left_join(
    reg_results %>% select(scenario, snp_id, slope, p.value,
                           fdr_p, r_squared, significant, direction),
    by = c("scenario", "snp_id")
  )

# Label data: position at the last (most recent) cohort
label_df <- plot_df %>%
  filter(snp_id %in% label_ids) %>%
  group_by(snp_id, scenario) %>%
  filter(mean_year == max(mean_year, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(label = sprintf("%s\nβ = %.4f/yr\nFDR = %.3f", snp_id, slope, fdr_p))

# Shared theme
theme_traj <- theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9,  colour = "grey35"),
    plot.caption     = element_text(size = 8,  colour = "grey50", hjust = 0),
    axis.text.x      = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey93", linewidth = 0.4),
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey96"),
    strip.text       = element_text(face = "bold"),
    plot.margin      = margin(10, 100, 10, 10)   # extra right margin for labels
  )

# ── Per-scenario trajectory plot ──────────────────────────────────────────────
make_scenario_plot <- function(scen) {
  d  <- plot_df  %>% filter(scenario == scen) 
  ld <- label_df %>% filter(scenario == scen)
  
  # Unique colour per significant SNP using Spectral palette
  sig_ids <- unique(d %>% filter(significant) %>% pull(snp_id))
  n_col   <- max(length(sig_ids), 1)
  sig_pal <- setNames(
    colorRampPalette(brewer.pal(min(n_col + 2, 11), "Spectral"))(n_col),
    sig_ids
  )
  
  n_sig_plot <- n_distinct(sig_ids)
  n_tot_plot <- n_distinct(d$snp_id)
  
  ggplot(d, aes(x = mean_year, y = alt_freq, group = snp_id)) +
    
    # ── Background: all non-significant SNPs in grey ──────────────────────────
    geom_line(
      data    = d %>% filter(!significant | is.na(significant)),
      colour  = "grey83", linewidth = 0.35, alpha = 0.6
    ) +
    geom_point(
      data    = d %>% filter(!significant | is.na(significant)),
      colour  = "grey78", size = 1.0, alpha = 0.5
    ) +
    
    # ── Foreground: significant SNPs, each with a distinct colour ─────────────
    geom_line(
      data      = d %>% filter(significant),
      aes(colour = snp_id), linewidth = 1.0
    ) +
    geom_point(
      data      = d %>% filter(significant),
      aes(colour = snp_id), size = 2.5
    ) +
    
    # ── Labels at last cohort ─────────────────────────────────────────────────
    geom_label_repel(
      data          = ld,
      aes(label = label, colour = snp_id),
      size          = 2.4,
      label.size    = 0.15,
      label.padding = unit(0.12, "lines"),
      nudge_x       = 2,
      direction     = "y",
      hjust         = 0,
      show.legend   = FALSE
    ) +
    
    # ── MAF = 0.5 reference line ──────────────────────────────────────────────
    geom_hline(yintercept = 0.5, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    
    # ── Per-cohort sample size below x-axis ──────────────────────────────────
    geom_text(
      data        = distinct(d, cohort, mean_year, cohort_size),
      aes(x = mean_year, y = -0.10, label = paste0("n=", cohort_size)),
      inherit.aes = FALSE, size = 2.8, colour = "grey40"
    ) +
    
    scale_colour_manual(values = sig_pal, guide = "none") +
    scale_x_continuous(
      name   = "Mean cohort release year",
      breaks = cohort_summary$mean_year,
      labels = round(cohort_summary$mean_year)
    ) +
    scale_y_continuous(
      name   = "Alternative allele frequency",
      limits = c(-0.13, 1.0),
      breaks = seq(0, 1, 0.2),
      labels = percent_format(accuracy = 1)
    ) +
    labs(
      title    = sprintf("Allele frequency trajectories — %s phase",
                         tools::toTitleCase(scen)),
      subtitle = sprintf(
        paste0("Top %d XP-CLR windows | ",
               "%d / %d candidate SNPs with significant trends ",
               "(FDR < %.2f, |β| \u2265 %.3f/yr)"),
        CONFIG$top_n_windows, n_sig_plot, n_tot_plot,
        CONFIG$pvalue_threshold, CONFIG$min_slope
      ),
      caption  = paste0(
        "Grey lines: non-significant SNPs.  ",
        "Coloured lines: FDR-significant linear trends.\n",
        "\u03b2 = change in alt allele frequency per year.  ",
        "n = number of lines in each decadal cohort."
      )
    ) +
    theme_traj
}

for (scen in detected_scenarios) {
  p   <- make_scenario_plot(scen)
  out <- file.path(CONFIG$output_dir,
                   paste0("plot_trajectories_", scen, ".pdf"))
  ggsave(out, p,
         width  = CONFIG$figure_width,
         height = CONFIG$figure_height,
         dpi    = CONFIG$dpi)
  message("  Saved: ", basename(out))
}

# ── Combined faceted plot (both scenarios side by side) ───────────────────────
p_combined <- ggplot(plot_df,
                     aes(x = mean_year, y = alt_freq, group = snp_id)) +
  
  geom_line(
    data   = plot_df %>% filter(!significant | is.na(significant)),
    colour = "grey84", linewidth = 0.3, alpha = 0.55
  ) +
  geom_line(
    data = plot_df %>% filter(significant),
    aes(colour = direction), linewidth = 0.9, alpha = 0.8
  ) +
  geom_point(
    data = plot_df %>% filter(significant),
    aes(colour = direction), size = 1.8, alpha = 0.8
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             colour = "grey50", linewidth = 0.35) +
  
  scale_colour_manual(
    name   = "Trend direction",
    values = DIRECTION_PALETTE[c("increasing", "decreasing")],
    labels = c("Increasing (alt allele favoured by selection)",
               "Decreasing (ref allele favoured by selection)")
  ) +
  scale_x_continuous(
    name   = "Mean cohort release year",
    breaks = cohort_summary$mean_year,
    labels = round(cohort_summary$mean_year)
  ) +
  scale_y_continuous(
    name   = "Alternative allele frequency",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = percent_format(accuracy = 1)
  ) +
  facet_wrap(
    ~ scenario, ncol = 2,
    labeller = as_labeller(c(adaptation = "Adaptation phase",
                             breeding   = "Breeding phase"))
  ) +
  labs(
    title    = "Temporal allele frequency trajectories — XP-CLR top sweep candidates",
    subtitle = sprintf(
      "Top %d windows per scenario | Coloured: FDR-significant trends (FDR < %.2f)",
      CONFIG$top_n_windows, CONFIG$pvalue_threshold
    ),
    caption  = paste0(
      "Each line = one candidate SNP.  Grey: non-significant.  ",
      "Coloured by direction of allele frequency change over breeding decades."
    )
  ) +
  theme_traj +
  theme(legend.position = "bottom")

ggsave(
  file.path(CONFIG$output_dir, "plot_trajectories_combined.pdf"),
  p_combined,
  width  = CONFIG$figure_width,
  height = CONFIG$figure_height * 0.75,
  dpi    = CONFIG$dpi
)
message("  Saved: plot_trajectories_combined.pdf")


# =============================================================================
# ── 8. CONSOLE SUMMARY ────────────────────────────────────────────────────────
# =============================================================================

cat("\n", strrep("=", 68), "\n")
cat("ANALYSIS COMPLETE\n")
cat(strrep("=", 68), "\n\n")
cat("Output directory: ", CONFIG$output_dir, "\n\n")

cat("Files written:\n")
cat("  allele_freq_by_cohort.csv\n")
cat("  regression_results.csv\n")
for (s in detected_scenarios)
  cat("  plot_trajectories_", s, ".pdf\n", sep = "")
cat("  plot_trajectories_combined.pdf\n\n")

cat("Cohort summary:\n")
print(cohort_summary)

cat("\nSignificant trajectories (FDR < ", CONFIG$pvalue_threshold,
    ", |slope| >= ", CONFIG$min_slope, "):\n", sep = "")

if (n_sig > 0) {
  reg_results %>%
    filter(significant) %>%
    count(scenario, direction, name = "n_snps") %>%
    print()
  
  cat("\nTop 5 strongest trajectories (by |slope|):\n")
  reg_results %>%
    filter(significant) %>%
    arrange(desc(abs(slope))) %>%
    select(scenario, snp_id, chr, xpclr_score, z_score,
           slope, fdr_p, r_squared, direction) %>%
    head(5) %>%
    print()
} else {
  cat("  None detected with current thresholds.\n")
  cat("  Suggestions:\n")
  cat("  - Increase top_n_windows (currently ", CONFIG$top_n_windows, ")\n")
  cat("  - Relax pvalue_threshold (currently ", CONFIG$pvalue_threshold, ")\n")
  cat("  - Lower min_slope (currently ", CONFIG$min_slope, ")\n")
  cat("  - Check cohort sizes — cohorts with < 5 lines reduce statistical power\n")
}

cat("\n")

# =============================================================================
# ── 8. COMBINED PUBLICATION FIGURE (MANHATTAN + TRAJECTORIES) ────────────────
# =============================================================================
cat("\nGenerating Combined 3-Panel Publication Figure...\n")

if(!require(patchwork)) install.packages("patchwork")
if(!require(ggrepel)) install.packages("ggrepel")
library(patchwork)
library(ggrepel)
library(scales)

# --- 1. Define Gene Labels ---
labels_adapt <- data.frame(
  chr_num = c(2,4,2,4),
  bp_mid = c(52.657,106.326,416.548,288.766), # Replace with your actual physical positions in Mb!
  xpclr_val = c(300,435,405,240),             # The height of the peak to attach the label to
  gene_name = c("Phytochrome","Meristem maintenance/development","Structural chrom maintenance","GIGANTEA")
)

labels_breed <- data.frame(
  chr_num = c(3, 6, 6),
  bp_mid = c(110.047, 292.989, 375.836), 
  xpclr_val = c(200, 235, 240),
  gene_name = c("Cellulose synthase", "Actin-depolymerizing factor","Meristem maintenance/development")
)

# --- 2. Helper Function: Calculate Cumulative X-Axis ---
# This uses your exact math to turn 7 chromosomes into one continuous X-axis
prep_cum_data <- function(df_sub, threshold_quantile = 0.99) {
  
  # FIXED: Strip everything before and including "Chr" so "Lcu.1GRN.Chr1" becomes just "1"
  df_sub <- df_sub %>% mutate(chr_num = as.numeric(gsub(".*Chr", "", chr)))
  
  # Calculate cumulative offsets
  chr_lens <- df_sub %>% group_by(chr_num) %>% summarize(chr_max = max(window_start, na.rm = TRUE), .groups = "drop")
  chr_lens$offset <- cumsum(as.numeric(chr_lens$chr_max)) - chr_lens$chr_max
  
  df_sub <- df_sub %>% 
    left_join(chr_lens, by = "chr_num") %>% 
    mutate(bp_cum = window_start + offset)
  
  # Calculate center of each chromosome for the X-axis labels
  axis_df <- df_sub %>% group_by(chr_num) %>% summarize(center = mean(bp_cum, na.rm = TRUE), .groups = "drop")
  
  thresh_val <- quantile(df_sub$xpclr_score, threshold_quantile, na.rm = TRUE)
  
  return(list(data = df_sub, axis = axis_df, threshold = thresh_val))
}

# Process the data already loaded in your environment (xpclr_all)
data_adapt <- prep_cum_data(xpclr_all %>% filter(scenario == "adaptation"))
data_breed <- prep_cum_data(xpclr_all %>% filter(scenario == "breeding"))

# Color Palette
chr_colors <- c("1"="#1b9e77", "2"="#d95f02", "3"="#7570b3", "4"="#e7298a", "5"="#66a61e", "6"="#e6ab02", "7"="#a6761d")

# --- 3. Helper Function: Plot Manhattan/Lollipop ---
make_lollipop_plot <- function(prep_data, show_x_labels = FALSE, labels_df = NULL) {
  df <- prep_data$data
  ax <- prep_data$axis
  thr <- prep_data$threshold
  
  p <- ggplot(df, aes(x = bp_cum, y = xpclr_score, color = as.character(chr_num))) +
    geom_segment(aes(xend = bp_cum, yend = 0), alpha = 0.8, linewidth = 0.4) +
    geom_point(size = 1.2, alpha = 0.9) +
    geom_hline(yintercept = thr, color = "firebrick", linetype = "dashed", linewidth = 0.8) +
    annotate("text", x = min(df$bp_cum) - (max(df$bp_cum)*0.01), y = thr, 
             label = round(thr, 2), color = "firebrick", hjust = 1, vjust = -0.5, size = 3.5, fontface = "bold") +
    scale_color_manual(values = chr_colors) +
    scale_x_continuous(breaks = ax$center, labels = paste("LG", ax$chr_num, sep=""), expand = c(0.02, 0.02)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) + # Extra room for ggrepel labels
    coord_cartesian(clip = "off") + 
    labs(x = NULL, y = "XP-CLR") +
    
    # 10pt Standard Theme
    theme_classic(base_size = 10) +
    theme(
      legend.position = "none",
      axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
      axis.text.y = element_text(color = "black"),
      axis.text.x = if(show_x_labels) element_text(face = "bold", color = "black") else element_blank(),
      axis.ticks.x = if(show_x_labels) element_line(color = "black") else element_blank(),
      axis.line.x = if(show_x_labels) element_line(color = "black") else element_blank(),
      plot.margin = margin(t = 5, r = 10, b = if(show_x_labels) 10 else 2, l = 15)
    )
  
  # Add ggrepel annotations
  if (!is.null(labels_df)) {
    # Dynamically apply the offset to the label coordinates
    offset_map <- df %>% distinct(chr_num, offset)
    # Convert Mb back to raw bp to match the window_start scale
    labels_mapped <- labels_df %>%
      mutate(bp_raw = bp_mid * 1e6) %>%
      left_join(offset_map, by = "chr_num") %>%
      mutate(bp_cum = bp_raw + offset)
    
    p <- p + geom_label_repel(
      data = labels_mapped,
      aes(x = bp_cum, y = xpclr_val, label = gene_name),
      inherit.aes = FALSE,
      color = "black", fill = "white", fontface = "bold", size = 3.5,
      box.padding = 1.2, point.padding = 0.5, segment.color = "black", segment.size = 0.6,
      min.segment.length = 0,
      nudge_y = max(df$xpclr_score) * 0.15 
    )
  }
  return(p)
}

# --- 4. Generate the Individual Plots ---
p_a <- make_lollipop_plot(data_adapt, show_x_labels = FALSE, labels_df = labels_adapt)
p_b <- make_lollipop_plot(data_breed, show_x_labels = TRUE, labels_df = labels_breed)

# Trajectory Plot (From previous step)
p_c <- ggplot(plot_df, aes(x = mean_year, y = alt_freq, group = snp_id)) +
  geom_line(data = plot_df %>% 
              filter(!significant | is.na(significant)), 
            colour = "grey84", 
            linewidth = 0.3, 
            alpha = 0.55) +
  geom_line(data = plot_df 
            %>% filter(significant), 
            aes(colour = direction), 
            linewidth = 0.9, 
            alpha = 0.8) +
  geom_point(data = plot_df %>% 
               filter(significant), 
             aes(colour = direction), 
             size = 1.8, 
             alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50", linewidth = 0.35) +
  scale_colour_manual(values = DIRECTION_PALETTE[c("increasing", "decreasing")], labels = c("alt allele favoured by selection", "ref allele favoured by selection)")) +
  scale_x_continuous(breaks = cohort_summary$mean_year, labels = round(cohort_summary$mean_year)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = percent_format(accuracy = 1)) +
  facet_wrap(~ scenario, ncol = 2, labeller = as_labeller(c(adaptation = "Adaptation phase", breeding = "Breeding phase"))) +
  labs(x = "Mean cohort release year", y = "Alternative allele frequency", color = NULL) +
  
  # 10pt Standard Theme
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey96", color = "black"),
    strip.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 15, r = 5, b = 5, l = 5)
  )

# --- 5. Assemble and Save the Final Plot ---
final_plot <- p_a / p_b / p_c + 
  plot_layout(heights = c(1, 1, 1.3)) +
  plot_annotation(tag_levels = 'a', tag_suffix = ')') & 
  theme(plot.tag = element_text(face = 'bold', size = 10))

# Save strictly to A4 width specifications
ggsave(file.path(CONFIG$output_dir, "Figure_XPCLR_Trajectories_Combined.png"), 
       plot = final_plot, 
       width = 17, 
       height = 20, 
       units = "cm", 
       dpi = 600, 
       bg = "white")

cat("Done! Full-page publication figure saved to: Figure_XPCLR_Trajectories_Combined.png\n")
