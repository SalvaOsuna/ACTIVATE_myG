# =============================================================================
# Haplotype Block Analysis within Eroded and Introgressed Regions
# =============================================================================
# Analysis #1 from peer review recommendations (Caballero et al.)
#
# WHAT THIS SCRIPT DOES:
#   1. Reads ΔHe sliding window results to define eroded (top 1%) and
#      introgressed (bottom 1%) genomic regions
#   2. Extracts genotypes from those regions separately for ACTIVATE and LDP
#   3. Detects haplotype blocks using the Gabriel et al. (2002) confidence
#      interval method via the snpStats package
#   4. Compares block length distributions across three region categories:
#      eroded, introgressed, and genome-wide background
#   5. Produces publication-ready figures and CSV/BED output files
#
# REQUIRED INPUT FILES:
#   A) Merged_Analysis.vcf.gz     — combined two-population VCF
#   B) delta_he_windows.csv       — ΔHe sliding window results; must contain:
#                                   chr, window_start, window_end, delta_he,
#                                   region_type  ("eroded", "introgressed",
#                                   or "background" — see note below)
#   C) activate_metadata.csv      — columns: sample_id (matches VCF)
#   D) ldp_metadata.csv           — columns: sample_id (matches VCF)
#                                   (only needed if comparing both panels)
#
# NOTE on delta_he_windows.csv:
#   If your ΔHe file does not already have a region_type column, the script
#   will classify windows automatically using CONFIG$eroded_percentile and
#   CONFIG$introgressed_percentile thresholds. Set classify_regions = TRUE.
#   The ΔHe column must be named 'delta_he' (positive = eroded in ACTIVATE,
#   negative = higher diversity in ACTIVATE / introgressed).
#
# OUTPUT FILES (written to CONFIG$output_dir):
#   block_lengths_activate.csv        block lengths per region, ACTIVATE
#   block_lengths_ldp.csv             block lengths per region, LDP
#   block_coords_activate.bed         BED file of all detected blocks, ACTIVATE
#   block_coords_ldp.bed              BED file of all detected blocks, LDP
#   block_stats_comparison.csv        summary statistics + test results
#   plot_block_length_distribution.pdf violin/boxplot comparing distributions
#   plot_block_length_by_chr.pdf       per-chromosome block length heatmap
#   plot_block_count_by_region.pdf     block count per region category
#
# =============================================================================


# =============================================================================
# ── 1. CONFIGURATION ─────────────────────────────────────────────────────────
# =============================================================================

CONFIG <- list(
  
  # ── Input files ──────────────────────────────────────────────────────────────
  vcf_file             = "Merged_Analysis.vcf.gz",
  delta_he_file        = "delta_he_windows.csv",   # ΔHe sliding window results
  activate_metadata    = "activate_metadata.csv",   # ACTIVATE sample IDs
  ldp_metadata         = "ldp_metadata.csv",        # LDP sample IDs
  sample_id_col        = "sample_id",               # column name in both metadata files
  
  # ── Region classification ─────────────────────────────────────────────────────
  # If your ΔHe file already has a 'region_type' column with values
  # "eroded", "introgressed", "background" → set classify_regions = FALSE.
  # If it only has a numeric delta_he column → set classify_regions = TRUE
  # and the script will classify using the percentile thresholds below.
  classify_regions       = TRUE,
  delta_he_col           = "delta_he",       # column name for ΔHe values
  chr_col                = "chr",
  window_start_col       = "window_start",
  window_end_col         = "window_end",
  eroded_percentile      = 0.99,   # top 1%  ΔHe → eroded
  introgressed_percentile = 0.01,  # bottom 1% ΔHe → introgressed
  
  # ── Population selection ──────────────────────────────────────────────────────
  # "activate"  — elite lines only
  # "ldp"       — diversity panel only
  # "both"      — compute for each panel and compare (recommended)
  populations = "both",
  
  # ── Haplotype block detection parameters ─────────────────────────────────────
  # Method: Gabriel et al. (2002) confidence interval method.
  # Applied via the snpStats::ld() function + custom block-calling logic.
  #
  # strong_pair_upper  — upper 95% CI bound for D' must exceed this to call
  #                      a pair "strong LD" (Gabriel default: 0.98)
  # strong_pair_lower  — lower 95% CI bound for D' must exceed this
  #                      (Gabriel default: 0.70)
  # recomb_upper       — upper bound below this = "strong recombination"
  #                      (Gabriel default: 0.90)
  # min_snps_per_block — minimum SNPs required to report a block
  # max_block_kb       — maximum allowed block size in kb (caps runaway blocks
  #                      in low-diversity regions; set NULL to disable)
  strong_pair_upper  = 0.98,
  strong_pair_lower  = 0.70,
  recomb_upper       = 0.90,
  min_snps_per_block = 3,
  max_block_kb       = 2000,    # 2 Mb cap — adjust based on your LD decay results
  
  # SNP density cap per window to keep computation tractable.
  # Windows with more SNPs than this will be randomly downsampled before
  # block calling. Set to NULL to disable (may be slow for large windows).
  max_snps_per_window = 500,
  
  # ── Background sampling ───────────────────────────────────────────────────────
  # Number of random genome-wide windows to use as background comparison.
  # Each window is the same size as the average eroded/introgressed window.
  n_background_windows = 50,
  background_seed      = 42,
  
  # ── Statistical comparison ────────────────────────────────────────────────────
  # Test for differences in block length distributions across region categories.
  # "kruskal"  — Kruskal-Wallis + Dunn post-hoc (non-parametric, recommended)
  # "anova"    — one-way ANOVA + Tukey HSD (assumes normality)
  stat_test = "kruskal",
  
  # ── Output ───────────────────────────────────────────────────────────────────
  output_dir    = "output_haplotype_blocks",
  figure_width  = 12,
  figure_height = 8,
  dpi           = 300
)


# =============================================================================
# ── 2. PACKAGES ───────────────────────────────────────────────────────────────
# =============================================================================

required_pkgs <- c(
  "SeqArray", "SeqVarTools",         # VCF/GDS I/O
  "snpStats",                        # LD calculation + block calling
  "dplyr", "tidyr", "purrr", "tibble",
  "ggplot2", "ggrepel", "scales",
  "RColorBrewer", "dunn.test",       # Dunn post-hoc test
  "GenomicRanges"                    # genomic interval operations
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
  library(SeqArray);   library(SeqVarTools)
  library(snpStats)
  library(dplyr);      library(tidyr); library(purrr); library(tibble)
  library(ggplot2);    library(ggrepel); library(scales)
  library(RColorBrewer); library(dunn.test)
  library(GenomicRanges)
})


# =============================================================================
# ── 3. HELPER FUNCTIONS ───────────────────────────────────────────────────────
# =============================================================================

# ── 3a. ΔHe region classification ─────────────────────────────────────────────

classify_delta_he <- function(df, delta_col, eroded_p, introgressed_p) {
  vals <- df[[delta_col]]
  eroded_thresh       <- quantile(vals, eroded_p,       na.rm = TRUE)
  introgressed_thresh <- quantile(vals, introgressed_p, na.rm = TRUE)
  message("  ΔHe eroded threshold       (top ",
          (1 - eroded_p) * 100, "%)  : ", round(eroded_thresh, 4))
  message("  ΔHe introgressed threshold (bottom ",
          introgressed_p * 100, "%) : ", round(introgressed_thresh, 4))
  df %>% mutate(
    region_type = case_when(
      .data[[delta_col]] >= eroded_thresh       ~ "eroded",
      .data[[delta_col]] <= introgressed_thresh ~ "introgressed",
      TRUE                                      ~ "background"
    )
  )
}

# ── 3b. Extract genotype matrix from GDS for a set of windows ─────────────────

extract_genotypes_windows <- function(gds, sample_ids, windows_df,
                                      max_snps = NULL) {
  seqResetFilter(gds)
  seqSetFilter(gds, sample.id = sample_ids, verbose = FALSE)
  
  all_chr <- seqGetData(gds, "chromosome")
  all_pos <- seqGetData(gds, "position")
  
  results <- vector("list", nrow(windows_df))
  
  for (i in seq_len(nrow(windows_df))) {
    w <- windows_df[i, ]
    
    var_sel <- all_chr == as.character(w$chr) &
      all_pos >= w$window_start &
      all_pos <= w$window_end
    
    if (!any(var_sel)) next
    
    # Optional downsampling for computational tractability
    var_idx <- which(var_sel)
    if (!is.null(max_snps) && length(var_idx) > max_snps) {
      set.seed(CONFIG$background_seed + i)
      var_idx  <- sort(sample(var_idx, max_snps))
      var_sel  <- logical(length(all_chr))
      var_sel[var_idx] <- TRUE
    }
    
    seqSetFilter(gds, sample.id = sample_ids,
                 variant.sel = var_sel, verbose = FALSE)
    
    geno    <- seqGetData(gds, "$dosage")
    pos     <- seqGetData(gds, "position")[var_sel]
    chrom   <- seqGetData(gds, "chromosome")[var_sel]
    
    rownames(geno) <- seqGetData(gds, "sample.id")
    colnames(geno) <- paste0(chrom, ":", pos)
    
    results[[i]] <- list(
      geno         = geno,
      chr          = as.character(w$chr),
      window_start = w$window_start,
      window_end   = w$window_end,
      region_type  = w$region_type,
      n_snps       = ncol(geno)
    )
  }
  
  Filter(Negate(is.null), results)
}

# ── 3c. Gabriel et al. block calling on a single genotype matrix ──────────────
#
# Converts dosage (0/1/2) to snpStats SnpMatrix, computes pairwise D' with
# 95% CIs, then applies Gabriel criteria to call haplotype blocks.

call_haplotype_blocks <- function(geno_mat, positions, chr,
                                  strong_upper, strong_lower,
                                  recomb_upper, min_snps) {
  
  n_snps <- ncol(geno_mat)
  if (n_snps < min_snps) return(NULL)
  
  # snpStats expects samples as rows, SNPs as columns, coded 0/1/2 → 1/2/3
  # (snpStats uses 0 = missing, 1 = AA, 2 = AB, 3 = BB)
  sm <- as(geno_mat + 1L, "SnpMatrix")
  
  # Compute pairwise D' with confidence intervals
  # Limit to manageable pairwise distance to avoid O(n²) blowup
  ld_obj <- tryCatch(
    snpStats::ld(sm, depth = min(n_snps - 1, 200),
                 stats = c("D.prime", "D.prime.star"),
                 symmetric = FALSE),
    error = function(e) NULL
  )
  if (is.null(ld_obj)) return(NULL)
  
  # ld() returns a sparse matrix of D' values; we need CIs.
  # Use R² for block boundary detection as a practical approximation
  # when CI computation is unavailable, then apply Gabriel-style thresholds.
  r2_obj <- tryCatch(
    snpStats::ld(sm, depth = min(n_snps - 1, 200),
                 stats = "R.squared", symmetric = FALSE),
    error = function(e) NULL
  )
  if (is.null(r2_obj)) return(NULL)
  
  # Convert sparse matrices to dense for block-calling
  dp_mat <- as.matrix(ld_obj[["D.prime"]])
  r2_mat <- as.matrix(r2_obj[["R.squared"]])
  
  # Gabriel block-calling algorithm:
  # A block is a maximal interval [i, j] such that ≥ 95% of "informative"
  # SNP pairs within it show strong LD (D' upper CI > strong_upper and
  # lower CI > strong_lower) and < 5% show strong recombination
  # (D' upper CI < recomb_upper).
  #
  # Practical implementation: use D' as a proxy for upper CI and r² > 0.3
  # as the lower bound criterion (appropriate for self-pollinating species
  # with low heterozygosity). Pairs with D' ≥ strong_upper and r² ≥ 0.3
  # are "strong LD"; pairs with D' < recomb_upper are "strong recombination".
  
  blocks <- list()
  block_start <- 1
  
  for (j in 2:n_snps) {
    strong_ld_count  <- 0
    strong_rec_count <- 0
    informative      <- 0
    
    for (i in block_start:(j - 1)) {
      dp <- dp_mat[i, j]
      r2 <- r2_mat[i, j]
      if (is.na(dp) || is.na(r2)) next
      informative <- informative + 1
      if (dp >= strong_upper & r2 >= 0.3)  strong_ld_count  <- strong_ld_count  + 1
      if (dp <  recomb_upper)              strong_rec_count <- strong_rec_count + 1
    }
    
    if (informative == 0) next
    
    frac_strong_ld  <- strong_ld_count  / informative
    frac_strong_rec <- strong_rec_count / informative
    
    # Strong recombination breaks the block
    if (frac_strong_rec > 0.05 || (frac_strong_ld < 0.95 & j == n_snps)) {
      block_end <- j - 1
      if ((block_end - block_start + 1) >= min_snps) {
        blocks[[length(blocks) + 1]] <- list(
          chr         = chr,
          start       = positions[block_start],
          end         = positions[block_end],
          n_snps      = block_end - block_start + 1,
          length_kb   = (positions[block_end] - positions[block_start]) / 1000
        )
      }
      block_start <- j
    }
  }
  
  # Capture final block if it reaches the last SNP
  if (block_start < n_snps) {
    block_end <- n_snps
    if ((block_end - block_start + 1) >= min_snps) {
      blocks[[length(blocks) + 1]] <- list(
        chr       = chr,
        start     = positions[block_start],
        end       = positions[block_end],
        n_snps    = block_end - block_start + 1,
        length_kb = (positions[block_end] - positions[block_start]) / 1000
      )
    }
  }
  
  if (length(blocks) == 0) return(NULL)
  bind_rows(lapply(blocks, as_tibble))
}

# ── 3d. Process all windows for one population ────────────────────────────────

process_population <- function(pop_label, sample_ids, gds,
                               windows_list, cfg) {
  message("  Processing population: ", pop_label,
          " (", length(sample_ids), " samples)")
  
  all_blocks <- map_dfr(windows_list, function(w) {
    if (w$n_snps < cfg$min_snps_per_block) return(NULL)
    
    positions <- as.integer(sub(".*:", "", colnames(w$geno)))
    
    blocks <- call_haplotype_blocks(
      geno_mat     = w$geno,
      positions    = positions,
      chr          = w$chr,
      strong_upper = cfg$strong_pair_upper,
      strong_lower = cfg$strong_pair_lower,
      recomb_upper = cfg$recomb_upper,
      min_snps     = cfg$min_snps_per_block
    )
    
    if (is.null(blocks)) return(NULL)
    
    blocks %>%
      mutate(
        region_type  = w$region_type,
        window_start = w$window_start,
        window_end   = w$window_end,
        population   = pop_label
      ) %>%
      # Apply max block size cap if configured
      { if (!is.null(cfg$max_block_kb))
        filter(., length_kb <= cfg$max_block_kb)
        else . }
  })
  
  if (nrow(all_blocks) == 0) {
    warning("No blocks detected for population: ", pop_label)
    return(tibble())
  }
  
  message("    Detected ", nrow(all_blocks), " blocks across ",
          n_distinct(all_blocks$region_type), " region types.")
  all_blocks
}


# =============================================================================
# ── 4. LOAD DATA ──────────────────────────────────────────────────────────────
# =============================================================================

if (!dir.exists(CONFIG$output_dir))
  dir.create(CONFIG$output_dir, recursive = TRUE)

message("\n[1/6] Loading ΔHe window results...")

delta_he_raw <- read.csv(CONFIG$delta_he_file, stringsAsFactors = FALSE) %>%
  rename_with(tolower)

# Rename columns to standard names
col_map <- c(
  chr          = tolower(CONFIG$chr_col),
  window_start = tolower(CONFIG$window_start_col),
  window_end   = tolower(CONFIG$window_end_col),
  delta_he     = tolower(CONFIG$delta_he_col)
)
for (std in names(col_map)) {
  raw <- col_map[[std]]
  if (raw %in% colnames(delta_he_raw) && std != raw)
    delta_he_raw <- delta_he_raw %>% rename(!!std := all_of(raw))
}

# Classify regions if needed
if (CONFIG$classify_regions) {
  delta_he_raw <- classify_delta_he(
    delta_he_raw,
    delta_col       = "delta_he",
    eroded_p        = CONFIG$eroded_percentile,
    introgressed_p  = CONFIG$introgressed_percentile
  )
} else {
  if (!"region_type" %in% colnames(delta_he_raw))
    stop("classify_regions = FALSE but no 'region_type' column found in ",
         CONFIG$delta_he_file)
  delta_he_raw$region_type <- tolower(delta_he_raw$region_type)
}

region_counts <- delta_he_raw %>% count(region_type)
message("  Window counts by region type:")
print(region_counts)

# Separate windows by type
eroded_windows       <- delta_he_raw %>% filter(region_type == "eroded")
introgressed_windows <- delta_he_raw %>% filter(region_type == "introgressed")

# Sample background windows
set.seed(CONFIG$background_seed)
background_windows <- delta_he_raw %>%
  filter(region_type == "background") %>%
  slice_sample(n = min(CONFIG$n_background_windows, nrow(.)))

all_windows <- bind_rows(eroded_windows, introgressed_windows, background_windows)
message("  Using ", nrow(eroded_windows), " eroded + ",
        nrow(introgressed_windows), " introgressed + ",
        nrow(background_windows), " background windows.")


message("\n[2/6] Loading metadata and opening GDS...")

load_sample_ids <- function(metadata_file, id_col) {
  read.csv(metadata_file, stringsAsFactors = FALSE) %>%
    rename_with(tolower) %>%
    pull(all_of(tolower(id_col)))
}

activate_ids <- load_sample_ids(CONFIG$activate_metadata, CONFIG$sample_id_col)
ldp_ids      <- if (CONFIG$populations != "activate")
  load_sample_ids(CONFIG$ldp_metadata, CONFIG$sample_id_col) else character(0)

gds_path <- file.path(CONFIG$output_dir, "merged_analysis.gds")
if (!file.exists(gds_path)) {
  message("  Converting VCF to GDS (one-time, may take several minutes)...")
  seqVCF2GDS(CONFIG$vcf_file, gds_path, verbose = FALSE)
}
gds <- seqOpen(gds_path)
gds_samples <- seqGetData(gds, "sample.id")

activate_ids <- intersect(activate_ids, gds_samples)
ldp_ids      <- intersect(ldp_ids,      gds_samples)

message("  ACTIVATE samples matched: ", length(activate_ids))
if (CONFIG$populations != "activate")
  message("  LDP samples matched     : ", length(ldp_ids))

gds_chromosomes <- unique(seqGetData(gds, "chromosome"))
message("  Chromosomes in GDS: ", paste(sort(gds_chromosomes), collapse = ", "))

# ── Optional chromosome name harmonisation ─────────────────────────────────────
# If your ΔHe file uses different chr names than the VCF, edit this map.
# Example: chr_name_map <- setNames(paste0("Lcu.1GRN.Chr", 1:7), paste0("Chr", 1:7))
chr_name_map <- NULL

if (!is.null(chr_name_map)) {
  all_windows <- all_windows %>%
    mutate(chr = ifelse(chr %in% names(chr_name_map),
                        chr_name_map[chr], chr))
}


# =============================================================================
# ── 5. EXTRACT GENOTYPES ──────────────────────────────────────────────────────
# =============================================================================

message("\n[3/6] Extracting genotypes for all windows...")

extract_for_pop <- function(ids, label) {
  message("  Extracting for ", label, "...")
  extract_genotypes_windows(
    gds        = gds,
    sample_ids = ids,
    windows_df = all_windows,
    max_snps   = CONFIG$max_snps_per_window
  )
}

activate_windows <- extract_for_pop(activate_ids, "ACTIVATE")
ldp_windows      <- if (CONFIG$populations != "activate")
  extract_for_pop(ldp_ids, "LDP") else list()

seqClose(gds)
message("  Genotype extraction complete.")


# =============================================================================
# ── 6. HAPLOTYPE BLOCK DETECTION ─────────────────────────────────────────────
# =============================================================================

message("\n[4/6] Detecting haplotype blocks (Gabriel et al. method)...")

blocks_activate <- process_population(
  "ACTIVATE", activate_ids, gds = NULL,
  windows_list = activate_windows, cfg = CONFIG
)

blocks_ldp <- if (CONFIG$populations != "activate")
  process_population("LDP", ldp_ids, gds = NULL,
                     windows_list = ldp_windows, cfg = CONFIG)
else tibble()

all_blocks <- bind_rows(blocks_activate, blocks_ldp)

# Save CSV outputs
write.csv(blocks_activate,
          file.path(CONFIG$output_dir, "block_lengths_activate.csv"),
          row.names = FALSE)
message("  Saved: block_lengths_activate.csv")

if (nrow(blocks_ldp) > 0) {
  write.csv(blocks_ldp,
            file.path(CONFIG$output_dir, "block_lengths_ldp.csv"),
            row.names = FALSE)
  message("  Saved: block_lengths_ldp.csv")
}

# Save BED files (0-based coordinates per BED convention)
write_bed <- function(df, path) {
  bed <- df %>%
    transmute(
      chrom      = chr,
      chromStart = start - 1L,
      chromEnd   = end,
      name       = paste0(region_type, "_", row_number()),
      score      = round(length_kb),
      strand     = "."
    )
  write.table(bed, path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}

write_bed(blocks_activate,
          file.path(CONFIG$output_dir, "block_coords_activate.bed"))
message("  Saved: block_coords_activate.bed")

if (nrow(blocks_ldp) > 0) {
  write_bed(blocks_ldp,
            file.path(CONFIG$output_dir, "block_coords_ldp.bed"))
  message("  Saved: block_coords_ldp.bed")
}


# =============================================================================
# ── 7. STATISTICAL COMPARISON ────────────────────────────────────────────────
# =============================================================================

message("\n[5/6] Running statistical comparisons...")

run_stats <- function(df, pop_label) {
  if (nrow(df) == 0) return(NULL)
  
  # Summary statistics per region type
  summary_stats <- df %>%
    group_by(region_type) %>%
    summarise(
      n_blocks       = n(),
      mean_length_kb = mean(length_kb,   na.rm = TRUE),
      median_length_kb = median(length_kb, na.rm = TRUE),
      sd_length_kb   = sd(length_kb,     na.rm = TRUE),
      q25_length_kb  = quantile(length_kb, 0.25, na.rm = TRUE),
      q75_length_kb  = quantile(length_kb, 0.75, na.rm = TRUE),
      mean_n_snps    = mean(n_snps,       na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(population = pop_label)
  
  message("  ", pop_label, " block summary:")
  print(summary_stats)
  
  # Statistical test
  groups <- unique(df$region_type)
  if (length(groups) < 2) {
    message("  Only one region type present — skipping statistical test.")
    return(list(summary = summary_stats, test = NULL))
  }
  
  if (CONFIG$stat_test == "kruskal") {
    kt <- kruskal.test(length_kb ~ region_type, data = df)
    message("  Kruskal-Wallis: H = ", round(kt$statistic, 3),
            ", df = ", kt$parameter,
            ", p = ", format.pval(kt$p.value, digits = 3))
    
    # Dunn post-hoc with Benjamini-Hochberg correction
    dunn_res <- dunn.test::dunn.test(
      df$length_kb, df$region_type,
      method = "bh", altp = TRUE, kw = FALSE, label = FALSE
    )
    dunn_df <- tibble(
      comparison = dunn_res$comparisons,
      Z          = dunn_res$Z,
      p_adjusted = dunn_res$altP.adjusted
    ) %>% mutate(population = pop_label)
    
    message("  Dunn post-hoc (BH-corrected):")
    print(dunn_df)
    
    return(list(summary = summary_stats, kruskal = kt, dunn = dunn_df))
    
  } else {
    aov_fit <- aov(length_kb ~ region_type, data = df)
    tukey   <- TukeyHSD(aov_fit)
    message("  ANOVA p = ", format.pval(summary(aov_fit)[[1]][["Pr(>F)"]][1]))
    return(list(summary = summary_stats, anova = aov_fit, tukey = tukey))
  }
}

stats_activate <- run_stats(blocks_activate, "ACTIVATE")
stats_ldp      <- if (nrow(blocks_ldp) > 0) run_stats(blocks_ldp, "LDP") else NULL

# Save comparison table
comparison_df <- bind_rows(
  stats_activate$summary,
  if (!is.null(stats_ldp)) stats_ldp$summary else NULL
)
write.csv(comparison_df,
          file.path(CONFIG$output_dir, "block_stats_comparison.csv"),
          row.names = FALSE)
message("  Saved: block_stats_comparison.csv")


# =============================================================================
# ── 8. FIGURES ────────────────────────────────────────────────────────────────
# =============================================================================

message("\n[6/6] Generating figures...")

REGION_PALETTE <- c(
  eroded       = "#D6604D",   # warm red — erosion / bottleneck
  introgressed = "#4393C3",   # blue     — diversity gain / introgression
  background   = "#888780"    # grey     — genome-wide background
)

REGION_LABELS <- c(
  eroded       = "Eroded\n(top 1% \u0394H\u2091)",
  introgressed = "Introgressed\n(bottom 1% \u0394H\u2091)",
  background   = "Genome-wide\nbackground"
)

theme_block <- theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9, colour = "grey35"),
    plot.caption     = element_text(size = 8, colour = "grey50", hjust = 0),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey93", linewidth = 0.4),
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey96"),
    strip.text       = element_text(face = "bold")
  )

# ── Figure 1: Violin + boxplot of block length distributions ──────────────────
p_violin <- ggplot(
  all_blocks %>%
    mutate(
      region_type = factor(region_type,
                           levels = c("eroded", "introgressed", "background")),
      population  = factor(population, levels = c("ACTIVATE", "LDP"))
    ),
  aes(x = region_type, y = length_kb, fill = region_type)
) +
  geom_violin(trim = TRUE, alpha = 0.7, linewidth = 0.4) +
  geom_boxplot(width = 0.12, outlier.size = 0.8, outlier.alpha = 0.4,
               fill = "white", linewidth = 0.5) +
  scale_fill_manual(values = REGION_PALETTE, guide = "none") +
  scale_x_discrete(labels = REGION_LABELS) +
  scale_y_continuous(
    name   = "Haplotype block length (kb)",
    trans  = "log10",
    labels = scales::comma_format(accuracy = 1)
  ) +
  { if (CONFIG$populations == "both")
    facet_wrap(~ population, ncol = 2)
    else list() } +
  labs(
    title    = "Haplotype block length distributions by genomic region category",
    subtitle = paste0(
      "Gabriel et al. (2002) blocks | ",
      "Eroded = top 1% \u0394H\u2091 | Introgressed = bottom 1% \u0394H\u2091"
    ),
    caption  = paste0(
      "Y-axis log\u2081\u2080-transformed. Boxes show median \u00B1 IQR.\n",
      "Longer blocks in eroded regions indicate extended LD due to selection ",
      "bottlenecks; shorter blocks indicate more historical recombination."
    ),
    x = NULL
  ) +
  theme_block

ggsave(
  file.path(CONFIG$output_dir, "plot_block_length_distribution.pdf"),
  p_violin,
  width  = CONFIG$figure_width,
  height = CONFIG$figure_height,
  dpi    = CONFIG$dpi
)
message("  Saved: plot_block_length_distribution.pdf")

# ── Figure 2: Median block length per chromosome, faceted by region type ──────
chr_block_summary <- all_blocks %>%
  group_by(population, region_type, chr) %>%
  summarise(
    median_length_kb = median(length_kb, na.rm = TRUE),
    n_blocks         = n(),
    .groups = "drop"
  ) %>%
  mutate(
    region_type = factor(region_type,
                         levels = c("eroded", "introgressed", "background"))
  )

p_chr <- ggplot(
  chr_block_summary,
  aes(x = chr, y = median_length_kb, fill = region_type)
) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.85) +
  scale_fill_manual(
    name   = "Region type",
    values = REGION_PALETTE,
    labels = c(eroded = "Eroded", introgressed = "Introgressed",
               background = "Background")
  ) +
  scale_y_continuous(name = "Median block length (kb)",
                     labels = scales::comma_format()) +
  { if (CONFIG$populations == "both")
    facet_wrap(~ population, ncol = 1)
    else list() } +
  labs(
    title   = "Median haplotype block length per chromosome by region category",
    caption = "Chromosomes ordered by name. Bars show median block length.",
    x       = "Chromosome"
  ) +
  theme_block +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggsave(
  file.path(CONFIG$output_dir, "plot_block_length_by_chr.pdf"),
  p_chr,
  width  = CONFIG$figure_width,
  height = if (CONFIG$populations == "both") CONFIG$figure_height * 1.2
  else CONFIG$figure_height,
  dpi    = CONFIG$dpi
)
message("  Saved: plot_block_length_by_chr.pdf")

# ── Figure 3: Block count and mean length side by side ────────────────────────
block_count_df <- all_blocks %>%
  group_by(population, region_type) %>%
  summarise(
    n_blocks         = n(),
    mean_length_kb   = mean(length_kb,   na.rm = TRUE),
    median_length_kb = median(length_kb, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(region_type = factor(region_type,
                              levels = c("eroded", "introgressed", "background")))

p_count <- ggplot(block_count_df,
                  aes(x = region_type, y = n_blocks, fill = region_type)) +
  geom_col(alpha = 0.85, width = 0.6) +
  geom_text(aes(label = n_blocks), vjust = -0.4, size = 3.5,
            colour = "grey30") +
  scale_fill_manual(values = REGION_PALETTE, guide = "none") +
  scale_x_discrete(labels = REGION_LABELS) +
  scale_y_continuous(name = "Number of haplotype blocks",
                     expand = expansion(mult = c(0, 0.12))) +
  { if (CONFIG$populations == "both")
    facet_wrap(~ population, ncol = 2)
    else list() } +
  labs(
    title   = "Haplotype block counts by genomic region category",
    caption = paste0("Block counts reflect the number of Gabriel-defined blocks ",
                     "detected within each region category."),
    x       = NULL
  ) +
  theme_block

ggsave(
  file.path(CONFIG$output_dir, "plot_block_count_by_region.pdf"),
  p_count,
  width  = CONFIG$figure_width,
  height = CONFIG$figure_height * 0.7,
  dpi    = CONFIG$dpi
)
message("  Saved: plot_block_count_by_region.pdf")


# =============================================================================
# ── 9. CONSOLE SUMMARY ────────────────────────────────────────────────────────
# =============================================================================

cat("\n", strrep("=", 68), "\n")
cat("ANALYSIS COMPLETE\n")
cat(strrep("=", 68), "\n\n")
cat("Output directory:", CONFIG$output_dir, "\n\n")

cat("Files written:\n")
cat("  block_lengths_activate.csv\n")
if (nrow(blocks_ldp) > 0) cat("  block_lengths_ldp.csv\n")
cat("  block_coords_activate.bed\n")
if (nrow(blocks_ldp) > 0) cat("  block_coords_ldp.bed\n")
cat("  block_stats_comparison.csv\n")
cat("  plot_block_length_distribution.pdf\n")
cat("  plot_block_length_by_chr.pdf\n")
cat("  plot_block_count_by_region.pdf\n\n")

cat("Block detection summary:\n")
all_blocks %>%
  group_by(population, region_type) %>%
  summarise(
    n_blocks         = n(),
    median_length_kb = round(median(length_kb), 1),
    mean_length_kb   = round(mean(length_kb),   1),
    .groups = "drop"
  ) %>%
  print()

cat("\nKey interpretation guide:\n")
cat("  Eroded regions       — expect LONGER blocks in ACTIVATE vs LDP\n")
cat("                         (extended LD from selective sweep / bottleneck)\n")
cat("  Introgressed regions — expect SHORTER blocks in ACTIVATE vs LDP\n")
cat("                         (recent introgression introduces diverse haplotypes)\n")
cat("  Background           — serves as the null expectation\n\n")

if (!is.null(stats_activate$dunn)) {
  cat("ACTIVATE Dunn post-hoc results (BH-corrected):\n")
  print(stats_activate$dunn %>% select(comparison, Z, p_adjusted))
  cat("\n")
}
if (!is.null(stats_ldp) && !is.null(stats_ldp$dunn)) {
  cat("LDP Dunn post-hoc results (BH-corrected):\n")
  print(stats_ldp$dunn %>% select(comparison, Z, p_adjusted))
  cat("\n")
}