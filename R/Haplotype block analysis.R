# =============================================================================
# Haplotype Block Analysis — Eroded vs Introgressed Genomic Regions
# =============================================================================
# Analysis #1 from peer review recommendations (Caballero et al.)
#
# WHAT THIS SCRIPT DOES:
#   1. Splits Merged_Analysis.vcf.gz into ACTIVATE and LDP sub-VCFs
#   2. Converts each to PLINK binary format (.bed/.bim/.fam)
#   3. Calls haplotype blocks in each panel using the Gabriel et al. (2002)
#      algorithm (PLINK --blocks)
#   4. Classifies blocks as falling in: eroded regions (top 1% ΔHe),
#      introgressed regions (bottom 1% ΔHe), or background (neither)
#   5. Compares block lengths across region types with Mann-Whitney U tests
#   6. Produces publication-ready ggplot2 figures and CSV outputs
#
# REQUIRED INPUT FILES:
#   A) Merged_Analysis.vcf.gz                      — merged two-population VCF
#   B) data/ACT187_samples.txt                     — plain text file, one ACTIVATE
#                                                    sample ID per line, no header
#   C) data/LDP324_samples.txt                     — plain text file, one LDP
#                                                    sample ID per line, no header
#   D) DeltaHe_Erosion_SignificantRegions.csv       — pre-filtered eroded windows;
#                                                    must contain at minimum:
#                                                    chr, window_start, window_end
#   E) DeltaHe_Introgression_SignificantRegions.csv — pre-filtered introgressed
#                                                    windows; same columns required
#
# REQUIRED EXTERNAL TOOL:
#   PLINK v1.9  — must be installed and accessible on PATH
#   Install: https://www.cog-genomics.org/plink/
#   Check:   system(plink_bin, " --version")
#
# OUTPUT FILES (written to CONFIG$output_dir):
#   activate_blocks.blocks.det     raw PLINK block output, ACTIVATE panel
#   ldp_blocks.blocks.det          raw PLINK block output, LDP panel
#   blocks_annotated_activate.csv  blocks with region classification, ACTIVATE
#   blocks_annotated_ldp.csv       blocks with region classification, LDP
#   mannwhitney_results.csv        test statistics for all comparisons
#   plot_block_lengths_violin.pdf  violin + boxplot by region type and panel
#   plot_block_lengths_ecdf.pdf    empirical CDF curves per region type
#   plot_block_density_genome.pdf  genome-wide block density track
#
# =============================================================================


# =============================================================================
# ── 1. CONFIGURATION ─────────────────────────────────────────────────────────
# =============================================================================

CONFIG <- list(
  
  # ── Input files ──────────────────────────────────────────────────────────────
  vcf_file              = "data/Merged_Analysis.vcf.gz",
  
  # Plain-text sample ID files — one ID per line, no header, read with readLines
  activate_ids_file     = "data/ACT187_samples.txt",
  ldp_ids_file          = "data/LDP324_samples.txt",
  
  # Pre-filtered ΔHe significant region files — no thresholding needed,
  # these already contain only the top/bottom 1% windows from your analysis.
  # Required columns: chr, window_start, window_end  (additional cols ignored)
  delta_he_erosion_file        = "Results/DeltaHe_Erosion_SignificantRegions.csv",
  delta_he_introgression_file  = "Results/DeltaHe_Introgression_SignificantRegions.csv",
  
  # ── PLINK block-calling parameters ─────────────────────────────────────────
  # Gabriel et al. 2002 algorithm (PLINK default --blocks)
  # Adjust these if your LD decay results suggest different window sizes
  plink_blocks_max_kb   = 3000,   # maximum block span in kb (3 Mb default)
  plink_maf             = 0.05,   # MAF filter for PLINK (matches your VCF QC)
  plink_geno            = 0.10,   # max missingness per SNP
  plink_mind            = 0.10,   # max missingness per sample
  
  # Chromosome name handling:
  # PLINK expects numeric or "chrN" chromosome names.
  # If your VCF uses "Lcu.1GRN.Chr1" format, set recode_chr = TRUE and
  # provide a mapping from VCF names to PLINK-compatible names.
  recode_chr = TRUE,
  chr_map    = setNames(
    as.character(1:7),                                    # PLINK names: 1–7
    paste0("Lcu.1GRN.Chr", 1:7)                          # VCF names
  ),
  
  # ── Statistical test ─────────────────────────────────────────────────────────
  # Minimum number of blocks required in a region category to run the test
  min_blocks_for_test = 10,
  
  # ── Output ──────────────────────────────────────────────────────────────────
  output_dir    = "output_haplotype_blocks",
  figure_width  = 12,
  figure_height = 8,
  dpi           = 300
)


# =============================================================================
# ── 2. PACKAGES ───────────────────────────────────────────────────────────────
# =============================================================================

required_pkgs <- c(
  "vcfR",                              # pure-R VCF reading and subsetting (no bcftools)
  "dplyr", "tidyr", "purrr", "tibble",
  "GenomicRanges", "IRanges",          # interval overlap: blocks ↔ ΔHe windows
  "ggplot2", "scales",
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
  library(vcfR)
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(GenomicRanges); library(IRanges)
  library(ggplot2); library(scales)
  library(RColorBrewer); library(broom)
})


# =============================================================================
# ── 3. HELPER FUNCTIONS ───────────────────────────────────────────────────────
# =============================================================================

# ── 3a. Run a system command and stop on failure ──────────────────────────────
run_cmd <- function(cmd, step = "") {
  message("  $ ", cmd)
  ret <- system(cmd)
  if (ret != 0)
    stop("Command failed (exit code ", ret, ")",
         if (nchar(step) > 0) paste0(" at step: ", step) else "")
  invisible(ret)
}

# ── 3b. Parse PLINK .blocks.det file ─────────────────────────────────────────
# PLINK --blocks det output columns:
#   CHR   BP1   BP2   KB   NSNPS   SNPS
parse_plink_blocks <- function(path, panel_name) {
  if (!file.exists(path))
    stop("PLINK block file not found: ", path,
         "\n  Check that PLINK ran successfully and output_dir is correct.")
  
  read.table(path, header = TRUE, stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    rename_with(tolower) %>%
    dplyr::rename(
      chr          = chr,
      block_start  = bp1,
      block_end    = bp2,
      block_kb     = kb,
      n_snps_block = nsnps
    ) %>%
    mutate(
      chr        = as.character(chr),
      block_len_bp = as.integer(block_end - block_start),
      panel      = panel_name
    ) %>%
    select(panel, chr, block_start, block_end, block_len_bp, block_kb, n_snps_block)
}

# ── 3c. Split VCF and recode chromosome names using vcfR (pure R, no bcftools)
# Reads the full merged VCF once, subsets to the requested sample IDs,
# optionally recodes chromosome names for PLINK compatibility, and writes
# a plain (uncompressed) VCF that PLINK can read directly on Windows.
#
# NOTE: vcfR loads the full VCF into memory. For very large VCFs (>5 GB)
# this may be slow. The result is cached — if the output VCF already exists
# the function returns immediately without re-reading the input.
split_vcf_r <- function(vcf_in, sample_ids, chr_map, output_dir, tag) {
  out_vcf <- file.path(output_dir, paste0(tag, "_panel.vcf"))
  
  if (file.exists(out_vcf)) {
    message("  Panel VCF for ", tag, " already exists — skipping split.")
    return(out_vcf)
  }
  
  message("  Reading merged VCF into R (this may take a few minutes)...")
  vcf_obj <- read.vcfR(vcf_in, verbose = FALSE)
  
  # Subset to requested samples
  gt_cols    <- colnames(vcf_obj@gt)          # first col is "FORMAT"
  sample_col <- gt_cols[-1]                   # actual sample names
  keep_cols  <- sample_col[sample_col %in% sample_ids]
  
  missing <- setdiff(sample_ids, sample_col)
  if (length(missing) > 0)
    warning(length(missing), " sample IDs from ", tag,
            " not found in VCF and will be skipped.
",
            "  First few missing: ", paste(head(missing, 5), collapse = ", "))
  
  if (length(keep_cols) == 0)
    stop("No matching samples found in VCF for panel ", tag,
         ".
  Check that sample IDs in ", CONFIG[[paste0(tolower(tag), "_ids_file")]],
         " match VCF column names exactly.")
  
  # Subset genotype matrix to FORMAT + matching samples
  vcf_obj@gt <- vcf_obj@gt[, c("FORMAT", keep_cols), drop = FALSE]
  
  # Recode chromosome names for PLINK (Lcu.1GRN.Chr1 → 1)
  if (CONFIG$recode_chr && length(chr_map) > 0) {
    message("  Recoding chromosome names for PLINK compatibility...")
    old_chrs <- vcf_obj@fix[, "CHROM"]
    new_chrs <- ifelse(old_chrs %in% names(chr_map),
                       chr_map[old_chrs], old_chrs)
    vcf_obj@fix[, "CHROM"] <- new_chrs
  }
  
  message("  Writing panel VCF for ", tag,
          " (", length(keep_cols), " samples, ",
          nrow(vcf_obj@fix), " variants)...")
  write.vcf(vcf_obj, file = out_vcf)
  message("  Saved: ", out_vcf)
  out_vcf
}

# ── 3d. Classify a block data frame by ΔHe region type ───────────────────────
# Uses GenomicRanges findOverlaps to intersect block coordinates with
# eroded and introgressed window coordinates.
classify_blocks <- function(blocks_df, eroded_gr, introgressed_gr) {
  
  # HELPER: Strip all text/punctuation from chromosome names (e.g., "Chr1" -> "1")
  # This guarantees PLINK output ("1") perfectly matches your CSVs ("Chr1" or "Lcu...Chr1")
  extract_chr_num <- function(chr_vec) {
    gsub("[^0-9]", "", as.character(chr_vec))
  }
  
  # 1. Build GRanges for blocks using purely numeric chromosome names
  blocks_gr <- GRanges(
    seqnames = extract_chr_num(blocks_df$chr),
    ranges   = IRanges(start = blocks_df$block_start,
                       end   = blocks_df$block_end)
  )
  
  # 2. Safely standardize the target GRanges chromosome names to pure numbers
  eroded_gr_norm <- eroded_gr
  seqlevels(eroded_gr_norm) <- extract_chr_num(seqlevels(eroded_gr_norm))
  
  introgressed_gr_norm <- introgressed_gr
  seqlevels(introgressed_gr_norm) <- extract_chr_num(seqlevels(introgressed_gr_norm))
  
  # 3. Calculate overlaps on the standardized coordinate system
  eroded_hits        <- countOverlaps(blocks_gr, eroded_gr_norm)        > 0
  introgressed_hits  <- countOverlaps(blocks_gr, introgressed_gr_norm)  > 0
  
  # 4. Apply classification back to the original dataframe
  blocks_df %>%
    mutate(
      region_type = case_when(
        eroded_hits & introgressed_hits ~ "both",       # rare edge case
        eroded_hits                     ~ "eroded",
        introgressed_hits               ~ "introgressed",
        TRUE                            ~ "background"
      )
    )
}

# ── 3e. Mann-Whitney U test wrapper ──────────────────────────────────────────
run_mwu <- function(df, group_col, value_col,
                    group_a, group_b, min_n = 10) {
  a <- df[[value_col]][df[[group_col]] == group_a]
  b <- df[[value_col]][df[[group_col]] == group_b]
  if (length(a) < min_n || length(b) < min_n) {
    return(tibble(
      comparison = paste0(group_a, " vs ", group_b),
      n_a = length(a), n_b = length(b),
      statistic = NA_real_, p_value = NA_real_,
      median_a = median(a), median_b = median(b),
      note = paste0("Skipped: n < ", min_n)
    ))
  }
  wt <- wilcox.test(a, b, exact = FALSE, correct = TRUE)
  tibble(
    comparison = paste0(group_a, " vs ", group_b),
    n_a        = length(a),
    n_b        = length(b),
    statistic  = wt$statistic,
    p_value    = wt$p.value,
    median_a   = median(a, na.rm = TRUE),
    median_b   = median(b, na.rm = TRUE),
    note       = ""
  )
}


# =============================================================================
# ── 4. SETUP ──────────────────────────────────────────────────────────────────
# =============================================================================

if (!dir.exists(CONFIG$output_dir))
  dir.create(CONFIG$output_dir, recursive = TRUE)

# ── Locate PLINK genetics tool (Windows) ─────────────────────────────────────
# IMPORTANT: On Windows, PuTTY also installs a plink.exe (an SSH client).
# If Sys.which("plink") returns the PuTTY path, the script will fail with
# "unknown option --vcf". You must point directly to the PLINK genetics tool.
#
# HOW TO SET THE CORRECT PATH:
#   1. Open Command Prompt and run: where plink
#      This lists every plink.exe on your PATH.
#   2. Find the one that is NOT in a PuTTY folder.
#      It should be something like:
#        C:\Users\YourName\plink_win64\plink.exe
#        C:\Tools\plink\plink.exe
#   3. Set it here using forward slashes:
#
options(plink.path = "C:/Users/Salva/Downloads/plink_win64_20250819/plink.exe")
#
# If PLINK genetics is not installed yet, download it from:
#   https://www.cog-genomics.org/plink/  (use the Windows 64-bit zip)
# Extract the zip and note the folder — then set the path above.

is_putty_plink <- function(path) {
  # Run with --version and check if output looks like PuTTY (SSH) or PLINK genetics
  ver <- tryCatch(
    system(paste0('"', path, '" --version'), intern = TRUE, ignore.stderr = TRUE),
    error = function(e) character(0)
  )
  # PuTTY plink prints SSH usage; PLINK genetics prints "PLINK v1.9..."
  length(ver) == 0 || !any(grepl("PLINK v", ver, ignore.case = TRUE))
}

plink_bin <- local({
  # Check for explicit path set via options() first
  opt <- getOption("plink.path", default = NULL)
  if (!is.null(opt) && nzchar(opt)) {
    if (!file.exists(opt))
      stop("PLINK path set in options('plink.path') does not exist:\n  ", opt)
    if (is_putty_plink(opt))
      stop(
        "The plink.exe at the path you set does not appear to be PLINK genetics:\n",
        "  ", opt, "\n",
        "  It may be PuTTY's plink.exe (an SSH client).\n",
        "  Run `where plink` in Command Prompt to find all plink.exe files,\n",
        "  then set options(plink.path = ...) to the PLINK genetics one."
      )
    return(opt)
  }
  
  # Search PATH, skipping any PuTTY installation
  all_plinks <- tryCatch(
    system("where plink", intern = TRUE, ignore.stderr = TRUE),
    error = function(e) character(0)
  )
  # Also try plink.exe explicitly
  all_plinks <- unique(c(all_plinks,
                         Sys.which("plink"), Sys.which("plink.exe")))
  all_plinks <- all_plinks[nzchar(all_plinks) & file.exists(all_plinks)]
  
  if (length(all_plinks) == 0)
    stop(
      "No plink.exe found on PATH.\n",
      "  Download PLINK v1.9 from: https://www.cog-genomics.org/plink/\n",
      "  Then set: options(plink.path = 'C:/path/to/plink.exe')"
    )
  
  # Filter out PuTTY plink
  genetics_plinks <- all_plinks[!sapply(all_plinks, is_putty_plink)]
  
  if (length(genetics_plinks) == 0) {
    putty_paths <- paste(all_plinks, collapse = "\n    ")
    stop(
      "Found plink.exe on PATH, but it appears to be PuTTY (SSH), not PLINK genetics.\n",
      "  PuTTY path(s) found:\n    ", putty_paths, "\n\n",
      "  Fix: set the path to your PLINK genetics executable explicitly:\n",
      "    options(plink.path = 'C:/Users/YourName/plink_win64/plink.exe')\n",
      "  Download PLINK genetics from: https://www.cog-genomics.org/plink/"
    )
  }
  
  if (length(genetics_plinks) > 1)
    message("  Multiple PLINK genetics executables found — using first:\n    ",
            paste(genetics_plinks, collapse = "\n    "))
  
  genetics_plinks[1]
})

plink_ver <- tryCatch(
  system(paste0('"', plink_bin, '" --version'), intern = TRUE,
         ignore.stderr = TRUE),
  error = function(e) "unknown version"
)
message("PLINK (genetics) found: ", plink_bin)
message("PLINK version         : ", plink_ver[1])
message("VCF subsetting        : using vcfR (pure R — no bcftools required)")


# =============================================================================
# ── 5. LOAD METADATA AND SPLIT VCF BY PANEL ──────────────────────────────────
# =============================================================================

message("\n[1/6] Loading sample IDs and splitting VCF by panel...")

# Load sample IDs — plain text, one ID per line, no header
activate_ids <- readLines(CONFIG$activate_ids_file)
ldp_ids      <- readLines(CONFIG$ldp_ids_file)

# Strip any blank lines or whitespace (common in files saved on Windows)
activate_ids <- trimws(activate_ids[nzchar(trimws(activate_ids))])
ldp_ids      <- trimws(ldp_ids[nzchar(trimws(ldp_ids))])

message("  ACTIVATE samples: ", length(activate_ids))
message("  LDP samples     : ", length(ldp_ids))

# Split the merged VCF into per-panel VCFs using vcfR (pure R, no bcftools).
# Output is a plain uncompressed VCF — PLINK on Windows reads these fine.
act_vcf <- split_vcf_r(CONFIG$vcf_file, activate_ids,
                       CONFIG$chr_map, CONFIG$output_dir, "activate")
ldp_vcf <- split_vcf_r(CONFIG$vcf_file, ldp_ids,
                       CONFIG$chr_map, CONFIG$output_dir, "ldp")


# =============================================================================
# ── 6. RECODE CHROMOSOME NAMES AND CONVERT TO PLINK FORMAT ───────────────────
# =============================================================================

message("\n[2/6] Preparing PLINK binary files...")

panels <- list(
  list(vcf = act_vcf, tag = "activate"),
  list(vcf = ldp_vcf, tag = "ldp")
)

for (p in panels) {
  bed_prefix <- file.path(CONFIG$output_dir, paste0(p$tag, "_plink"))
  
  if (!file.exists(paste0(bed_prefix, ".bed"))) {
    message("  Converting ", p$tag, " VCF to PLINK binary format...")
    # Chr names are already recoded to integers by split_vcf_r, so no
    # additional renaming step is needed here.
    # Quoting the paths handles Windows directory paths with spaces.
    run_cmd(
      paste0(
        '"', plink_bin, '"',
        " --vcf ", '"', p$vcf, '"',
        " --double-id",          # use sample ID as both FID and IID
        " --chr-set 7",          # lentil has 7 chromosomes
        " --allow-extra-chr",    # NEW: Prevents crashing on non-integer chromosome names
        " --set-missing-var-ids @:#", # NEW: Assigns unique IDs to SNPs missing RSIDs
        " --vcf-half-call m",    # NEW: Safely handles any weird half-missing genotypes
        " --maf ",  CONFIG$plink_maf,
        " --geno ", CONFIG$plink_geno,
        " --mind ", CONFIG$plink_mind,
        " --make-bed",
        " --out ", '"', bed_prefix, '"'
      ),
      step = paste("VCF to BED:", p$tag)
    )
  } else {
    message("  PLINK BED for ", p$tag, " already exists — skipping conversion.")
  }
}


# =============================================================================
# ── 7. HAPLOTYPE BLOCK CALLING (GABRIEL ET AL. 2002) ─────────────────────────
# =============================================================================

message("\n[3/6] Calling haplotype blocks with PLINK (Gabriel et al. 2002)...")

for (p in panels) {
  bed_prefix   <- file.path(CONFIG$output_dir, paste0(p$tag, "_plink"))
  block_prefix <- file.path(CONFIG$output_dir, paste0(p$tag, "_blocks"))
  block_file   <- paste0(block_prefix, ".blocks.det")
  
  if (!file.exists(block_file)) {
    message("  Calling blocks for ", p$tag, "...")
    run_cmd(
      paste0(
        '"', plink_bin, '"',
        " --bfile ",      '"', bed_prefix,   '"',
        " --chr-set 7",
        " --allow-extra-chr",      # NEW: Required to read the string-based chr names in the .bim file
        " --blocks no-pheno-req",
        " --blocks-max-kb ", CONFIG$plink_blocks_max_kb,
        " --out ",        '"', block_prefix, '"'
      ),
      step = paste("block calling:", p$tag)
    )
  } else {
    message("  Block file for ", p$tag, " already exists — skipping.")
  }
  
  if (!file.exists(block_file))
    stop("Expected block file not found after PLINK run: ", block_file,
         "\n  Check PLINK log: ", block_prefix, ".log")
}

# Parse block files
message("  Parsing block output files...")
blocks_activate <- parse_plink_blocks(
  file.path(CONFIG$output_dir, "activate_blocks.blocks.det"), "ACTIVATE")
blocks_ldp      <- parse_plink_blocks(
  file.path(CONFIG$output_dir, "ldp_blocks.blocks.det"),      "LDP")

message("  Blocks detected — ACTIVATE: ", nrow(blocks_activate),
        "  |  LDP: ", nrow(blocks_ldp))

# =============================================================================
# ── 8. LOAD ΔHe WINDOWS AND CLASSIFY REGIONS ─────────────────────────────────
# =============================================================================

message("\n[4/6] Loading pre-filtered ΔHe regions and classifying blocks...")

# Helper: read a SignificantRegions CSV and normalise column names flexibly
read_delta_he <- function(path, region_label) {
  if (!file.exists(path))
    stop("ΔHe file not found: ", path)
  df <- read.csv(path, stringsAsFactors = FALSE) %>%
    rename_with(tolower)
  
  # Identify chr/position columns regardless of exact naming convention
  chr_col   <- intersect(c("chr", "chrom", "chromosome"),   colnames(df))[1]
  start_col <- intersect(c("window_start", "start", "bp1"), colnames(df))[1]
  end_col   <- intersect(c("window_end", "stop", "end", "bp2"), colnames(df))[1]
  
  if (any(is.na(c(chr_col, start_col, end_col))))
    stop("Cannot identify chr/start/end columns in ", path,
         "\n  Found: ", paste(colnames(df), collapse = ", "),
         "\n  Expected one of: chr/chrom/chromosome, ",
         "window_start/start/bp1, window_end/stop/end/bp2")
  
  df %>%
    transmute(
      # NEW: Clean the chromosome column. 
      # gsub(".*Chr", "", ...) deletes everything up to and including "Chr",
      # turning "Lcu.1GRN.Chr4" directly into "4" so it matches PLINK perfectly.
      chr          = gsub(".*Chr", "", as.character(.data[[chr_col]])),
      window_start = as.integer(.data[[start_col]]),
      window_end   = as.integer(.data[[end_col]])
    )
}

eroded_windows       <- read_delta_he(CONFIG$delta_he_erosion_file,       "eroded")
introgressed_windows <- read_delta_he(CONFIG$delta_he_introgression_file, "introgressed")

message("  Eroded windows loaded       : ", nrow(eroded_windows))
message("  Introgressed windows loaded : ", nrow(introgressed_windows))
message("  Using pre-filtered SignificantRegions files — no percentile thresholding applied.")

# Build GRanges objects for interval overlap
to_gr <- function(df) {
  GRanges(seqnames = df$chr,
          ranges   = IRanges(start = df$window_start, end = df$window_end))
}

eroded_gr        <- to_gr(eroded_windows)
introgressed_gr  <- to_gr(introgressed_windows)

# Classify blocks from each panel
blocks_activate_cl <- classify_blocks(blocks_activate, eroded_gr, introgressed_gr)
blocks_ldp_cl      <- classify_blocks(blocks_ldp,      eroded_gr, introgressed_gr)

# Combine and save
blocks_all <- bind_rows(blocks_activate_cl, blocks_ldp_cl)

write.csv(blocks_activate_cl,
          file.path(CONFIG$output_dir, "blocks_annotated_activate.csv"),
          row.names = FALSE)
write.csv(blocks_ldp_cl,
          file.path(CONFIG$output_dir, "blocks_annotated_ldp.csv"),
          row.names = FALSE)

message("  Region classification summary:")
blocks_all %>%
  count(panel, region_type) %>%
  tidyr::pivot_wider(names_from = region_type, values_from = n, values_fill = 0) %>%
  print()


# =============================================================================
# ── 9. STATISTICAL TESTS ──────────────────────────────────────────────────────
# =============================================================================

message("\n[5/6] Running Mann-Whitney U tests...")

# Comparisons to run per panel:
#   eroded vs background
#   introgressed vs background
#   eroded vs introgressed
comparisons <- list(
  c("eroded",        "background"),
  c("introgressed",  "background"),
  c("eroded",        "introgressed")
)

mwu_results <- map_dfr(c("ACTIVATE", "LDP"), function(pan) {
  df <- blocks_all %>% filter(panel == pan)
  map_dfr(comparisons, function(comp) {
    run_mwu(df,
            group_col = "region_type",
            value_col = "block_len_bp",
            group_a   = comp[1],
            group_b   = comp[2],
            min_n     = CONFIG$min_blocks_for_test) %>%
      mutate(panel = pan, .before = 1)
  })
}) %>%
  mutate(
    fdr_p = p.adjust(p_value, method = "BH"),
    significant = fdr_p < 0.05
  )

write.csv(mwu_results,
          file.path(CONFIG$output_dir, "mannwhitney_results.csv"),
          row.names = FALSE)
message("  Saved: mannwhitney_results.csv")
message("  Significant comparisons (FDR < 0.05):")
print(mwu_results %>%
        filter(significant) %>%
        select(panel, comparison, n_a, n_b,
               median_a, median_b, p_value, fdr_p))


# =============================================================================
# ── 10. FIGURES ───────────────────────────────────────────────────────────────
# =============================================================================

message("\n[6/6] Generating figures...")

# ── Colour scheme ─────────────────────────────────────────────────────────────
REGION_COLOURS <- c(
  eroded        = "#D6604D",   # warm red
  introgressed  = "#4393C3",   # blue
  background    = "grey72"
)

REGION_LABELS <- c(
  eroded        = "Eroded\n(top 1% \u0394He)",
  introgressed  = "Introgressed\n(bottom 1% \u0394He)",
  background    = "Background"
)

PANEL_LABELS <- c(ACTIVATE = "ACTIVATE panel", LDP = "LDP")

# Ordered factor for consistent x-axis ordering
blocks_all <- blocks_all %>%
  mutate(
    region_type = factor(region_type,
                         levels = c("eroded", "introgressed", "background")),
    panel       = factor(panel, levels = c("ACTIVATE", "LDP"))
  )

shared_theme <- theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9, colour = "grey35"),
    plot.caption     = element_text(size = 8, colour = "grey50", hjust = 0),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey93", linewidth = 0.4),
    strip.background = element_rect(fill = "grey96"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "none"
  )

# ── Figure 1: Violin + boxplot of block lengths by region and panel ───────────
# Add significance brackets from MWU results
sig_brackets <- mwu_results %>%
  filter(significant) %>%
  select(panel, comparison, fdr_p) %>%
  mutate(
    group1 = sub(" vs .*", "", comparison),
    group2 = sub(".* vs ", "", comparison),
    label  = ifelse(fdr_p < 0.001, "***",
                    ifelse(fdr_p < 0.01, "**", "*"))
  )

p_violin <- ggplot(blocks_all,
                   aes(x = region_type, y = block_len_bp / 1000,
                       fill = region_type, colour = region_type)) +
  
  geom_violin(alpha = 0.35, linewidth = 0.4, trim = TRUE) +
  geom_boxplot(width = 0.18, alpha = 0.85, linewidth = 0.5,
               outlier.size = 0.6, outlier.alpha = 0.4,
               colour = "grey30", fill = "white") +
  
  # Median label
  stat_summary(fun = median, geom = "text",
               aes(label = round(after_stat(y), 1)),
               vjust = -0.5, size = 2.8,
               colour = "grey20") +
  
  scale_fill_manual(values   = REGION_COLOURS) +
  scale_colour_manual(values = REGION_COLOURS) +
  scale_x_discrete(labels = REGION_LABELS) +
  scale_y_continuous(
    name   = "Haplotype block length (kb)",
    labels = comma_format(accuracy = 1),
    trans  = "log10"                          # log scale — block lengths are skewed
  ) +
  facet_wrap(~ panel, ncol = 2,
             labeller = as_labeller(PANEL_LABELS)) +
  labs(
    title    = "Haplotype block lengths by genomic region type",
    subtitle = paste0(
      "Gabriel et al. (2002) algorithm | Max block size: ",
      CONFIG$plink_blocks_max_kb, " kb | ",
      "Significance: Mann-Whitney U, BH-corrected (* FDR<0.05, ** FDR<0.01, *** FDR<0.001)"
    ),
    caption  = paste0(
      "Eroded regions: top 1% \u0394He (He_LDP \u2212 He_ACTIVATE).  ",
      "Introgressed regions: bottom 1% \u0394He.\n",
      "Y-axis on log\u2081\u2080 scale. White box = interquartile range; ",
      "number = median block length (kb)."
    ),
    x = NULL
  ) +
  shared_theme

ggsave(file.path(CONFIG$output_dir, "plot_block_lengths_violin.pdf"),
       p_violin,
       width = CONFIG$figure_width, height = CONFIG$figure_height,
       dpi   = CONFIG$dpi)
message("  Saved: plot_block_lengths_violin.pdf")

# ── Figure 2: Empirical CDF curves ────────────────────────────────────────────
p_ecdf <- ggplot(blocks_all,
                 aes(x = block_len_bp / 1000,
                     colour = region_type, linetype = region_type)) +
  
  stat_ecdf(linewidth = 0.8, geom = "step") +
  
  scale_colour_manual(
    name   = "Region type",
    values = REGION_COLOURS,
    labels = REGION_LABELS
  ) +
  scale_linetype_manual(
    name   = "Region type",
    values = c(eroded = "solid", introgressed = "dashed", background = "dotted"),
    labels = REGION_LABELS
  ) +
  scale_x_continuous(
    name   = "Haplotype block length (kb)",
    trans  = "log10",
    labels = comma_format(accuracy = 1)
  ) +
  scale_y_continuous(
    name   = "Cumulative proportion",
    labels = percent_format(accuracy = 1)
  ) +
  facet_wrap(~ panel, ncol = 2,
             labeller = as_labeller(PANEL_LABELS)) +
  labs(
    title   = "Empirical CDF of haplotype block lengths by region type",
    caption = paste0(
      "Curves shifted to the right indicate longer blocks (stronger/older LD).\n",
      "Eroded regions with longer blocks than background indicate LD extension ",
      "due to selection-driven haplotype fixation."
    )
  ) +
  shared_theme +
  theme(legend.position = "bottom",
        legend.title    = element_text(size = 9),
        legend.text     = element_text(size = 8))

ggsave(file.path(CONFIG$output_dir, "plot_block_lengths_ecdf.pdf"),
       p_ecdf,
       width = CONFIG$figure_width, height = CONFIG$figure_height * 0.75,
       dpi   = CONFIG$dpi)
message("  Saved: plot_block_lengths_ecdf.pdf")

# ── Figure 3: Genome-wide block density track ─────────────────────────────────
# Bins the genome into 1 Mb windows and counts blocks per bin per panel,
# coloured by dominant region type in that bin.
bin_size <- 1e6   # 1 Mb bins

density_df <- blocks_all %>%
  mutate(bin = floor(block_start / bin_size) * bin_size) %>%
  group_by(panel, chr, bin, region_type) %>%
  summarise(n_blocks = n(), .groups = "drop") %>%
  group_by(panel, chr, bin) %>%
  mutate(prop = n_blocks / sum(n_blocks)) %>%
  ungroup() %>%
  mutate(chr_num = as.integer(gsub("[^0-9]", "", chr)))  # numeric chr for ordering

p_density <- ggplot(density_df,
                    aes(x = bin / 1e6, y = n_blocks, fill = region_type)) +
  
  geom_col(width = bin_size / 1e6, position = "stack", alpha = 0.85) +
  
  scale_fill_manual(
    name   = "Region type",
    values = REGION_COLOURS,
    labels = REGION_LABELS
  ) +
  scale_x_continuous(name = "Chromosomal position (Mb)",
                     labels = comma_format(accuracy = 1)) +
  scale_y_continuous(name = "Number of haplotype blocks per Mb") +
  facet_grid(panel ~ chr_num,
             scales = "free_x", space = "free_x",
             labeller = labeller(
               panel   = PANEL_LABELS,
               chr_num = function(x) paste0("Chr", x)
             )) +
  labs(
    title   = "Genome-wide haplotype block density",
    subtitle = "Blocks per 1 Mb window, coloured by region classification",
    caption  = paste0(
      "Orange = blocks overlapping eroded regions (top 1% \u0394He).  ",
      "Blue = introgressed regions (bottom 1% \u0394He).  ",
      "Grey = background."
    )
  ) +
  shared_theme +
  theme(
    legend.position  = "bottom",
    axis.text.x      = element_text(size = 7, angle = 45, hjust = 1),
    panel.spacing.x  = unit(0.15, "lines")
  )

ggsave(file.path(CONFIG$output_dir, "plot_block_density_genome.pdf"),
       p_density,
       width  = CONFIG$figure_width * 1.4,
       height = CONFIG$figure_height,
       dpi    = CONFIG$dpi)
message("  Saved: plot_block_density_genome.pdf")


# =============================================================================
# ── 11. CONSOLE SUMMARY ───────────────────────────────────────────────────────
# =============================================================================

cat("\n", strrep("=", 68), "\n")
cat("ANALYSIS COMPLETE\n")
cat(strrep("=", 68), "\n\n")
cat("Output directory:", CONFIG$output_dir, "\n\n")

cat("Files written:\n")
cat("  blocks_annotated_activate.csv\n")
cat("  blocks_annotated_ldp.csv\n")
cat("  mannwhitney_results.csv\n")
cat("  plot_block_lengths_violin.pdf\n")
cat("  plot_block_lengths_ecdf.pdf\n")
cat("  plot_block_density_genome.pdf\n\n")

cat("Block count summary:\n")
blocks_all %>%
  count(panel, region_type) %>%
  tidyr::pivot_wider(names_from = region_type,
                     values_from = n, values_fill = 0) %>%
  print()

cat("\nMedian block lengths (kb) by panel and region:\n")
blocks_all %>%
  group_by(panel, region_type) %>%
  summarise(
    n          = n(),
    median_kb  = round(median(block_kb, na.rm = TRUE), 2),
    mean_kb    = round(mean(block_kb,   na.rm = TRUE), 2),
    max_kb     = round(max(block_kb,    na.rm = TRUE), 2),
    .groups    = "drop"
  ) %>%
  print()

cat("\nMann-Whitney U results:\n")
print(mwu_results %>%
        select(panel, comparison, n_a, n_b,
               median_a, median_b, p_value, fdr_p, significant))
cat("\n")
