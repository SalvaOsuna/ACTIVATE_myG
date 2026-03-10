#Set reticulate

# Load reticulate
library(reticulate)

# --- 1. Set Up the Isolated Python Environment ---
env_name <- "xpclr_env"

# Tell reticulate to explicitly use this environment
use_virtualenv(env_name, required = TRUE)

# Extract the exact absolute path to the Python executable
py_bin <- py_exe()

# Locate the xpclr executable in the exact same folder as the Python executable
xpclr_exe <- file.path(dirname(py_bin), "xpclr")
if (.Platform$OS.type == "windows") {
  xpclr_exe <- paste0(xpclr_exe, ".exe")
}

cat("Environment ready! Using XP-CLR executable:", xpclr_exe, "\n")

# --- 2. Define the XP-CLR Python Execution Function ---
run_python_xpclr <- function(vcf_file, ref_file, obj_file, out_prefix) {
  
  # Ensure the output directory exists
  if (!dir.exists("Results")) dir.create("Results")
  
  # Loop through the 7 lentil chromosomes
  for (i in 1:7) {
    chr_name <- paste0("Lcu.1GRN.Chr", i)
    out_file <- paste0("Results/", out_prefix, "_Chr", i)
    
    cat(sprintf("\nRunning %s on %s...\n", out_prefix, chr_name))
    
    # Build the command-line arguments (Removed "-m xpclr")
    args <- c(
      "--out", out_file,
      "--format", "vcf",
      "--input", vcf_file,
      "--samplesA", ref_file,
      "--samplesB", obj_file,
      "--chr", chr_name,
      "--size", "50000",
      "--step", "10000",
      "--phased", "False"
    )
    
    # Execute the command directly using the xpclr executable
    system2(xpclr_exe, args)
  }
}

# --- 3. Execute Both Scenarios ---
# Ensure this points to your exported VCF
vcf_data <- "data/Merged_Analysis.vcf.gz"

# SCENARIO 1: The Adaptation Scan
cat("\n=== Starting Scenario 1: Adaptation ===\n")
run_python_xpclr(
  vcf_file   = vcf_data,
  ref_file   = "data/LDP227_Non_Temperate.txt",
  obj_file   = "data/ACT187_samples.txt",
  out_prefix = "Py_Scenario1_Adaptation"
)

# SCENARIO 2: The Breeding Scan
cat("\n=== Starting Scenario 2: Breeding ===\n")
run_python_xpclr(
  vcf_file   = vcf_data,
  ref_file   = "data/LDP97_Temperate.txt",
  obj_file   = "data/ACT187_samples.txt",
  out_prefix = "Py_Scenario2_Breeding"
)

reticulate::use_virtualenv("xpclr_env", required = TRUE)
reticulate::py_install("numpy==1.26.4", envname = "xpclr_env", pip = TRUE)
