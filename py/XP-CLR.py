import os
import sys
import subprocess
import numpy as np

# --- 1. The Level 2 Auto-Patch ---
# We are upgrading the patch to explicitly catch and sanitize empty strings ('') 
# that scikit-allel improperly scraped from the VCF.
print("Applying advanced patch to handle empty strings in the VCF...")
patched = False
for path in sys.path:
    methods_file = os.path.join(path, "xpclr", "methods.py")
    if os.path.exists(methods_file):
        with open(methods_file, "r") as f:
            content = f.read()
        
        # We look for either the original line OR the one we patched in the previous step
        old_line_1 = "distance = np.abs(dq - dq.mean())"
        old_line_2 = "distance = np.abs(dq.astype(float) - dq.astype(float).mean())"
        
        # New robust line that sanitizes empty strings before doing math
        new_line = "dq_clean = np.array([float(x) if str(x).strip() not in ['', 'None'] else 0.0 for x in dq]); distance = np.abs(dq_clean - dq_clean.mean())"
        
        if old_line_2 in content:
            content = content.replace(old_line_2, new_line)
            with open(methods_file, "w") as f:
                f.write(content)
            print("  -> Patch upgraded! Empty strings are now sanitized.")
            patched = True
            break
        elif old_line_1 in content:
            content = content.replace(old_line_1, new_line)
            with open(methods_file, "w") as f:
                f.write(content)
            print("  -> Patch applied! Empty strings are now sanitized.")
            patched = True
            break

if not patched:
    print("  -> Could not locate methods.py or patch already applied.")


# --- 2. Setup paths and directories ---
scripts_dir = os.path.dirname(sys.executable)
xpclr_script = os.path.join(scripts_dir, "xpclr-script.py")
if not os.path.exists(xpclr_script):
    xpclr_script = os.path.join(scripts_dir, "xpclr") # Fallback

os.makedirs("Results", exist_ok=True)
vcf_data = "data/Merged_Analysis.vcf.gz"


# --- 3. Define and Execute ---
def run_xpclr_scenario(ref_file, obj_file, out_prefix):
    print(f"\n{'='*50}")
    print(f"Starting {out_prefix}")
    print(f"{'='*50}")
    
    for i in range(1, 8):
        chr_name = f"Lcu.1GRN.Chr{i}"
        out_file = f"Results/{out_prefix}_Chr{i}"
        
        print(f"  Processing {chr_name}...")
        
        cmd = [
            sys.executable,
            xpclr_script,
            "--out", out_file,
            "--format", "vcf",
            "--input", vcf_data,
            "--samplesA", ref_file,
            "--samplesB", obj_file,
            "--chr", chr_name,
            "--size", "50000",
            "--step", "10000"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            print(f"    -> Success!")
        else:
            print(f"    -> ERROR! Here is what went wrong:")
            print(result.stderr)
            break 

# Execute Scenario 2: Breeding
run_xpclr_scenario(
    ref_file="data/LDP97_Temperate.txt", 
    obj_file="data/ACT187_samples.txt", 
    out_prefix="Py_Scenario2_Breeding"
)

# Execute Scenario 1: Adaptation
run_xpclr_scenario(
    ref_file="data/LDP227_Non_Temperate.txt", 
    obj_file="data/ACT187_samples.txt", 
    out_prefix="Py_Scenario1_Adaptation"
)


