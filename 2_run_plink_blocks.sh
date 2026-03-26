#!/bin/bash

# Configuration
# Assuming PLINK is in the PATH. If not, replace 'plink' with the full path (e.g. /home/salvador/bin/plink)
PLINK_CMD="plink"
BASE_DIR="c:/Users/Salva/OneDrive - University of Saskatchewan/UsasK/github/ACTIVATE_myG"
RESULTS_DIR="${BASE_DIR}/Results"

cd "$RESULTS_DIR" || exit 1

echo "========================================================="
echo "        Haplotype Blocking for ACTIVATE Population       "
echo "========================================================="
# 1. Convert VCF to Bed/Bim/Fam
$PLINK_CMD --vcf activate_merged_snps.vcf.gz --make-bed --out activate_merged_snps_plink --double-id --allow-extra-chr
# 2. Run Blocks
$PLINK_CMD --bfile activate_merged_snps_plink --blocks no-pheno-req --out activate_blocks --allow-extra-chr

echo ""
echo "========================================================="
echo "          Haplotype Blocking for LDP Population          "
echo "========================================================="
# 1. Convert VCF to Bed/Bim/Fam
$PLINK_CMD --vcf ldp_merged_snps.vcf.gz --make-bed --out ldp_merged_snps_plink --double-id --allow-extra-chr
# 2. Run Blocks
$PLINK_CMD --bfile ldp_merged_snps_plink --blocks no-pheno-req --out ldp_blocks --allow-extra-chr

echo ""
echo "========================================================="
echo "Done! The files activate_blocks.blocks.det and "
echo "ldp_blocks.blocks.det are ready for downstream analysis."
echo "========================================================="
