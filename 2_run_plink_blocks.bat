@echo off
set PLINK_CMD="C:\Users\Salva\Downloads\plink_win64_20250819\plink.exe"
set RESULTS_DIR="c:\Users\Salva\OneDrive - University of Saskatchewan\UsasK\github\ACTIVATE_myG\Results"

cd /d %RESULTS_DIR%

echo =========================================================
echo         Haplotype Blocking for ACTIVATE Population       
echo =========================================================
%PLINK_CMD% --vcf activate_merged_snps.vcf --make-bed --out activate_merged_snps_plink --double-id --allow-extra-chr
%PLINK_CMD% --bfile activate_merged_snps_plink --blocks no-pheno-req --out activate_blocks --allow-extra-chr

echo.
echo =========================================================
echo           Haplotype Blocking for LDP Population          
echo =========================================================
%PLINK_CMD% --vcf ldp_merged_snps.vcf --make-bed --out ldp_merged_snps_plink --double-id --allow-extra-chr
%PLINK_CMD% --bfile ldp_merged_snps_plink --blocks no-pheno-req --out ldp_blocks --allow-extra-chr

echo.
echo =========================================================
echo Done! The files activate_blocks.blocks.det and 
echo ldp_blocks.blocks.det are ready for downstream analysis.
echo =========================================================
pause
