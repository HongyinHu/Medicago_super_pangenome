#!/usr/bin/env bash
set -euo pipefail

RUN=path/to/project/N_5.pod_coil/08_SHP_multiallelic_candidate_20260717
GWAS=path/to/project/N_5.pod_coil/03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714
PLINK=path/to/home/anaconda3/envs/biosofeware/bin/plink
BFILE=$GWAS/01_structure/plink/coil144.SNP.analysis
MODEL=$GWAS/01_structure/model/SNP.model.tsv
PHENO=$RUN/inputs/coil143_binary_plink.pheno
ASSOC=$RUN/uncorrected/SNP_unadjusted.assoc
OUT=${ASSOC%.assoc}

mkdir -p "$RUN"/{inputs,uncorrected,logs,summary}
exec > >(tee -a "$RUN/logs/uncorrected_snp_assoc.log") 2>&1

[[ -x "$PLINK" ]]
[[ -s "$BFILE.bed" && -s "$BFILE.bim" && -s "$BFILE.fam" ]]
[[ -s "$MODEL" ]]

awk -F '\t' 'NR > 1 {if ($2 == 1) p = 2; else if ($2 == 0) p = 1; else next; print $1, $1, p}' "$MODEL" > "$PHENO"
[[ $(wc -l < "$PHENO") -eq 143 ]]

"$PLINK" \
    --bfile "$BFILE" \
    --pheno "$PHENO" \
    --assoc \
    --allow-no-sex \
    --out "$OUT"

[[ -s "$ASSOC" ]]
printf 'samples\t%s\nvariants\t%s\nassociation_rows\t%s\nmodel\tPLINK_allelic_no_covariates_no_GRM\n' \
    "$(wc -l < "$PHENO")" "$(wc -l < "$BFILE.bim")" "$(( $(wc -l < "$ASSOC") - 1 ))" > "$RUN/summary/uncorrected_snp_assoc_audit.tsv"
echo "[$(date '+%F %T %Z')] unadjusted SNP allelic association complete"
