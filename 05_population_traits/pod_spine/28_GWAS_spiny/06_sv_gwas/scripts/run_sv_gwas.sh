#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/N_4.pod_spiny/28_GWAS_spiny
STAGE05="$ROOT/05_sv_genotyping"
STAGE="$ROOT/06_sv_gwas"
STRUCT="$ROOT/03_population_structure"
PYTHON=path/to/home/anaconda3/envs/panpop/bin/python
GEMMA=path/to/home/anaconda3/envs/biosofeware/bin/gemma

QC_VCF="$STAGE05/cohort/medicago143.smoove.square.gq20.dp5.cr80.maf05.mac5.bnd1.vcf.gz"
SAMPLE_ORDER="$STAGE05/summary/sv_analysis_samples.order.txt"
GENO="$STAGE/inputs/medicago143.SV.gemma.geno.txt"
ANNO="$STAGE/inputs/medicago143.SV.gemma.anno.txt"
MAPPING="$STAGE/inputs/medicago143.SV.marker_map.tsv.gz"
INPUT_SUMMARY="$STAGE/summary/sv_gemma_input_counts.tsv"
KIN="$STRUCT/results/gemma/output/SNP_kinship.cXX.txt"
COVAR="$STRUCT/summary/covariates.pc5.txt"
PHENO="$STRUCT/summary/phenotype.0_1.txt"
ASSOC="$STAGE/results/gemma/output/SV.pc5.assoc.txt"

mkdir -p "$STAGE/inputs" "$STAGE/results/gemma" "$STAGE/results/summary" "$STAGE/logs" "$STAGE/summary"
test -s "$STAGE05/summary/sv_cohort_qc.done"
test -s "$QC_VCF"
test -s "$QC_VCF.csi"
test -s "$SAMPLE_ORDER"
test -s "$KIN"
test -s "$COVAR"
test -s "$PHENO"
[[ "$(wc -l < "$SAMPLE_ORDER")" -eq 143 ]]
[[ "$(wc -l < "$COVAR")" -eq 143 ]]
[[ "$(wc -l < "$PHENO")" -eq 143 ]]
[[ "$(wc -l < "$KIN")" -eq 143 ]]

if [[ ! -s "$GENO" || ! -s "$ANNO" || ! -s "$MAPPING" ]]; then
    rm -f "$GENO" "$ANNO" "$MAPPING" "$INPUT_SUMMARY"
    "$PYTHON" "$STAGE/scripts/vcf_to_bimbam.py" \
        --vcf "$QC_VCF" \
        --sample-order "$SAMPLE_ORDER" \
        --geno "$GENO" \
        --anno "$ANNO" \
        --mapping "$MAPPING" \
        --summary "$INPUT_SUMMARY"
fi

variants=$(awk -F '\t' '$1=="variants" {print $2}' "$INPUT_SUMMARY")
[[ -n "$variants" && "$variants" -gt 0 ]]
[[ "$(wc -l < "$GENO")" -eq "$variants" ]]
[[ "$(wc -l < "$ANNO")" -eq "$variants" ]]

if [[ ! -s "$ASSOC" ]]; then
    (
        cd "$STAGE/results/gemma"
        "$GEMMA" \
            -g "$GENO" \
            -a "$ANNO" \
            -p "$PHENO" \
            -k "$KIN" \
            -c "$COVAR" \
            -miss 0.20 \
            -maf 0.05 \
            -lmm 4 \
            -o SV.pc5
    )
fi
test -s "$ASSOC"

"$PYTHON" "$STAGE/scripts/summarize_sv_gemma.py" \
    --assoc "$ASSOC" \
    --mapping "$MAPPING" \
    --outdir "$STAGE/results/summary"

{
    printf 'completed\t%s\n' "$(date '+%F %T %Z')"
    printf 'samples\t143\n'
    printf 'phenotype\t0=spineless,1=spiny\n'
    printf 'model\tGEMMA_LMM_Wald_primary\n'
    printf 'kinship\tSNP_LD_pruned_centered\n'
    printf 'covariates\tintercept_plus_PC1-PC5\n'
    printf 'genotype_mask\tGQ<20_or_DP<5_to_missing\n'
    printf 'site_filter\tcall_rate>=0.80,MAF>=0.05,MAC>=5\n'
    printf 'BND_representation\tone_record_per_reciprocal_pair\n'
    printf 'association_file\t%s\n' "$ASSOC"
} > "$STAGE/summary/sv_gwas.done"

echo "[$(date '+%F %T')] SV GWAS completed"
