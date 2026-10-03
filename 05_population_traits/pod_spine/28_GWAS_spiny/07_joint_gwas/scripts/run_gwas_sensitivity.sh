#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/N_4.pod_spiny/28_GWAS_spiny
STAGE="$ROOT/07_joint_gwas"
STRUCT="$ROOT/03_population_structure"
SVSTAGE="$ROOT/06_sv_gwas"
GEMMA=path/to/home/anaconda3/envs/biosofeware/bin/gemma
PYTHON=path/to/home/anaconda3/envs/panpop/bin/python
OUT="$STAGE/sensitivity/gemma"
PHENO="$STRUCT/summary/phenotype.0_1.txt"
KIN="$STRUCT/results/gemma/output/SNP_kinship.cXX.txt"
PC5="$STRUCT/summary/covariates.pc5.txt"
PC3="$STAGE/sensitivity/covariates.pc3.txt"
SV_GENO="$SVSTAGE/inputs/medicago143.SV.gemma.geno.txt"
SV_ANNO="$SVSTAGE/inputs/medicago143.SV.gemma.anno.txt"

mkdir -p "$OUT" "$STAGE/sensitivity" "$STAGE/logs/slurm" "$STAGE/summary"
test -s "$ROOT/04_snp_indel_gwas/summary/snp_indel_gwas.done"
test -s "$SVSTAGE/summary/sv_gwas.done"
awk '{print $1, $2, $3, $4}' "$PC5" > "$PC3"
[[ "$(wc -l < "$PC3")" -eq 143 ]]

run_plink() {
    local variant=$1
    local model=$2
    local covar=${3:-}
    local assoc="$OUT/output/${variant}.${model}.assoc.txt"
    [[ -s "$assoc" ]] && return 0
    mkdir -p "$OUT"
    (
        cd "$OUT"
        cmd=("$GEMMA" -bfile "$STRUCT/results/plink/analysis.${variant}" -k "$KIN" -p "$PHENO" -lmm 4 -o "${variant}.${model}")
        [[ -n "$covar" ]] && cmd+=( -c "$covar" )
        "${cmd[@]}"
    )
}

run_sv() {
    local model=$1
    local covar=${2:-}
    local assoc="$OUT/output/SV.${model}.assoc.txt"
    [[ -s "$assoc" ]] && return 0
    (
        cd "$OUT"
        cmd=("$GEMMA" -g "$SV_GENO" -a "$SV_ANNO" -k "$KIN" -p "$PHENO" -miss 0.20 -maf 0.05 -lmm 4 -o "SV.${model}")
        [[ -n "$covar" ]] && cmd+=( -c "$covar" )
        "${cmd[@]}"
    )
}

for variant in SNP INDEL; do
    run_plink "$variant" Konly
    run_plink "$variant" K_PC3 "$PC3"
done
run_sv Konly
run_sv K_PC3 "$PC3"

"$PYTHON" "$STAGE/scripts/summarize_sensitivity.py" \
    --stage-root "$ROOT" \
    --sensitivity-dir "$OUT/output" \
    --out "$STAGE/results/tables/model_sensitivity.tsv"

test -s "$STAGE/results/tables/model_sensitivity.tsv"
printf 'completed\t%s\nmodels	Konly,K_PC3,K_PC5_primary\n' "$(date '+%F %T %Z')" \
    > "$STAGE/summary/sensitivity.done"
echo "[$(date '+%F %T')] GWAS sensitivity analysis completed"
