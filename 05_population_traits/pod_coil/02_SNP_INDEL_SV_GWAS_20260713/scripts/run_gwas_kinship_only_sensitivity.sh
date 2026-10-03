#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/N_5.pod_coil/02_SNP_INDEL_SV_GWAS_20260713
GEMMA=path/to/home/anaconda3/envs/biosofeware/bin/gemma
PYTHON=path/to/home/anaconda3/envs/panpop/bin/python
KIN="$ROOT/01_structure/gemma/output/coil_SNP_kinship.cXX.txt"
PHENO="$ROOT/01_structure/summary/phenotype.0_1.txt"
OUT="$ROOT/02_snp_indel_gwas/kinship_only"
LOG="$ROOT/logs/run_gwas_kinship_only_sensitivity.log"

mkdir -p "$OUT/gemma" "$OUT/results"
exec > >(tee -a "$LOG") 2>&1

echo "[$(date '+%F %T %Z')] kinship-only sensitivity GWAS start on $(hostname)"
test -s "$KIN"
test -s "$PHENO"

run_one() {
    local label=$1
    local result="$OUT/gemma/output/${label}.coil.kinship_only.assoc.txt"
    if [[ ! -s "$result" ]]; then
        (
            cd "$OUT/gemma"
            "$GEMMA" -bfile "$ROOT/01_structure/plink/coil.$label" \
                -k "$KIN" -p "$PHENO" -lmm 4 -o "${label}.coil.kinship_only"
        )
    fi
    test -s "$result"
}

run_one SNP &
pid_snp=$!
run_one INDEL &
pid_indel=$!
wait "$pid_snp"
wait "$pid_indel"

"$PYTHON" "$ROOT/scripts/summarize_gwas.py" \
    --assoc "SNP=$OUT/gemma/output/SNP.coil.kinship_only.assoc.txt" \
    --assoc "INDEL=$OUT/gemma/output/INDEL.coil.kinship_only.assoc.txt" \
    --outdir "$OUT/results" --prefix coil_SNP_INDEL_kinship_only

{
    printf 'completed\t%s\n' "$(date '+%F %T %Z')"
    printf 'samples\t125\n'
    printf 'model\tGEMMA_LMM_Wald_kinship_only\n'
    printf 'kinship\tLD_pruned_SNP_centered\n'
    printf 'covariates\tintercept_only\n'
} > "$OUT/kinship_only.done"

echo "[$(date '+%F %T %Z')] kinship-only sensitivity GWAS completed"
