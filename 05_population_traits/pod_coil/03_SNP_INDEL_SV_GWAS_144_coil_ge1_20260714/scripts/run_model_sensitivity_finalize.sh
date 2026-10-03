#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/N_5.pod_coil/03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714
SENS=$ROOT/05_model_sensitivity_20260715
PYTHON=path/to/home/anaconda3/envs/panpop/bin/python
export MPLCONFIGDIR="$SENS/tmp/matplotlib"
mkdir -p "$SENS"/{diagnostics,results,summary,tmp/matplotlib,logs}
LOG="$SENS/logs/finalize_model_sensitivity.log"
exec > >(tee -a "$LOG") 2>&1

for label in SNP INDEL; do
  for chrom in $(seq 1 8); do
    test -s "$SENS/summary/tasks/${label}.chr${chrom}.done"
  done
done
test -s "$SENS/summary/tasks/SV.done"

echo "[$(date '+%F %T %Z')] baseline calibration and structure diagnostics"
"$PYTHON" "$ROOT/scripts/diagnose_model_calibration.py" \
  --root "$ROOT" --outdir "$SENS/diagnostics"
"$PYTHON" "$ROOT/scripts/make_structure_diagnostics.py" \
  --root "$ROOT" --outdir "$SENS/diagnostics"

echo "[$(date '+%F %T %Z')] A/B/C association summary"
"$PYTHON" "$ROOT/scripts/summarize_model_sensitivity.py" \
  --root "$ROOT" --outdir "$SENS/results"

for required in \
  "$SENS/results/ABC_model_comparison.tsv" \
  "$SENS/results/ABC_lambda_by_MAF.tsv" \
  "$SENS/results/ABC_candidate_peak_stability.tsv" \
  "$SENS/results/ABC_model_comparison_QQ.png" \
  "$SENS/results/Model_A_SNP_INDEL_SV_Manhattan.png" \
  "$SENS/results/Model_B_SNP_INDEL_SV_Manhattan.png" \
  "$SENS/results/Model_C_SNP_INDEL_SV_Manhattan.png" \
  "$SENS/diagnostics/PCA_trait_structure.png" \
  "$SENS/diagnostics/GRM_heatmap.png" \
  "$SENS/diagnostics/GRM_pairwise_distribution.png"; do
  test -s "$required"
done

{
  printf 'completed\t%s\n' "$(date '+%F %T %Z')"
  printf 'baseline_model\tA: GMMAT PC1-PC5 + SNP GRM\n'
  printf 'sensitivity_model_B\tPLINK logistic PC1-PC3, no GRM\n'
  printf 'sensitivity_model_C\tGMMAT PC1-PC2 + SNP GRM\n'
  printf 'selection_rule\tCalibration, trait-PC confounding, MAF-stratified lambda, and peak stability; not lambda alone\n'
} > "$SENS/summary/model_sensitivity.done"
echo "[$(date '+%F %T %Z')] model sensitivity finalize complete"
