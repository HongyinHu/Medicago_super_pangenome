#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/N_5.pod_coil/03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714
OUT=$ROOT/06_maf_sensitivity_modelA_20260715
PYTHON=path/to/home/anaconda3/envs/panpop/bin/python
RSCRIPT=path/to/home/anaconda3/envs/biosofeware/bin/Rscript
mkdir -p "$OUT/diagnostics" "$OUT/logs" "$OUT/tmp/matplotlib"
export MPLCONFIGDIR="$OUT/tmp/matplotlib"
LOG="$OUT/logs/diagnostics.log"
exec > >(tee -a "$LOG") 2>&1

"$PYTHON" "$ROOT/scripts/diagnose_model_calibration.py" --root "$ROOT" --outdir "$OUT/diagnostics"
"$PYTHON" "$ROOT/scripts/make_structure_diagnostics.py" --root "$ROOT" --outdir "$OUT/diagnostics"
gapit=$("$RSCRIPT" -e 'cat(requireNamespace("GAPIT", quietly=TRUE))' 2>/dev/null || true)
printf 'software\tstatus\treason\nFarmCPU_BLINK_GAPIT\t%s\t%s\n' \
  "${gapit:-FALSE}" "GAPIT package was not installed in the validated analysis R environment; no alternative software substituted." > "$OUT/diagnostics/farmcpu_blink_environment.tsv"
for required in "$OUT/diagnostics/PCA_trait_structure.png" "$OUT/diagnostics/GRM_heatmap.png" "$OUT/diagnostics/GRM_pairwise_distribution.png" "$OUT/diagnostics/phenotype_pc_association.tsv"; do
  test -s "$required"
done
printf 'completed\t%s\n' "$(date '+%F %T %Z')" > "$OUT/summary/diagnostics.done"
