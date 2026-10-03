#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/N_5.pod_coil/03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714
OUT=$ROOT/06_maf_sensitivity_modelA_20260715
TAG=maf005
mkdir -p "$OUT/$TAG/raw/SNP" "$OUT/$TAG/raw/INDEL" "$OUT/$TAG/raw/SV" "$OUT/summary/tasks" "$OUT/logs"
LOG="$OUT/logs/${TAG}.baseline_reuse.log"
exec > >(tee -a "$LOG") 2>&1

for label in SNP INDEL; do
  for chr in $(seq 1 8); do
    source="$ROOT/02_snp_indel_gwas/gmmat/${label}.chr${chr}.coil_ge1.pc5.gmmat.score.tsv"
    target="$OUT/$TAG/raw/$label/${label}.chr${chr}.Model_A_${TAG}.gmmat.score.tsv"
    test -s "$source"
    ln -sfn "$source" "$target"
    test -s "$target"
    done_file="$OUT/summary/tasks/${TAG}.${label}.chr${chr}.done"
    printf 'baseline_reused\t%s\nmaf_cutoff\t0.05\nscore_file\t%s\n' "$(date '+%F %T %Z')" "$source" > "$done_file"
  done
done

source=$(find "$ROOT/03_sv_gwas" -type f -name 'SV*.pc5.gmmat.score.tsv' | head -n 1)
test -n "$source"
target="$OUT/$TAG/raw/SV/SV.Model_A_${TAG}.gmmat.score.tsv"
ln -sfn "$source" "$target"
test -s "$target"
printf 'baseline_reused\t%s\nmaf_cutoff\t0.05\nscore_file\t%s\n' "$(date '+%F %T %Z')" "$source" > "$OUT/summary/tasks/${TAG}.SV.done"
printf 'completed\t%s\nmethod\tReused exact Model A baseline, which already used MAF >= 0.05\n' "$(date '+%F %T %Z')" > "$OUT/summary/${TAG}.prepared.done"
