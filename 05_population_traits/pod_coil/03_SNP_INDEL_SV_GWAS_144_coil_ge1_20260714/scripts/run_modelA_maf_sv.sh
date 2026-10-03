#!/usr/bin/env bash
set -euo pipefail

MAF=${1:?Usage: run_modelA_maf_sv.sh maf}
awk -v maf="$MAF" 'BEGIN { exit !(maf > 0 && maf < 0.5) }'
MAF_TAG=$(awk -v maf="$MAF" 'BEGIN { printf "maf%03d", int(maf * 100 + 0.5) }')

ROOT=path/to/project/N_5.pod_coil/03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714
OUT=$ROOT/06_maf_sensitivity_modelA_20260715
PLINK=path/to/home/anaconda3/envs/biosofeware/bin/plink
RSCRIPT=path/to/home/anaconda3/envs/biosofeware/bin/Rscript
export GWAS_ROOT="$ROOT"
export R_LIBS_USER="$ROOT/software/Rlib"
export TMPDIR="$OUT/tmp/$MAF_TAG/SV"

mkdir -p "$TMPDIR" "$OUT/$MAF_TAG/bfiles/SV" "$OUT/$MAF_TAG/raw/SV" "$OUT/logs" "$OUT/summary/tasks"
LOG="$OUT/logs/${MAF_TAG}.SV.log"
exec > >(tee -a "$LOG") 2>&1

test -s "$ROOT/summary/SV_gwas.done"
SOURCE="$ROOT/03_sv_gwas/plink/coil144.SV.GQ20.analysis"
MODEL="$ROOT/03_sv_gwas/model/SV.model.tsv"
GRM="$ROOT/01_structure/plink/coil144.SNP.grm"
FILTERED="$OUT/$MAF_TAG/bfiles/SV/SV.${MAF_TAG}"
RESULT="$OUT/$MAF_TAG/raw/SV/SV.Model_A_${MAF_TAG}.gmmat.score.tsv"
DONE="$OUT/summary/tasks/${MAF_TAG}.SV.done"
test -s "$SOURCE.bed"

echo "[$(date '+%F %T %Z')] SV, Model A, MAF >= $MAF on $(hostname)"
if [[ ! -s "$FILTERED.bed" || ! -s "$FILTERED.bim" || ! -s "$FILTERED.fam" ]]; then
  "$PLINK" --bfile "$SOURCE" --chr-set 8 no-xy --allow-no-sex --maf "$MAF" --make-bed --out "$FILTERED"
fi
test -s "$FILTERED.bed"
test -s "$FILTERED.bim"
test -s "$FILTERED.fam"
MARKERS=$(wc -l < "$FILTERED.bim")
if [[ "$MARKERS" -eq 0 ]]; then
  echo "No SV markers remain after MAF filtering" >&2
  exit 1
fi

if [[ ! -s "$RESULT" ]]; then
  "$RSCRIPT" --vanilla "$ROOT/scripts/run_gmmat_gwas_maf.R" \
    "$FILTERED" "$MODEL" "$GRM.rel" "$GRM.rel.id" "$RESULT" \
    "SV_Model_A_${MAF_TAG}" 5 "${SLURM_CPUS_PER_TASK:-2}" "$MAF"
fi
test -s "$RESULT"
awk -v min_maf="$MAF" 'NR > 1 { af=$8 + 0; maf=(af <= 0.5 ? af : 1-af); if (maf + 1e-9 < min_maf) { print "Below requested MAF:", $2, maf > "/dev/stderr"; exit 1 } }' "$RESULT"

tmp_done="$DONE.tmp.$$"
{
  printf 'completed\t%s\n' "$(date '+%F %T %Z')"
  printf 'model\tA: GMMAT binomial mixed model; PC1-PC5; SNP GRM\n'
  printf 'variant_class\tSV\n'
  printf 'maf_cutoff\t%s\n' "$MAF"
  printf 'samples\t%s\n' "$(wc -l < "$FILTERED.fam")"
  printf 'markers_after_plink_maf_filter\t%s\n' "$MARKERS"
  printf 'score_file\t%s\n' "$RESULT"
} > "$tmp_done"
mv "$tmp_done" "$DONE"
echo "[$(date '+%F %T %Z')] complete"
