#!/usr/bin/env bash
set -euo pipefail

STEP="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="${PY:-python3}"

REFA_VCF="$STEP/results/RefA/panSV/panSV.RefA.genotyped.18species.vcf.gz"
REFB_VCF="$STEP/results/RefB/panSV/panSV.RefB.genotyped.18species.vcf.gz"
REFA_PAV="$STEP/results/RefA/panSV/panSV.RefA.PAV.matrix.tsv"
REFB_PAV="$STEP/results/RefB/panSV/panSV.RefB.PAV.matrix.tsv"
PAF="$STEP/results/meta-panSV/RefA_to_RefB.minimap2.paf"
OUT_DIR="$STEP/results/meta-panSV/strict"
LOG_DIR="$STEP/logs"

mkdir -p "$OUT_DIR" "$LOG_DIR" "$STEP/summary"

for f in "$REFA_VCF" "$REFB_VCF" "$REFA_PAV" "$REFB_PAV" "$PAF"; do
  if [[ ! -s "$f" ]]; then
    echo "missing required input: $f" >&2
    exit 1
  fi
done

"$PY" "$STEP/scripts/cross_ref_match_strict.py" \
  --refa-vcf "$REFA_VCF" \
  --refb-vcf "$REFB_VCF" \
  --refa-pav "$REFA_PAV" \
  --refb-pav "$REFB_PAV" \
  --paf "$PAF" \
  --out-dir "$OUT_DIR" \
  --distance "${STRICT_DISTANCE:-1000}" \
  --min-ro "${STRICT_MIN_RO:-0.5}" \
  --max-len-ratio "${STRICT_MAX_LEN_RATIO:-2.0}" \
  --min-pav-jaccard "${STRICT_MIN_PAV_JACCARD:-0.30}" \
  --pav-min-presence "${STRICT_PAV_MIN_PRESENCE:-2}" \
  > "$LOG_DIR/strict_meta_panSV.out" \
  2> "$LOG_DIR/strict_meta_panSV.err"

touch "$STEP/summary/strict_meta_panSV.done"
