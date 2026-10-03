#!/usr/bin/env bash
set -euo pipefail

STEP="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="${PY:-python3}"

REFA_VCF="$STEP/results/RefA/panSV/panSV.RefA.genotyped.18species.vcf.gz"
REFB_VCF="$STEP/results/RefB/panSV/panSV.RefB.genotyped.18species.vcf.gz"
REFA_PAV="$STEP/results/RefA/panSV/panSV.RefA.PAV.matrix.tsv"
REFB_PAV="$STEP/results/RefB/panSV/panSV.RefB.PAV.matrix.tsv"
PAF="$STEP/results/meta-panSV/RefA_to_RefB.minimap2.paf"
OUT_DIR="$STEP/results/meta-panSV/publication"
LOG_DIR="$STEP/logs"

mkdir -p "$OUT_DIR" "$LOG_DIR" "$STEP/summary"

for f in "$REFA_VCF" "$REFB_VCF" "$REFA_PAV" "$REFB_PAV" "$PAF"; do
  if [[ ! -s "$f" ]]; then
    echo "missing required input: $f" >&2
    exit 1
  fi
done

"$PY" "$STEP/scripts/cross_ref_match_publication.py" \
  --refa-vcf "$REFA_VCF" \
  --refb-vcf "$REFB_VCF" \
  --refa-pav "$REFA_PAV" \
  --refb-pav "$REFB_PAV" \
  --paf "$PAF" \
  --out-dir "$OUT_DIR" \
  --distance "${PUB_DISTANCE:-500}" \
  --min-ro "${PUB_MIN_RO:-0.7}" \
  --max-len-ratio "${PUB_MAX_LEN_RATIO:-1.5}" \
  --ins-distance "${PUB_INS_DISTANCE:-500}" \
  --ins-max-len-ratio "${PUB_INS_MAX_LEN_RATIO:-1.5}" \
  --ins-min-jaccard "${PUB_INS_MIN_JACCARD:-0.35}" \
  --ins-min-containment "${PUB_INS_MIN_CONTAINMENT:-0.70}" \
  --min-pav-jaccard "${PUB_MIN_PAV_JACCARD:-0.50}" \
  --pav-min-presence "${PUB_PAV_MIN_PRESENCE:-2}" \
  --min-mapq "${PUB_MIN_MAPQ:-20}" \
  --max-dv "${PUB_MAX_DV:-0.12}" \
  --min-synteny-block "${PUB_MIN_SYNTENY_BLOCK:-10000}" \
  --max-synteny-blocks "${PUB_MAX_SYNTENY_BLOCKS:-1}" \
  > "$LOG_DIR/publication_meta_panSV.out" \
  2> "$LOG_DIR/publication_meta_panSV.err"

touch "$STEP/summary/publication_meta_panSV.done"
