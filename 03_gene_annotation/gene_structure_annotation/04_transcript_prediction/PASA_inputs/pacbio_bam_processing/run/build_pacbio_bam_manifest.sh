#!/usr/bin/env bash
set -euo pipefail

BASE=path/to/project/N_1.coding_gene_anno
RNA="$BASE/00_inputs/my_species_RNA"
OUT="$BASE/04_transcript_prediction/PASA_inputs/pacbio_bam_processing"

mkdir -p "$OUT"/{run,logs,manifests,work,locks,status,final}

MAN="$OUT/manifests/pacbio_bam_manifest.tsv"
TMP="$MAN.tmp"
COMMON_PRIMER="$OUT/manifests/common_NEB_Clontech_primer.fa"
COMMON_SOURCE="$RNA/genome_410/iso-seq/primer.fa"

if [ -e "$COMMON_SOURCE" ]; then
  cp -L "$COMMON_SOURCE" "$COMMON_PRIMER"
fi

printf 'sample\tinput_bam\treal_bam\tsize_bytes\tpbi\tprimer\treal_primer\tstate\tnote\n' > "$TMP"

declare -A first_by_key

find "$RNA" -mindepth 3 -maxdepth 3 \( -type f -o -type l \) -path '*/iso-seq/*.bam' | sort | while read -r bam; do
  sample=$(echo "$bam" | awk -F/ '{for (i=1; i<=NF; i++) if ($i ~ /^genome_/) print $i}' | tail -1)
  real_bam=$(readlink -f "$bam" 2>/dev/null || echo "$bam")
  size_bytes=$(stat -Lc '%s' "$bam" 2>/dev/null || echo 0)
  pbi=no
  [ -e "$bam.pbi" ] && pbi=yes

  primer=""
  if [ -e "$RNA/$sample/iso-seq/primer.fa" ]; then
    primer="$RNA/$sample/iso-seq/primer.fa"
  elif [ "${ALLOW_COMMON_PRIMER_FOR_MISSING:-0}" = "1" ] && [ -s "$COMMON_PRIMER" ]; then
    primer="$COMMON_PRIMER"
  fi
  real_primer=""
  if [ -n "$primer" ]; then
    real_primer=$(readlink -f "$primer" 2>/dev/null || echo "$primer")
  fi

  state=missing_primer
  note=needs_primer_for_lima_refine
  if [ -n "$real_primer" ]; then
    key="${real_bam}|${real_primer}"
    if [ -n "${first_by_key[$key]:-}" ]; then
      state=duplicate_ready
      note="duplicate_of=${first_by_key[$key]}"
    else
      state=ready_full
      if [ "$primer" = "$COMMON_PRIMER" ] && [ ! -e "$RNA/$sample/iso-seq/primer.fa" ]; then
        note=ccs_lima_refine_cluster2_common_primer_inferred_by_motif_scan
      else
        note=ccs_lima_refine_cluster2
      fi
      first_by_key[$key]="$sample"
    fi
  fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$sample" "$bam" "$real_bam" "$size_bytes" "$pbi" "$primer" "$real_primer" "$state" "$note" >> "$TMP"
done

mv "$TMP" "$MAN"

awk -F'\t' 'NR > 1 && $8 == "ready_full" {print $1}' "$MAN" > "$OUT/manifests/ready_full.samples"
awk -F'\t' 'NR > 1 && $8 == "duplicate_ready" {print $1"\t"$9}' "$MAN" > "$OUT/manifests/duplicate_ready.samples"
awk -F'\t' 'NR > 1 && $8 == "missing_primer" {print $1}' "$MAN" > "$OUT/manifests/missing_primer.samples"

echo "manifest=$MAN"
awk -F'\t' 'NR > 1 {count[$8]++} END {for (s in count) print s, count[s]}' "$MAN" | sort
