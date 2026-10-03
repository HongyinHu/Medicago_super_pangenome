#!/usr/bin/env bash
set -euo pipefail
BASE=path/to/project/N_1.coding_gene_anno
PASA=$BASE/04_transcript_prediction/PASA_inputs
OUT=$PASA/combined
MAN=$PASA/manifests/pasa_transcripts_manifest.tsv
mkdir -p "$OUT" "$PASA/manifests"
printf 'sample\tcombined_fasta\tn_total\tn_illumina\tn_isoseq_existing\tstatus\tsources\n' > "$MAN"
tail -n +2 "$BASE/01_genome_versions/genome_manifest.tsv" | while IFS=$'\t' read -r sample soft unmask rest; do
  [ -n "$sample" ] || continue
  sources=()
  illum="$PASA/stringtie/$sample/${sample}.illumina_stringtie.transcripts.fa"
  iso="$PASA/isoseq_existing_fasta/$sample/${sample}.isoseq_existing.transcripts.fa"
  [ -s "$illum" ] && sources+=("$illum")
  [ -s "$iso" ] && sources+=("$iso")
  out="$OUT/${sample}.pasa_transcripts.fa"
  if [ "${#sources[@]}" -eq 0 ]; then
    printf '%s\t%s\t0\t0\t0\tmissing\t\n' "$sample" "$out" >> "$MAN"
    continue
  fi
  tmp="$out.tmp"
  : > "$tmp"
  for f in "${sources[@]}"; do cat "$f" >> "$tmp"; done
  awk 'BEGIN{RS=">"; ORS=""} NR>1{split($0,a,"\n"); h=a[1]; seq=""; for(i=2;i<=length(a);i++) seq=seq a[i]; gsub(/[[:space:]]/,"",seq); if(length(seq)>=100 && !seen[seq]++){print ">"h"\n"; for(i=1;i<=length(seq);i+=80) print substr(seq,i,80)"\n"}}' "$tmp" > "$out"
  rm -f "$tmp"
  ntot=$(grep -c '^>' "$out" || true)
  nillum=0; niso=0
  [ -s "$illum" ] && nillum=$(grep -c '^>' "$illum" || true)
  [ -s "$iso" ] && niso=$(grep -c '^>' "$iso" || true)
  printf '%s\t%s\t%s\t%s\t%s\tok\t%s\n' "$sample" "$out" "$ntot" "$nillum" "$niso" "$(IFS=,; echo "${sources[*]}")" >> "$MAN"
done
