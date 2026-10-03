#!/usr/bin/env bash
set -euo pipefail
BASE=path/to/project/N_1.coding_gene_anno
PASA=$BASE/04_transcript_prediction/PASA_inputs
RNA=$BASE/00_inputs/my_species_RNA
OUT=$PASA/isoseq_existing_fasta
MAN=$PASA/manifests/isoseq_existing_fasta_manifest.tsv
PENDING=$PASA/manifests/isoseq_raw_pending.tsv
mkdir -p "$OUT" "$PASA/manifests"
# Rebuild this directory from curated transcript FASTA products only.
find "$OUT" -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} +
printf 'sample\tinput_fasta\toutput_fasta\tn_seq\tstatus\n' > "$MAN"
printf 'sample\tinput_path\treason\n' > "$PENDING"
clean_fasta(){
  sample=$1; infile=$2; outfile=$3
  case "$infile" in *.gz) gzip -cd "$infile";; *) cat "$infile";; esac | \
  awk -v sample="$sample" 'BEGIN{n=0; RS=">"; ORS=""}
    NR>1{split($0,a,"\n"); h=a[1]; gsub(/[^A-Za-z0-9_.:-]/,"_",h); seq=""; for(i=2;i<=length(a);i++) seq=seq a[i]; gsub(/[[:space:]]/,"",seq); gsub(/[^A-Za-z]/,"",seq); if(length(seq)>=100){n++; print ">"sample"|isoseq_existing|"n"|"h"\n"; for(i=1;i<=length(seq);i+=80) print toupper(substr(seq,i,80))"\n"}}' > "$outfile.tmp"
  mv "$outfile.tmp" "$outfile"
}
for d in "$RNA"/genome_*; do
  [ -d "$d" ] || continue
  sample=$(basename "$d")
  iso="$d/iso-seq"
  [ -d "$iso" ] || continue
  # PASA-ready or near-ready transcript FASTA only. Skip raw all.SMRT.RNA.fa and ONT/PacBio raw fastq.
  mapfile -t fas < <(find "$iso" -maxdepth 1 \( -type f -o -type l \) \( -iname 'unigenes_isoform.fa' -o -iname '*transcript*.fa' -o -iname '*isoform*.fa' \) 2>/dev/null | grep -Evi 'primer|barcode|adapter' | sort)
  if [ "${#fas[@]}" -gt 0 ]; then
    outdir="$OUT/$sample"; mkdir -p "$outdir"
    tmpcat="$outdir/${sample}.isoseq_existing.transcripts.fa.tmpcat"
    : > "$tmpcat"
    for fa in "${fas[@]}"; do
      one="$outdir/$(basename "$fa" | sed 's/\.gz$//').clean.fa"
      clean_fasta "$sample" "$fa" "$one"
      cat "$one" >> "$tmpcat"
    done
    final="$outdir/${sample}.isoseq_existing.transcripts.fa"
    awk 'BEGIN{RS=">"; ORS=""} NR>1{split($0,a,"\n"); h=a[1]; seq=""; for(i=2;i<=length(a);i++) seq=seq a[i]; gsub(/[[:space:]]/,"",seq); if(length(seq)>=100 && !seen[seq]++){print ">"h"\n"; for(i=1;i<=length(seq);i+=80) print substr(seq,i,80)"\n"}}' "$tmpcat" > "$final"
    rm -f "$tmpcat"
    nseq=$(grep -c '^>' "$final" || true)
    printf '%s\t%s\t%s\t%s\tok\n' "$sample" "$(IFS=,; echo "${fas[*]}")" "$final" "$nseq" >> "$MAN"
  fi
  find "$iso" -maxdepth 1 \( -type f -o -type l \) \( -iname '*.bam' -o -iname '*.fq.gz' -o -iname '*.fastq.gz' -o -iname 'all.SMRT.RNA.fa' -o -iname '*all.pacbio.rna.fa' \) 2>/dev/null | while read -r raw; do
    printf '%s\t%s\t%s\n' "$sample" "$raw" "raw_or_unverified_long_read_not_direct_pasa_fasta" >> "$PENDING"
  done
done
