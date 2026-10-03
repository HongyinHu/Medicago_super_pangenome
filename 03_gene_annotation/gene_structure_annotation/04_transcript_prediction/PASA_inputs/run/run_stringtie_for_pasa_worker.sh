#!/usr/bin/env bash
set -euo pipefail
THREADS=${THREADS:-16}
BASE=path/to/project/N_1.coding_gene_anno
PASA=$BASE/04_transcript_prediction/PASA_inputs
OUT=$PASA/stringtie
LOCK=$PASA/locks
LOG=$PASA/logs
mkdir -p "$OUT" "$LOCK" "$LOG"
source path/to/home/anaconda3/etc/profile.d/conda.sh
conda activate braker3_env
host=$(hostname)
exec >> "$LOG/stringtie_for_pasa_worker_${host}_$(date +%Y%m%d_%H%M%S).log" 2>&1
echo "[$(date)] START host=$host threads=$THREADS"
tail -n +2 "$BASE/01_genome_versions/genome_manifest.tsv" | while IFS=$'\t' read -r sample soft unmask rest; do
  [ -n "$sample" ] || continue
  outdir="$OUT/$sample"
  done="$outdir/stringtie_for_pasa.done"
  fail="$outdir/stringtie_for_pasa.failed"
  lockdir="$LOCK/${sample}.stringtie.lock"
  [ -s "$done" ] && continue
  mapfile -t bams < <(find "$BASE/04_transcript_prediction/illumina_align/$sample" -maxdepth 1 -type f -name '*.hisat2.sorted.bam' -size +0c 2>/dev/null | sort)
  [ "${#bams[@]}" -gt 0 ] || continue
  if ! mkdir "$lockdir" 2>/dev/null; then continue; fi
  cleanup(){ rm -rf "$lockdir"; }
  trap cleanup EXIT
  mkdir -p "$outdir/individual"
  rm -f "$fail"
  echo "[$(date)] STRINGTIE sample=$sample n_bam=${#bams[@]} genome=$unmask"
  gtf_list="$outdir/gtf.list"
  : > "$gtf_list"
  ok=1
  for bam in "${bams[@]}"; do
    lib=$(basename "$bam" .hisat2.sorted.bam)
    gtf="$outdir/individual/${lib}.stringtie.gtf"
    if [ ! -s "$gtf" ]; then
      stringtie -p "$THREADS" -l "${sample}_${lib}" -o "$gtf" "$bam" || ok=0
    fi
    [ -s "$gtf" ] && echo "$gtf" >> "$gtf_list"
  done
  if [ "$ok" -ne 1 ] || [ ! -s "$gtf_list" ]; then
    echo -e "failed\tstringtie_individual\t$(date)\t$host" > "$fail"
    cleanup; trap - EXIT; continue
  fi
  n_gtf=$(wc -l < "$gtf_list")
  final_gtf="$outdir/${sample}.stringtie_for_pasa.gtf"
  if [ "$n_gtf" -gt 1 ]; then
    stringtie --merge -p "$THREADS" -o "$final_gtf" "$gtf_list" || { echo -e "failed\tstringtie_merge\t$(date)\t$host" > "$fail"; cleanup; trap - EXIT; continue; }
  else
    cp "$(cat "$gtf_list")" "$final_gtf"
  fi
  raw_fa="$outdir/${sample}.illumina_stringtie.raw.fa"
  clean_fa="$outdir/${sample}.illumina_stringtie.transcripts.fa"
  gffread "$final_gtf" -g "$unmask" -w "$raw_fa"
  awk -v sample="$sample" 'BEGIN{n=0; RS=">"; ORS=""} NR>1{split($0,a,"\n"); h=a[1]; gsub(/[^A-Za-z0-9_.:-]/,"_",h); seq=""; for(i=2;i<=length(a);i++) seq=seq a[i]; gsub(/[[:space:]]/,"",seq); if(length(seq)>=100){n++; print ">"sample"|illumina_stringtie|"n"|"h"\n"; for(i=1;i<=length(seq);i+=80) print toupper(substr(seq,i,80))"\n"}}' "$raw_fa" > "$clean_fa"
  nseq=$(grep -c '^>' "$clean_fa" || true)
  echo -e "ok\t$(date)\t$host\t$THREADS\tn_bam=${#bams[@]}\tn_transcripts=$nseq\t$clean_fa" > "$done"
  echo "[$(date)] DONE sample=$sample n_transcripts=$nseq"
  cleanup
  trap - EXIT
done
echo "[$(date)] FINISH host=$host"
