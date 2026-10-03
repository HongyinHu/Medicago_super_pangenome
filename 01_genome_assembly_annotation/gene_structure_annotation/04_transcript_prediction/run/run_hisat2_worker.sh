#!/usr/bin/env bash
set -euo pipefail
THREADS=${THREADS:-16}
BASE=path/to/project/N_1.coding_gene_anno
WORK=$BASE/04_transcript_prediction
GENOMES=$BASE/01_genome_versions/unmask
IDXDIR=$WORK/indices/hisat2
OUTROOT=$WORK/illumina_align
LOCKROOT=$WORK/locks
LOGROOT=$WORK/logs
mkdir -p "$IDXDIR" "$OUTROOT" "$LOCKROOT" "$LOGROOT"
source path/to/home/anaconda3/etc/profile.d/conda.sh
conda activate trans_analysis
host=$(hostname)
log="$LOGROOT/hisat2_worker_${host}_$(date +%Y%m%d_%H%M%S).log"
exec >> "$log" 2>&1
echo "[$(date)] START host=$host threads=$THREADS"
tail -n +2 "$WORK/rnaseq_pairs.tsv" | while IFS=$'\t' read -r sample pair r1 r2; do
  [ -n "$sample" ] || continue
  safe_pair=$(echo "$pair" | sed 's/[^A-Za-z0-9_.-]/_/g')
  outdir="$OUTROOT/$sample"
  done="$outdir/${safe_pair}.done"
  bam="$outdir/${safe_pair}.hisat2.sorted.bam"
  fail="$outdir/${safe_pair}.failed"
  lock="$LOCKROOT/${sample}__${safe_pair}.lock"
  [ -s "$done" ] && continue
  if ! mkdir "$lock" 2>/dev/null; then
    continue
  fi
  cleanup(){ rm -rf "$lock"; }
  trap cleanup EXIT
  mkdir -p "$outdir"
  rm -f "$fail"
  genome="$GENOMES/${sample}.unmasked.fa"
  prefix="$IDXDIR/${sample}"
  if [ ! -s "$genome" ]; then
    echo "[$(date)] MISSING genome $genome"; echo missing_genome > "$fail"; cleanup; trap - EXIT; continue
  fi
  if [ ! -s "${prefix}.1.ht2" ] && [ ! -s "${prefix}.1.ht2l" ]; then
    ilock="$LOCKROOT/${sample}.hisat2_index.lock"
    if mkdir "$ilock" 2>/dev/null; then
      echo "[$(date)] BUILD_INDEX $sample"
      if hisat2-build -p "$THREADS" "$genome" "$prefix"; then
        rm -rf "$ilock"
      else
        rm -rf "$ilock"
        echo index_failed > "$fail"
        echo "[$(date)] INDEX_FAILED $sample"
        cleanup; trap - EXIT; continue
      fi
    else
      echo "[$(date)] WAIT_INDEX $sample"
      for i in $(seq 1 720); do
        if [ -s "${prefix}.1.ht2" ] || [ -s "${prefix}.1.ht2l" ]; then break; fi
        sleep 30
      done
    fi
  fi
  if [ ! -s "${prefix}.1.ht2" ] && [ ! -s "${prefix}.1.ht2l" ]; then
    echo "[$(date)] INDEX_FAILED_OR_TIMEOUT $sample"; echo index_failed_or_timeout > "$fail"; cleanup; trap - EXIT; continue
  fi
  echo "[$(date)] ALIGN sample=$sample pair=$pair r1=$r1 r2=$r2"
  set +e
  hisat2 --dta -p "$THREADS" -x "$prefix" -1 "$r1" -2 "$r2" 2> "$outdir/${safe_pair}.hisat2.log" | samtools sort -@ 4 -o "$bam.tmp" -
  status=("${PIPESTATUS[@]}")
  set -e
  rc=${status[0]:-99}
  sort_rc=${status[1]:-99}
  if [ "$rc" -eq 0 ] && [ "$sort_rc" -eq 0 ] && [ -s "$bam.tmp" ]; then
    mv "$bam.tmp" "$bam"
    samtools index -@ 4 "$bam"
    echo -e "ok\t$(date)\t$host\t$THREADS\t$bam" > "$done"
    echo "[$(date)] DONE $sample $pair"
  else
    rm -f "$bam.tmp"
    echo -e "failed\thisat2_rc=$rc\tsort_rc=$sort_rc\t$(date)" > "$fail"
    echo "[$(date)] FAILED $sample $pair rc=$rc sort=$sort_rc"
  fi
  cleanup
  trap - EXIT
done
echo "[$(date)] FINISH host=$host"
