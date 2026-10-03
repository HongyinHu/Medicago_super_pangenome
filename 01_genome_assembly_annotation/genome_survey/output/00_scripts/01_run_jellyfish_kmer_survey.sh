#!/usr/bin/env bash
set -euo pipefail
export PATH="path/to/home/anaconda3/bin:$PATH"
command -v jellyfish >/dev/null || { echo "ERROR: jellyfish not found in PATH" >&2; exit 127; }
BASE=path/to/project/N_2.DNA_survy
DATA="$BASE/data"
OUT="$BASE/output"
SAMPLE_TABLE="$OUT/01_sample_table/samples.tsv"
STATUS="$OUT/03_estimates/pipeline_status.tsv"
K_LIST="${K_LIST:-21}"
THREADS="${THREADS:-32}"
JF_SIZE="${JF_SIZE:-8G}"
KEEP_JF="${KEEP_JF:-0}"
mkdir -p "$OUT/01_sample_table" "$OUT/02_jellyfish" "$OUT/03_estimates" "$OUT/logs"

printf "sample\tR1\tR2\n" > "$SAMPLE_TABLE"
for r1 in "$DATA"/*_R1.fq.gz; do
  [ -e "$r1" ] || continue
  sample=$(basename "$r1" _R1.fq.gz)
  r2="$DATA/${sample}_R2.fq.gz"
  if [ -e "$r2" ]; then
    printf "%s\t%s\t%s\n" "$sample" "$r1" "$r2" >> "$SAMPLE_TABLE"
  else
    printf "Missing R2 for %s\n" "$sample" >&2
  fi
done

printf "time\tsample\tk\tstatus\tnote\n" > "$STATUS"
while IFS=$'\t' read -r sample r1 r2; do
  [ "$sample" = "sample" ] && continue
  for k in $K_LIST; do
    kdir="$OUT/02_jellyfish/k${k}"
    mkdir -p "$kdir"
    jf="$kdir/${sample}.k${k}.jf"
    hist="$kdir/${sample}.k${k}.histo"
    log="$OUT/logs/${sample}.k${k}.log"
    if [ -s "$hist" ]; then
      printf "%s\t%s\t%s\tskipped\texisting_histo\n" "$(date '+%F %T')" "$sample" "$k" >> "$STATUS"
      continue
    fi
    printf "%s\t%s\t%s\tstart\tthreads=%s size=%s\n" "$(date '+%F %T')" "$sample" "$k" "$THREADS" "$JF_SIZE" >> "$STATUS"
    {
      echo "[$(date '+%F %T')] START sample=$sample k=$k"
      echo "R1=$r1"
      echo "R2=$r2"
      echo "THREADS=$THREADS JF_SIZE=$JF_SIZE"
      jellyfish count -C -m "$k" -s "$JF_SIZE" -t "$THREADS" -o "$jf" <(gzip -cd "$r1" "$r2")
      jellyfish histo -t "$THREADS" "$jf" > "$hist"
      if [ "$KEEP_JF" = "0" ]; then
        rm -f "$jf"
      fi
      echo "[$(date '+%F %T')] DONE sample=$sample k=$k hist=$hist"
    } > "$log" 2>&1
    printf "%s\t%s\t%s\tdone\t%s\n" "$(date '+%F %T')" "$sample" "$k" "$hist" >> "$STATUS"
    python3 "$OUT/00_scripts/02_estimate_genome_size.py" "$OUT/03_estimates/genome_size_estimates.tsv" "$OUT/02_jellyfish"/k*/*.histo || true
    Rscript "$OUT/00_scripts/03_plot_histograms.R" "$OUT/03_estimates/histogram_png" "$OUT/02_jellyfish"/k*/*.histo >/dev/null 2>&1 || true
    Rscript "$OUT/00_scripts/05_plot_histograms_pdf.R" "$OUT/03_estimates/histogram_pdf" "$OUT/02_jellyfish"/k*/*.histo >/dev/null 2>&1 || true
  done
done < "$SAMPLE_TABLE"
python3 "$OUT/00_scripts/02_estimate_genome_size.py" "$OUT/03_estimates/genome_size_estimates.tsv" "$OUT/02_jellyfish"/k*/*.histo
Rscript "$OUT/00_scripts/03_plot_histograms.R" "$OUT/03_estimates/histogram_png" "$OUT/02_jellyfish"/k*/*.histo >/dev/null 2>&1 || true
    Rscript "$OUT/00_scripts/05_plot_histograms_pdf.R" "$OUT/03_estimates/histogram_pdf" "$OUT/02_jellyfish"/k*/*.histo >/dev/null 2>&1 || true
printf "%s\tALL\t%s\tfinished\tall_histograms_processed\n" "$(date '+%F %T')" "$K_LIST" >> "$STATUS"
