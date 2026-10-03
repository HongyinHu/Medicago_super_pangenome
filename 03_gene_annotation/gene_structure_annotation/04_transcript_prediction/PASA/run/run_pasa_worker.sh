#!/usr/bin/env bash
set -euo pipefail

THREADS=${THREADS:-12}
BASE=path/to/project/N_1.coding_gene_anno
PASA_IN=$BASE/04_transcript_prediction/PASA_inputs
PASA_RUN=$BASE/04_transcript_prediction/PASA
LOCKROOT=$PASA_RUN/locks
LOGROOT=$PASA_RUN/logs
TEMPLATE=path/to/home/anaconda3/envs/pasa_env/opt/pasa-2.5.3/pasa_conf/pasa.alignAssembly.Template.txt

mkdir -p "$PASA_RUN" "$LOCKROOT" "$LOGROOT"
source path/to/home/anaconda3/etc/profile.d/conda.sh
conda activate pasa_env

host=$(hostname)
log="$LOGROOT/pasa_worker_${host}_$(date +%Y%m%d_%H%M%S).log"
exec >> "$log" 2>&1

echo "[$(date)] START PASA worker host=$host threads=$THREADS"

tail -n +2 "$BASE/01_genome_versions/genome_manifest.tsv" | while IFS=$'\t' read -r sample soft unmask rest; do
  [ -n "$sample" ] || continue
  trans="$PASA_IN/combined/${sample}.pasa_transcripts.fa"
  outdir="$PASA_RUN/$sample"
  done="$outdir/pasa.done"
  fail="$outdir/pasa.failed"
  lock="$LOCKROOT/${sample}.lock"

  [ -s "$done" ] && continue
  [ -s "$trans" ] || continue
  [ -s "$unmask" ] || continue
  if ! mkdir "$lock" 2>/dev/null; then
    continue
  fi
  cleanup(){ rm -rf "$lock"; }
  trap cleanup EXIT

  mkdir -p "$outdir"
  rm -f "$fail"
  db="$outdir/${sample}.sqlite"
  conf="$outdir/${sample}.alignAssembly.config"
  if [ -s "$TEMPLATE" ]; then
    sed \
      -e "s#<__DATABASE__>#$db#g" \
      -e 's#<__MIN_PERCENT_ALIGNED__>#90#g' \
      -e 's#<__MIN_AVG_PER_ID__>#95#g' \
      "$TEMPLATE" > "$conf"
  else
    {
      echo "DATABASE=$db"
      echo "validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=90"
      echo "validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=95"
      echo "subcluster_builder.dbi:-m=50"
    } > "$conf"
  fi

  nseq=$(grep -c '^>' "$trans" 2>/dev/null || echo 0)
  echo "[$(date)] PASA sample=$sample transcripts=$nseq genome=$unmask outdir=$outdir"
  rm -f "$db"
  set +e
  (
    cd "$outdir"
    Launch_PASA_pipeline.pl \
      -c "$conf" -C -R \
      -g "$unmask" \
      -t "$trans" \
      --ALIGNERS gmap,blat \
      --CPU "$THREADS"
  ) > "$outdir/pasa.stdout.log" 2> "$outdir/pasa.stderr.log"
  rc=$?
  set -e

  if [ "$rc" -eq 0 ] && find "$outdir" -maxdepth 2 -type f \( -name '*pasa_assemblies*.gff3' -o -name '*assemblies*.fasta' -o -name '*.pasa_assemblies.gff3' \) -size +0c | grep -q .; then
    echo -e "ok\t$(date)\t$host\t$THREADS\ttranscripts=$nseq" > "$done"
    echo "[$(date)] DONE PASA $sample"
  else
    echo -e "failed\trc=$rc\t$(date)\t$host\ttranscripts=$nseq" > "$fail"
    echo "[$(date)] FAILED PASA $sample rc=$rc"
    tail -60 "$outdir/pasa.stderr.log" || true
    tail -60 "$outdir/pasa.stdout.log" || true
  fi
  cleanup
  trap - EXIT
done

echo "[$(date)] FINISH PASA worker host=$host"
