#!/usr/bin/env bash
set -euo pipefail

sample=${1:?sample required}
threads=${THREADS:-12}
BASE=path/to/project/N_1.coding_gene_anno
PASA_IN=$BASE/04_transcript_prediction/PASA_inputs
PASA_RUN=$BASE/04_transcript_prediction/PASA
LOCKROOT=$PASA_RUN/locks
TEMPLATE=path/to/home/anaconda3/envs/pasa_env/opt/pasa-2.5.3/pasa_conf/pasa.alignAssembly.Template.txt
host=$(hostname)
stamp=$(date +%Y%m%d_%H%M%S)
logroot=$BASE/00_logs_recovery
mkdir -p "$LOCKROOT" "$logroot"
log=$logroot/pasa_one_${sample}_${stamp}_${host}.log
exec > >(tee -a "$log") 2>&1

cd "$BASE"
echo "[$(date)] START PASA sample=$sample host=$host threads=$threads"

trans="$PASA_IN/combined/${sample}.pasa_transcripts.fa"
unmask=$(awk -F'\t' -v s="$sample" 'NR>1 && $1==s{print $3; exit}' "$BASE/01_genome_versions/genome_manifest.tsv")
outdir="$PASA_RUN/$sample"
done_file="$outdir/pasa.done"
fail_file="$outdir/pasa.failed"
lock="$LOCKROOT/${sample}.lock"
archive_root="$PASA_RUN/rerun_archives"

[ -s "$trans" ] || { echo "missing transcript fasta: $trans"; exit 2; }
[ -s "$unmask" ] || { echo "missing genome fasta: $unmask"; exit 3; }
rm -f "$trans.fai"
if [ -s "$done_file" ]; then
  echo "already done: $done_file"
  exit 0
fi
if ! mkdir "$lock" 2>/dev/null; then
  echo "lock exists: $lock"
  exit 0
fi
{
  echo "host=$host"
  echo "pid=$$"
  echo "started=$(date '+%F %T %Z')"
  echo "script=$0"
} > "$lock/owner"
cleanup(){ rm -rf "$lock"; }
trap cleanup EXIT

mkdir -p "$archive_root"
if [ -d "$outdir" ]; then
  mv "$outdir" "$archive_root/${sample}_before_pasa_one_${stamp}"
  echo "archived old PASA dir to $archive_root/${sample}_before_pasa_one_${stamp}"
fi
mkdir -p "$outdir"

source path/to/home/anaconda3/etc/profile.d/conda.sh
conda activate pasa_env

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
echo "[$(date)] Launch_PASA_pipeline sample=$sample transcripts=$nseq genome=$unmask"
set +e
(
  cd "$outdir"
  Launch_PASA_pipeline.pl \
    -c "$conf" -C -R \
    -g "$unmask" \
    -t "$trans" \
    --ALIGNERS gmap,blat \
    --CPU "$threads"
) > "$outdir/pasa.stdout.log" 2> "$outdir/pasa.stderr.log"
rc=$?
set -e

if [ "$rc" -eq 0 ] && find "$outdir" -maxdepth 2 -type f \( -name '*pasa_assemblies*.gff3' -o -name '*assemblies*.fasta' -o -name '*.pasa_assemblies.gff3' \) -size +0c | grep -q .; then
  echo -e "ok\t$(date)\t$host\t$threads\ttranscripts=$nseq" > "$done_file"
  rm -f "$fail_file"
  echo "[$(date)] DONE PASA $sample"
else
  echo -e "failed\trc=$rc\t$(date)\t$host\ttranscripts=$nseq" > "$fail_file"
  echo "[$(date)] FAILED PASA $sample rc=$rc"
  tail -80 "$outdir/pasa.stderr.log" || true
  tail -80 "$outdir/pasa.stdout.log" || true
  exit 10
fi

echo "[$(date)] FINISH PASA sample=$sample log=$log"
