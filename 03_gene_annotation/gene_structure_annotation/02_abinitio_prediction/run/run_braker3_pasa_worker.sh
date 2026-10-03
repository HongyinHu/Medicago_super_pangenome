#!/usr/bin/env bash
set -euo pipefail

THREADS=${THREADS:-16}
BASE=path/to/project/N_1.coding_gene_anno
AB=$BASE/02_abinitio_prediction
PASA_RUN=$BASE/04_transcript_prediction/PASA
PASA_IN=$BASE/04_transcript_prediction/PASA_inputs
OUTROOT=$AB/BRAKER3_PASA_RNA_protein
LOCKROOT=$AB/locks_pasa_braker
LOGROOT=$AB/logs
PROT=$AB/inputs/related_species.clean.pep.fa

mkdir -p "$OUTROOT" "$LOCKROOT" "$LOGROOT"
source path/to/home/anaconda3/etc/profile.d/conda.sh
conda activate braker3_env

export HOME=path/to/home
export AUGUSTUS_CONFIG_PATH=$AB/augustus_config
export GENEMARK_PATH=path/to/home/sofeware/gmes_linux_64_4
export PATH=$GENEMARK_PATH:$PATH
export PERL5LIB=path/to/home/anaconda3/envs/braker3_env/lib/perl5/site_perl:path/to/home/anaconda3/envs/braker3_env/lib/perl5/vendor_perl:path/to/home/anaconda3/envs/braker3_env/lib/perl5/core_perl${PERL5LIB:+:$PERL5LIB}

host=$(hostname)
log="$LOGROOT/braker3_pasa_worker_${host}_$(date +%Y%m%d_%H%M%S).log"
exec >> "$log" 2>&1
echo "[$(date)] START BRAKER3_PASA worker host=$host threads=$THREADS"

if [ ! -s "$PROT" ]; then
  echo "[$(date)] build related protein set"
  bash "$AB/run/build_related_protein_for_braker3.sh"
fi

find_pasa_transcripts() {
  local sample="$1"
  local outdir="$PASA_RUN/$sample"
  find "$outdir" -maxdepth 2 -type f \( -name '*assemblies*.fasta' -o -name '*pasa*.fasta' -o -name '*.assemblies.fasta' \) -size +0c 2>/dev/null | sort | head -1
}

tail -n +2 "$BASE/01_genome_versions/genome_manifest.tsv" | while IFS=$'\t' read -r sample soft unmask rest; do
  [ -n "$sample" ] || continue
  pasa_done="$PASA_RUN/$sample/pasa.done"
  [ -s "$pasa_done" ] || continue
  [ -s "$soft" ] || continue

  outdir="$OUTROOT/$sample"
  done="$outdir/braker3_pasa.done"
  fail="$outdir/braker3_pasa.failed"
  lock="$LOCKROOT/${sample}.lock"
  [ -s "$done" ] && continue
  if ! mkdir "$lock" 2>/dev/null; then
    continue
  fi
  cleanup(){ rm -rf "$lock"; }
  trap cleanup EXIT

  mkdir -p "$outdir/transcript_bam"
  rm -f "$fail"

  trans=$(find_pasa_transcripts "$sample")
  if [ -z "$trans" ] || [ ! -s "$trans" ]; then
    trans="$PASA_IN/combined/${sample}.pasa_transcripts.fa"
  fi
  if [ ! -s "$trans" ]; then
    echo -e "failed\tmissing_pasa_transcripts\t$(date)\t$host" > "$fail"
    cleanup; trap - EXIT; continue
  fi

  bam="$outdir/transcript_bam/${sample}.pasa_transcripts.minimap2.sorted.bam"
  if [ ! -s "$bam" ]; then
    echo "[$(date)] minimap2 transcript BAM sample=$sample trans=$trans"
    minimap2 -t "$THREADS" -ax splice:hq --secondary=no "$soft" "$trans" \
      | samtools sort -@ "$THREADS" -o "$bam"
    samtools index "$bam"
  fi

  species="medpan_${sample}_pasa"
  echo "[$(date)] BRAKER3 sample=$sample species=$species genome=$soft bam=$bam prot=$PROT"
  set +e
  braker.pl \
    --genome="$soft" \
    --species="$species" \
    --prot_seq="$PROT" \
    --bam="$bam" \
    --softmasking \
    --gff3 \
    --threads="$THREADS" \
    --workingdir="$outdir" \
    > "$outdir/braker3_pasa.stdout.log" 2> "$outdir/braker3_pasa.stderr.log"
  rc=$?
  set -e

  if [ "$rc" -eq 0 ] && { [ -s "$outdir/braker.gff3" ] || [ -s "$outdir/augustus.hints.gff3" ] || [ -s "$outdir/augustus.ab_initio.gff3" ]; }; then
    echo -e "ok\t$(date)\t$host\t$THREADS\tbam=$bam\tprot=$PROT\ttranscripts=$trans" > "$done"
    echo "[$(date)] DONE BRAKER3 $sample"
  else
    echo -e "failed\trc=$rc\t$(date)\t$host\tbam=$bam\tprot=$PROT\ttranscripts=$trans" > "$fail"
    echo "[$(date)] FAILED BRAKER3 $sample rc=$rc"
    tail -80 "$outdir/braker3_pasa.stderr.log" || true
  fi
  cleanup
  trap - EXIT
done

echo "[$(date)] FINISH BRAKER3_PASA worker host=$host"
