#!/usr/bin/env bash
set -euo pipefail
THREADS=${THREADS:-24}
SAMPLE=genome_Msa1
BASE=path/to/project/N_1.coding_gene_anno
AB3=$BASE/02_abinitio_prediction/BRAKER3_RNA_protein
PROT=$AB3/inputs/related_species.all_refs.gffread_completed.sanitized.no_stop.pep.fa
SOFT=$BASE/01_genome_versions/softmask/${SAMPLE}.softmasked.fa
BAMDIR=$BASE/04_transcript_prediction/illumina_align/$SAMPLE
BAMS=$(printf '%s,' \
  "$BAMDIR/jingA.hisat2.sorted.bam" \
  "$BAMDIR/jingB.hisat2.sorted.bam" \
  "$BAMDIR/jingC.hisat2.sorted.bam" \
  "$BAMDIR/yeA.hisat2.sorted.bam" \
  "$BAMDIR/yeB.hisat2.sorted.bam" \
  "$BAMDIR/yeC.hisat2.sorted.bam" | sed 's/,$//')
LOCKROOT=$AB3/locks
LOGROOT=$AB3/logs
OUTDIR=$AB3/$SAMPLE
mkdir -p "$LOCKROOT" "$LOGROOT" "$OUTDIR"
LOG=$LOGROOT/${SAMPLE}_rna_protein_one_$(hostname)_$(date +%Y%m%d_%H%M%S).log
exec >> "$LOG" 2>&1
log(){ echo "[$(date '+%F %T %Z')] $*"; }
fail(){ echo -e "failed\t$*\t$(date '+%F %T %Z')\t$(hostname)" > "$OUTDIR/braker3_rna_protein.failed"; log "FAILED: $*"; exit 1; }
log "START sample=$SAMPLE host=$(hostname) threads=$THREADS"
[ -s "$PROT" ] || fail "missing protein FASTA $PROT"
[ -s "$SOFT" ] || fail "missing softmask FASTA $SOFT"
for b in ${BAMS//,/ }; do [ -s "$b" ] || fail "missing BAM $b"; done
if [ -s "$OUTDIR/braker3_rna_protein.done" ]; then log "already done"; exit 0; fi
if [ -s "$OUTDIR/braker3_rna_protein.failed" ]; then
  mkdir -p "$AB3/failed_attempts/${SAMPLE}_before_one_rerun_$(date +%Y%m%d_%H%M%S)"
  cp -a "$OUTDIR/braker3_rna_protein.failed" "$OUTDIR/braker3_rna_protein.stderr.log" "$OUTDIR/braker3_rna_protein.stdout.log" "$AB3/failed_attempts/${SAMPLE}_before_one_rerun_$(date +%Y%m%d_%H%M%S)/" 2>/dev/null || true
  rm -f "$OUTDIR/braker3_rna_protein.failed"
fi
LOCK=$LOCKROOT/${SAMPLE}.lock
if ! mkdir "$LOCK" 2>/dev/null; then log "locked: $LOCK"; [ -f "$LOCK/owner" ] && cat "$LOCK/owner"; exit 0; fi
cleanup(){ rm -rf "$LOCK"; }
trap cleanup EXIT
{
  echo "host=$(hostname)"
  echo "pid=$$"
  echo "started=$(date '+%F %T %Z')"
  echo "script=$0"
  echo "threads=$THREADS"
  echo "bam=$BAMS"
} > "$LOCK/owner"
source path/to/home/anaconda3/etc/profile.d/conda.sh
conda activate braker3_env
export AUGUSTUS_CONFIG_PATH=$AB3/augustus_config
export GENEMARK_PATH=path/to/home/sofeware/GeneMark-ETP-main/bin
export PROTHINT_PATH=path/to/home/sofeware/GeneMark-ETP-main/bin/gmes/ProtHint/bin
export PATH=$GENEMARK_PATH:$PROTHINT_PATH:$GENEMARK_PATH/gmes:$PATH
export PERL5LIB=path/to/home/anaconda3/envs/braker3_env/lib/perl5/site_perl:path/to/home/anaconda3/envs/braker3_env/lib/perl5/vendor_perl:path/to/home/anaconda3/envs/braker3_env/lib/perl5/core_perl${PERL5LIB:+:$PERL5LIB}
species="medpan_${SAMPLE}_braker3rp"
rm -rf "$AB3/augustus_config/species/$species"
log "RUN braker sample=$SAMPLE species=$species bams=$BAMS prot=$PROT soft=$SOFT"
set +e
braker.pl --gff3 --softmasking --threads "$THREADS" --prot_seq="$PROT" --genome="$SOFT" --bam="$BAMS" --species="$species" --AUGUSTUS_ab_initio --workingdir="$OUTDIR" > "$OUTDIR/braker3_rna_protein.stdout.log" 2> "$OUTDIR/braker3_rna_protein.stderr.log"
rc=$?
set -e
if [ "$rc" -eq 0 ] && { [ -s "$OUTDIR/braker.gff3" ] || [ -s "$OUTDIR/Augustus/augustus.hints.gff3" ] || [ -s "$OUTDIR/Augustus/augustus.ab_initio.gff3" ]; }; then
  echo -e "ok\t$(date '+%F %T %Z')\t$(hostname)\t$THREADS\t6\t$BAMS" > "$OUTDIR/braker3_rna_protein.done"
  log "DONE $SAMPLE"
else
  echo -e "failed\trc=$rc\t$(date '+%F %T %Z')\t$(hostname)\t6\t$BAMS" > "$OUTDIR/braker3_rna_protein.failed"
  log "FAILED $SAMPLE rc=$rc"
  tail -120 "$OUTDIR/braker3_rna_protein.stderr.log" || true
  exit "$rc"
fi
