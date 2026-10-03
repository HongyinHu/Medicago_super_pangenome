#!/usr/bin/env bash
set -euo pipefail
THREADS=${THREADS:-30}
BASE=path/to/project/N_1.coding_gene_anno
AB3=$BASE/02_abinitio_prediction/BRAKER3_RNA_protein
MAN=$AB3/inputs/braker3_rna_protein_bam_manifest.tsv
PROT=$AB3/inputs/related_species.all_refs.gffread_completed.sanitized.pep.fa
LOCKROOT=$AB3/locks
LOGROOT=$AB3/logs
mkdir -p "$LOCKROOT" "$LOGROOT"
source path/to/home/anaconda3/etc/profile.d/conda.sh
conda activate braker3_env
export AUGUSTUS_CONFIG_PATH=$AB3/augustus_config
export GENEMARK_PATH=path/to/home/sofeware/gmes_linux_64_4
export PROTHINT_PATH=path/to/home/sofeware/gmes_linux_64_4/ProtHint/bin
export PATH=$PROTHINT_PATH:$GENEMARK_PATH:$PATH
export PERL5LIB=path/to/home/anaconda3/envs/braker3_env/lib/perl5/site_perl:path/to/home/anaconda3/envs/braker3_env/lib/perl5/vendor_perl:path/to/home/anaconda3/envs/braker3_env/lib/perl5/core_perl${PERL5LIB:+:$PERL5LIB}
host=$(hostname)
log="$LOGROOT/braker3_rna_protein_worker_${host}_$(date +%Y%m%d_%H%M%S).log"
exec >> "$log" 2>&1
echo "[$(date)] START host=$host threads=$THREADS prot=$PROT manifest=$MAN"
tail -n +2 "$MAN" | while IFS=$'\t' read -r sample soft bam_csv n_bam status; do
  [ -n "$sample" ] || continue
  [ "$status" = ready ] || { echo "[$(date)] SKIP sample=$sample status=$status n_bam=$n_bam"; continue; }
  outdir="$AB3/$sample"
  done="$outdir/braker3_rna_protein.done"
  fail="$outdir/braker3_rna_protein.failed"
  lock="$LOCKROOT/${sample}.lock"
  [ -s "$done" ] && continue
  if ! mkdir "$lock" 2>/dev/null; then continue; fi
  cleanup(){ rm -rf "$lock"; }
  trap cleanup EXIT
  mkdir -p "$outdir"
  rm -f "$fail"
  species="medpan_${sample}_braker3rp"
  rm -rf "$AB3/augustus_config/species/$species"
  echo "[$(date)] BRAKER3_RNA_PROTEIN sample=$sample species=$species n_bam=$n_bam bam=$bam_csv"
  set +e
  braker.pl --gff3 --threads "$THREADS" --prot_seq="$PROT" --genome="$soft" --bam="$bam_csv" --species="$species" --AUGUSTUS_ab_initio --workingdir="$outdir" > "$outdir/braker3_rna_protein.stdout.log" 2> "$outdir/braker3_rna_protein.stderr.log"
  rc=$?
  set -e
  if [ "$rc" -eq 0 ] && { [ -s "$outdir/braker.gff3" ] || [ -s "$outdir/augustus.hints.gff3" ] || [ -s "$outdir/augustus.ab_initio.gff3" ]; }; then
    echo -e "ok\t$(date)\t$host\t$THREADS\t$n_bam\t$bam_csv" > "$done"
    echo "[$(date)] DONE $sample"
  else
    echo -e "failed\trc=$rc\t$(date)\t$host\t$n_bam\t$bam_csv" > "$fail"
    echo "[$(date)] FAILED $sample rc=$rc"
    tail -80 "$outdir/braker3_rna_protein.stderr.log" || true
  fi
  cleanup
  trap - EXIT
done
echo "[$(date)] FINISH host=$host"
