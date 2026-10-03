#!/usr/bin/env bash
set -euo pipefail
THREADS=${THREADS:-16}
BASE=path/to/project/N_1.coding_gene_anno
AB3=$BASE/02_abinitio_prediction/BRAKER3_protein
PROT=$AB3/inputs/related_species.all_refs.gffread_completed.pep.fa
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
log="$LOGROOT/braker3_protein_worker_${host}_$(date +%Y%m%d_%H%M%S).log"
exec >> "$log" 2>&1
echo "[$(date)] START host=$host threads=$THREADS prot=$PROT prothint=$PROTHINT_PATH"
awk 'NR>1{print $1"\t"$2}' "$BASE/01_genome_versions/genome_manifest.tsv" | while IFS=$'\t' read -r sample soft; do
  [ -n "$sample" ] || continue
  outdir="$AB3/$sample"
  done="$outdir/braker3_protein.done"
  fail="$outdir/braker3_protein.failed"
  lock="$LOCKROOT/${sample}.lock"
  [ -s "$done" ] && continue
  if ! mkdir "$lock" 2>/dev/null; then continue; fi
  cleanup(){ rm -rf "$lock"; }
  trap cleanup EXIT
  mkdir -p "$outdir"
  rm -f "$fail"
  species="medpan_${sample}_braker3p"
  rm -rf "$AB3/augustus_config/species/$species"
  echo "[$(date)] BRAKER3_PROTEIN sample=$sample species=$species"
  set +e
  braker.pl --genome="$soft" --prot_seq="$PROT" --species="$species" --gff3 --AUGUSTUS_ab_initio --threads="$THREADS" --workingdir="$outdir" > "$outdir/braker3_protein.stdout.log" 2> "$outdir/braker3_protein.stderr.log"
  rc=$?
  set -e
  if [ "$rc" -eq 0 ] && { [ -s "$outdir/braker.gff3" ] || [ -s "$outdir/augustus.hints.gff3" ] || [ -s "$outdir/augustus.ab_initio.gff3" ]; }; then
    echo -e "ok\t$(date)\t$host\t$THREADS" > "$done"
    echo "[$(date)] DONE $sample"
  else
    echo -e "failed\trc=$rc\t$(date)\t$host" > "$fail"
    echo "[$(date)] FAILED $sample rc=$rc"
    tail -60 "$outdir/braker3_protein.stderr.log" || true
  fi
  cleanup
  trap - EXIT
done
echo "[$(date)] FINISH host=$host"
