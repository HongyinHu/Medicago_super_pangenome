#!/usr/bin/env bash
set -euo pipefail
THREADS=${THREADS:-16}
BASE=path/to/project/N_1.coding_gene_anno
ABP=$BASE/02_abinitio_prediction/BRAKER3_protein
RNAMAN=$BASE/02_abinitio_prediction/BRAKER3_RNA_protein/inputs/braker3_rna_protein_bam_manifest.tsv
PROT=${PROT:-$BASE/02_abinitio_prediction/BRAKER3_RNA_protein/inputs/related_species.all_refs.gffread_completed.sanitized.no_stop.pep.fa}
[ -s "$PROT" ] || PROT=$ABP/inputs/related_species.all_refs.gffread_completed.pep.fa
LOCKROOT=$ABP/locks
LOGROOT=$ABP/logs
mkdir -p "$LOCKROOT" "$LOGROOT"
source path/to/home/anaconda3/etc/profile.d/conda.sh
conda activate braker3_env
export AUGUSTUS_CONFIG_PATH=$ABP/augustus_config
export GENEMARK_PATH=path/to/home/sofeware/gmes_linux_64_4
export PROTHINT_PATH=path/to/home/sofeware/gmes_linux_64_4/ProtHint/bin
export PATH=$PROTHINT_PATH:$GENEMARK_PATH:$PATH
export PERL5LIB=path/to/home/anaconda3/envs/braker3_env/lib/perl5/site_perl:path/to/home/anaconda3/envs/braker3_env/lib/perl5/vendor_perl:path/to/home/anaconda3/envs/braker3_env/lib/perl5/core_perl${PERL5LIB:+:$PERL5LIB}
host=$(hostname)
log="$LOGROOT/braker3_protein_only_worker_${host}_$(date +%Y%m%d_%H%M%S).log"
exec >> "$log" 2>&1
echo "[$(date)] START protein-only host=$host threads=$THREADS prot=$PROT"
if [ ! -s "$PROT" ]; then echo "[$(date)] ERROR missing protein FASTA $PROT"; exit 2; fi
awk 'BEGIN{FS="\t"} NR>1 && $5=="ready"{ready[$1]=1} END{for (s in ready) print s}' "$RNAMAN" > "$ABP/inputs/rna_ready_samples.current.txt"
awk 'NR>1{print $1"\t"$2}' "$BASE/01_genome_versions/genome_manifest.tsv" | while IFS=$'\t' read -r sample soft; do
  [ -n "${sample:-}" ] || continue
  if grep -qx "$sample" "$ABP/inputs/rna_ready_samples.current.txt"; then echo "[$(date)] SKIP RNA_ready $sample"; continue; fi
  outdir="$ABP/$sample"
  done_file="$outdir/braker3_protein.done"
  fail_file="$outdir/braker3_protein.failed"
  lock="$LOCKROOT/${sample}.lock"
  [ -s "$done_file" ] && { echo "[$(date)] SKIP done $sample"; continue; }
  [ -s "$fail_file" ] && { echo "[$(date)] SKIP failed $sample"; continue; }
  [ -s "$soft" ] || { mkdir -p "$outdir"; echo -e "failed\tmissing_softmask\t$(date)\t$host" > "$fail_file"; continue; }
  if ! mkdir "$lock" 2>/dev/null; then echo "[$(date)] SKIP locked $sample"; continue; fi
  echo "$host $$ $(date '+%F %T %Z') braker3_protein_only threads=$THREADS" > "$lock/owner"
  cleanup(){ rm -rf "$lock"; }
  trap cleanup EXIT
  mkdir -p "$outdir"
  species="medpan_${sample}_braker3p"
  rm -rf "$ABP/augustus_config/species/$species"
  echo "[$(date)] RUN protein-only sample=$sample species=$species"
  set +e
  braker.pl --genome="$soft" --prot_seq="$PROT" --species="$species" --gff3 --softmasking --AUGUSTUS_ab_initio --threads="$THREADS" --workingdir="$outdir" > "$outdir/braker3_protein.stdout.log" 2> "$outdir/braker3_protein.stderr.log"
  rc=$?
  set -e
  if [ "$rc" -eq 0 ] && { [ -s "$outdir/braker.gff3" ] || [ -s "$outdir/augustus.hints.gff3" ] || [ -s "$outdir/augustus.ab_initio.gff3" ]; }; then
    echo -e "ok\t$(date)\t$host\t$THREADS\tprotein_only" > "$done_file"
    echo "[$(date)] DONE $sample"
  else
    echo -e "failed\trc=$rc\t$(date)\t$host\tprotein_only" > "$fail_file"
    echo "[$(date)] FAILED $sample rc=$rc"
    tail -120 "$outdir/braker3_protein.stderr.log" || true
  fi
  cleanup
  trap - EXIT
done
echo "[$(date)] FINISH protein-only host=$host"
