#!/usr/bin/env bash
# Proteins of the x = 8 genomes with inter-chromosomal changes vs M. truncatula R108
# (same R108 annotation as ED5a), for the dot-plot evidence.  Protein files are
# the ones behind the karyotype runs (in.gff1 of 04_wgdi_p0.2, .gff -> .pep).
set -euo pipefail
R=path/to/project/37.karyotype_reconstruction
OUT=$R/output_ED5a_redo_20260929
DOT=$OUT/07_karyotype_evolution/dotplot
DMD=path/to/home/anaconda3/envs/blast/bin/diamond
cd $DOT
[ -e R108_old.pep.dmnd ] || $DMD makedb --in $R/data/genome_R108/old/genome_R108.pep -d R108_old.pep --quiet
run() {
  lab=$1
  gff=$(readlink -f $OUT/04_wgdi_p0.2/$lab/in.gff1)
  pep=${gff%.gff}.pep
  n1=$(grep -c '>' $pep); n2=$(wc -l < $gff)
  $DMD blastp --threads 16 --db R108_old.pep --query $pep --out ${lab}_vs_R108.blastp.txt \
    --outfmt 6 --sensitive --max-target-seqs 10 --evalue 1e-5 --quiet
  echo "$lab pep=$n1 gff=$n2 hits=$(wc -l < ${lab}_vs_R108.blastp.txt)"
}
export -f run; export OUT DMD
for lab in Malbus Mlan Mrut Marc Mrad Medg Mlup Msuf Msec Msat_zm4; do echo $lab; done | \
  xargs -P 10 -I{} bash -c 'run {}'
