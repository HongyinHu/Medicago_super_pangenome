#!/usr/bin/env bash
# M. praecox proteins vs M. truncatula R108 proteins (same R108 annotation as ED5a,
# data/genome_R108/old) for the dot-plot evidence of the x = 7 karyotype.
set -euo pipefail
R=path/to/project/37.karyotype_reconstruction
OUT=$R/output_ED5a_redo_20260929/07_karyotype_evolution/dotplot
DMD=path/to/home/anaconda3/envs/blast/bin/diamond
mkdir -p $OUT && cd $OUT
$DMD makedb --in $R/data/genome_R108/old/genome_R108.pep -d R108_old.pep --quiet
$DMD blastp --threads 32 --db R108_old.pep --query $R/data/genome_410/genome_410.pep \
  --out Mpra_vs_R108.blastp.txt --outfmt 6 --sensitive --max-target-seqs 10 --evalue 1e-5 --quiet
wc -l Mpra_vs_R108.blastp.txt
