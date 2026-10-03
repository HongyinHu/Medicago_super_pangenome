#!/usr/bin/env bash
# M. polymorpha (T2T, genome_nan annotation used in the karyotype runs) vs
# M. truncatula R108 proteins (same R108 annotation as ED5a) for the dot plots.
set -euo pipefail
R=path/to/project/37.karyotype_reconstruction
OUT=$R/output_ED5a_redo_20260929/07_karyotype_evolution/dotplot
DMD=path/to/home/anaconda3/envs/blast/bin/diamond
cd $OUT
[ -e R108_old.pep.dmnd ] || $DMD makedb --in $R/data/genome_R108/old/genome_R108.pep -d R108_old.pep --quiet
cmp <(cut -f2 $R/data/genome_nan/genome_nan.gff | sort) \
    <(cut -f2 $R/output_ED5a_redo_20260929/04_wgdi_p0.2/Mpol/in.gff1 | sort) && echo "Mpol gene IDs match karyotype run"
$DMD blastp --threads 32 --db R108_old.pep --query $R/data/genome_nan/genome_nan.pep \
  --out Mpol_vs_R108.blastp.txt --outfmt 6 --sensitive --max-target-seqs 10 --evalue 1e-5 --quiet
wc -l Mpol_vs_R108.blastp.txt
