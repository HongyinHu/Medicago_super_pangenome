#!/usr/bin/env bash
set -euo pipefail

# Test workflow for sequence-level A17 vs R108 centromere synteny.
# Run on lz10.

source ~/.bashrc 2>/dev/null || true
conda activate biosofeware

BASE=path/to/project/10.centromere_analysis
OUT=${BASE}/output_synthen/02_cent_synteny/genome_A17_genome_R108
PREV=${BASE}/output_synthen/A17_R108_synteny_lz10_run

mkdir -p "${OUT}"
cd "${OUT}"

ln -sf ../../../output_all/00_genome/genome_A17.fa genome_A17.fa
ln -sf ../../../output_all/00_genome/genome_R108.fa genome_R108.fa

cp "${PREV}/genome_A17.CENH3_vs_Input.q20.10k.log2.bedgraph" .
cp "${PREV}/genome_R108.CENH3_vs_Input.q20.10k.log2.bedgraph" .
cp "${PREV}/genome_A17.pericentromere_6Mb.bed" .
cp "${PREV}/genome_R108.pericentromere_6Mb.bed" .
cp "${PREV}/chr_pairs1-4.tsv" .
cp "${PREV}/chr_pairs5-8.tsv" .

if [ ! -s genome_A17_vs_genome_R108.asm5.paf ]; then
  minimap2 -x asm5 -t 24 -c --eqx genome_R108.fa genome_A17.fa > genome_A17_vs_genome_R108.asm5.paf
fi

ls -lh genome_A17_vs_genome_R108.asm5.paf
wc -l genome_A17_vs_genome_R108.asm5.paf
