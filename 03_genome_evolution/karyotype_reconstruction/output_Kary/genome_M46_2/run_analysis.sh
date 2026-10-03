#!/usr/bin/env bash
set -euo pipefail
ROOT=path/to/project/37.karyotype_reconstruction
OUT=$ROOT/output_Kary/genome_M46_2
cd "$OUT"
sed 's/genome_M46/genome_M46_2/g' ../genome_M46/total.conf > total.conf
path/to/home/anaconda3/envs/blast/bin/diamond makedb --in aak.pep -d aak.pep
path/to/home/anaconda3/envs/blast/bin/diamond blastp --threads 10 --db aak.pep --query genome_M46_2.pep --out all.blastp.txt --outfmt 6 --sensitive --max-target-seqs 10 --evalue 1e-5
WGDI=path/to/home/anaconda3/envs/wgdi/bin/wgdi
"$WGDI" -d total.conf
"$WGDI" -icl total.conf
"$WGDI" -ks total.conf
"$WGDI" -bi total.conf
"$WGDI" -c total.conf
"$WGDI" -bk total.conf
"$WGDI" -km total.conf
"$WGDI" -k total.conf