#!/usr/bin/env bash
set -euo pipefail
BASE=path/to/project/6.genome_quality_assess/output/X_genome_M46_2
export PATH=path/to/home/anaconda3/envs/biosofeware/bin:$PATH
MERYL=path/to/home/anaconda3/envs/biosofeware/bin/meryl
export MERQURY=path/to/home/sofeware/merqury
R1="$BASE/data/survey_R1.fq.gz"
R2="$BASE/data/survey_R2.fq.gz"
FA="$BASE/data/M46_genome_Chr_reorder.genome.ctg.fa"
cd "$BASE/assembly_evaluate/Merqury"
"$MERYL" k=19 threads=8 memory=16 count output read_R1.meryl "$R1"
"$MERYL" k=19 threads=8 memory=16 count output read_R2.meryl "$R2"
"$MERYL" union-sum output read-db.meryl read_R1.meryl read_R2.meryl
"$MERQURY/merqury.sh" read-db.meryl "$FA" out
