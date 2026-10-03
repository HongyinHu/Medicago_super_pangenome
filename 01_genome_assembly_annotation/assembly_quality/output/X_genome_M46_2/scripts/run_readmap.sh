#!/usr/bin/env bash
set -euo pipefail
BASE=path/to/project/6.genome_quality_assess/output/X_genome_M46_2
BWA=path/to/software/bwa/bwa-0.7.17/bwa
SAMTOOLS=path/to/home/anaconda3/envs/biosofeware/bin/samtools
FA="$BASE/assembly_evaluate/ReadMap/M46_genome.fa"
R1="$BASE/data/survey_R1.fq.gz"
R2="$BASE/data/survey_R2.fq.gz"
cp -n "$BASE/data/M46_genome_Chr_reorder.genome.ctg.fa" "$BASE/assembly_evaluate/ReadMap/M46_genome.fa"
cd "$BASE/assembly_evaluate/ReadMap"
"$BWA" index -a bwtsw "$FA"
"$BWA" mem -t 20 "$FA" "$R1" "$R2" 2> bwa_mem.log | "$SAMTOOLS" view -b - | "$SAMTOOLS" flagstat - > flagstat.txt
