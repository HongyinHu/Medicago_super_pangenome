#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/10.centromere_analysis/output_key
CONDA=path/to/home/anaconda3/bin/conda

cd "$ROOT"
echo "GENOMES"
ls -lh data/00_genome/genome_*.fa data/00_genome/genome_*.fa.fai 2>/dev/null || true

echo "REPEAT_INPUTS"
find -L data/03_repeat -maxdepth 3 -type f | grep -Ei '(TElib|TEanno|TRASH_arrays|monomers|consensus|\.lib|\.fa$|\.gff3$)' | head -300 || true

echo "CHIPSEQ_READS"
find -L data/02_chipseq -maxdepth 3 -type f | grep -Ei '(clean\.R[12]\.fq\.gz|genome_.*_(CENH3|Input)_R[12]\.fq\.gz)' | sort | head -200 || true

echo "TOOLS_GENOME_REPEAT"
"$CONDA" run -n genome_repeat bash -lc 'for x in RepeatMasker trf EDTA.pl RepeatModeler; do printf "%s\t" "$x"; command -v "$x" || true; done'

echo "TOOLS_BLAST"
"$CONDA" run -n blast bash -lc 'for x in blastn makeblastdb; do printf "%s\t" "$x"; command -v "$x" || true; done'

echo "TOOLS_MUMMER"
"$CONDA" run -n mummer_env bash -lc 'for x in nucmer show-coords delta-filter; do printf "%s\t" "$x"; command -v "$x" || true; done'

echo "TOOLS_MINIMAP2"
"$CONDA" run -n minimap2_env bash -lc 'for x in minimap2; do printf "%s\t" "$x"; command -v "$x" || true; done'

echo "TOOLS_BIOSOFEWARE"
"$CONDA" run -n biosofeware bash -lc 'for x in samtools bedtools python; do printf "%s\t" "$x"; command -v "$x" || true; done'
