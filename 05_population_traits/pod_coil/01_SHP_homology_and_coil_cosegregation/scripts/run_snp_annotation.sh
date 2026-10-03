#!/usr/bin/env bash
set -euo pipefail
ROOT=path/to/project
RUN="$ROOT/N_5.pod_coil/01_SHP_homology_and_coil_cosegregation"
REF="$ROOT/N_4.pod_spiny/00_data/3.two_ref/genome_Msa.fa"
GFF="$ROOT/N_4.pod_spiny/00_data/3.two_ref/genome_Msa.gff"
SAMTOOLS=path/to/home/anaconda3/envs/biosofeware/bin/samtools
mkdir -p "$RUN/reference"
$SAMTOOLS faidx "$REF" Chr3:53194386-53205316 > "$RUN/reference/Chr14120.plusminus2kb.fa"
$SAMTOOLS faidx "$REF" Chr4:91171307-91180034 > "$RUN/reference/Chr24229.plusminus2kb.fa"
python3 "$RUN/scripts/annotate_gene_snps.py" --association "$RUN/tables/Chr14120.population_SNP_coil_association.tsv" \
  --gff "$GFF" --gene Chr14120 --region-fasta "$RUN/reference/Chr14120.plusminus2kb.fa" \
  --output "$RUN/tables/Chr14120.population_SNP_coil_association.annotated.tsv"
python3 "$RUN/scripts/annotate_gene_snps.py" --association "$RUN/tables/Chr24229.population_SNP_coil_association.tsv" \
  --gff "$GFF" --gene Chr24229 --region-fasta "$RUN/reference/Chr24229.plusminus2kb.fa" \
  --output "$RUN/tables/Chr24229.population_SNP_coil_association.annotated.tsv"
touch "$RUN/SNP_annotation.done"
