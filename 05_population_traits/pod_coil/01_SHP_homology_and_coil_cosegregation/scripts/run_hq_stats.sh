#!/usr/bin/env bash
set -euo pipefail
ROOT=path/to/project
RUN="$ROOT/N_5.pod_coil/01_SHP_homology_and_coil_cosegregation"
DATA18="$ROOT/N_4.pod_spiny/00_data"
ENV=path/to/home/anaconda3/envs/biosofeware/bin
python3 "$RUN/scripts/build_rbh_tables.py" --work "$RUN" --species "$RUN/species.list" --stage final
for gene in Chr14120.1 Chr24229.1; do
  "$ENV/mafft" --auto "$RUN/rbh/strict_HQ_${gene}.pep" > "$RUN/alignments/strict_HQ_${gene}.pep.aln.fa" 2> "$RUN/logs/HQ_${gene}.pep.mafft.log"
  "$ENV/mafft" --auto "$RUN/rbh/strict_HQ_${gene}.cds" > "$RUN/alignments/strict_HQ_${gene}.cds.aln.fa" 2> "$RUN/logs/HQ_${gene}.cds.mafft.log"
  python3 "$RUN/scripts/association.py" cross --alignment "$RUN/alignments/strict_HQ_${gene}.pep.aln.fa" \
    --phenotype "$DATA18/5.phenotye_info/phenotype.tsv" --output "$RUN/tables/${gene}.HQ_protein_site_coil_association.tsv"
  python3 "$RUN/scripts/association.py" cross --alignment "$RUN/alignments/strict_HQ_${gene}.cds.aln.fa" \
    --phenotype "$DATA18/5.phenotye_info/phenotype.tsv" --output "$RUN/tables/${gene}.HQ_cds_site_coil_association.tsv"
done
touch "$RUN/HQ_stats.done"
