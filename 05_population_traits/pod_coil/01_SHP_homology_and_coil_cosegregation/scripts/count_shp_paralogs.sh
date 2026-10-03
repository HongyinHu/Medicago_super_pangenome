#!/usr/bin/env bash
set -euo pipefail

RUN=path/to/project/N_5.pod_coil/01_SHP_homology_and_coil_cosegregation
OUT="$RUN/copy_number"
ENV=path/to/home/anaconda3/envs/biosofeware/bin
PYTHON=python3

mkdir -p "$OUT"

# Use the published Medicago SHP and Arabidopsis SHP1/SHP2 proteins to recover
# full-length candidate homologs. AG/STK-class references are retained in the
# phylogeny so that closely related paralogs are not counted as SHP copies.
$PYTHON "$RUN/scripts/fasta_tools.py" \
  --fasta "$RUN/rbh/Arabidopsis_MADS_refs.pep" \
  --ids AT3G58780.1,AT2G42830.1 \
  --output "$OUT/AtSHP1_SHP2.pep.fa"
cat "$RUN/rbh/JX308825.MtruSHP.pep.fa" "$OUT/AtSHP1_SHP2.pep.fa" \
  > "$OUT/SHP_selection_queries.pep.fa"

: > "$OUT/candidate_summary.tsv"
printf 'genome\tcandidate_id\n' >> "$OUT/candidate_summary.tsv"
: > "$OUT/Msa_R108_candidates.pep.fa"

for item in "Msa:genome_Msa" "R108:genome_R108"; do
  label=${item%%:*}
  species=${item#*:}
  "$ENV/blastp" -query "$OUT/SHP_selection_queries.pep.fa" \
    -db "$RUN/db/$species" -max_target_seqs 100 -evalue 1e-20 \
    -outfmt '6 qseqid sseqid pident length qcovs evalue bitscore' \
    -out "$OUT/${label}.SHP_queries_vs_proteome.tsv"

  awk '$3 >= 35 && $5 >= 70 && $6 <= 1e-20 {print $2}' \
    "$OUT/${label}.SHP_queries_vs_proteome.tsv" | sort -u \
    > "$OUT/${label}.candidate_ids.txt"

  awk -v g="$label" '{print g "\t" $0}' "$OUT/${label}.candidate_ids.txt" \
    >> "$OUT/candidate_summary.tsv"

  $PYTHON "$RUN/scripts/fasta_tools.py" \
    --fasta "$RUN/proteomes/${species}.pep" \
    --ids "$OUT/${label}.candidate_ids.txt" \
    --output "$OUT/${label}.candidates.pep.fa"
  awk -v p="$label" '/^>/{sub(/^>/, ">" p "|")} {print}' \
    "$OUT/${label}.candidates.pep.fa" >> "$OUT/Msa_R108_candidates.pep.fa"
done

cat "$OUT/Msa_R108_candidates.pep.fa" \
    "$RUN/rbh/Arabidopsis_MADS_refs.pep" \
    "$RUN/rbh/JX308825.MtruSHP.pep.fa" \
    > "$OUT/SHP_AG_STK_candidates_and_refs.pep.fa"

"$ENV/mafft" --localpair --maxiterate 1000 \
  "$OUT/SHP_AG_STK_candidates_and_refs.pep.fa" \
  > "$OUT/SHP_AG_STK_candidates_and_refs.pep.aln.fa" \
  2> "$OUT/mafft.log"

"$ENV/iqtree2" -s "$OUT/SHP_AG_STK_candidates_and_refs.pep.aln.fa" \
  -m MFP -B 1000 -alrt 1000 -T 4 \
  --prefix "$OUT/SHP_AG_STK_ML" \
  > "$OUT/iqtree.stdout.log" 2>&1

cat "$RUN/rbh/Arabidopsis_MADS_refs.pep" \
    "$RUN/rbh/JX308825.MtruSHP.pep.fa" > "$OUT/reference_MADS.pep.fa"
"$ENV/makeblastdb" -in "$OUT/reference_MADS.pep.fa" -dbtype prot \
  -out "$OUT/reference_MADS" > "$OUT/makeblastdb.log"
"$ENV/blastp" -query "$OUT/Msa_R108_candidates.pep.fa" \
  -db "$OUT/reference_MADS" -max_target_seqs 8 -evalue 1e-10 \
  -outfmt '6 qseqid sseqid pident length qcovs evalue bitscore' \
  -out "$OUT/candidates_vs_reference_MADS.tsv"
awk 'BEGIN{OFS="\t"; print "candidate","best_reference","pident","qcovs","evalue","bitscore"} \
     !seen[$1]++ {print $1,$2,$3,$5,$6,$7}' \
  "$OUT/candidates_vs_reference_MADS.tsv" \
  > "$OUT/candidates_best_reference.tsv"

touch "$OUT/SHP_copy_number_analysis.done"
