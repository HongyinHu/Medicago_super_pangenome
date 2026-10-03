#!/usr/bin/env bash
set -euo pipefail
ROOT=path/to/project/37.karyotype_reconstruction
OUT=$ROOT/output_Kary/genome_M46_2
cd "$OUT"
ln -sfn ../..path/to/data aak.gff
ln -sfn ../..path/to/data aak.lens
ln -sfn ../..path/to/data aak.cds
ln -sfn ../..path/to/data aak.pep
ln -sfn ../genome_M46/aak.ancestor.txt aak.ancestor.txt
path/to/home/anaconda3/envs/biopython/bin/python "$ROOT/my_run/s1.get_gff_lens.py" \
  "$ROOT/data/genome_M46_2/genome_M46_2.genome.fa" \
  "$ROOT/data/genome_M46_2/finally.genome.anno.gff" \
  "$ROOT/data/genome_M46_2/finally.genome.anno.cds" \
  "$ROOT/data/genome_M46_2/finally.genome.anno.pep" \
  genome_M46_2
cat genome_M46_2.cds aak.cds > all.cds
cat genome_M46_2.pep aak.pep > all.pep
test "$(wc -l < genome_M46_2.lens)" -eq 7
test -s genome_M46_2.gff
test -s all.cds
test -s all.pep
printf 'prepared_chromosomes=%s\n' "$(wc -l < genome_M46_2.lens)"
printf 'M46_2_gene_rows=%s\n' "$(wc -l < genome_M46_2.gff)"
printf 'M46_2_cds_records=%s\n' "$(grep -c '^>' genome_M46_2.cds)"
printf 'M46_2_pep_records=%s\n' "$(grep -c '^>' genome_M46_2.pep)"