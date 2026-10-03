#!/usr/bin/env bash
# Map T2T M. sativa subsp. caerulea genes (CENH3 coordinate system) onto the AMK gene
# order (pan-genome Msa annotation) by protein identity, for chromosomes 3, 5 and 6.
set -euo pipefail
R=path/to/project/37.karyotype_reconstruction
A=path/to/project/9.T2T_gene_anno_new/output/4.evm_combind/genome_Msa/test1
OUT=$R/output_ED5a_redo_20260929/07_karyotype_evolution/cen5_check
DMD=path/to/home/anaconda3/envs/blast/bin/diamond
PY=path/to/home/anaconda3/envs/wgdi/bin/python
cd $OUT
for c in 3 5 6; do
  awk -F'\t' -v C="Chr$c" '$1==C && $3=="gene"{split($9,a,"[=;]"); print a[2]"\t"$4"\t"$5}' $A/genome_Msa.anno.gff > t2t_chr${c}_genes.tsv
  $PY - $c <<'EOF'
import sys
c = sys.argv[1]
R = "path/to/project/37.karyotype_reconstruction"
ids = set()
for l in open("t2t_chr%s_genes.tsv" % c):
    g = l.split("\t")[0]; ids.update([g, g + ".1"])
def sub(src, keep_ids, dst):
    keep, out = False, open(dst, "w")
    for line in open(src):
        if line.startswith(">"):
            keep = line[1:].split()[0] in keep_ids
        if keep: out.write(line)
sub("path/to/project/9.T2T_gene_anno_new/output/4.evm_combind/genome_Msa/test1/genome_Msa.anno.pep",
    ids, "t2t_chr%s.pep" % c)
amk = set(l.split("\t")[1] for l in open(R + "path/to/data") if l.split("\t")[0] == c)
sub(R + "path/to/data", amk, "aak%s.pep" % c)
EOF
  $DMD makedb --in aak${c}.pep -d aak${c} --quiet
  $DMD blastp --threads 16 --db aak${c} --query t2t_chr${c}.pep --out t2t_vs_aak${c}.tsv --outfmt 6 \
    --max-target-seqs 1 --evalue 1e-10 --id 98 --query-cover 90 --subject-cover 90 --quiet
  echo "Chr$c: $(wc -l < t2t_chr${c}_genes.tsv) genes, $(wc -l < t2t_vs_aak${c}.tsv) mapped"
done
