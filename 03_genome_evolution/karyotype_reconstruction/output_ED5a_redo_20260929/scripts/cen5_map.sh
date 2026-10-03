#!/usr/bin/env bash
# Map the T2T M. sativa subsp. caerulea annotation (CENH3 coordinates) onto the AMK
# gene order (pan-genome Msa annotation) by protein identity, Chr5 only.
set -euo pipefail
R=path/to/project/37.karyotype_reconstruction
T=path/to/project/9.T2T_gene_anno_new/output/4.evm_combind/genome_Msa/test1
OUT=$R/output_ED5a_redo_20260929/07_karyotype_evolution/cen5_check
DMD=path/to/home/anaconda3/envs/blast/bin/diamond
mkdir -p $OUT && cd $OUT
PEP=$T/genome_Msa.anno.pep
awk -F'\t' '$1=="Chr5" && $3=="gene"{split($9,a,"[=;]"); print a[2]"\t"$4"\t"$5}' $T/genome_Msa.anno.gff > t2t_chr5_genes.tsv
wc -l t2t_chr5_genes.tsv
# proteins of T2T Chr5 genes (first isoform)
cut -f1 t2t_chr5_genes.tsv | awk '{print $1".1"; print $1}' > ids.txt
path/to/home/anaconda3/envs/wgdi/bin/python - "$PEP" <<'EOF'
import sys
ids = set(l.strip() for l in open("ids.txt"))
keep, out = False, open("t2t_chr5.pep", "w")
for line in open(sys.argv[1]):
    if line.startswith(">"):
        keep = line[1:].split()[0] in ids
    if keep:
        out.write(line)
EOF
grep -c '>' t2t_chr5.pep
awk -F'\t' '$1==5' $R/data/Kary_aak/aak.gff | cut -f2 > aak5.ids
path/to/home/anaconda3/envs/wgdi/bin/python - <<'EOF'
ids = set(l.strip() for l in open("aak5.ids"))
keep, out = False, open("aak5.pep", "w")
for line in open("path/to/project/37.karyotype_reconstruction/data/Kary_aak/aak.pep"):
    if line.startswith(">"):
        keep = line[1:].split()[0] in ids
    if keep:
        out.write(line)
EOF
$DMD makedb --in aak5.pep -d aak5 --quiet
$DMD blastp --threads 16 --db aak5 --query t2t_chr5.pep --out t2t_vs_aak5.tsv --outfmt 6 \
  --max-target-seqs 1 --evalue 1e-10 --id 98 --query-cover 90 --subject-cover 90 --quiet
wc -l t2t_vs_aak5.tsv
