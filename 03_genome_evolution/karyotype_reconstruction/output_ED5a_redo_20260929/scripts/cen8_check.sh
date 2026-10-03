#!/usr/bin/env bash
# Is the AMK8 breakpoint of M. praecox (AMK8 3659|3661) at the ancestral CEN8?
# Same approach as CEN5: T2T Msa Chr8 genes -> AMK8 by protein identity; R108 CEN8 flanks.
set -euo pipefail
R=path/to/project/37.karyotype_reconstruction
T=path/to/project
A=$T/9.T2T_gene_anno_new/output/4.evm_combind/genome_Msa/test1
OUT=$R/output_ED5a_redo_20260929/07_karyotype_evolution/cen5_check
DMD=path/to/home/anaconda3/envs/blast/bin/diamond
PY=path/to/home/anaconda3/envs/wgdi/bin/python
cd $OUT
awk -F'\t' '$1=="Chr8" && $3=="gene"{split($9,a,"[=;]"); print a[2]"\t"$4"\t"$5}' $A/genome_Msa.anno.gff > t2t_chr8_genes.tsv
$PY - <<'EOF'
ids = set()
for l in open("t2t_chr8_genes.tsv"):
    g = l.split("\t")[0]; ids.update([g, g + ".1"])
keep, out = False, open("t2t_chr8.pep", "w")
for line in open("path/to/project/9.T2T_gene_anno_new/output/4.evm_combind/genome_Msa/test1/genome_Msa.anno.pep"):
    if line.startswith(">"):
        keep = line[1:].split()[0] in ids
    if keep: out.write(line)
aak8 = set(l.split("\t")[1] for l in open("path/to/project/37.karyotype_reconstruction/data/Kary_aak/aak.gff") if l.split("\t")[0] == "8")
keep, out = False, open("aak8.pep", "w")
for line in open("path/to/project/37.karyotype_reconstruction/data/Kary_aak/aak.pep"):
    if line.startswith(">"):
        keep = line[1:].split()[0] in aak8
    if keep: out.write(line)
EOF
$DMD makedb --in aak8.pep -d aak8 --quiet
$DMD blastp --threads 16 --db aak8 --query t2t_chr8.pep --out t2t_vs_aak8.tsv --outfmt 6 \
  --max-target-seqs 1 --evalue 1e-10 --id 98 --query-cover 90 --subject-cover 90 --quiet
$PY - <<'EOF'
import pandas as pd
R = "path/to/project/37.karyotype_reconstruction"
T = "path/to/project"
t = pd.read_csv("t2t_chr8_genes.tsv", sep="\t", header=None, names=["gene", "s", "e"])
h = pd.read_csv("t2t_vs_aak8.tsv", sep="\t", header=None).sort_values(11, ascending=False).drop_duplicates(0)
h["gene"] = h[0].str.replace(r"\.\d+$", "", regex=True)
aak = pd.read_csv(R + "path/to/data", sep="\t", header=None, usecols=range(6),
                  names=["chr", "id", "s", "e", "str", "order"]).set_index("id")
h["amk"] = aak.loc[h[1], "order"].values
m = t.merge(h[["gene", "amk"]], on="gene").sort_values("s")
bed = pd.read_csv(T + "/10.centromere_analysis/output_synthen/06_functional_centromere_sequence_composition/"
                  "all_species_TableS11_style/regions/genome_Msa.functional_centromere.TableS10.ChrOnly.bed",
                  sep="\t", header=None)
cs, ce = bed[bed[0] == "Chr8"][[1, 2]].values[0]
L = m[m.e < cs].tail(3); Rr = m[m.s > ce].head(3)
print("Msa CEN8 (T2T): %.2f-%.2f Mb" % (cs / 1e6, ce / 1e6))
for _, r in L.iterrows(): print("  left  flank %s %.2f Mb -> AMK8 %d" % (r.gene, r.s / 1e6, r.amk))
for _, r in Rr.iterrows(): print("  right flank %s %.2f Mb -> AMK8 %d" % (r.gene, r.s / 1e6, r.amk))
for a in (3659, 3661):
    k = (m.amk - a).abs().idxmin()
    print("  AMK8 %d (Mpra breakpoint) ~ T2T Msa %s %.2f Mb" % (a, m.loc[k, "gene"], m.loc[k, "s"] / 1e6))
# R108 CEN8 flanks -> AMK8
cen = pd.read_csv(T + "/10.centromere_analysis/output_synthen/03_cent_TE/genome_R108/CENH3_functional_centromere/"
                  "genome_R108.CENH3.functional_centromere.final.bed", sep="\t", header=None)
c8 = cen[cen[0] == "Chr8"]; rs, re_ = c8[1].min(), c8[2].max()
rg = pd.read_csv(R + "path/to/data", sep="\t", header=None, usecols=range(6),
                 names=["chr", "gene", "s", "e", "str", "order"])
r8 = rg[rg.chr.astype(str) == "8"].sort_values("s")
ra = pd.read_csv(R + "/output_Kary/genome_R108/all.blastp.txt", sep="\t", header=None)
ra = ra.sort_values([0, 11], ascending=[True, False]).drop_duplicates(0).set_index(0)
print("R108 CEN8: %.2f-%.2f Mb" % (rs / 1e6, re_ / 1e6))
for lab, sub in (("left", r8[r8.e < rs].tail(4)), ("right", r8[r8.s > re_].head(4))):
    for _, r in sub.iterrows():
        a = ra.loc[r.gene, 1] if r.gene in ra.index else None
        if a in aak.index:
            print("  %s flank %s %.2f Mb -> AMK%s %d" % (lab, r.gene, r.s / 1e6, aak.loc[a, "chr"], aak.loc[a, "order"]))
# Mpra genes at the AMK8 breakpoint -> R108
bl = pd.read_csv(R + "/output_ED5a_redo_20260929/07_karyotype_evolution/dotplot/Mpra_vs_R108.blastp.txt",
                 sep="\t", header=None).sort_values([0, 11], ascending=[True, False]).drop_duplicates(0).set_index(0)
b = pd.read_csv(R + "/output_ED5a_redo_20260929/04_wgdi_p0.2/Mpra/Mpra_aak.ortholog_blocks.1to1.csv")
b = b[b.chr2.astype(str) == "8"]
pg = pd.read_csv(R + "/output_ED5a_redo_20260929/04_wgdi_p0.2/Mpra/in.gff1", sep="\t", header=None,
                 usecols=range(6), names=["chr", "gene", "s", "e", "str", "order"])
pg["chr"] = pg.chr.astype(str); pg = pg.set_index(["chr", "order"])
rgi = rg.set_index("gene")
for _, r in b.iterrows():
    for q, a in zip(map(int, r.block1.split("_")), map(int, r.block2.split("_"))):
        if 3640 <= a <= 3680:
            g = pg.loc[(str(r.chr1), q), "gene"]
            s = bl.loc[g, 1] if g in bl.index else None
            pos = "R108 %s Chr%s %.2f Mb" % (s, rgi.loc[s, "chr"], rgi.loc[s, "s"] / 1e6) if s in rgi.index else "-"
            print("  AMK8 %d <- Mpra Chr%s %s %.2f Mb -> %s" % (a, r.chr1, g, pg.loc[(str(r.chr1), q), "s"] / 1e6, pos))
EOF
