#!/usr/bin/env python3
"""Do the AMK5 breakpoints of M. praecox Chr5 (and M. polymorpha) fall at CEN5?

1. T2T M. sativa subsp. caerulea Chr5 genes -> AMK5 gene index (protein identity,
   cen5_check/t2t_vs_aak5.tsv); CEN5 anchors Msa29724 / Msa29749.
2. AMK5 breakpoints of Mpra Chr5 and Mpol from the one-to-one ortholog blocks
   (04_wgdi_p0.2), with the flanking species genes.
3. The same Mpra flanking genes located on R108 Chr5 (Mpra_vs_R108 best hits)
   relative to the R108 CENH3-defined CEN5.
"""
import os

import numpy as np
import pandas as pd

ROOT = "path/to/project/37.karyotype_reconstruction"
OUT = os.path.join(ROOT, "output_ED5a_redo_20260929")
CK = os.path.join(OUT, "07_karyotype_evolution", "cen5_check")
DOT = os.path.join(OUT, "07_karyotype_evolution", "dotplot")
T2T = "path/to/project"
MSA_CEN5_USER = (45.70e6, 48.65e6)      # core CEN5 given by the user (T2T Msa)
R108_CEN_BED = (T2T + "/10.centromere_analysis/output_synthen/03_cent_TE/genome_R108/"
                "CENH3_functional_centromere/genome_R108.CENH3.functional_centromere.final.bed")
MSA_CEN_BED = (T2T + "/10.centromere_analysis/output_synthen/06_functional_centromere_sequence_"
               "composition/all_species_TableS11_style/regions/genome_Msa.functional_centromere."
               "TableS10.ChrOnly.bed")

# ---------- 1. T2T Msa Chr5 -> AMK5
t2t = pd.read_csv(os.path.join(CK, "t2t_chr5_genes.tsv"), sep="\t", header=None,
                  names=["gene", "start", "end"])
hit = pd.read_csv(os.path.join(CK, "t2t_vs_aak5.tsv"), sep="\t", header=None)
hit = hit.sort_values(11, ascending=False).drop_duplicates(0)
hit["gene"] = hit[0].str.replace(r"\.\d+$", "", regex=True)
aak = pd.read_csv(os.path.join(ROOT, "data/Kary_aak/aak.gff"), sep="\t", header=None,
                  usecols=range(7), names=["chr", "id", "s", "e", "str", "order", "msa"])
aak5 = aak[aak.chr.astype(str) == "5"].set_index("id")
hit["amk_order"] = aak5.loc[hit[1], "order"].values
m = t2t.merge(hit[["gene", "amk_order", 2]], on="gene", how="left").rename(columns={2: "pid"})
m.to_csv(os.path.join(CK, "t2t_msa_chr5_to_AMK5.tsv"), sep="\t", index=False)
print("T2T Msa Chr5 genes: {}, mapped to AMK5: {}".format(len(m), m.amk_order.notnull().sum()))
for g in ("Msa29724", "Msa29749"):
    r = m[m.gene == g].iloc[0]
    print("  {} at {:.2f} Mb -> AMK5 gene {}".format(g, r.start / 1e6, r.amk_order))
mm = m.dropna(subset=["amk_order"])
inside = mm[(mm.start >= MSA_CEN5_USER[0]) & (mm.end <= MSA_CEN5_USER[1])]
print("  mapped genes inside 45.70-48.65 Mb: {} (AMK5 {}-{})".format(
    len(inside), inside.amk_order.min() if len(inside) else "-", inside.amk_order.max() if len(inside) else "-"))
bed = pd.read_csv(MSA_CEN_BED, sep="\t", header=None)
print("  Msa CENH3 bed Chr5:", bed[bed[0] == "Chr5"].values.tolist())


def t2t_pos(amk_idx):
    """T2T Msa coordinate of the AMK5 gene (nearest mapped)."""
    k = (mm.amk_order - amk_idx).abs().idxmin()
    return mm.loc[k, "amk_order"], mm.loc[k, "gene"], mm.loc[k, "start"] / 1e6


# ---------- 2. AMK5 breakpoints in Mpra / Mpol
def pairs(lab, chr1, chr2="5"):
    b = pd.read_csv(os.path.join(OUT, "04_wgdi_p0.2", lab, lab + "_aak.ortholog_blocks.1to1.csv"))
    b = b[(b.chr1.astype(str) == chr1) & (b.chr2.astype(str) == chr2)]
    rows = []
    for _, r in b.iterrows():
        for q, a in zip(r.block1.split("_"), r.block2.split("_")):
            rows.append((int(q), int(a), r.id))
    return pd.DataFrame(rows, columns=["q", "amk", "block"]).sort_values("q")


def qgenes(lab):
    g = pd.read_csv(os.path.join(OUT, "04_wgdi_p0.2", lab, "in.gff1"), sep="\t", header=None,
                    usecols=range(6), names=["chr", "gene", "s", "e", "str", "order"])
    g["chr"] = g["chr"].astype(str)
    return g


def report(lab, chr1, bps):
    p = pairs(lab, chr1)
    g = qgenes(lab)
    g = g[g.chr == chr1].set_index("order")
    print("\n## {} Chr{}: AMK5 breakpoints".format(lab, chr1))
    out = []
    for a in bps:
        r = p.iloc[(p.amk - a).abs().argsort().iloc[0]]
        gene = g.loc[r.q, "gene"]
        o, tg, mbp = t2t_pos(r.amk)
        print("  AMK5 gene {:>5} <- {} Chr{} gene {:>5} ({}; {:.2f} Mb) | T2T Msa: {} {:.2f} Mb".format(
            int(r.amk), lab, chr1, int(r.q), gene, g.loc[r.q, "s"] / 1e6, tg, mbp))
        out.append((a, gene, int(r.q)))
    return out


mpra = report("Mpra", "5", [2030, 2033, 4303, 4304])
# where do the CEN5 anchors (AMK5 2711 / 2731) sit on Mpra Chr5?
pp = pairs("Mpra", "5")
gg = qgenes("Mpra")
gg = gg[gg.chr == "5"].set_index("order")
print("  AMK5 2690-2850 on Mpra Chr5 (CEN5 anchors 2711/2731):")
seg = pp[pp.amk.between(2690, 2850)]
for _, r in seg.iterrows():
    print("    AMK5 {:>5} <- Mpra Chr5 gene {:>5} {:.2f} Mb (block {})".format(
        r.amk, r.q, gg.loc[r.q, "s"] / 1e6, r.block))
mpol3 = report("Mpol", "3", [2731])
mpol5 = report("Mpol", "5", [2711])

# ---------- 3. Mpra flanking genes on R108 Chr5
rg = pd.read_csv(os.path.join(ROOT, "data/genome_R108/old/genome_R108.gff"), sep="\t", header=None,
                 usecols=range(6), names=["chr", "gene", "s", "e", "str", "order"]).set_index("gene")
bl = pd.read_csv(os.path.join(DOT, "Mpra_vs_R108.blastp.txt"), sep="\t", header=None)
bl = bl.sort_values([0, 11], ascending=[True, False]).drop_duplicates(0).set_index(0)
cen = pd.read_csv(R108_CEN_BED, sep="\t", header=None)
c5 = cen[cen[0] == "Chr5"]
print("\n## R108 CENH3 functional centromere Chr5:", c5.values.tolist())
# R108 genes flanking its CEN5 -> AMK5 index (R108 vs AMK hits of the ED5a run)
r5 = rg[rg.chr.astype(str) == "5"].sort_values("s")
cs, ce = c5[1].min(), c5[2].max()
left = r5[r5.e < cs].tail(3)
right = r5[r5.s > ce].head(3)
ra = pd.read_csv(os.path.join(ROOT, "output_Kary/genome_R108/all.blastp.txt"), sep="\t", header=None)
ra = ra.sort_values([0, 11], ascending=[True, False]).drop_duplicates(0).set_index(0)
for lab_, sub in (("left", left), ("right", right)):
    for gid, r in sub.iterrows():
        a = ra.loc[gid, 1] if gid in ra.index else None
        ai = aak.set_index("id").loc[a] if a in set(aak.id) else None
        print("  R108 CEN5 {} flank {} {:.2f} Mb -> AMK{} gene {}".format(
            lab_, gid, r.s / 1e6, ai.chr if ai is not None else "-", ai.order if ai is not None else "-"))
pg = qgenes("Mpra")
pg = pg[pg.chr == "5"].set_index("order")
for a, gene, q in mpra:
    # the breakpoint gene and its 4 neighbours on the same side (same AMK5 block)
    side = range(q - 4, q + 1) if a in (2030, 4304) else range(q, q + 5)
    for qq in side:
        gg = pg.loc[qq, "gene"] if qq in pg.index else None
        if gg in bl.index:
            s = bl.loc[gg, 1]
            if s in rg.index:
                print("  Mpra Chr5 {} (AMK5~{}) -> R108 {} Chr{} {:.2f} Mb".format(
                    gg, a, s, rg.loc[s, "chr"], rg.loc[s, "s"] / 1e6))
