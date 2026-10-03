#!/usr/bin/env python3
"""Supplementary table: inter-AMK breakpoints of M. polymorpha located on two
T2T references (M. sativa subsp. caerulea and M. truncatula R108) relative to
their CENH3-defined centromeres.

Flanking genes of each junction come from the one-to-one ortholog blocks
(04_wgdi_p0.2/Mpol).  Msa orthologues: AMK gene index -> T2T Msa gene by protein
identity (cen5_check/t2t_vs_aak{3,5,6}.tsv; the AMK uses the pan-genome Msa
annotation, the CENH3 data the T2T one).  R108 orthologues: best DIAMOND hit of
the Mpol gene on the R108 chromosome homologous to the AMK chromosome (the R108
annotation used here and the CENH3 data share the same T2T assembly).
"""
import os

import pandas as pd

ROOT = "path/to/project/37.karyotype_reconstruction"
OUT = os.path.join(ROOT, "output_ED5a_redo_20260929")
EVO = os.path.join(OUT, "07_karyotype_evolution")
CK = os.path.join(EVO, "cen5_check")
T2T = "path/to/project/10.centromere_analysis/output_synthen"
MSA_BED = T2T + ("/06_functional_centromere_sequence_composition/all_species_TableS11_style/"
                 "regions/genome_Msa.functional_centromere.TableS10.ChrOnly.bed")
MPO_BED = T2T + ("/06_functional_centromere_sequence_composition/all_species_TableS11_style/"
                 "regions/genome_Mpo.functional_centromere.TableS10.ChrOnly.bed")
R108_BED = T2T + "/03_cent_TE/genome_R108/CENH3_functional_centromere/genome_R108.CENH3.functional_centromere.final.bed"
LAB = "Mpol"


def wgdi_gff(path):
    g = pd.read_csv(path, sep="\t", header=None)
    g = g.iloc[:, :7] if g.shape[1] >= 7 else g
    g.columns = ["chr", "gene", "s", "e", "str", "order", "orig"][:g.shape[1]]
    g["chr"] = g["chr"].astype(str)
    if "orig" not in g:
        g["orig"] = g["gene"]
    return g


def bed(path):
    b = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2])
    return {str(r[0]).replace("Chr", ""): (r[1], r[2]) for _, r in b.iterrows()}


def rel(pos, cen):
    s, e = cen
    if pos < s:
        return "%.2f Mb left of CEN" % ((s - pos) / 1e6)
    if pos > e:
        return "%.2f Mb right of CEN" % ((pos - e) / 1e6)
    return "within CEN"


# ---- junctions (same rule as dotplot_evidence.py: segments >= 100 genes)
seg = pd.read_csv(os.path.join(EVO, "segments_min30.tsv"), sep="\t")
seg["chr"] = seg["chr"].astype(str)
seg = seg[(seg.species == LAB) & (seg.genes >= 100)]
J = []
for ch, g in seg.groupby("chr"):
    m = []
    for r in g.sort_values("start").to_dict("records"):
        if m and m[-1]["AMK"] == r["AMK"]:
            m[-1]["last"] = r
        else:
            r["first"] = r["last"] = r
            m.append(r)
    for a, b in zip(m, m[1:]):
        la, fb = a["last"], b["first"]
        J.append((ch, la["AMK"], la["amk_end"] if la["ori"] == "+" else la["amk_start"],
                  fb["AMK"], fb["amk_start"] if fb["ori"] == "+" else fb["amk_end"]))

# ---- one-to-one ortholog pairs of Mpol
blk = pd.read_csv(os.path.join(OUT, "04_wgdi_p0.2", LAB, LAB + "_aak.ortholog_blocks.1to1.csv"))
P = []
for _, r in blk.iterrows():
    for q, a in zip(r.block1.split("_"), r.block2.split("_")):
        P.append((str(r.chr1), int(q), int(r.chr2), int(a)))
P = pd.DataFrame(P, columns=["chr", "q", "amk_chr", "amk"])
qg = wgdi_gff(os.path.join(OUT, "04_wgdi_p0.2", LAB, "in.gff1"))
qgi = qg.set_index(["chr", "order"])

# ---- Msa T2T mapping
msa = {}
for c in (3, 5, 6):
    t = pd.read_csv(os.path.join(CK, "t2t_chr%d_genes.tsv" % c), sep="\t", header=None, names=["gene", "s", "e"])
    h = pd.read_csv(os.path.join(CK, "t2t_vs_aak%d.tsv" % c), sep="\t", header=None)
    h = h.sort_values(11, ascending=False).drop_duplicates(0)
    h["gene"] = h[0].str.replace(r"\.\d+$", "", regex=True)
    aak = wgdi_gff(os.path.join(ROOT, "data/Kary_aak/aak.gff")).set_index("gene")
    h["amk"] = aak.loc[h[1], "order"].values
    msa[c] = t.merge(h[["gene", "amk"]], on="gene")


def msa_orth(c, idx):
    """T2T Msa gene for an AMK index; when several T2T genes map to the same AMK gene
    (near-identical copies), keep the one positioned among its AMK neighbours."""
    m = msa[c]
    hit = m[m.amk == idx]
    if len(hit) == 0:
        hit = m.iloc[[(m.amk - idx).abs().argmin()]]
    if len(hit) > 1:
        nb = m[(m.amk - idx).abs().between(1, 15)]
        ref = nb.s.median()
        hit = hit.iloc[[(hit.s - ref).abs().argmin()]]
    r = hit.iloc[0]
    return r.gene, r.s, int(r.amk)


# ---- R108: syntenic orthologue of the AMK gene, from the R108-AMK one-to-one blocks
rg = wgdi_gff(os.path.join(ROOT, "data/genome_R108/old/genome_R108.gff"))
rgi = rg.set_index(["chr", "order"])
rb = pd.read_csv(os.path.join(OUT, "04_wgdi_p0.2", "Mtru_R108", "Mtru_R108_aak.ortholog_blocks.1to1.csv"))
RP = []
for _, r in rb.iterrows():
    for q, a in zip(r.block1.split("_"), r.block2.split("_")):
        RP.append((str(r.chr1), int(q), int(r.chr2), int(a)))
RP = pd.DataFrame(RP, columns=["chr", "q", "amk_chr", "amk"])


def r108_orth(c, idx):
    p = RP[RP.amk_chr == c]
    r = p.iloc[(p.amk - idx).abs().argsort().iloc[0]]
    g = rgi.loc[(r.chr, r.q)]
    return g.orig, g.s, r.chr, int(r.amk)

cen_msa, cen_r108, cen_mpo = bed(MSA_BED), bed(R108_BED), bed(MPO_BED)
# CEN5 of M. sativa subsp. caerulea: core CENH3 domain used in the figure source data
# (TableS10 lists a broader 40.70-48.65 Mb interval)
cen_msa["5"] = (45.70e6, 48.65e6)
rows = []
for k, (ch, a_chr, a_bp, b_chr, b_bp) in enumerate(J, 1):
    for side, amk_chr, amk_idx in (("left", a_chr, a_bp), ("right", b_chr, b_bp)):
        pp = P[(P.chr == ch) & (P.amk_chr == amk_chr)]
        r = pp.iloc[(pp.amk - amk_idx).abs().argsort().iloc[0]]
        qgene = qgi.loc[(ch, r.q)]
        mg, mpos, midx = msa_orth(amk_chr, int(r.amk))
        rid, rpos, rchr, ridx = r108_orth(amk_chr, int(r.amk))
        rtxt = rel(rpos, cen_r108[rchr])
        rows.append({
            "Junction": "J%d" % k,
            "M. polymorpha chromosome": "Chr" + ch,
            "AMK junction": "AMK%d|AMK%d" % (a_chr, b_chr),
            "Side": side,
            "M. polymorpha gene": qgene.orig,
            "M. polymorpha position (Mb)": round(qgene.s / 1e6, 2),
            "AMK chromosome": "AMK%d" % amk_chr,
            "AMK gene index": int(r.amk),
            "M. sativa subsp. caerulea orthologue": mg,
            "M. sativa subsp. caerulea position (Mb)": "Chr%d:%.2f" % (amk_chr, mpos / 1e6),
            "M. sativa subsp. caerulea CENH3 domain (Mb)": "CEN%d %.2f-%.2f" % ((amk_chr,) + tuple(x / 1e6 for x in cen_msa[str(amk_chr)])),
            "Relative to M. sativa CEN": rel(mpos, cen_msa[str(amk_chr)]),
            "M. truncatula R108 orthologue": rid,
            "M. truncatula R108 position (Mb)": "Chr%s:%.2f" % (rchr, rpos / 1e6),
            "M. truncatula R108 CENH3 domain (Mb)": "CEN%s %.2f-%.2f" % ((rchr,) + tuple(x / 1e6 for x in cen_r108[rchr])),
            "Relative to R108 CEN": rtxt,
        })
T = pd.DataFrame(rows)
dst = os.path.join(EVO, "SupTable_Mpol_breakpoints_Msa_R108")
T.to_csv(dst + ".tsv", sep="\t", index=False)
try:
    T.to_excel(dst + ".xlsx", index=False)
except Exception as exc:
    print("xlsx not written:", exc)
pd.set_option("display.width", 300)
pd.set_option("display.max_columns", 30)
print(T.to_string(index=False))
print("\nM. polymorpha CENH3 domains (TableS10):", {k: v for k, v in cen_mpo.items() if k in ("3", "5")})
