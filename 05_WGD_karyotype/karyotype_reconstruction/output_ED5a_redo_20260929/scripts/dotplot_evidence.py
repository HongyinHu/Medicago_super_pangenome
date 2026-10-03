#!/usr/bin/env python3
"""Dot-plot evidence for a karyotype change (layout after the AZAK/TLI/Aal figure).

Rows    : rearranged chromosomes of the focal species, with an AMK colour bar.
Columns : the AMK chromosomes involved | the same chromosomes of an x = 8
          reference species.
Dots    : protein hits (diamond, e <= 1e-5, bitscore >= 100) in gene order;
          best hit red, second hit blue, 3rd-5th grey (WGDI dotplot convention).
Green dashed lines: junctions between AMK segments on the focal chromosomes
(07_karyotype_evolution/inter_AMK_junctions.tsv); cyan dashed lines: the
corresponding breakpoints on the AMK chromosomes.
Usage: dotplot_evidence.py Mpra
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("RUN", "04_wgdi_p0.2")
import plot_ed5a as P  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = "path/to/project/37.karyotype_reconstruction"
OUT = os.path.join(ROOT, "output_ED5a_redo_20260929")
EVO = os.path.join(OUT, "07_karyotype_evolution")
DOT = os.path.join(EVO, "dotplot")

CASES = {
    "Mpra": dict(
        name="Medicago praecox", x=7, rows=["5", "6"], amk=[5, 6, 8],
        ref_label=r"$\mathit{M.\ truncatula}$ R108", ref_key="Mtru_R108", ref_cols=["5", "6", "8"],
        amk_blast=os.path.join(ROOT, "output_Kary/genome_410/all.blastp.txt"),
        ref_blast=os.path.join(DOT, "Mpra_vs_R108.blastp.txt"),
        ref_gff=os.path.join(ROOT, "data/genome_R108/old/genome_R108.gff")),
    "Mpol": dict(
        name="Medicago polymorpha", x=7, rows=["3", "5"], amk=[3, 5, 6],
        ref_label=r"$\mathit{M.\ truncatula}$ R108", ref_key="Mtru_R108", ref_cols=["3", "5", "6"],
        amk_blast=os.path.join(ROOT, "output_Kary/genome_Mpo_rerun_20260618/all.blastp.txt"),
        ref_blast=os.path.join(DOT, "Mpol_vs_R108.blastp.txt"),
        ref_gff=os.path.join(ROOT, "data/genome_R108/old/genome_R108.gff")),
    # all seven chromosomes are AMK mosaics; the R108 columns would repeat the AMK
    # columns (R108 is fully collinear with the AMK), so only the AMK is shown
    "Mfis": dict(
        name="Medicago fischeriana", x=7, rows=[str(i) for i in range(1, 8)],
        amk=list(range(1, 9)), ref_cols=[],
        amk_blast=os.path.join(ROOT, "output_Kary/genome_M46_2/all.blastp.txt")),
}
# x = 8 genomes with inter-chromosomal changes: (species chromosomes, AMK involved,
# minimum segment size in genes used to call junctions; 30 for the small translocations)
X8 = {
    "Malbus": (["6", "8"], [5, 6], 100),
    "Mlan": (["1", "6"], [1, 6], 100),
    "Mrut": (["5", "7"], [5, 7], 30),
    "Marc": (["5", "7"], [5, 7], 30),
    "Mrad": (["2", "6"], [2, 6], 100),
    "Medg": (["1", "3"], [1, 3], 100),
    "Mlup": (["5", "6"], [5, 6], 100),
    "Msuf": (["1", "6", "8"], [1, 6, 8], 100),
    "Msec": (["2", "6", "7", "8"], [2, 6, 7, 8], 100),
    "Msat_zm4": (["2", "7"], [2, 7], 30),
}
R108_GFF = os.path.join(ROOT, "data/genome_R108/old/genome_R108.gff")


def x8_case(lab):
    rows, amk, mg = X8[lab]
    sp = pd.read_csv(os.path.join(OUT, "species_table.tsv"), sep="\t").set_index("label")
    man = pd.read_csv(os.path.join(OUT, "04_wgdi_p0.2", "run_manifest.tsv"), sep="\t").set_index("label")
    return dict(name=sp.loc[lab, "full_name"], x=int(sp.loc[lab, "basic_x"]), rows=rows, amk=amk,
                min_genes=mg, ref_label=r"$\mathit{M.\ truncatula}$ R108", ref_key="Mtru_R108",
                ref_cols=[str(a) for a in amk], amk_blast=man.loc[lab, "blast"],
                ref_blast=os.path.join(DOT, lab + "_vs_R108.blastp.txt"), ref_gff=R108_GFF)


def junctions(lab, min_genes):
    """Inter-AMK junctions from segments_min30.tsv (same rule as junctions.py)."""
    seg = pd.read_csv(os.path.join(EVO, "segments_min30.tsv"), sep="\t")
    seg["chr"] = seg["chr"].astype(str)
    seg = seg[(seg.species == lab) & (seg.genes >= min_genes)]
    out = []
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
            out.append(dict(species=lab, chr=ch, amk_a=la["AMK"],
                            bp_a=la["amk_end"] if la["ori"] == "+" else la["amk_start"],
                            amk_b=fb["AMK"],
                            bp_b=fb["amk_start"] if fb["ori"] == "+" else fb["amk_end"],
                            species_pos=int(la["end"])))
    return pd.DataFrame(out, columns=["species", "chr", "amk_a", "bp_a", "amk_b", "bp_b",
                                      "species_pos"])


RANK_COL = {1: "#D62728", 2: "#1F5FBF"}
OTHER = "#B8B8B8"
JCOL, BCOL = "#1FA34A", "#17A2B8"


def read_gff(path):
    g = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 5], names=["chr", "gene", "order"])
    g["chr"] = g["chr"].astype(str)
    return g.set_index("gene")


def read_hits(path, qgff, sgff, rows, cols):
    b = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 10, 11], names=["q", "s", "e", "bit"])
    b = b[(b.e <= 1e-5) & (b.bit >= 100)]
    b = b.sort_values(["q", "bit"], ascending=[True, False])
    b["rank"] = b.groupby("q").cumcount() + 1
    b = b[b["rank"] <= 5]
    b = b[b.q.isin(qgff.index) & b.s.isin(sgff.index)]
    b["qchr"] = qgff.loc[b.q, "chr"].values
    b["qord"] = qgff.loc[b.q, "order"].values
    b["schr"] = sgff.loc[b.s, "chr"].values
    b["sord"] = sgff.loc[b.s, "order"].values
    return b[b.qchr.isin(rows) & b.schr.isin(cols)]


def main(lab):
    c = CASES[lab] if lab in CASES else x8_case(lab)
    d = os.path.join(OUT, "04_wgdi_p0.2", lab)
    qgff = read_gff(os.path.join(d, "in.gff1"))
    agff = read_gff(os.path.join(d, "in.gff2"))
    amk_cols = [str(a) for a in c["amk"]]
    ha = read_hits(c["amk_blast"], qgff, agff, c["rows"], amk_cols)
    qlen = qgff.groupby("chr")["order"].max()
    alen = agff.groupby("chr")["order"].max()
    lens = P.read_lens(os.path.join(d, "in.lens1"))
    segs = P.fill_gaps(P.read_km(os.path.join(d, "km_result.txt")), lens)
    if c["ref_cols"]:
        rgff = read_gff(c["ref_gff"])
        hr = read_hits(c["ref_blast"], qgff, rgff, c["rows"], c["ref_cols"])
        rlen = rgff.groupby("chr")["order"].max()
        rd = os.path.join(OUT, "04_wgdi_p0.2", c["ref_key"])
        rsegs = P.fill_gaps(P.read_km(os.path.join(rd, "km_result.txt")),
                            P.read_lens(os.path.join(rd, "in.lens1")))
    else:
        hr, rlen, rsegs = ha.iloc[0:0], {}, None
    J = junctions(lab, c.get("min_genes", 100))
    J = J[J.chr.astype(str).isin(c["rows"])]

    # ---- geometry (inches): 1 inch = K genes
    K = 7000.0
    gap, bar = 0.06, 0.09
    colw = [alen[a] / K for a in amk_cols] + [rlen[r] / K for r in c["ref_cols"]]
    rowh = [qlen[r] / K for r in c["rows"]]
    left, top_pad, mid = 0.75, 0.55, 0.22
    if not c["ref_cols"]:
        mid = 0.0
    W = max(left + sum(colw) + gap * (len(colw) - 1) + mid + 0.55, left + 6.0)
    H = top_pad + sum(rowh) + gap * (len(rowh) - 1) + 0.45
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")

    xs, x = [], left
    for i, w in enumerate(colw):
        if i == len(amk_cols):
            x += mid
        xs.append(x)
        x += w + gap
    ys, y = [], H - top_pad
    for h in rowh:
        ys.append(y)            # top of the row
        y -= h + gap

    # column headers
    for i, cx in enumerate(xs):
        if i < len(amk_cols):
            a = int(amk_cols[i])
            ax.add_patch(Rectangle((cx, H - top_pad + 0.03), colw[i], bar,
                                   fc=P.AMK_COLOR["AMK%d" % a], ec="none"))
            ax.text(cx + colw[i] / 2, H - top_pad + 0.03 + bar / 2, str(a), ha="center",
                    va="center", fontsize=7, color="white", fontweight="bold")
        else:
            rc = c["ref_cols"][i - len(amk_cols)]
            g = rsegs[rsegs.chr == rc]
            tot = rlen[rc]
            for _, s in g.iterrows():
                ax.add_patch(Rectangle((cx + colw[i] * (s.drawn_start - 1) / tot, H - top_pad + 0.03),
                                       colw[i] * (s.drawn_end - s.drawn_start + 1) / tot, bar,
                                       fc=P.AMK_COLOR[s.AMK], ec="none"))
            ax.text(cx + colw[i] / 2, H - top_pad + 0.03 + bar / 2, "Chr" + rc, ha="center",
                    va="center", fontsize=6, color="white", fontweight="bold")
    ax.text(xs[0] + (xs[len(amk_cols) - 1] + colw[len(amk_cols) - 1] - xs[0]) / 2, H - 0.12,
            r"AMK (= $\mathit{M.\ sativa}$ subsp. $\mathit{caerulea}$, T2T)", ha="center",
            va="top", fontsize=7.5)
    if c["ref_cols"]:
        ax.text(xs[len(amk_cols)] + (xs[-1] + colw[-1] - xs[len(amk_cols)]) / 2, H - 0.12,
                c["ref_label"] + " (x = 8)", ha="center", va="top", fontsize=7.5)

    # row colour bars and labels
    for j, r in enumerate(c["rows"]):
        g = segs[segs.chr == r]
        tot = qlen[r]
        for _, s in g.iterrows():
            ax.add_patch(Rectangle((left - 0.03 - bar, ys[j] - rowh[j] * s.drawn_end / tot), bar,
                                   rowh[j] * (s.drawn_end - s.drawn_start + 1) / tot,
                                   fc=P.AMK_COLOR[s.AMK], ec="none"))
        ax.text(left - 0.06 - bar, ys[j] - rowh[j] / 2, "Chr" + r, rotation=90, ha="right",
                va="center", fontsize=7, fontweight="bold")
    ax.text(0.06, H - top_pad - sum(rowh) / 2, c["name"] + "  (x = %d)" % c["x"], rotation=90,
            ha="left", va="center", fontsize=7.5, style="italic")

    # panels with dots
    for j, r in enumerate(c["rows"]):
        for i, cx in enumerate(xs):
            if i < len(amk_cols):
                t, hits, tl = amk_cols[i], ha, alen[amk_cols[i]]
            else:
                t = c["ref_cols"][i - len(amk_cols)]
                hits, tl = hr, rlen[t]
            ax.add_patch(Rectangle((cx, ys[j] - rowh[j]), colw[i], rowh[j], fc="none",
                                   ec="#555555", lw=0.6, zorder=3))
            sub = hits[(hits.qchr == r) & (hits.schr == t)]
            for rk in (5, 4, 3, 2, 1):
                s = sub[sub["rank"] == rk]
                ax.scatter(cx + colw[i] * s.sord / tl, ys[j] - rowh[j] * s.qord / qlen[r],
                           s=0.25 if rk > 1 else 0.35, c=RANK_COL.get(rk, OTHER), lw=0,
                           rasterized=True, zorder=2)

    # junctions: green across the row, cyan on the AMK columns
    for _, jn in J.iterrows():
        j = c["rows"].index(str(jn.chr))
        yj = ys[j] - rowh[j] * (jn.species_pos + 0.5) / qlen[str(jn.chr)]
        ax.plot([left, xs[-1] + colw[-1]], [yj, yj], color=JCOL, lw=0.6, ls=(0, (3, 2)), zorder=4)
        for amk_id, bp in ((jn.amk_a, jn.bp_a), (jn.amk_b, jn.bp_b)):
            if str(amk_id) in amk_cols:
                i = amk_cols.index(str(amk_id))
                xb = xs[i] + colw[i] * bp / alen[str(amk_id)]
                ax.plot([xb, xb], [ys[j], ys[j] - rowh[j]], color=BCOL, lw=0.5, ls=(0, (2, 2)),
                        zorder=4)
                ax.add_patch(plt.Circle((xb, yj), 0.07, fc="none", ec="#F28E2B", lw=0.9,
                                        ls=(0, (2, 1)), zorder=5))

    # legend
    ly = 0.22
    for k, (lbl, col) in enumerate([("best hit", RANK_COL[1]), ("second hit", RANK_COL[2]),
                                    ("3rd-5th hit", OTHER)]):
        ax.scatter([left + k * 1.0], [ly], s=10, c=col)
        ax.text(left + k * 1.0 + 0.07, ly, lbl, fontsize=6, va="center")
    ax.plot([left + 2.7, left + 2.95], [ly, ly], color=JCOL, lw=0.8, ls=(0, (3, 2)))
    ax.text(left + 3.0, ly, "junction between AMK segments", fontsize=6, va="center")
    ax.plot([left + 4.55, left + 4.55], [ly - 0.07, ly + 0.07], color=BCOL, lw=0.8, ls=(0, (2, 2)))
    ax.text(left + 4.62, ly, "AMK breakpoint", fontsize=6, va="center")

    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(DOT, "dotplot_evidence_%s.%s" % (lab, ext)),
                    dpi=400 if ext == "png" else 300)
    print("saved", len(ha), "AMK hits", len(hr), "ref hits", len(J), "junctions")


if __name__ == "__main__":
    main(sys.argv[1])
