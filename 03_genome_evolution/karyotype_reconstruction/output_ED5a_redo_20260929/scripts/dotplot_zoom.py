#!/usr/bin/env python3
"""Zoomed dot plots for small inter-chromosomal translocations.

2 x 2 panels: rows = (host chromosome window around the inserted segment,
donor chromosome window where the segment would sit ancestrally);
columns = (AMK window of the host region, AMK window of the moved segment).
A segment that has moved (not been copied) shows a collinear stretch in the
host row and a gap in the donor row.  AMK gene order = M. sativa subsp.
caerulea T2T.  Usage: dotplot_zoom.py Mrut Marc Msat_zm4
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("RUN", "04_wgdi_p0.2")
import plot_ed5a as P  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = "path/to/project/37.karyotype_reconstruction"
OUT = os.path.join(ROOT, "output_ED5a_redo_20260929")
EVO = os.path.join(OUT, "07_karyotype_evolution")
DOT = os.path.join(EVO, "dotplot")
# species: (host chromosome, donor AMK)
SMALL = {"Mrut": ("5", 7), "Marc": ("5", 7), "Msat_zm4": ("2", 7)}
PAD_HOST, PAD_AMK = 300, 250
RANK_COL = {1: "#D62728", 2: "#1F5FBF"}
OTHER = "#B8B8B8"
JCOL, BCOL = "#1FA34A", "#17A2B8"
K = 260.0          # genes per inch


def read_gff(path):
    g = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2, 5], names=["chr", "gene", "bp", "order"])
    g["chr"] = g["chr"].astype(str)
    return g.set_index("gene")


def hits(path, qg, sg):
    b = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 10, 11], names=["q", "s", "e", "bit"])
    b = b[(b.e <= 1e-5) & (b.bit >= 100)].sort_values(["q", "bit"], ascending=[True, False])
    b["rank"] = b.groupby("q").cumcount() + 1
    b = b[(b["rank"] <= 5) & b.q.isin(qg.index) & b.s.isin(sg.index)]
    for side, g in (("q", qg), ("s", sg)):
        b[side + "chr"] = g.loc[b[side], "chr"].values
        b[side + "ord"] = g.loc[b[side], "order"].values
    return b


def mb(g, ch, order):
    r = g[(g.chr == ch) & (g.order == int(order))]
    return r.bp.iloc[0] / 1e6 if len(r) else float("nan")


def main(lab):
    host, donor = SMALL[lab]
    d = os.path.join(OUT, "04_wgdi_p0.2", lab)
    man = pd.read_csv(os.path.join(OUT, "04_wgdi_p0.2", "run_manifest.tsv"), sep="\t").set_index("label")
    sp = pd.read_csv(os.path.join(OUT, "species_table.tsv"), sep="\t").set_index("label")
    qg, ag = read_gff(os.path.join(d, "in.gff1")), read_gff(os.path.join(d, "in.gff2"))
    H = hits(man.loc[lab, "blast"], qg, ag)
    seg = pd.read_csv(os.path.join(EVO, "segments_min30.tsv"), sep="\t")
    seg["chr"] = seg["chr"].astype(str)
    ins = seg[(seg.species == lab) & (seg.chr == host) & (seg.AMK == donor)].iloc[0]
    s0, s1, a0, a1 = int(ins.start), int(ins.end), int(ins.amk_start), int(ins.amk_end)
    host_amk = int(seg[(seg.species == lab) & (seg.chr == host) & (seg.AMK != donor)]
                   .sort_values("genes").iloc[-1].AMK)
    qlen = qg.groupby("chr")["order"].max()
    alen = ag.groupby("chr")["order"].max()

    # windows
    hw = (max(1, s0 - PAD_HOST), min(qlen[host], s1 + PAD_HOST))
    b1 = H[(H["rank"] == 1)]
    flank = b1[(b1.qchr == host) & (b1.qord.between(*hw)) & (~b1.qord.between(s0, s1)) &
               (b1.schr == str(host_amk))]
    ahw = (int(np.percentile(flank.sord, 2)) - 50, int(np.percentile(flank.sord, 98)) + 50)
    adw = (max(1, a0 - PAD_AMK), min(alen[str(donor)], a1 + PAD_AMK))
    don = b1[(b1.schr == str(donor)) & b1.sord.between(*adw) & (b1.qchr != host)]
    dchr = don.qchr.mode().iloc[0]
    don = don[don.qchr == dchr]
    # flanks of the ancestral position: AMK genes just before / after the moved interval
    lf = don[don.sord.between(a0 - 150, a0 - 1)].qord
    rf = don[don.sord.between(a1 + 1, a1 + 150)].qord
    lm = int(lf.median()) if len(lf) else int(rf.median())
    rm = int(rf.median()) if len(rf) else lm
    flank_note = ""
    if abs(rm - lm) <= 600:
        dw = (max(1, min(lm, rm) - 250), min(qlen[dchr], max(lm, rm) + 250))
    else:        # flanks separated by other rearrangements: show the left flank only
        dw = (max(1, lm - 250), min(qlen[dchr], lm + 250))
        flank_note = "; right flank (AMK%d %d+) lies elsewhere on Chr%s (%.2f Mb)" % (
            donor, a1 + 1, dchr, mb(qg, dchr, rm))
    kept = b1[(b1.qchr == dchr) & b1.sord.between(a0, a1) & (b1.schr == str(donor))]
    moved = b1[(b1.qchr == host) & b1.qord.between(s0, s1) & (b1.schr == str(donor))]
    n_amk = a1 - a0 + 1

    rows = [(host, hw), (dchr, dw)]
    cols = [(str(host_amk), ahw), (str(donor), adw)]
    colw = [(c[1][1] - c[1][0]) / K for c in cols]
    rowh = [(r[1][1] - r[1][0]) / K for r in rows]
    left, top, gap = 1.0, 0.75, 0.12
    Wf = left + sum(colw) + gap + 0.4
    Hf = top + sum(rowh) + gap + 0.95
    fig = plt.figure(figsize=(Wf, Hf))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, Wf)
    ax.set_ylim(0, Hf)
    ax.axis("off")
    xs = [left, left + colw[0] + gap]
    ys = [Hf - top, Hf - top - rowh[0] - gap]
    for i, (c, w) in enumerate(cols):
        ax.add_patch(Rectangle((xs[i], Hf - top + 0.05), colw[i], 0.1,
                               fc=P.AMK_COLOR["AMK" + c], ec="none"))
        ax.text(xs[i] + colw[i] / 2, Hf - top + 0.2, "AMK%s: genes %d-%d" % (c, w[0], w[1]),
                ha="center", va="bottom", fontsize=6.5)
    for j, (r, w) in enumerate(rows):
        ax.text(left - 0.08, ys[j] - rowh[j] / 2,
                "Chr%s\n%.2f-%.2f Mb" % (r, mb(qg, r, w[0]), mb(qg, r, w[1])),
                ha="right", va="center", fontsize=6.5, linespacing=1.3)
        for i, (c, cw) in enumerate(cols):
            x0, y0 = xs[i], ys[j]
            ax.add_patch(Rectangle((x0, y0 - rowh[j]), colw[i], rowh[j], fc="none", ec="#555555",
                                   lw=0.6, zorder=3))
            sub = H[(H.qchr == r) & (H.schr == c) & H.qord.between(*w) & H.sord.between(*cw)]
            for rk in (5, 4, 3, 2, 1):
                s = sub[sub["rank"] == rk]
                ax.scatter(x0 + (s.sord - cw[0]) / K, y0 - (s.qord - w[0]) / K,
                           s=3 if rk == 1 else 1.5, c=RANK_COL.get(rk, OTHER), lw=0, zorder=2)
    # insertion boundaries on the host row, moved AMK interval on the donor column
    for q in (s0 - 0.5, s1 + 0.5):
        yq = ys[0] - (q - hw[0]) / K
        ax.plot([left, xs[1] + colw[1]], [yq, yq], color=JCOL, lw=0.7, ls=(0, (3, 2)), zorder=4)
    for a in (a0 - 0.5, a1 + 0.5):
        xa = xs[1] + (a - adw[0]) / K
        ax.plot([xa, xa], [ys[0], ys[1] - rowh[1]], color=BCOL, lw=0.7, ls=(0, (2, 2)), zorder=4)
    ax.add_patch(Rectangle((xs[1] + (a0 - adw[0]) / K, ys[1] - rowh[1]), n_amk / K, rowh[1],
                           fc="#F28E2B", alpha=0.10, ec="none", zorder=1))
    name = sp.loc[lab, "full_name"]
    ax.text(0.08, Hf - 0.12, name, fontsize=8, style="italic", va="top")
    note = ("AMK%d genes %d-%d: %d best hits on Chr%s (genes %d-%d, %.2f-%.2f Mb), "
            "%d left on Chr%s" % (donor, a0, a1, len(moved), host, s0, s1,
                                  mb(qg, host, s0), mb(qg, host, s1), len(kept), dchr)) + flank_note
    ax.text(left, 0.55, note, fontsize=6.3, va="top")
    ax.text(left, 0.35, "red: best hit; blue: second hit; grey: 3rd-5th hit; green: boundaries of "
            "the inserted segment; cyan/orange: its AMK interval", fontsize=5.8, color="#555555",
            va="top")
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(DOT, "dotplot_zoom_%s.%s" % (lab, ext)), dpi=400)
    print(lab, note)


if __name__ == "__main__":
    for lab in sys.argv[1:]:
        main(lab)
