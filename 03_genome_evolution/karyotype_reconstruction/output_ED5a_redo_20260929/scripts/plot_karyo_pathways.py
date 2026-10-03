#!/usr/bin/env python3
"""Panel b: step-by-step inter-chromosomal rearrangement pathways from the AMK.

Inputs are AMK chromosomes (numbered), final products are the species'
chromosomes drawn from the WGDI -km data (Chr#), intermediates (I1, I2, I3) are
the most parsimonious schematic states consistent with the breakpoints in
segments_min30.tsv.  Intra-chromosomal inversions and orientation are not
modelled.  M. fischeriana: event order cannot be resolved from the data, so each
chromosome is shown with the AMK pieces (>= 100 genes) it is composed of.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("RUN", "04_wgdi_p0.2")
import plot_ed5a as P  # noqa: E402
import plot_karyo_evolution as K  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402
import pandas as pd  # noqa: E402

INK, MUTED, CL = P.INK, P.MUTED, K.CL
OPC = "#2E7D32"
AMK_N = {1: 6320, 2: 5612, 3: 6726, 4: 6736, 5: 5544, 6: 5318, 7: 6466, 8: 6214}
SCALE = 0.40 / 6736.0      # inches per gene
BW = 0.058                 # bar width
PLUS = 0.12                # gap holding a "+"
ROW_H = 0.66


def amk(a):
    return ("AMK", [(a, 1, AMK_N[a])], str(a))


def inter(name, pieces):
    return ("I", pieces, name)


class Chroms:
    def __init__(self, data):
        self.data = data

    def get(self, sp, ch):
        lens, segs = self.data[sp]
        g = segs[segs["chr"] == ch].sort_values("drawn_start")
        return ("CHR", [(int(r.AMK[3:]), int(r.drawn_start), int(r.drawn_end))
                        for _, r in g.iterrows()], "Chr" + ch)


def draw_bar(ax, x, ytop, pieces, scale=SCALE, w=BW):
    tot = sum(e - s + 1 for _, s, e in pieces)
    if tot * scale < 0.035:            # keep very short pieces visible
        scale = 0.035 / tot
    h = tot * scale
    style = "round,pad=0,rounding_size={}".format(w * 0.45)
    shape = FancyBboxPatch((x, ytop - h), w, h, boxstyle=style, fc=P.BACKBONE, ec="none", zorder=4)
    ax.add_patch(shape)
    y = ytop
    for a, s, e in pieces:
        hh = (e - s + 1) * scale
        p = Rectangle((x, y - hh), w, hh, fc=P.AMK_COLOR["AMK%d" % a], ec="white", lw=0.2, zorder=5)
        ax.add_patch(p)
        p.set_clip_path(shape)
        y -= hh
    ax.add_patch(FancyBboxPatch((x, ytop - h), w, h, boxstyle=style, fc="none", ec=P.EDGE,
                                lw=0.25, zorder=6))
    return h


def main_colour(item):
    kind, pieces, _ = item
    a = max(pieces, key=lambda p: p[2] - p[1])[0]
    return P.AMK_COLOR["AMK%d" % a]


def draw_group(ax, x, ytop, items, scale=SCALE):
    for k, it in enumerate(items):
        h = draw_bar(ax, x, ytop, it[1], scale)
        kind = it[0]
        ax.text(x + BW / 2, ytop - h - 0.02, it[2], ha="center", va="top",
                fontsize=4.6 if kind != "CHR" else 4.4,
                color=MUTED if kind == "AMK" else (INK if kind == "CHR" else "#8A8A88"),
                fontweight="bold" if kind == "CHR" else "normal")
        x += BW
        if k < len(items) - 1:
            ax.text(x + PLUS / 2, ytop - 0.14, "+", ha="center", va="center", fontsize=6.5, color=INK)
            x += PLUS
    return x


def draw_op(ax, x, ytop, op, cols):
    y = ytop - 0.15
    if op == "RCT":
        ax.annotate("", xy=(x + 0.12, y), xytext=(x + 0.02, y),
                    arrowprops=dict(arrowstyle="-|>", lw=0.5, color=INK, mutation_scale=4))
        cx, s = x + 0.22, 0.075
        ax.plot([cx - s, cx + s], [y + s * 1.2, y - s * 1.2], color=cols[0], lw=2.6,
                solid_capstyle="round", zorder=5)
        ax.plot([cx - s, cx + s], [y - s * 1.2, y + s * 1.2], color=cols[1], lw=2.6,
                solid_capstyle="round", zorder=5)
        ax.text(cx, y - s * 1.2 - 0.04, "RCT", ha="center", va="top", fontsize=4.6, color=OPC,
                fontweight="bold")
        ax.annotate("", xy=(x + 0.42, y), xytext=(x + 0.32, y),
                    arrowprops=dict(arrowstyle="-|>", lw=0.5, color=INK, mutation_scale=4))
        return x + 0.46
    ax.annotate("", xy=(x + 0.34, y), xytext=(x + 0.03, y),
                arrowprops=dict(arrowstyle="-|>", lw=0.6, color=INK, mutation_scale=5))
    ax.text(x + 0.185, y + 0.025, op, ha="center", va="bottom", fontsize=4.6, color=OPC,
            fontweight="bold")
    return x + 0.38


def draw_step(ax, x, ytop, ins, op, outs):
    x = draw_group(ax, x, ytop, ins)
    x = draw_op(ax, x + 0.02, ytop, op, [main_colour(ins[0]), main_colour(ins[1])])
    return draw_group(ax, x + 0.02, ytop, outs)


def box(ax, x, y, w, h, title, clade, subtitle=None):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.04",
                                fc="white", ec=CL[clade], lw=0.8, zorder=1))
    ax.text(x + 0.07, y + h - 0.05, title, ha="left", va="top", fontsize=5.6, style="italic",
            color=INK, linespacing=1.1)
    if subtitle:
        ax.text(x + w - 0.07, y + h - 0.05, subtitle, ha="right", va="top", fontsize=5,
                color=HILITE if "x = 7" in subtitle else MUTED,
                fontweight="bold" if "x = 7" in subtitle else "normal")


HILITE = P.HILITE


def mfis_pieces(seg, ch, min_genes=100):
    g = seg[(seg.species == "Mfis") & (seg.chr.astype(str) == ch) & (seg.genes >= min_genes)]
    g = g.sort_values("start")
    out = []
    for _, r in g.iterrows():
        if out and out[-1][0] == r.AMK:
            out[-1][1].append((int(r.AMK), int(r.amk_start), int(r.amk_end)))
        else:
            out.append((int(r.AMK), [(int(r.AMK), int(r.amk_start), int(r.amk_end))]))
    return out


def main():
    sp, data = K.load_species()
    C = Chroms(data)
    seg = pd.read_csv(os.path.join(K.DST, "segments_min30.tsv"), sep="\t")
    fig = plt.figure(figsize=(10.6, 4.05))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 10.6)
    ax.set_ylim(0, 4.05)
    ax.set_aspect("equal")
    ax.axis("off")

    # ---------------- row 1: single-step events
    W1, H1, Y1 = 1.6, 0.9, 3.05
    singles = [
        ("Melilotus albus", "out", "x = 8", [amk(5), amk(6)], "RCT",
         [C.get("Malbus", "6"), C.get("Malbus", "8")]),
        ("Medicago lanigera", "III", "x = 8", [amk(1), amk(6)], "RCT",
         [C.get("Mlan", "1"), C.get("Mlan", "6")]),
        ("M. ruthenica + M. archiducis-nicolai", "III", "x = 8", [amk(5), amk(7)], "Tr",
         [C.get("Mrut", "5"), C.get("Mrut", "7")]),
        ("Medicago radiata", "II", "x = 8", [amk(2), amk(6)], "RCT",
         [C.get("Mrad", "2"), C.get("Mrad", "6")]),
        ("Medicago edgeworthii", "II", "x = 8", [amk(1), amk(3)], "RCT",
         [C.get("Medg", "1"), C.get("Medg", "3")]),
        ("Medicago lupulina", "I", "x = 8", [amk(5), amk(6)], "RCT",
         [C.get("Mlup", "5"), C.get("Mlup", "6")]),
    ]
    x = 0.12
    for title, clade, sub, ins, op, outs in singles:
        if title.startswith("M. ruthenica"):
            title = "Medicago ruthenica +\nM. archiducis-nicolai"
        box(ax, x, Y1, W1, H1, title, clade, sub)
        draw_step(ax, x + 0.1, Y1 + H1 - 0.24, ins, op, outs)
        x += W1 + 0.12
    ax.text(0.0, 4.05, "b", fontsize=9, fontweight="bold", va="top")

    # ---------------- row 2: multi-step events
    Y2T = 2.9                      # top of row-2 boxes
    W2 = 1.62

    def multi(x, title, clade, sub, steps):
        h = 0.3 + ROW_H * len(steps)
        box(ax, x, Y2T - h, W2, h, title, clade, sub)
        yt = Y2T - 0.28
        for k, (ins, op, outs) in enumerate(steps):
            ax.add_patch(plt.Circle((x + 0.1, yt - 0.12), 0.045, fc="none", ec=MUTED, lw=0.5))
            ax.text(x + 0.1, yt - 0.122, str(k + 1), fontsize=4.8, color=MUTED,
                    va="center", ha="center")
            draw_step(ax, x + 0.2, yt, ins, op, outs)
            yt -= ROW_H
        return x + W2 + 0.13

    x = 0.12
    # M. suffruticosa: RCT(AMK1/AMK8) then RCT(I1/AMK6)
    I1 = inter("I1", [(8, 4, 1414), (1, 409, 6319)])
    x = multi(x, "Medicago suffruticosa", "I", "x = 8", [
        ([amk(1), amk(8)], "RCT", [C.get("Msuf", "8"), I1]),
        ([I1, amk(6)], "RCT", [C.get("Msuf", "1"), C.get("Msuf", "6")]),
    ])
    # M. secundiflora: two independent RCTs
    x = multi(x, "Medicago secundiflora", "I", "x = 8", [
        ([amk(2), amk(6)], "RCT", [C.get("Msec", "2"), C.get("Msec", "6")]),
        ([amk(7), amk(8)], "RCT", [C.get("Msec", "7"), C.get("Msec", "8")]),
    ])
    # M. polymorpha: RCT(AMK5/AMK6) -> RCT(I1/AMK3) -> EEJ(I3 + I2)
    I1p = inter("I1", [(5, 1, 2711), (6, 1, 2717)])
    I2p = inter("I2", [(6, 3468, 5299), (5, 2731, 5537)])
    I3p = inter("I3", [(3, 1781, 6722)])
    x = multi(x, "Medicago polymorpha", "I", "x = 7", [
        ([amk(5), amk(6)], "RCT", [I1p, I2p]),
        ([I1p, amk(3)], "RCT", [C.get("Mpol", "5"), I3p]),
        ([I3p, I2p], "EEJ", [C.get("Mpol", "3")]),
    ])
    # M. praecox: RCT(AMK6/AMK8) -> NCF(I1 into AMK5)
    I1r = inter("I1", [(8, 4, 3659), (6, 1, 501)])
    x = multi(x, "Medicago praecox", "I", "x = 7", [
        ([amk(6), amk(8)], "RCT", [C.get("Mpra", "6"), I1r]),
        ([I1r, amk(5)], "NCF", [C.get("Mpra", "5")]),
    ])

    # ---------------- M. fischeriana composition
    fw, fh = 3.35, 2.46
    fx, fy = x, Y2T - fh
    box(ax, fx, fy, fw, fh, "Medicago fischeriana", "II", "x = 7")
    ax.text(fx + 0.07, fy + fh - 0.17, "≥7 inter-chromosomal events (≥6 RCT + 1 fusion); "
            "order not resolvable - AMK pieces (≥100 genes) of each chromosome",
            fontsize=4.4, color=MUTED, va="top")
    sc = SCALE * 0.72
    for i in range(7):
        ch = str(i + 1)
        col, row = (0, i) if i < 4 else (1, i - 4)
        cx0 = fx + 0.12 + col * 1.65
        yt = fy + fh - 0.36 - row * 0.52
        pieces = mfis_pieces(seg, ch)
        xx = cx0
        for k, (a, pcs) in enumerate(pieces):
            draw_bar(ax, xx, yt, pcs, sc)
            ax.text(xx + BW / 2, yt - sum(e - s + 1 for _, s, e in pcs) * sc - 0.015, str(a),
                    ha="center", va="top", fontsize=4.2, color=MUTED)
            xx += BW
            if k < len(pieces) - 1:
                ax.text(xx + 0.05, yt - 0.1, "+", ha="center", va="center", fontsize=5.5)
                xx += 0.1
        ax.annotate("", xy=(xx + 0.22, yt - 0.1), xytext=(xx + 0.04, yt - 0.1),
                    arrowprops=dict(arrowstyle="-|>", lw=0.5, color=INK, mutation_scale=4))
        c = C.get("Mfis", ch)
        h = draw_bar(ax, xx + 0.27, yt, c[1], sc)
        ax.text(xx + 0.27 + BW / 2, yt - h - 0.015, c[2], ha="center", va="top", fontsize=4.4,
                fontweight="bold", color=INK)

    # ---------------- note
    ax.text(0.12, 0.08, "Numbers: AMK chromosomes; Chr#: chromosomes of the species (drawn from "
            "data); I1-I3: inferred intermediates (schematic, most parsimonious). "
            "Orientation and intra-chromosomal inversions are not modelled.",
            fontsize=4.8, color=MUTED, va="bottom")

    for ext in ("pdf", "svg", "png"):
        fig.savefig(os.path.join(K.DST, "karyotype_pathways." + ext), dpi=500 if ext == "png" else None)
    print("saved")


if __name__ == "__main__":
    main()
