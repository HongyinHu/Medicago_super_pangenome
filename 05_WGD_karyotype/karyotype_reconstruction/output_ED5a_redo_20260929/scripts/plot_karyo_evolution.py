#!/usr/bin/env python3
"""Karyotype evolution of Medicago from the AMK, drawn on the species tree
(layout after the Brassicaceae ACBK figure).

Tree topology and node ages: species tree of Fig. 1b (Ma).
Tip karyotypes: WGDI -km on Ks-filtered, >=20-gene, one-to-one ortholog blocks
(04_wgdi_p0.2), gaps split at midpoints for display (as ED5a v2 filled).
Event boxes: inter-chromosomal rearrangements inferred from the segment table
(07_karyotype_evolution/segments_min30.tsv); each box shows the AMK chromosomes
involved -> event type / minimum number -> the resulting chromosomes of the
species.  No inter-chromosomal rearrangement is shared by two or more species
except the small AMK7 -> AMK5 translocation of M. ruthenica + M. archiducis-nicolai,
so every internal node keeps the AMK.
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("RUN", "04_wgdi_p0.2")
import plot_ed5a as P  # noqa: E402  (fonts, palette, readers, draw_karyotype, fill_gaps)

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402
import pandas as pd  # noqa: E402

OUT = P.OUT
W = P.W
DST = os.path.join(OUT, "07_karyotype_evolution")

INK, MUTED, HILITE = P.INK, P.MUTED, P.HILITE
CL = {"III": "#B23A55", "II": "#E08A2C", "I": "#3F8F7A", "stem": "#7B5EA7", "out": "#333333"}

# ---- species tree (Fig. 1b); tips left -> right
TIPS = ["Malbus", "Mlan", "Marc", "Mrut", "Mrad", "Medg", "Mfis", "Mlup", "Msuf", "Mcar",
        "Msec", "Morb", "Mpol", "Mpra", "Mtru_R108", "Mmar", "Mcre", "Msat_zm4", "Msat_cae"]
# node: (children, age Ma, clade colour key)
NODES = {
    "n_rm": (["Marc", "Mrut"], 1.15, "III"),
    "cIII": (["Mlan", "n_rm"], 5.05, "III"),
    "n_ef": (["Medg", "Mfis"], 4.22, "II"),
    "cII": (["Mrad", "n_ef"], 4.93, "II"),
    "n_ls": (["Mlup", "Msuf"], 4.68, "I"),
    "n_sc": (["Msat_zm4", "Msat_cae"], 1.09, "I"),
    "n_cr": (["Mcre", "n_sc"], 2.56, "I"),
    "n_ma": (["Mmar", "n_cr"], 3.01, "I"),
    "n_pt": (["Mpra", "Mtru_R108"], 3.62, "I"),
    "n_40": (["n_pt", "n_ma"], 4.03, "I"),
    "n_po": (["Mpol", "n_40"], 4.28, "I"),
    "n_or": (["Morb", "n_po"], 4.52, "I"),
    "n_se": (["Msec", "n_or"], 4.78, "I"),
    "n_ca": (["Mcar", "n_se"], 5.13, "I"),
    "cI": (["n_ls", "n_ca"], 5.50, "I"),
    "n_617": (["cII", "cI"], 6.17, "stem"),
    "crown": (["cIII", "n_617"], 6.63, "stem"),
}
TIP_CLADE = {t: "III" for t in ["Mlan", "Marc", "Mrut"]}
TIP_CLADE.update({t: "II" for t in ["Mrad", "Medg", "Mfis"]})
TIP_CLADE.update({t: "I" for t in TIPS[7:]})
TIP_CLADE["Malbus"] = "out"

# ---- inter-chromosomal events (curated from segments_min30.tsv)
# (species or node, [AMK inputs], label, [species chromosomes out], box centre y, note)
EVENTS = [
    ("Malbus", [5, 6], "RCT", ["6", "8"], 4.05),
    ("Mlan", [1, 6], "RCT", ["1", "6"], 3.55),
    ("n_rm", [5, 7], "Tr\n~2 Mb", ["5"], 4.05),
    ("Mrad", [2, 6], "RCT", ["2", "6"], 4.25),
    ("Medg", [1, 3], "RCT", ["1", "3"], 3.05),
    ("Mfis", list(range(1, 9)), "≥7 events\n≥6 RCT\n+ fusion", [], 3.95),
    ("Mlup", [5, 6], "RCT", ["5", "6"], 3.05),
    ("Msuf", [1, 6, 8], "2×RCT", ["1", "6", "8"], 3.95),
    ("Msec", [2, 6], "RCT", ["2", "6"], 4.05),
    ("Msec", [7, 8], "RCT", ["7", "8"], 3.0),
    ("Mpol", [3, 5, 6], "2×RCT\n+ EEJ", ["3", "5"], 3.55),
    ("Mpra", [5, 6, 8], "RCT\n+ NCF", ["5", "6"], 2.72),
]
BOX_WIDTH = {"Mfis": 0.56}

FIG_W, FIG_H = 10.6, 7.7
X0, DX = 0.62, 0.52
Y_T0 = 2.12                     # y of age 0
K = (6.25 - Y_T0) / 6.63        # inches per Myr
ROW_A_TOP, ROW_B_TOP = 2.02, 1.18
TBW, TBH = 0.46, 0.34           # tip box bar area
BOXW = 0.46                     # event box width


def yt(age):
    return Y_T0 + K * age


def tip_xy(i):
    return X0 + i * DX, (ROW_A_TOP if i % 2 == 0 else ROW_B_TOP)


def load_species():
    sp = pd.read_csv(os.path.join(OUT, "species_table.tsv"), sep="\t").set_index("label")
    data = {}
    for lab in TIPS:
        d = os.path.join(W, lab)
        lens = P.read_lens(os.path.join(d, "in.lens1"))
        segs = P.fill_gaps(P.read_km(os.path.join(d, "km_result.txt")), lens)
        data[lab] = (lens, segs)
    return sp, data


def amk_lens():
    return P.read_lens(os.path.join(W, "Mfis", "in.lens2"))


def bar(ax, x, ytop, w, h, segs_rows, total):
    """One vertical chromosome bar with coloured segments; segs_rows = [(start,end,AMK)]."""
    style = "round,pad=0,rounding_size={}".format(w * 0.45)
    shape = FancyBboxPatch((x, ytop - h), w, h, boxstyle=style, fc=P.BACKBONE, ec="none", zorder=4)
    ax.add_patch(shape)
    for s, e, a in segs_rows:
        y1 = ytop - h * (s - 1) / float(total)
        y0 = ytop - h * e / float(total)
        p = Rectangle((x, y0), w, y1 - y0, fc=P.AMK_COLOR["AMK%d" % a], ec="white", lw=0.2, zorder=5)
        ax.add_patch(p)
        p.set_clip_path(shape)
    ax.add_patch(FancyBboxPatch((x, ytop - h), w, h, boxstyle=style, fc="none", ec=P.EDGE,
                                lw=0.25, zorder=6))


def chrom_rows(segs, ch):
    g = segs[segs["chr"] == ch]
    return [(r.drawn_start, r.drawn_end, int(r.AMK[3:])) for _, r in g.iterrows()]


def event_box(ax, cx, cy, ev, data, alen, clade):
    lab, ins, text, outs, _ = ev
    BOXW = BOX_WIDTH.get(lab, 0.46)
    bw, hmax = 0.045, 0.22
    n_in = len(ins)
    compact = n_in > 3
    lines = text.count("\n") + 1
    h = (0.06 + hmax + 0.05 + 0.085 * lines + 0.05 + (hmax if outs else 0) + 0.06) if not compact \
        else (0.06 + 0.13 + 0.05 + 0.085 * lines + 0.06)
    x0, y0 = cx - BOXW / 2, cy - h / 2
    ax.add_patch(FancyBboxPatch((x0, y0), BOXW, h, boxstyle="round,pad=0,rounding_size=0.03",
                                fc="white", ec=CL[clade], lw=0.7, zorder=3))
    y = y0 + h - 0.06
    amax = alen["genes"].max()
    if compact:
        slot = (BOXW - 0.06) / n_in
        for k, a in enumerate(ins):
            bar(ax, x0 + 0.03 + k * slot + (slot - 0.03) / 2, y, 0.03,
                0.13 * alen.loc[a - 1, "genes"] / amax, [(1, alen.loc[a - 1, "genes"], a)],
                alen.loc[a - 1, "genes"])
        y -= 0.13 + 0.05
    else:
        gap = 0.07
        tot = n_in * bw + (n_in - 1) * gap
        xs = cx - tot / 2
        for k, a in enumerate(ins):
            n = alen.loc[a - 1, "genes"]
            bar(ax, xs + k * (bw + gap), y, bw, hmax * n / amax, [(1, n, a)], n)
            ax.text(xs + k * (bw + gap) + bw / 2, y - hmax * n / amax - 0.012, str(a),
                    ha="center", va="top", fontsize=4.2, color=MUTED)
            if k < n_in - 1:
                ax.text(xs + k * (bw + gap) + bw + gap / 2, y - hmax * 0.45, "+", ha="center",
                        va="center", fontsize=6, color=INK)
        y -= hmax + 0.05
    ax.annotate("", xy=(cx, y - 0.085 * lines - 0.02), xytext=(cx, y + 0.01),
                arrowprops=dict(arrowstyle="-|>", lw=0.5, color=INK, shrinkA=0, shrinkB=0,
                                mutation_scale=4), zorder=6)
    ax.text(cx + 0.035, y - 0.085 * lines / 2 - 0.005, text, ha="left", va="center", fontsize=4.6,
            color="#2E7D32", fontweight="bold", linespacing=1.0, zorder=7)
    y -= 0.085 * lines + 0.05
    if outs:
        lens, segs = data[lab] if lab in data else data["Mrut"]
        ln = lens.set_index("chr")["genes"]
        smax = max(ln[c] for c in outs)
        gap = 0.07
        tot = len(outs) * bw + (len(outs) - 1) * gap
        xs = cx - tot / 2
        for k, c in enumerate(outs):
            bar(ax, xs + k * (bw + gap), y, bw, hmax * ln[c] / smax, chrom_rows(segs, c), ln[c])
    return y0, y0 + h


def main():
    sp, data = load_species()
    alen = amk_lens().reset_index(drop=True)
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, FIG_W)
    ax.set_ylim(0, FIG_H)
    ax.set_aspect("equal")
    ax.axis("off")

    # ---- node coordinates
    pos = {t: (tip_xy(i)[0], None) for i, t in enumerate(TIPS)}

    def place(n):
        if n in pos and pos[n][1] is None and n in TIPS:
            return pos[n][0]
        ch, age, _ = NODES[n]
        xs = [place(c) for c in ch]
        pos[n] = ((xs[0] + xs[-1]) / 2, yt(age))
        return pos[n][0]
    place("crown")

    def col(n):
        return CL[NODES[n][2]] if n in NODES else CL[TIP_CLADE[n]]

    # ---- branches
    for n, (ch, age, ck) in NODES.items():
        nx, ny = pos[n]
        xs = [pos[c][0] for c in ch]
        ax.plot([min(xs), max(xs)], [ny, ny], color=CL[ck] if ck != "stem" else CL["stem"],
                lw=1.1, solid_capstyle="butt", zorder=1)
        for c in ch:
            cx = pos[c][0]
            if c in NODES:
                cy = pos[c][1]
            else:
                cy = tip_xy(TIPS.index(c))[1]
            ax.plot([cx, cx], [ny, cy], color=col(c), lw=1.1, zorder=1)
        ax.add_patch(plt.Circle((nx, ny), 0.028, fc=INK, ec="none", zorder=2))
        if n not in ("crown",):
            ax.text(nx + 0.035, ny + 0.03, "{:.2f}".format(age), fontsize=4.2, color=MUTED,
                    ha="left", va="bottom", zorder=2)

    # ---- AMK root box and Melilotus outgroup (dashed, not to scale)
    cx, cy = pos["crown"]
    bw_, bh_ = 1.05, 0.72
    bx, by = cx - bw_ / 2, cy + 0.14
    root_y = by + bh_ + 0.16
    ax.plot([cx, cx], [cy, root_y], color=INK, lw=1.1, zorder=1)
    mx = tip_xy(0)[0]
    ax.plot([mx, cx], [root_y, root_y], color=INK, lw=1.1, zorder=1)
    ax.plot([mx, mx], [root_y, tip_xy(0)[1]], color=CL["out"], lw=1.1, ls=(0, (3, 2)), zorder=1)
    ax.add_patch(plt.Circle((cx, root_y), 0.03, fc=INK, ec="none", zorder=2))
    # axis break: the 18.41 Ma root is not to scale
    yb = by + bh_ + 0.08
    for d in (-0.025, 0.025):
        ax.plot([cx - 0.05, cx + 0.05], [yb + d - 0.02, yb + d + 0.02], color=INK, lw=0.6, zorder=4)
        ax.add_patch(Rectangle((cx - 0.04, yb - 0.018), 0.08, 0.036, fc="white", ec="none", zorder=3))
    ax.text(cx + 0.04, root_y + 0.03, "18.41", fontsize=4.2, color=MUTED, ha="left", va="bottom")
    ax.text(cx + 0.04, cy - 0.03, "6.63", fontsize=4.2, color=MUTED, ha="left", va="top")
    # AMK box on the Medicago stem, just above the crown node
    ax.add_patch(FancyBboxPatch((bx, by), bw_, bh_, boxstyle="round,pad=0,rounding_size=0.04",
                                fc="white", ec=INK, lw=0.7, zorder=3))
    amax = alen["genes"].max()
    for k in range(8):
        n = alen.loc[k, "genes"]
        xk = bx + 0.1 + k * 0.11
        bar(ax, xk, by + bh_ - 0.07, 0.06, 0.36 * n / amax, [(1, n, k + 1)], n)
        ax.text(xk + 0.03, by + bh_ - 0.07 - 0.36 - 0.02, str(k + 1), fontsize=5, ha="center",
                va="top", color=INK)
    ax.text(bx + bw_ / 2, by + 0.05, "AMK (x = 8)", fontsize=6.5, fontweight="bold",
            ha="center", va="bottom", color=INK)
    ax.text(bx - 0.06, by + bh_ / 2, "Medicago", fontsize=6, style="italic", ha="right",
            va="center", color=INK)

    # clade labels
    for key, node, dx in (("III", "cIII", -0.05), ("II", "cII", -0.05), ("I", "cI", -0.05)):
        x, y = pos[node]
        ax.text(x + dx, y + 0.06, "Clade-" + key, fontsize=6.2, fontweight="bold",
                color=CL[key], ha="right", va="bottom")

    # ---- tips
    for i, lab in enumerate(TIPS):
        x, ytop = tip_xy(i)
        lens, segs = data[lab]
        s = sp.loc[lab]
        x7 = int(s.basic_x) == 7
        pad = 0.035
        ax.add_patch(FancyBboxPatch((x - TBW / 2 - pad, ytop - TBH - pad), TBW + 2 * pad, TBH + 2 * pad,
                                    boxstyle="round,pad=0,rounding_size=0.03", fc="white",
                                    ec=HILITE if x7 else "#9A9A98", lw=0.7 if x7 else 0.5,
                                    ls=(0, (2.5, 1.5)) if x7 else "-", zorder=3))
        P.draw_karyotype(ax, x - TBW / 2, ytop - TBH, TBW, TBH, lens, segs)
        lines = P.label_lines(s.full_name)
        yy = ytop - TBH - pad - 0.03
        for t, it in lines:
            ax.text(x, yy, t, ha="center", va="top", fontsize=5, style="italic" if it else "normal",
                    color=INK)
            yy -= 0.078
        ax.text(x, yy, "x = {}".format(int(s.basic_x)), ha="center", va="top", fontsize=5,
                color=HILITE if x7 else MUTED, fontweight="bold" if x7 else "normal")

    # ---- event boxes
    for ev in EVENTS:
        lab = ev[0]
        if lab in NODES:
            ch = NODES[lab][0]
            ex = pos[lab][0]
            # sit on the stem between clade node and this node
            clade = NODES[lab][2]
            # draw on the stem: x of the node
        else:
            ex = pos[lab][0]
            clade = TIP_CLADE[lab]
        event_box(ax, ex, ev[4], ev, data, alen, clade if clade != "out" else "out")

    # ---- time axis
    ax.plot([0.18, 0.18], [yt(0), yt(6)], color=MUTED, lw=0.5)
    for t in range(0, 7):
        ax.plot([0.15, 0.18], [yt(t), yt(t)], color=MUTED, lw=0.5)
        ax.text(0.13, yt(t), str(t), fontsize=4.8, ha="right", va="center", color=MUTED)
    ax.text(0.06, yt(6.4), "Ma", fontsize=5, ha="center", va="bottom", color=MUTED)

    # ---- legend
    lx, ly = FIG_W - 2.55, FIG_H - 0.3
    for k, (t, d) in enumerate([
            ("RCT", "reciprocal translocation"),
            ("EEJ", "end-to-end joining (fusion)"),
            ("NCF", "nested chromosome fusion"),
            ("Tr", "small translocation (shared)")]):
        ax.text(lx, ly - k * 0.12, t, fontsize=5.2, color="#2E7D32", fontweight="bold", va="top")
        ax.text(lx + 0.27, ly - k * 0.12, d, fontsize=5.2, color=INK, va="top")
    ax.text(lx, ly - 4 * 0.12 - 0.04, "All internal nodes retain the AMK (x = 8);\n"
            "intra-chromosomal inversions are not shown.\nNode ages in Ma (Fig. 1b).",
            fontsize=4.8, color=MUTED, va="top", linespacing=1.3)

    # ---- event table with segment evidence (segments >= 100 genes)
    seg = pd.read_csv(os.path.join(DST, "segments_min30.tsv"), sep="\t")
    seg["chr"] = seg["chr"].astype(str)
    rows = []
    for lab, ins, text, outs, _ in EVENTS:
        spc = "Mrut+Marc" if lab == "n_rm" else lab
        who = ["Mrut", "Marc"] if lab == "n_rm" else [lab]
        chroms = outs if outs else sorted(seg[seg.species == lab].chr.unique(), key=int)
        ev = []
        for w in who:
            for c in chroms:
                g = seg[(seg.species == w) & (seg.chr == c) & ((seg.genes >= 100) | (lab == "n_rm"))]
                ev.append("{} chr{}: ".format(w, c) + " | ".join(
                    "AMK{}[{}-{}]{} {}g {}Mb".format(r.AMK, r.amk_start, r.amk_end, r.ori, r.genes, r.Mb)
                    for _, r in g.iterrows()))
        rows.append({"branch": spc, "event": text.replace("\n", " "),
                     "AMK_involved": ",".join("AMK%d" % a for a in ins),
                     "resulting_chromosomes": ",".join(chroms), "evidence": " || ".join(ev)})
    pd.DataFrame(rows).to_csv(os.path.join(DST, "karyotype_events.tsv"), sep="\t", index=False)

    os.makedirs(DST, exist_ok=True)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(os.path.join(DST, "karyotype_evolution." + ext),
                    dpi=500 if ext == "png" else None)
    print("saved")


if __name__ == "__main__":
    main()
