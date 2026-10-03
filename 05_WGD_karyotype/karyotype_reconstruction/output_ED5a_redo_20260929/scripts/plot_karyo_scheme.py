#!/usr/bin/env python3
"""Karyotype-evolution scheme (Brassicaceae ACBK style, not a dated binary tree).

Only nodes that carry karyotype information are resolved: the AMK, the three
clades of Fig. 1b (Clade-III, Clade-II + Clade-I) and the M. ruthenica +
M. archiducis-nicolai node (shared small AMK7 -> AMK5 translocation).  Species
whose karyotype equals that of their ancestor hang from a dashed comb; branches
carrying an inter-chromosomal rearrangement are solid and bear an event box.
Branch lengths are not to scale.  Events/tip data as in plot_karyo_evolution.py.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("RUN", "04_wgdi_p0.2")
import plot_ed5a as P  # noqa: E402
import plot_karyo_evolution as K  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

INK, MUTED, HILITE, CL = P.INK, P.MUTED, P.HILITE, K.CL
TIPS = K.TIPS
FIG_W, FIG_H = 10.6, 6.35
COMB_Y = 4.25                  # dashed clade combs
TOP_Y = 4.85                   # AMK -> Clade-III / (II + I) split
II_I_Y = 4.55                  # Clade-II / Clade-I split
CLADE_TIPS = {"III": ["Mlan", "Marc", "Mrut"], "II": ["Mrad", "Medg", "Mfis"],
              "I": TIPS[7:]}
# event box centre y (staggered so that neighbouring boxes never touch)
BOX_Y = {("Malbus", 0): 3.45, ("Mlan", 0): 2.72, ("n_rm", 0): 3.72, ("Mrad", 0): 3.5,
         ("Medg", 0): 2.72, ("Mfis", 0): 3.55, ("Mlup", 0): 2.72, ("Msuf", 0): 3.5,
         ("Msec", 0): 3.5, ("Msec", 1): 2.66, ("Mpol", 0): 3.5, ("Mpra", 0): 2.72}
RM_NODE_Y = 3.2                # split below the shared Tr box


def main():
    sp, data = K.load_species()
    alen = K.amk_lens().reset_index(drop=True)
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, FIG_W)
    ax.set_ylim(0, FIG_H)
    ax.set_aspect("equal")
    ax.axis("off")
    tx = {t: K.tip_xy(i)[0] for i, t in enumerate(TIPS)}
    ty = {t: K.tip_xy(i)[1] for i, t in enumerate(TIPS)}
    changed = {e[0] for e in K.EVENTS}

    # ---- clade combs
    cx_clade = {}
    for key, tips in CLADE_TIPS.items():
        xs = [tx[t] for t in tips]
        if key == "III":
            xs = [tx["Mlan"], (tx["Marc"] + tx["Mrut"]) / 2]
        cx_clade[key] = (min(xs) + max(xs)) / 2
        if key == "II":   # centre would coincide with the M. edgeworthii branch
            cx_clade[key] = (tx["Mrad"] + tx["Medg"]) / 2
        ax.add_patch(plt.Circle((cx_clade[key], COMB_Y), 0.03, fc=INK, ec="none", zorder=2))
        ax.plot([min(xs), max(xs)], [COMB_Y, COMB_Y], color=CL[key], lw=1.1, ls=(0, (3, 2)))
        lx_, ha_ = (cx_clade[key] - 0.05, "right") if key == "II" else (min(xs) - 0.05, "left")
        ax.text(lx_, COMB_Y + 0.05, "Clade-" + key, color=CL[key], fontsize=6.5,
                fontweight="bold", ha=ha_, va="bottom")
    # ---- backbone: AMK -> III | (II + I)
    x_iii_i = (cx_clade["II"] + cx_clade["I"]) / 2
    stem = CL["stem"]
    ax.plot([cx_clade["III"], x_iii_i], [TOP_Y, TOP_Y], color=stem, lw=1.1)
    ax.plot([cx_clade["III"], cx_clade["III"]], [TOP_Y, COMB_Y], color=CL["III"], lw=1.1)
    ax.plot([x_iii_i, x_iii_i], [TOP_Y, II_I_Y], color=stem, lw=1.1)
    ax.plot([cx_clade["II"], cx_clade["I"]], [II_I_Y, II_I_Y], color=stem, lw=1.1)
    for key in ("II", "I"):
        ax.plot([cx_clade[key]] * 2, [II_I_Y, COMB_Y], color=CL[key], lw=1.1)
    ax.add_patch(plt.Circle((x_iii_i, II_I_Y), 0.03, fc=INK, ec="none", zorder=2))
    amk_x_bar = (cx_clade["III"] + x_iii_i) / 2
    ax.add_patch(plt.Circle((amk_x_bar, TOP_Y), 0.03, fc=INK, ec="none", zorder=2))

    # ---- AMK box on the Medicago stem and Melilotus outgroup
    bw_, bh_ = 1.05, 0.72
    bx, by = amk_x_bar - bw_ / 2, TOP_Y + 0.14
    root_y = by + bh_ + 0.16
    ax.plot([amk_x_bar] * 2, [TOP_Y, root_y], color=INK, lw=1.1, zorder=1)
    ax.add_patch(FancyBboxPatch((bx, by), bw_, bh_, boxstyle="round,pad=0,rounding_size=0.04",
                                fc="white", ec=INK, lw=0.7, zorder=3))
    amax = alen["genes"].max()
    for k in range(8):
        n = alen.loc[k, "genes"]
        xk = bx + 0.1 + k * 0.11
        K.bar(ax, xk, by + bh_ - 0.07, 0.06, 0.36 * n / amax, [(1, n, k + 1)], n)
        ax.text(xk + 0.03, by + bh_ - 0.45, str(k + 1), fontsize=5, ha="center", va="top", color=INK)
    ax.text(bx + bw_ / 2, by + 0.05, "AMK (x = 8)", fontsize=6.5, fontweight="bold",
            ha="center", va="bottom", color=INK)
    ax.text(bx - 0.06, by + bh_ / 2, "Medicago", fontsize=6, style="italic", ha="right",
            va="center", color=INK)
    mx = tx["Malbus"]
    ax.plot([mx, amk_x_bar], [root_y, root_y], color=INK, lw=1.1)
    ax.add_patch(plt.Circle((amk_x_bar, root_y), 0.03, fc=INK, ec="none", zorder=2))
    ax.plot([mx, mx], [root_y, ty["Malbus"]], color=CL["out"], lw=1.1, ls=(0, (3, 2)), zorder=1)

    # ---- species branches
    for key, tips in CLADE_TIPS.items():
        for t in tips:
            if t in ("Marc", "Mrut"):
                continue
            solid = t in changed
            ax.plot([tx[t]] * 2, [COMB_Y, ty[t]], color=CL[key], lw=1.1,
                    ls="-" if solid else (0, (3, 2)), zorder=1)
    # M. ruthenica + M. archiducis-nicolai: solid stem with the shared Tr, then dashed tips
    rx = (tx["Marc"] + tx["Mrut"]) / 2
    ax.plot([rx, rx], [COMB_Y, RM_NODE_Y], color=CL["III"], lw=1.1, zorder=1)
    ax.plot([tx["Marc"], tx["Mrut"]], [RM_NODE_Y, RM_NODE_Y], color=CL["III"], lw=1.1,
            ls=(0, (3, 2)))
    ax.add_patch(plt.Circle((rx, RM_NODE_Y), 0.03, fc=INK, ec="none", zorder=2))
    for t in ("Marc", "Mrut"):
        ax.plot([tx[t]] * 2, [RM_NODE_Y, ty[t]], color=CL["III"], lw=1.1, ls=(0, (3, 2)), zorder=1)

    # ---- event boxes
    seen = {}
    for ev in K.EVENTS:
        lab = ev[0]
        k = seen.get(lab, 0)
        seen[lab] = k + 1
        x = rx if lab == "n_rm" else tx[lab]
        clade = "III" if lab == "n_rm" else K.TIP_CLADE[lab]
        K.event_box(ax, x, BOX_Y[(lab, k)], ev, data, alen, clade)

    # ---- tips
    for i, lab in enumerate(TIPS):
        x, ytop = K.tip_xy(i)
        lens, segs = data[lab]
        s = sp.loc[lab]
        x7 = int(s.basic_x) == 7
        pad = 0.035
        ax.add_patch(FancyBboxPatch((x - K.TBW / 2 - pad, ytop - K.TBH - pad), K.TBW + 2 * pad,
                                    K.TBH + 2 * pad, boxstyle="round,pad=0,rounding_size=0.03",
                                    fc="white", ec=HILITE if x7 else "#9A9A98",
                                    lw=0.7 if x7 else 0.5, ls=(0, (2.5, 1.5)) if x7 else "-",
                                    zorder=3))
        P.draw_karyotype(ax, x - K.TBW / 2, ytop - K.TBH, K.TBW, K.TBH, lens, segs)
        yy = ytop - K.TBH - pad - 0.03
        for t, it in P.label_lines(s.full_name):
            ax.text(x, yy, t, ha="center", va="top", fontsize=5, style="italic" if it else "normal",
                    color=INK)
            yy -= 0.078
        ax.text(x, yy, "x = {}".format(int(s.basic_x)), ha="center", va="top", fontsize=5,
                color=HILITE if x7 else MUTED, fontweight="bold" if x7 else "normal")

    ax.text(0.0, FIG_H - 0.05, "a", fontsize=9, fontweight="bold", va="top")

    # ---- legend
    lx, ly = FIG_W - 2.55, FIG_H - 0.25
    items = [("RCT", "reciprocal translocation"), ("EEJ", "end-to-end joining (fusion)"),
             ("NCF", "nested chromosome fusion"), ("Tr", "small translocation (shared)")]
    for k, (t, d) in enumerate(items):
        ax.text(lx, ly - k * 0.12, t, fontsize=5.2, color="#2E7D32", fontweight="bold", va="top")
        ax.text(lx + 0.27, ly - k * 0.12, d, fontsize=5.2, color=INK, va="top")
    y2 = ly - 4 * 0.12 - 0.06
    ax.plot([lx, lx + 0.2], [y2 - 0.035] * 2, color=INK, lw=1.0)
    ax.text(lx + 0.27, y2, "inter-chromosomal rearrangement on this branch", fontsize=5.2,
            color=INK, va="top")
    ax.plot([lx, lx + 0.2], [y2 - 0.155] * 2, color=INK, lw=1.0, ls=(0, (3, 2)))
    ax.text(lx + 0.27, y2 - 0.12, "karyotype as in the ancestor above", fontsize=5.2,
            color=INK, va="top")
    ax.text(lx, y2 - 0.28, "Topology as in Fig. 1b; branch lengths not to scale;\n"
            "intra-chromosomal inversions not shown.", fontsize=4.8, color=MUTED, va="top",
            linespacing=1.3)

    for ext in ("pdf", "svg", "png"):
        fig.savefig(os.path.join(K.DST, "karyotype_scheme." + ext), dpi=500 if ext == "png" else None)
    print("saved")


if __name__ == "__main__":
    main()
