#!/usr/bin/env python3
"""Redraw ED5a as a radial layout: AMK in the centre, 19 genomes around it.

Input : 02_wgdi/<label>/km_result.txt (WGDI -km on Ks-filtered, >=20-gene,
        genome-wide one-to-one orthologous blocks), in.lens1, in.gff1.
Output: 03_figure/ED5a_radial.{pdf,png,svg} and ED5a_source_data.tsv
Coordinates are WGDI gene-order indices; bar heights are gene counts scaled to
the largest chromosome of each species (same convention as the previous ED5a).
Species are placed clockwise from the top in the tip order of the ED5a tree.
"""
import glob
import itertools
import math
import os

import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
import pandas as pd

OUT = ("path/to/project/37.karyotype_reconstruction/"
       "output_ED5a_redo_20260929")
W = os.path.join(OUT, os.environ.get("RUN", "02_wgdi"))
FIG = os.path.join(OUT, os.environ.get("FIGDIR", "03_figure"))

for ttf in glob.glob(os.path.join(OUT, "fonts", "*.ttf")):
    font_manager.fontManager.addfont(ttf)
matplotlib.rcParams.update({
    "font.family": "Arial", "pdf.fonttype": 42, "svg.fonttype": "none",
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
})

# WGDI ancestor colour -> AMK id and display colour (validated palette, ED5a hue order)
AMK = {
    "royalblue": ("AMK1", "#3F5FA8"),
    "red": ("AMK2", "#E3262A"),
    "#99cc00": ("AMK3", "#7FB53A"),
    "deepskyblue": ("AMK4", "#29A3D8"),
    "#339966": ("AMK5", "#1E7B4F"),
    "#ffcc00": ("AMK6", "#E0A11B"),
    "fuchsia": ("AMK7", "#9B4FB0"),
    "#aa6e20": ("AMK8", "#A5561A"),
}
AMK_COLOR = {v[0]: v[1] for v in AMK.values()}
BACKBONE = "#E4E4E2"
EDGE = "#8A8A88"
INK = "#222222"
MUTED = "#6B6B6B"
HILITE = "#C0392B"

FIG_W, FIG_H = 7.2, 8.7          # inches (183 mm wide, Nature ED full width)
CX, CY = FIG_W / 2, FIG_H / 2
RX, RY = 3.08, 3.72              # ellipse carrying the panel centres
PW, PH = 0.70, 0.44              # panel bar-area width / height (inches)
FS = 5.5                         # label font size (pt)
LH = 0.088                       # label line height (inches)
PAD = 0.045                      # frame padding around bars


def read_lens(path):
    df = pd.read_csv(path, sep="\t", header=None, names=["chr", "bp", "genes"])
    df["chr"] = df["chr"].astype(str)
    return df


def read_gff(path):
    df = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2, 3, 5],
                     names=["chr", "gene", "start", "end", "order"])
    df["chr"] = df["chr"].astype(str)
    return df.set_index(["chr", "order"])


def read_km(path):
    df = pd.read_csv(path, sep="\t", header=None, names=["chr", "start", "end", "color", "cls"])
    df["chr"] = df["chr"].astype(str)
    df["AMK"] = [AMK[c.strip().lower()][0] for c in df["color"]]
    return df


FILL_GAPS = os.environ.get("FILL_GAPS", "none")   # "none" | "midpoint"


def fill_gaps(segs, lens):
    """For display only: split each unassigned interval between two adjacent
    segments at its midpoint; merge neighbours that carry the same AMK.
    Adds drawn_start/drawn_end; start/end keep the WGDI -km evidence coordinates."""
    out = []
    glen = lens.set_index("chr")["genes"]
    for ch, g in segs.groupby("chr", sort=False):
        g = g.sort_values("start").to_dict("records")
        for r in g:
            r["drawn_start"], r["drawn_end"] = r["start"], r["end"]
            r["evidence_genes"] = r["end"] - r["start"] + 1
        for a, b in zip(g, g[1:]):
            gap = b["start"] - a["end"] - 1
            if gap > 0:
                a["drawn_end"] = a["end"] + gap // 2
                b["drawn_start"] = a["drawn_end"] + 1
        g[0]["drawn_start"], g[-1]["drawn_end"] = 1, int(glen[ch])
        merged = [g[0]]
        for r in g[1:]:
            if r["AMK"] == merged[-1]["AMK"]:
                merged[-1]["drawn_end"] = r["drawn_end"]
                merged[-1]["end"] = r["end"]
                merged[-1]["evidence_genes"] += r["evidence_genes"]
            else:
                merged.append(r)
        out.extend(merged)
    return pd.DataFrame(out)


def draw_karyotype(ax, x0, y0, w, h, lens, segs, bar_frac=0.55):
    """Vertical chromosome bars (chr1 left), grey where no AMK segment is assigned."""
    slot = w / len(lens)
    bw = slot * bar_frac
    gmax = lens["genes"].max()
    for i, (_, c) in enumerate(lens.iterrows()):
        bx = x0 + i * slot + (slot - bw) / 2
        bh = h * c["genes"] / float(gmax)
        top = y0 + h
        style = "round,pad=0,rounding_size={}".format(bw * 0.45)
        shape = FancyBboxPatch((bx, top - bh), bw, bh, boxstyle=style, fc=BACKBONE, ec="none", zorder=2)
        ax.add_patch(shape)
        s0, s1 = ("drawn_start", "drawn_end") if "drawn_start" in segs else ("start", "end")
        for _, s in segs[segs["chr"] == c["chr"]].iterrows():
            y_hi = top - bh * (s[s0] - 1) / float(c["genes"])
            y_lo = top - bh * s[s1] / float(c["genes"])
            p = Rectangle((bx, y_lo), bw, y_hi - y_lo, fc=AMK_COLOR[s["AMK"]],
                          ec="white", lw=0.25, zorder=3)
            ax.add_patch(p)
            p.set_clip_path(shape)
        ax.add_patch(FancyBboxPatch((bx, top - bh), bw, bh, boxstyle=style,
                                    fc="none", ec=EDGE, lw=0.3, zorder=4))


def label_lines(full_name):
    """[(text, italic)] : genus / epithet / optional rank or accession line."""
    genus, rest = full_name.split(" ", 1)
    parts = rest.split(" ")
    if len(parts) == 2 and any(ch.isdigit() for ch in parts[1]):
        # accession names (ZM4, R108) share the epithet line, set roman
        return [(genus, True), (r"$\mathit{%s}$ %s" % (parts[0], parts[1]), False)]
    lines = [(genus, True), (parts[0], True)]
    if len(parts) > 1:
        extra = " ".join(p if p in ("subsp.", "var.", "cv.") or any(ch.isdigit() for ch in p)
                         else r"$\mathit{%s}$" % p for p in parts[1:])
        lines.append((extra, False))
    return lines


def main():
    os.makedirs(FIG, exist_ok=True)
    sp = pd.read_csv(os.path.join(OUT, "species_table.tsv"), sep="\t").sort_values("order")
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, FIG_W)
    ax.set_ylim(0, FIG_H)
    ax.set_aspect("equal")
    ax.axis("off")

    # ---- centre: AMK, doubles as the numbered colour key
    anc = pd.read_csv(os.path.join(W, "Mfis", "in.ancestor"), sep="\t", header=None,
                      names=["chr", "start", "end", "color", "cls"])
    anc["chr"] = anc["chr"].astype(str)
    anc["AMK"] = [AMK[c.strip().lower()][0] for c in anc["color"]]
    alens = read_lens(os.path.join(W, "Mfis", "in.lens2"))
    AR = 1.0
    aw, ah = 1.24, 0.74
    ax.add_patch(plt.Circle((CX, CY), AR, fc="#F6F6F4", ec="#D5D5D2", lw=0.6, zorder=1))
    draw_karyotype(ax, CX - aw / 2, CY - ah / 2 + 0.06, aw, ah, alens, anc, bar_frac=0.6)
    for i in range(8):
        ax.text(CX - aw / 2 + (i + 0.5) * aw / 8, CY - ah / 2 + 0.01, str(i + 1),
                ha="center", va="top", fontsize=6.5, color=INK)
    ax.text(CX, CY + ah / 2 + 0.16, "AMK", ha="center", va="bottom",
            fontsize=10, fontweight="bold", color=INK)
    ax.text(CX, CY - ah / 2 - 0.16, "Ancestral ", ha="right", va="top", fontsize=5.5, color=MUTED)
    ax.text(CX, CY - ah / 2 - 0.16, "Medicago", ha="left", va="top", fontsize=5.5,
            color=MUTED, style="italic")
    ax.text(CX, CY - ah / 2 - 0.25, "karyotype (x = 8)", ha="center", va="top",
            fontsize=5.5, color=MUTED)

    # ---- satellites
    src_rows, boxes = [], []
    n = len(sp)
    for k, (_, s) in enumerate(sp.iterrows()):
        theta = math.radians(90 - 360.0 * k / n)
        px, py = CX + RX * math.cos(theta), CY + RY * math.sin(theta)
        d = os.path.join(W, s.label)
        lens = read_lens(os.path.join(d, "in.lens1"))
        segs = read_km(os.path.join(d, "km_result.txt"))
        if FILL_GAPS == "midpoint":
            segs = fill_gaps(segs, lens)
        gff = read_gff(os.path.join(d, "in.gff1"))
        assert len(lens) == int(s.basic_x) == segs["chr"].nunique(), s.label
        x7 = int(s.basic_x) == 7
        lines = label_lines(s.full_name)
        label_h = (len(lines) + 1) * LH

        # bars sit at (px, py); labels go on the outer side so spokes never cross them
        fx0, fy0 = px - PW / 2 - PAD, py - PH / 2 - PAD
        fw, fh = PW + 2 * PAD, PH + 2 * PAD
        above = math.sin(theta) > 0.3
        if above:
            boxes.append((s.label, fx0, fy0, fx0 + fw, fy0 + fh + 0.03 + label_h))
        else:
            boxes.append((s.label, fx0, fy0 - 0.03 - label_h, fx0 + fw, fy0 + fh))

        # spoke: AMK disc edge -> nearest point of the frame (clipped at the frame)
        phi = math.atan2(py - CY, px - CX)
        sx, sy = CX + (AR + 0.02) * math.cos(phi), CY + (AR + 0.02) * math.sin(phi)
        tx = (fw / 2 + 0.03) / max(abs(math.cos(phi)), 1e-9)
        ty = (fh / 2 + 0.03) / max(abs(math.sin(phi)), 1e-9)
        t = min(tx, ty)
        ax.plot([sx, px - t * math.cos(phi)], [sy, py - t * math.sin(phi)],
                color="#C4C4C1", lw=0.5, zorder=0, solid_capstyle="round")

        ax.add_patch(FancyBboxPatch((fx0, fy0), fw, fh, boxstyle="round,pad=0,rounding_size=0.035",
                                    fc="white", ec=HILITE if x7 else "#CFCFCC",
                                    lw=0.8 if x7 else 0.5,
                                    ls=(0, (2.5, 1.5)) if x7 else "-", zorder=1))
        draw_karyotype(ax, px - PW / 2, py - PH / 2, PW, PH, lens, segs)

        y = fy0 + fh + 0.03 + label_h if above else fy0 - 0.03
        for text, italic in lines:
            ax.text(px, y, text, ha="center", va="top", fontsize=FS, color=INK,
                    style="italic" if italic else "normal")
            y -= LH
        ax.text(px, y, "x = {}".format(int(s.basic_x)), ha="center", va="top", fontsize=FS,
                color=HILITE if x7 else MUTED, fontweight="bold" if x7 else "normal")

        chrom_genes = lens.set_index("chr")["genes"]
        for _, g in segs.iterrows():
            a = gff.loc[(g["chr"], int(g["start"]))]
            b = gff.loc[(g["chr"], int(g["end"]))]
            src_rows.append({
                "order": int(s.order), "label": s.label, "species": s.full_name,
                "basic_x": int(s.basic_x), "analysis_dir": d,
                "chromosome": g["chr"], "AMK": g["AMK"], "display_color": AMK_COLOR[g["AMK"]],
                "start_gene_index": int(g["start"]), "end_gene_index": int(g["end"]),
                "segment_gene_count": int(g["end"] - g["start"] + 1),
                "chromosome_gene_count": int(chrom_genes[g["chr"]]),
                "start_gene_id": a["gene"], "start_bp": int(a["start"]),
                "end_gene_id": b["gene"], "end_bp": int(b["end"]),
            })
            if "drawn_start" in segs:
                ds, de = int(g["drawn_start"]), int(g["drawn_end"])
                src_rows[-1].update({
                    "drawn_start_gene_index": ds, "drawn_end_gene_index": de,
                    "evidence_gene_count": int(g["evidence_genes"]),
                    # genes painted for display only (no -km evidence)
                    "display_filled_genes": (de - ds + 1) - int(g["evidence_genes"]),
                })

    # ---- layout QC: no panel/label box may overlap another or leave the canvas
    problems = []
    for a, b in itertools.combinations(boxes, 2):
        if a[1] < b[3] and b[1] < a[3] and a[2] < b[4] and b[2] < a[4]:
            problems.append("overlap {} / {}".format(a[0], b[0]))
    for bx in boxes:
        if bx[1] < 0.05 or bx[2] < 0.05 or bx[3] > FIG_W - 0.05 or bx[4] > FIG_H - 0.05:
            problems.append("off-canvas {}".format(bx[0]))
    print("layout QC:", "PASS" if not problems else "; ".join(problems))

    for ext in ("pdf", "png", "svg"):
        fig.savefig(os.path.join(FIG, "ED5a_radial." + ext), dpi=600 if ext == "png" else None)
    pd.DataFrame(src_rows).to_csv(os.path.join(FIG, "ED5a_source_data.tsv"), sep="\t", index=False)
    print("segments:", len(src_rows))


if __name__ == "__main__":
    main()
