#!/usr/bin/env python3
import argparse
import os
import math

from reportlab.lib import colors
from reportlab.lib.pagesizes import landscape, A4
from reportlab.lib.units import mm
from reportlab.pdfgen import canvas


W, H = landscape(A4)


def hex_color(value):
    return colors.HexColor(value)


COL = {
    "ink": hex_color("#17212b"),
    "muted": hex_color("#5a6673"),
    "panel": hex_color("#f8fafc"),
    "line": hex_color("#8b96a3"),
    "spiny": hex_color("#b9d9ee"),
    "spineless": hex_color("#f6c99e"),
    "ref": hex_color("#ded8f4"),
    "sv": hex_color("#9d7bb8"),
    "sv2": hex_color("#f2b179"),
    "green": hex_color("#7bbf87"),
    "red": hex_color("#d98074"),
    "blue": hex_color("#6fa8cf"),
    "yellow": hex_color("#f5d77f"),
    "gray": hex_color("#e5e7eb"),
}


def set_font(c, name="Helvetica", size=8, color=None):
    c.setFont(name, size)
    if color is not None:
        c.setFillColor(color)


def panel(c, x, y, w, h, number, title):
    c.setFillColor(COL["panel"])
    c.setStrokeColor(hex_color("#c7cdd4"))
    c.roundRect(x, y, w, h, 5, fill=1, stroke=1)
    c.setFillColor(colors.white)
    c.setStrokeColor(COL["ink"])
    c.circle(x + 6 * mm, y + h - 6 * mm, 3.8 * mm, fill=1, stroke=1)
    set_font(c, "Helvetica-Bold", 9, COL["ink"])
    c.drawCentredString(x + 6 * mm, y + h - 8.5 * mm, str(number))
    set_font(c, "Helvetica-Bold", 10, COL["ink"])
    c.drawString(x + 11 * mm, y + h - 8.5 * mm, title)


def text(c, x, y, s, size=7.5, bold=False, color=None, center=False):
    set_font(c, "Helvetica-Bold" if bold else "Helvetica", size, color or COL["ink"])
    if center:
        c.drawCentredString(x, y, s)
    else:
        c.drawString(x, y, s)


def arrow(c, x1, y1, x2, y2, color=None, width=1.0):
    c.setStrokeColor(color or COL["ink"])
    c.setFillColor(color or COL["ink"])
    c.setLineWidth(width)
    c.line(x1, y1, x2, y2)
    ang = math.atan2(y2 - y1, x2 - x1)
    L = 4.2 * mm
    spread = 0.45
    p1 = (x2, y2)
    p2 = (x2 - L * math.cos(ang - spread), y2 - L * math.sin(ang - spread))
    p3 = (x2 - L * math.cos(ang + spread), y2 - L * math.sin(ang + spread))
    path = c.beginPath()
    path.moveTo(*p1)
    path.lineTo(*p2)
    path.lineTo(*p3)
    path.close()
    c.drawPath(path, fill=1, stroke=0)


def chromosome(c, x, y, w, label, segments=None):
    c.setStrokeColor(COL["ink"])
    c.setLineWidth(1.5)
    c.roundRect(x, y - 2 * mm, w, 4 * mm, 2 * mm, fill=0, stroke=1)
    c.setLineWidth(0.4)
    for i in range(18):
        xx = x + (i + 1) * w / 20
        c.line(xx, y - 5 * mm, xx, y - 2.5 * mm)
    if segments:
        for sx, sw, col in segments:
            c.setFillColor(col)
            c.rect(x + sx * w, y - 2 * mm, sw * w, 4 * mm, fill=1, stroke=0)
    text(c, x - 1 * mm, y + 4 * mm, label, 7, bold=True)


def doc_icon(c, x, y, fill=colors.white):
    c.setFillColor(fill)
    c.setStrokeColor(COL["ink"])
    c.roundRect(x, y, 9 * mm, 12 * mm, 1.5 * mm, fill=1, stroke=1)
    c.setStrokeColor(COL["muted"])
    for i in range(3):
        c.line(x + 2 * mm, y + (4 + 2 * i) * mm, x + 7 * mm, y + (4 + 2 * i) * mm)


def dna_icon(c, x, y, scale=1.0, color1=None, color2=None):
    color1 = color1 or COL["blue"]
    color2 = color2 or COL["sv"]
    c.setLineWidth(1)
    for i in range(10):
        yy = y + i * 1.25 * mm * scale
        x1 = x + math.sin(i * 0.8) * 2.4 * mm * scale
        x2 = x - math.sin(i * 0.8) * 2.4 * mm * scale + 7 * mm * scale
        c.setStrokeColor(color1 if i % 2 == 0 else color2)
        c.line(x1, yy, x2, yy + 1.0 * mm * scale)
    c.setStrokeColor(color1)
    pts1 = []
    pts2 = []
    for i in range(28):
        yy = y + i * 0.45 * mm * scale
        pts1.append((x + math.sin(i * 0.33) * 2.6 * mm * scale, yy))
        pts2.append((x + 7 * mm * scale - math.sin(i * 0.33) * 2.6 * mm * scale, yy))
    for pts, col in ((pts1, color1), (pts2, color2)):
        c.setStrokeColor(col)
        p = c.beginPath()
        p.moveTo(*pts[0])
        for px, py in pts[1:]:
            p.lineTo(px, py)
        c.drawPath(p, stroke=1, fill=0)


def sample_card(c, x, y, color):
    c.setFillColor(color)
    c.setStrokeColor(hex_color("#b8c2cc"))
    c.roundRect(x, y, 12 * mm, 16 * mm, 2 * mm, fill=1, stroke=1)
    dna_icon(c, x + 2.5 * mm, y + 2.5 * mm, scale=0.7)


def mini_bar(c, x, y, frac, color, label):
    c.setStrokeColor(COL["ink"])
    c.setFillColor(colors.white)
    c.roundRect(x, y, 24 * mm, 3.2 * mm, 1.5 * mm, fill=1, stroke=1)
    c.setFillColor(color)
    c.roundRect(x, y, 24 * mm * frac, 3.2 * mm, 1.5 * mm, fill=1, stroke=0)
    text(c, x + 26 * mm, y - 0.5 * mm, label, 6.4)


def funnel(c, x, y, w, h):
    c.setFillColor(hex_color("#f1f5f9"))
    c.setStrokeColor(COL["ink"])
    p = c.beginPath()
    p.moveTo(x, y + h)
    p.lineTo(x + w, y + h)
    p.lineTo(x + w * 0.62, y + h * 0.45)
    p.lineTo(x + w * 0.56, y)
    p.lineTo(x + w * 0.44, y)
    p.lineTo(x + w * 0.38, y + h * 0.45)
    p.close()
    c.drawPath(p, fill=1, stroke=1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    c = canvas.Canvas(args.out, pagesize=(W, H))
    c.setTitle("Phenotype-guided comparative SV cosegregation workflow")
    c.setFillColor(colors.white)
    c.rect(0, 0, W, H, fill=1, stroke=0)

    margin = 4 * mm
    gap = 3 * mm
    left_w = 150 * mm
    right_w = W - 2 * margin - left_w - gap
    top_h = 64 * mm
    bottom_h = H - 2 * margin - top_h - gap
    x1 = margin
    x2 = margin + left_w + gap
    y_top = H - margin - top_h
    y_bot = margin

    # Panel 1
    panel(c, x1, y_top, left_w, top_h, 1, "Input Data & Phenotype Groups")
    group_w = 48 * mm
    gy = y_top + 11 * mm
    c.setFillColor(COL["spiny"])
    c.roundRect(x1 + 6 * mm, gy, group_w, 44 * mm, 4, fill=1, stroke=0)
    c.setFillColor(COL["spineless"])
    c.roundRect(x1 + left_w - group_w - 6 * mm, gy, group_w, 44 * mm, 4, fill=1, stroke=0)
    text(c, x1 + 30 * mm, gy + 37 * mm, "Spiny group", 9, True, center=True)
    text(c, x1 + left_w - 30 * mm, gy + 37 * mm, "Spineless group", 9, True, center=True)
    for i in range(4):
        sample_card(c, x1 + 10 * mm + i * 9 * mm, gy + 20 * mm, colors.white)
        sample_card(c, x1 + left_w - group_w + 1 * mm + i * 9 * mm, gy + 20 * mm, colors.white)
    doc_icon(c, x1 + 14 * mm, gy + 6 * mm)
    doc_icon(c, x1 + 34 * mm, gy + 6 * mm)
    doc_icon(c, x1 + left_w - group_w + 5 * mm, gy + 6 * mm)
    doc_icon(c, x1 + left_w - group_w + 25 * mm, gy + 6 * mm)
    text(c, x1 + 18.5 * mm, gy + 1.5 * mm, "VCF", 7, True, center=True)
    text(c, x1 + 38.5 * mm, gy + 1.5 * mm, "BAM", 7, True, center=True)
    text(c, x1 + left_w - group_w + 9.5 * mm, gy + 1.5 * mm, "VCF", 7, True, center=True)
    text(c, x1 + left_w - group_w + 29.5 * mm, gy + 1.5 * mm, "BAM", 7, True, center=True)
    cx = x1 + left_w / 2 - 15 * mm
    c.setFillColor(colors.white)
    c.setStrokeColor(COL["ink"])
    c.roundRect(cx, gy + 9 * mm, 30 * mm, 36 * mm, 4, fill=1, stroke=1)
    dna_icon(c, cx + 9 * mm, gy + 18 * mm, scale=1.35, color1=COL["sv"], color2=COL["green"])
    text(c, cx + 15 * mm, gy + 5 * mm, "panSV / target\nintervals", 7, True, center=True)
    arrow(c, x1 + 54 * mm, gy + 28 * mm, cx - 1 * mm, gy + 28 * mm)
    arrow(c, x1 + left_w - 54 * mm, gy + 28 * mm, cx + 31 * mm, gy + 28 * mm)

    # Panel 2
    panel(c, x2, y_top, right_w, 39 * mm, 2, "SV Equivalence Rules")
    base_x = x2 + 28 * mm
    chromosome(c, base_x, y_top + 25 * mm, 92 * mm, "Reference", [(0.35, 0.20, COL["sv"]), (0.60, 0.18, COL["sv2"])])
    chromosome(c, base_x, y_top + 12 * mm, 92 * mm, "Sample", [(0.42, 0.30, COL["sv2"])])
    text(c, base_x + 46 * mm, y_top + 32 * mm, "reciprocal overlap >= 0.6", 7, True, center=True)
    arrow(c, base_x + 36 * mm, y_top + 30 * mm, base_x + 56 * mm, y_top + 30 * mm)
    text(c, base_x + 7 * mm, y_top + 33 * mm, "breakpoint\n<= 800 bp", 6, center=True)
    text(c, base_x + 85 * mm, y_top + 33 * mm, "length ratio\n0.5-2", 6, center=True)
    c.setFillColor(colors.white)
    c.setStrokeColor(COL["line"])
    c.roundRect(x2 + 6 * mm, y_top + 5 * mm, right_w - 12 * mm, 8 * mm, 3, fill=1, stroke=1)
    text(c, x2 + right_w / 2, y_top + 7.5 * mm, "Strict SV equivalence first; gene/intron-level rescue for shifted breakpoints", 7, True, center=True)

    # Panel 3
    panel(c, x1, y_bot, 115 * mm, bottom_h, 3, "Pass 1: Event-level PAV Screening")
    chromosome(c, x1 + 12 * mm, y_bot + bottom_h - 32 * mm, 80 * mm, "SV event", [(0.15, 0.18, COL["blue"]), (0.48, 0.16, COL["sv"]), (0.78, 0.14, COL["sv2"])])
    c.setFillColor(colors.white)
    c.setStrokeColor(COL["ink"])
    c.roundRect(x1 + 41 * mm, y_bot + bottom_h - 61 * mm, 32 * mm, 14 * mm, 3, fill=1, stroke=1)
    text(c, x1 + 57 * mm, y_bot + bottom_h - 53 * mm, "SV-by-SV\nPAV scan", 7.2, True, center=True)
    mini_bar(c, x1 + 12 * mm, y_bot + bottom_h - 76 * mm, 0.90, COL["blue"], "spiny 90%")
    mini_bar(c, x1 + 12 * mm, y_bot + bottom_h - 84 * mm, 0.10, COL["sv2"], "spineless 10%")
    c.setFillColor(colors.white)
    c.roundRect(x1 + 55 * mm, y_bot + 10 * mm, 48 * mm, 17 * mm, 3, fill=1, stroke=1)
    text(c, x1 + 79 * mm, y_bot + 21 * mm, "80% rule", 8, True, center=True)
    text(c, x1 + 79 * mm, y_bot + 14 * mm, "match rate >= 0.8\nin both groups", 6.4, center=True)
    text(c, x1 + 12 * mm, y_bot + 11 * mm, "Retained:\ngroup-enriched", 7, True, COL["green"])
    text(c, x1 + 12 * mm, y_bot + 4.5 * mm, "Discarded:\nno pattern", 7, True, COL["red"])

    # Panel 4
    panel(c, x1 + 118 * mm, y_bot, 58 * mm, bottom_h, 4, "Pass 2: Gene/Intron Rescue")
    funnel(c, x1 + 134 * mm, y_bot + 36 * mm, 27 * mm, 45 * mm)
    for i, (dx, dy, col) in enumerate([
        (3, 34, COL["spiny"]), (11, 31, COL["spineless"]), (18, 28, COL["spiny"]),
        (8, 22, COL["ref"]), (16, 18, COL["spineless"]), (12, 9, COL["yellow"]),
    ]):
        doc_icon(c, x1 + 134 * mm + dx * mm, y_bot + 36 * mm + dy * mm, fill=col)
    text(c, x1 + 147.5 * mm, y_bot + 27 * mm, "Aggregate all DELs\nwithin target intron", 7.1, True, center=True)
    c.setFillColor(colors.white)
    c.roundRect(x1 + 127 * mm, y_bot + 9 * mm, 42 * mm, 14 * mm, 3, fill=1, stroke=1)
    text(c, x1 + 148 * mm, y_bot + 17.5 * mm, "States", 7.5, True, center=True)
    text(c, x1 + 148 * mm, y_bot + 11.2 * mm, "intact | small indel\nlarge DEL | ambiguous", 5.8, center=True)

    # Panel 5
    panel(c, x2, y_bot + bottom_h - 54 * mm, right_w, 54 * mm, 5, "Coverage Validation")
    chromosome(c, x2 + 12 * mm, y_bot + bottom_h - 26 * mm, right_w - 24 * mm, "Candidate DEL", [(0.18, 0.18, COL["green"]), (0.38, 0.28, COL["sv"]), (0.72, 0.18, COL["green"])])
    text(c, x2 + 35 * mm, y_bot + bottom_h - 14 * mm, "left flank", 6.5, center=True)
    text(c, x2 + right_w / 2, y_bot + bottom_h - 14 * mm, "SV interval", 6.5, center=True)
    text(c, x2 + right_w - 35 * mm, y_bot + bottom_h - 14 * mm, "right flank", 6.5, center=True)
    for dx, label, col in [(18, "pass\nflanks high\ninside low", COL["green"]), (72, "remove\nlow flank\ncoverage", COL["red"])]:
        c.setFillColor(colors.white)
        c.roundRect(x2 + dx * mm, y_bot + bottom_h - 48 * mm, 45 * mm, 16 * mm, 3, fill=1, stroke=1)
        text(c, x2 + (dx + 22.5) * mm, y_bot + bottom_h - 38 * mm, label, 6.4, True, col, center=True)

    # Panel 6
    panel(c, x2, y_bot, right_w, bottom_h - 57 * mm, 6, "Final Output: High-confidence Candidate Set")
    for i in range(7):
        sample_card(c, x2 + 12 * mm + i * 12 * mm, y_bot + 42 * mm, COL["spiny"] if i % 2 == 0 else COL["spineless"])
    arrow(c, x2 + 99 * mm, y_bot + 50 * mm, x2 + 117 * mm, y_bot + 50 * mm)
    doc_icon(c, x2 + 121 * mm, y_bot + 42 * mm, fill=COL["ref"])
    doc_icon(c, x2 + 146 * mm, y_bot + 42 * mm, fill=colors.white)
    text(c, x2 + 125.5 * mm, y_bot + 37 * mm, "candidate.tsv", 6.3, True, center=True)
    text(c, x2 + 150.5 * mm, y_bot + 37 * mm, "QC + IGV\nmanifest", 6.3, True, center=True)
    c.setFillColor(colors.white)
    c.roundRect(x2 + 10 * mm, y_bot + 10 * mm, right_w - 20 * mm, 16 * mm, 3, fill=1, stroke=1)
    text(c, x2 + right_w / 2, y_bot + 20 * mm, "Candidate prioritization, not GWAS or causality by itself", 8, True, center=True)
    text(c, x2 + right_w / 2, y_bot + 13.5 * mm, "Requires PCR/amplicon, expression or splicing validation for causal claims", 6.7, COL["muted"], center=True)

    # Cross-panel arrows
    arrow(c, x1 + left_w * 0.50, y_top - 1 * mm, x1 + 55 * mm, y_bot + bottom_h - 18 * mm, width=1.2)
    arrow(c, x1 + 115 * mm, y_bot + bottom_h / 2, x1 + 118 * mm, y_bot + bottom_h / 2, width=1.2)
    arrow(c, x1 + 176 * mm, y_bot + bottom_h / 2, x2 - 1 * mm, y_bot + bottom_h - 28 * mm, width=1.2)
    arrow(c, x2 + right_w / 2, y_bot + bottom_h - 57 * mm, x2 + right_w / 2, y_bot + bottom_h - 68 * mm, width=1.2)

    c.showPage()
    c.save()


if __name__ == "__main__":
    main()
