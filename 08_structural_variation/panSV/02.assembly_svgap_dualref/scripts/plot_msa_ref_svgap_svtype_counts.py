#!/usr/bin/env python3
import argparse
import csv
import math
import os
from collections import Counter

FILES = {
    'DEL': 'All.DELs.50bplarge.bed.combined.sorted.txt',
    'INS': 'All.INTs.50bplarge.bed.combined.sorted.txt',
    'DEL_small': 'All.DELs.50bpsmall.bed.combined.sorted.txt',
    'INS_small': 'All.INTs.50bpsmall.bed.combined.sorted.txt',
    'CNV': 'All.CNVs.bed',
    'INV': 'All.INVs.bed',
    'TRL': 'All.TRLs.bed',
}

def data_line_count(path):
    n = 0
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            if line.strip():
                n += 1
    return n

def consecutive_unique_count(path):
    n = 0
    prev = None
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            if line != prev:
                n += 1
            prev = line
    return n

def consecutive_unique_by_col(path, col0):
    c = Counter()
    prev = None
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            if line == prev:
                prev = line
                continue
            parts = line.rstrip('\n').split('\t')
            val = parts[col0] if len(parts) > col0 else 'NA'
            c[val] += 1
            prev = line
    return c

def write_tsv(path, header, rows):
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(header)
        w.writerows(rows)

def format_count(n):
    return f'{n:,}'

def svg_bar(path, rows, title, subtitle, log_scale=False):
    # rows: list of (svtype, count)
    width, height = 920, 520
    left, right, top, bottom = 120, 175, 86, 72
    plot_w = width - left - right
    plot_h = height - top - bottom
    n = len(rows)
    gap = 18
    bar_h = min(46, (plot_h - gap*(n-1)) / n)
    max_v = max(v for _, v in rows) if rows else 1
    if log_scale:
        max_plot = math.log10(max_v + 1)
        xval = lambda v: math.log10(v + 1) / max_plot if max_plot else 0
        axis_label = 'log10(count + 1) scale; labels show exact counts'
    else:
        xval = lambda v: v / max_v if max_v else 0
        axis_label = 'linear count scale'
    colors = {
        'DEL': '#4E79A7',
        'INS': '#F28E2B',
        'INV': '#59A14F',
        'TRL': '#B07AA1',
        'CNV': '#E15759',
    }
    lines = []
    lines.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">')
    lines.append('<rect width="100%" height="100%" fill="white"/>')
    lines.append('<style>text{font-family:Arial, Helvetica, sans-serif; fill:#17202A;} .title{font-size:24px;font-weight:700;} .subtitle{font-size:14px;fill:#4B5563;} .axis{stroke:#111;stroke-width:1.6;} .grid{stroke:#E5E7EB;stroke-width:1;} .label{font-size:16px;font-weight:600;} .tick{font-size:12px;fill:#4B5563;} .value{font-size:14px;font-weight:600;}</style>')
    lines.append(f'<text class="title" x="{left}" y="36">{title}</text>')
    lines.append(f'<text class="subtitle" x="{left}" y="60">{subtitle}</text>')
    y0 = top + plot_h
    x0 = left
    # x-axis ticks
    if log_scale:
        ticks = [1, 10, 100, 1000, 10000, 100000, 1000000]
        ticks = [t for t in ticks if t <= max_v*1.08]
        for t in ticks:
            x = x0 + plot_w * xval(t)
            lines.append(f'<line class="grid" x1="{x:.1f}" y1="{top-8}" x2="{x:.1f}" y2="{y0}"/>')
            lines.append(f'<text class="tick" text-anchor="middle" x="{x:.1f}" y="{y0+24}">{format_count(t)}</text>')
    else:
        step = 100000
        ticks = list(range(0, int(math.ceil(max_v/step))*step + step, step))
        for t in ticks:
            if t > max_v*1.03: continue
            x = x0 + plot_w * (t / max_v if max_v else 0)
            lines.append(f'<line class="grid" x1="{x:.1f}" y1="{top-8}" x2="{x:.1f}" y2="{y0}"/>')
            lines.append(f'<text class="tick" text-anchor="middle" x="{x:.1f}" y="{y0+24}">{format_count(t)}</text>')
    lines.append(f'<line class="axis" x1="{x0}" y1="{y0}" x2="{x0+plot_w}" y2="{y0}"/>')
    lines.append(f'<line class="axis" x1="{x0}" y1="{top-8}" x2="{x0}" y2="{y0}"/>')
    for i, (sv, count) in enumerate(rows):
        y = top + i * (bar_h + gap)
        bw = max(1.5, plot_w * xval(count))
        fill = colors.get(sv, '#6B7280')
        lines.append(f'<text class="label" text-anchor="end" x="{left-16}" y="{y+bar_h/2+5:.1f}">{sv}</text>')
        lines.append(f'<rect x="{x0}" y="{y:.1f}" width="{bw:.1f}" height="{bar_h:.1f}" fill="{fill}"/>')
        value_x = x0 + bw + 10
        if value_x > width - right + 86:
            value_x = x0 + bw - 10
            anchor = 'end'
            color = 'white'
        else:
            anchor = 'start'
            color = '#17202A'
        lines.append(f'<text class="value" text-anchor="{anchor}" x="{value_x:.1f}" y="{y+bar_h/2+5:.1f}" fill="{color}">{format_count(count)}</text>')
    lines.append(f'<text class="subtitle" text-anchor="middle" x="{left+plot_w/2:.1f}" y="{height-18}">{axis_label}</text>')
    lines.append('</svg>')
    with open(path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--combined-dir', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--prefix', default='Msa_ref.SVGAP')
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    d = args.combined_dir

    counts = {
        'DEL': data_line_count(os.path.join(d, FILES['DEL'])),
        'INS': data_line_count(os.path.join(d, FILES['INS'])),
        'CNV': consecutive_unique_count(os.path.join(d, FILES['CNV'])),
        'INV': consecutive_unique_count(os.path.join(d, FILES['INV'])),
        'TRL': consecutive_unique_count(os.path.join(d, FILES['TRL'])),
    }
    small = {
        'DEL_small_lt50': data_line_count(os.path.join(d, FILES['DEL_small'])),
        'INS_small_lt50': data_line_count(os.path.join(d, FILES['INS_small'])),
    }
    cnv_sub = consecutive_unique_by_col(os.path.join(d, FILES['CNV']), 3)

    rows = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)
    count_rows = []
    source = {
        'DEL': FILES['DEL'],
        'INS': FILES['INS'],
        'CNV': FILES['CNV'],
        'INV': FILES['INV'],
        'TRL': FILES['TRL'],
    }
    rule = {
        'DEL': 'merged sorted rows excluding header; SVLEN >= 50 bp',
        'INS': 'merged sorted rows excluding header; SVLEN >= 50 bp',
        'CNV': 'consecutive duplicate rows collapsed from All.CNVs.bed',
        'INV': 'consecutive duplicate rows collapsed from All.INVs.bed',
        'TRL': 'consecutive duplicate rows collapsed from All.TRLs.bed',
    }
    total = sum(counts.values())
    for sv, c in rows:
        count_rows.append([sv, c, f'{100*c/total:.4f}', source[sv], rule[sv]])
    write_tsv(os.path.join(args.outdir, f'{args.prefix}.largeSV_type_counts.tsv'), ['SVTYPE','COUNT','PERCENT','SOURCE_FILE','COUNT_RULE'], count_rows)

    supp_rows = [
        ['DEL_large_ge50', counts['DEL'], FILES['DEL']],
        ['INS_large_ge50', counts['INS'], FILES['INS']],
        ['DEL_small_lt50', small['DEL_small_lt50'], FILES['DEL_small']],
        ['INS_small_lt50', small['INS_small_lt50'], FILES['INS_small']],
        ['DEL_total_with_lt50', counts['DEL'] + small['DEL_small_lt50'], 'large + small combined sorted files'],
        ['INS_total_with_lt50', counts['INS'] + small['INS_small_lt50'], 'large + small combined sorted files'],
    ]
    write_tsv(os.path.join(args.outdir, f'{args.prefix}.indel_size_supplement.tsv'), ['GROUP','COUNT','SOURCE'], supp_rows)
    write_tsv(os.path.join(args.outdir, f'{args.prefix}.CNV_subtype_counts.tsv'), ['CNV_SUBTYPE','COUNT'], sorted(cnv_sub.items(), key=lambda kv: kv[1], reverse=True))

    subtitle = f'genome_Msa reference, SVGAP Msa_ref CombinedSV, large SV only; total = {total:,}'
    svg_bar(os.path.join(args.outdir, f'{args.prefix}.largeSV_type_counts.sorted_bar.linear.svg'), rows, 'SVGAP genome_Msa SV counts by type', subtitle, log_scale=False)
    svg_bar(os.path.join(args.outdir, f'{args.prefix}.largeSV_type_counts.sorted_bar.log10.svg'), rows, 'SVGAP genome_Msa SV counts by type', subtitle, log_scale=True)
    with open(os.path.join(args.outdir, f'{args.prefix}.largeSV_type_counts.README.txt'), 'w') as fh:
        fh.write('Input directory: ' + d + '\n')
        fh.write('Main plot/table count only genome collinearity SVGAP Msa_ref CombinedSV results. Read-mapping results are not used.\n')
        fh.write('DEL and INS main counts use 50bplarge merged combined sorted files. The <50 bp DEL/INS counts are written only as a supplement.\n')
        fh.write('CNV/INV/TRL files contain adjacent duplicate rows; adjacent exact duplicates are collapsed before counting.\n')
        fh.write('Linear SVG preserves true magnitude; log10 SVG makes low-count CNV/INV/TRL visible.\n')
    print('DONE', args.outdir)
    print('SVTYPE\tCOUNT\tPERCENT')
    for sv, c in rows:
        print(f'{sv}\t{c}\t{100*c/total:.4f}')
    print('SUPPLEMENT')
    for r in supp_rows:
        print('\t'.join(map(str, r)))
    print('CNV_SUBTYPES')
    for k, v in sorted(cnv_sub.items(), key=lambda kv: kv[1], reverse=True):
        print(f'{k}\t{v}')

if __name__ == '__main__':
    main()
