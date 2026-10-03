#!/usr/bin/env python3
import argparse
import csv
import html
import os
import re
from collections import defaultdict, Counter

BIN = 100000
CATEGORIES4 = ["DNA", "LTR", "Other_TE", "Non_TE"]
CATEGORIES3 = ["DNA", "LTR", "Non_TE_or_other"]

DNA_PATTERNS = [
    "DNA", "MITE", "TIR", "Helitron", "DTA", "DTC", "DTH", "DTM", "DTT", "DTP", "DHH", "DXX",
    "Mutator", "hAT", "CACTA", "PIF", "Harbinger", "Tc1", "Mariner", "Pong", "Pogo", "CMC", "MULE",
]
LTR_PATTERNS = ["LTR", "Gypsy", "Copia"]
IGNORE_PATTERNS = ["Simple_repeat", "Low_complexity", "Satellite", "rRNA", "tRNA", "snRNA", "scRNA", "srpRNA"]

ATTR_RE = re.compile(r'(?:^|;)(?:Classification|classification)=([^;]+)')

def parse_args():
    ap = argparse.ArgumentParser(description="Classify Msa read-mapping PAVs by overlap with EDTA TE annotation and draw stacked percent bars.")
    ap.add_argument('--pav', required=True)
    ap.add_argument('--gff', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--prefix', default='Msa_read_mapping_only')
    ap.add_argument('--del-min-frac', type=float, default=0.5)
    ap.add_argument('--ins-window', type=int, default=50)
    return ap.parse_args()

def classify_te(feature, attrs):
    m = ATTR_RE.search(attrs)
    cls = m.group(1) if m else ''
    text = f'{feature};{cls};{attrs}'
    lower = text.lower()
    for pat in IGNORE_PATTERNS:
        if pat.lower() in lower:
            return None
    if any(pat.lower() in lower for pat in LTR_PATTERNS):
        return 'LTR'
    if any(pat.lower() in lower for pat in DNA_PATTERNS):
        return 'DNA'
    # LINE/SINE/non-LTR retrotransposons and unknown EDTA repeat classes are kept separate.
    if any(x in lower for x in ['line', 'sine', 'retrotransposon', 'repeat_region', 'repeat', 'helitron', 'mite', 'transposon']):
        return 'Other_TE'
    if cls and cls not in ('.', 'NA', 'Unknown'):
        return 'Other_TE'
    return None

def merge_intervals(intervals):
    intervals.sort()
    merged = []
    for s, e in intervals:
        if not merged or s > merged[-1][1] + 1:
            merged.append([s, e])
        elif e > merged[-1][1]:
            merged[-1][1] = e
    return [(s, e) for s, e in merged]

def load_te_index(gff):
    raw = defaultdict(lambda: defaultdict(list))
    total_records = 0
    used_records = 0
    with open(gff, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            total_records += 1
            chrom, feature, start_s, end_s, attrs = parts[0], parts[2], parts[3], parts[4], parts[8]
            cat = classify_te(feature, attrs)
            if not cat:
                continue
            try:
                start = int(start_s); end = int(end_s)
            except ValueError:
                continue
            if end < start:
                start, end = end, start
            raw[chrom][cat].append((start, end))
            used_records += 1
    merged = defaultdict(dict)
    index = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for chrom, bycat in raw.items():
        for cat, intervals in bycat.items():
            m = merge_intervals(intervals)
            merged[chrom][cat] = m
            for s, e in m:
                b0 = s // BIN
                b1 = e // BIN
                for b in range(b0, b1 + 1):
                    index[chrom][cat][b].append((s, e))
    stats = []
    for chrom in sorted(merged):
        for cat in sorted(merged[chrom]):
            bp = sum(e - s + 1 for s, e in merged[chrom][cat])
            stats.append((chrom, cat, len(merged[chrom][cat]), bp))
    return index, stats, total_records, used_records

def union_overlap_bp(intervals, qs, qe):
    ovs = []
    for s, e in intervals:
        if e < qs or s > qe:
            continue
        ovs.append((max(s, qs), min(e, qe)))
    if not ovs:
        return 0
    ovs.sort()
    total = 0
    cs, ce = ovs[0]
    for s, e in ovs[1:]:
        if s > ce + 1:
            total += ce - cs + 1
            cs, ce = s, e
        elif e > ce:
            ce = e
    total += ce - cs + 1
    return total

def query(index, chrom, qs, qe):
    overlaps = {cat: 0 for cat in ('DNA', 'LTR', 'Other_TE')}
    if chrom not in index:
        return overlaps
    b0 = qs // BIN
    b1 = qe // BIN
    for cat in ('DNA', 'LTR', 'Other_TE'):
        cand = []
        seen = set()
        bybin = index[chrom].get(cat, {})
        for b in range(b0, b1 + 1):
            for iv in bybin.get(b, []):
                if iv not in seen:
                    seen.add(iv)
                    cand.append(iv)
        overlaps[cat] = union_overlap_bp(cand, qs, qe)
    return overlaps

def assign_category(svtype, span_bp, overlaps, del_min_frac):
    best_cat, best_bp = max(overlaps.items(), key=lambda kv: kv[1])
    if svtype == 'DEL':
        frac = best_bp / span_bp if span_bp else 0.0
        if best_bp > 0 and frac >= del_min_frac:
            return best_cat, best_bp, frac
        return 'Non_TE', best_bp, frac
    # INS is a breakpoint event here; the window is only a proximity test.
    window_bp = span_bp if span_bp else 1
    frac = best_bp / window_bp if window_bp else 0.0
    if best_bp > 0:
        return best_cat, best_bp, frac
    return 'Non_TE', 0, 0.0

def write_tsv(path, rows, header):
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(header)
        w.writerows(rows)

def stacked_svg(path, percent_rows, categories, title, colors, legend_labels):
    width, height = 640, 520
    left, right, top, bottom = 92, 190, 40, 96
    plot_w = width - left - right
    plot_h = height - top - bottom
    bar_w = 88
    xs = {'DEL': left + plot_w * 0.30, 'INS': left + plot_w * 0.70}
    order = ['DEL', 'INS']
    label = {'DEL': 'Deletion', 'INS': 'Insertion'}
    data = {r['SVTYPE']: r for r in percent_rows}
    lines = []
    lines.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">')
    lines.append('<rect width="100%" height="100%" fill="white"/>')
    lines.append('<style>text{font-family:Arial, Helvetica, sans-serif; fill:#111;} .axis{stroke:#111;stroke-width:2;} .tick{stroke:#111;stroke-width:1.5;} .small{font-size:14px;} .label{font-size:20px;font-weight:600;} .legend{font-size:15px;} .xtick{font-size:17px;}</style>')
    y0 = top + plot_h
    x0 = left
    lines.append(f'<line class="axis" x1="{x0}" y1="{top}" x2="{x0}" y2="{y0}"/>')
    lines.append(f'<line class="axis" x1="{x0}" y1="{y0}" x2="{left+plot_w}" y2="{y0}"/>')
    for t in [0, 25, 50, 75, 100]:
        y = y0 - plot_h * t / 100.0
        lines.append(f'<line class="tick" x1="{x0-7}" y1="{y:.1f}" x2="{x0}" y2="{y:.1f}"/>')
        lines.append(f'<text class="small" text-anchor="end" x="{x0-12}" y="{y+5:.1f}">{t}</text>')
    lines.append(f'<text class="label" text-anchor="middle" transform="translate(28,{top+plot_h/2:.1f}) rotate(-90)">Percent (%)</text>')
    for sv in order:
        row = data.get(sv, {})
        x = xs[sv]
        curr = y0
        for cat in categories:
            pct = float(row.get(cat, 0.0))
            h = plot_h * pct / 100.0
            if h > 0:
                y = curr - h
                lines.append(f'<rect x="{x-bar_w/2:.1f}" y="{y:.1f}" width="{bar_w}" height="{h:.1f}" fill="{colors[cat]}"/>')
                curr = y
        lines.append(f'<text class="xtick" text-anchor="end" transform="translate({x+6:.1f},{y0+70}) rotate(-45)">{label[sv]}</text>')
    lx = left + plot_w + 42
    ly = top + 122
    lines.append(f'<text class="legend" x="{lx}" y="{ly-20}" font-weight="600">Type of PAVs</text>')
    for i, cat in enumerate(categories):
        yy = ly + i * 24
        lines.append(f'<rect x="{lx}" y="{yy-12}" width="14" height="14" fill="{colors[cat]}"/>')
        lines.append(f'<text class="legend" x="{lx+22}" y="{yy}">{html.escape(legend_labels.get(cat, cat))}</text>')
    if title:
        lines.append(f'<text x="{width-28}" y="30" font-size="28" font-family="Georgia, serif" font-weight="700" text-anchor="end">{html.escape(title)}</text>')
    lines.append('</svg>')
    with open(path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    idx, te_stats, total_gff, used_gff = load_te_index(args.gff)
    write_tsv(os.path.join(args.outdir, f'{args.prefix}.EDTA_TE_interval_stats.tsv'), te_stats, ['CHROM','TE_CATEGORY','MERGED_INTERVALS','MERGED_BP'])

    event_rows = []
    event_counts = {sv: Counter() for sv in ('DEL', 'INS')}
    presence_counts = {sv: Counter() for sv in ('DEL', 'INS')}
    n_rows = 0
    with open(args.pav, 'r', encoding='utf-8', errors='replace') as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader)
        sample_cols = header[6:]
        for row in reader:
            if len(row) < 6:
                continue
            svtype = row[4]
            if svtype not in ('DEL', 'INS'):
                continue
            n_rows += 1
            sv_id, chrom = row[0], row[1]
            try:
                pos = int(float(row[2])); end = int(float(row[3]))
            except ValueError:
                continue
            if end < pos:
                pos, end = end, pos
            if svtype == 'DEL':
                qs, qe = pos, end
                span = max(1, qe - qs + 1)
            else:
                qs = max(1, pos - args.ins_window)
                qe = pos + args.ins_window
                span = max(1, qe - qs + 1)
            overlaps = query(idx, chrom, qs, qe)
            cat, best_bp, frac = assign_category(svtype, span, overlaps, args.del_min_frac)
            presence = sum(1 for v in row[6:] if v == '1')
            event_counts[svtype][cat] += 1
            presence_counts[svtype][cat] += presence
            event_rows.append([
                sv_id, chrom, row[2], row[3], svtype, row[5], qs, qe, span, cat,
                overlaps['DNA'], overlaps['LTR'], overlaps['Other_TE'], best_bp, f'{frac:.6f}', presence
            ])
    assign_path = os.path.join(args.outdir, f'{args.prefix}.DEL_INS_TE_origin.event_assignments.tsv')
    write_tsv(assign_path, event_rows, ['SV_ID','CHROM','POS','END','SVTYPE','SVLEN','QUERY_START','QUERY_END','QUERY_BP','TE_CATEGORY','DNA_OVERLAP_BP','LTR_OVERLAP_BP','OTHER_TE_OVERLAP_BP','BEST_OVERLAP_BP','BEST_OVERLAP_FRACTION','PRESENCE_COUNT'])

    count_rows4 = []
    pct_rows4 = []
    count_rows3 = []
    pct_rows3 = []
    presence_rows4 = []
    for sv in ('DEL', 'INS'):
        total = sum(event_counts[sv].values())
        row4 = {'SVTYPE': sv}
        rowp4 = {'SVTYPE': sv}
        row3 = {'SVTYPE': sv}
        rowp3 = {'SVTYPE': sv}
        for cat in CATEGORIES4:
            row4[cat] = event_counts[sv][cat]
            rowp4[cat] = (100.0 * event_counts[sv][cat] / total) if total else 0.0
        row3['DNA'] = event_counts[sv]['DNA']
        row3['LTR'] = event_counts[sv]['LTR']
        row3['Non_TE_or_other'] = event_counts[sv]['Non_TE'] + event_counts[sv]['Other_TE']
        for cat in CATEGORIES3:
            rowp3[cat] = (100.0 * row3[cat] / total) if total else 0.0
        count_rows4.append([sv] + [row4[c] for c in CATEGORIES4] + [total])
        pct_rows4.append(rowp4)
        count_rows3.append([sv] + [row3[c] for c in CATEGORIES3] + [total])
        pct_rows3.append(rowp3)
        pres_total = sum(presence_counts[sv].values())
        presence_rows4.append([sv] + [presence_counts[sv][c] for c in CATEGORIES4] + [pres_total])

    write_tsv(os.path.join(args.outdir, f'{args.prefix}.DEL_INS_TE_origin.event_counts.4class.tsv'), count_rows4, ['SVTYPE'] + CATEGORIES4 + ['TOTAL'])
    write_tsv(os.path.join(args.outdir, f'{args.prefix}.DEL_INS_TE_origin.event_counts.3class.tsv'), count_rows3, ['SVTYPE'] + CATEGORIES3 + ['TOTAL'])
    write_tsv(os.path.join(args.outdir, f'{args.prefix}.DEL_INS_TE_origin.presence_weighted_counts.4class.tsv'), presence_rows4, ['SVTYPE'] + CATEGORIES4 + ['TOTAL_PRESENCE'])
    # percent TSVs
    pct4_rows = []
    pct3_rows = []
    for r in pct_rows4:
        pct4_rows.append([r['SVTYPE']] + [f'{r[c]:.4f}' for c in CATEGORIES4])
    for r in pct_rows3:
        pct3_rows.append([r['SVTYPE']] + [f'{r[c]:.4f}' for c in CATEGORIES3])
    write_tsv(os.path.join(args.outdir, f'{args.prefix}.DEL_INS_TE_origin.event_percent.4class.tsv'), pct4_rows, ['SVTYPE'] + CATEGORIES4)
    write_tsv(os.path.join(args.outdir, f'{args.prefix}.DEL_INS_TE_origin.event_percent.3class.tsv'), pct3_rows, ['SVTYPE'] + CATEGORIES3)

    colors4 = {'DNA':'#E59A3A', 'LTR':'#67C7D8', 'Other_TE':'#8E8E8E', 'Non_TE':'#33A08B'}
    colors3 = {'DNA':'#E59A3A', 'LTR':'#67C7D8', 'Non_TE_or_other':'#33A08B'}
    stacked_svg(os.path.join(args.outdir, f'{args.prefix}.DEL_INS_TE_origin.stacked_percent_bar.4class.svg'), pct_rows4, CATEGORIES4, '', colors4, {'DNA':'DNA', 'LTR':'LTR', 'Other_TE':'Other TE', 'Non_TE':'Non-TE'})
    stacked_svg(os.path.join(args.outdir, f'{args.prefix}.DEL_INS_TE_origin.stacked_percent_bar.3class.svg'), pct_rows3, CATEGORIES3, '', colors3, {'DNA':'DNA', 'LTR':'LTR', 'Non_TE_or_other':'Non-TE/other'})

    with open(os.path.join(args.outdir, f'{args.prefix}.DEL_INS_TE_origin.README.txt'), 'w') as fh:
        fh.write('Input PAV: ' + args.pav + '\n')
        fh.write('Input EDTA GFF: ' + args.gff + '\n')
        fh.write(f'GFF records scanned: {total_gff}; TE-like records indexed: {used_gff}\n')
        fh.write(f'PAV rows scanned for DEL/INS: {n_rows}\n')
        fh.write(f'DEL rule: TE if max merged TE-category overlap fraction >= {args.del_min_frac}\n')
        fh.write(f'INS rule: breakpoint +/- {args.ins_window} bp overlaps merged TE-category interval; this is breakpoint-proximity, not inserted-sequence RepeatMasker classification.\n')
        fh.write('Other_TE includes EDTA LINE/SINE/non-LTR/unknown repeat classes when classified as TE-like; simple_repeat/low_complexity/satellite are not treated as TE-caused.\n')
        fh.write('Primary plot files: stacked_percent_bar.4class.svg (scientific full classes) and stacked_percent_bar.3class.svg (DNA/LTR/non-DNA-LTR style).\n')
    print('DONE')
    print('outdir', args.outdir)
    print('event_counts_4class')
    for row in count_rows4:
        print('\t'.join(map(str, row)))
    print('event_percent_4class')
    for row in pct4_rows:
        print('\t'.join(map(str, row)))
    print('event_counts_3class')
    for row in count_rows3:
        print('\t'.join(map(str, row)))
    print('event_percent_3class')
    for row in pct3_rows:
        print('\t'.join(map(str, row)))

if __name__ == '__main__':
    main()
