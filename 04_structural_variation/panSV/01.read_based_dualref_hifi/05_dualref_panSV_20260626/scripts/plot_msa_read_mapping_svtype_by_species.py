#!/usr/bin/env python3
import csv
from pathlib import Path
from collections import Counter, defaultdict

STEP = Path('path/to/project/N_3.call_SV/01.read_based_dualref_hifi/05_dualref_panSV_20260626')
PAV = STEP / 'results/RefA/panSV/panSV.RefA.PAV.matrix.tsv'
PLOT_DIR = STEP / 'results/RefA/panSV/plots_read_mapping_only'
PLOT_DIR.mkdir(parents=True, exist_ok=True)

LONG_TSV = PLOT_DIR / 'Msa_read_mapping_only.SVTYPE_count_by_species.long.tsv'
WIDE_TSV = PLOT_DIR / 'Msa_read_mapping_only.SVTYPE_count_by_species.wide.tsv'
TOTAL_TSV = PLOT_DIR / 'Msa_read_mapping_only.SV_count_by_species.total.tsv'
SVG = PLOT_DIR / 'Msa_read_mapping_only.SVTYPE_count_by_species.stacked_bar.svg'
README = PLOT_DIR / 'Msa_read_mapping_only.SVTYPE_count_by_species.README.txt'

TYPE_ORDER = ['DEL', 'INS', 'DUP', 'INV', 'CNV', 'TRL']
COLORS = {
    'DEL': '#4C78A8',
    'INS': '#F58518',
    'DUP': '#54A24B',
    'INV': '#B279A2',
    'CNV': '#E45756',
    'TRL': '#72B7B2',
}

def esc(s):
    return str(s).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

def fmt_k(v):
    if v >= 1000000:
        return f'{v/1000000:.1f}M'
    if v >= 1000:
        return f'{v/1000:.0f}k'
    return str(v)

counts = defaultdict(Counter)
totals = Counter()
svtypes_seen = set()
rows = 0
with open(PAV, newline='', encoding='utf-8', errors='replace') as fh:
    reader = csv.DictReader(fh, delimiter='\t')
    samples = reader.fieldnames[6:]
    for row in reader:
        rows += 1
        svtype = row.get('SVTYPE', '.') or '.'
        svtypes_seen.add(svtype)
        for sample in samples:
            if row.get(sample) == '1':
                counts[sample][svtype] += 1
                totals[sample] += 1

svtypes = [x for x in TYPE_ORDER if x in svtypes_seen] + sorted(x for x in svtypes_seen if x not in TYPE_ORDER)
samples = sorted(samples, key=lambda sample: (-totals[sample], sample))

with open(LONG_TSV, 'w', newline='', encoding='utf-8') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['species', 'SVTYPE', 'present_count'])
    for sample in samples:
        for svtype in svtypes:
            w.writerow([sample, svtype, counts[sample][svtype]])

with open(WIDE_TSV, 'w', newline='', encoding='utf-8') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['species'] + svtypes + ['TOTAL'])
    for sample in samples:
        w.writerow([sample] + [counts[sample][svtype] for svtype in svtypes] + [totals[sample]])

with open(TOTAL_TSV, 'w', newline='', encoding='utf-8') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['species', 'total_present_read_mapping_SV'])
    for sample in samples:
        w.writerow([sample, totals[sample]])

width = 1600
height = 900
left = 120
right = 260
top = 70
bottom = 220
plot_w = width - left - right
plot_h = height - top - bottom
max_total = max(totals.values()) if totals else 1
step = 50000
if max_total > 500000:
    step = 100000
elif max_total < 100000:
    step = 20000
ymax = ((max_total + step - 1) // step) * step
bar_gap = 14
bar_w = max(18, (plot_w - bar_gap * (len(samples) - 1)) / max(1, len(samples)))

def y_pos(value):
    return top + plot_h - (value / ymax) * plot_h

parts = []
parts.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">')
parts.append('<rect width="100%" height="100%" fill="white"/>')
parts.append(f'<text x="{width/2}" y="32" text-anchor="middle" font-family="Arial" font-size="24" font-weight="700">genome_Msa read-mapping SV counts by species and SV type</text>')
parts.append(f'<text x="{width/2}" y="58" text-anchor="middle" font-family="Arial" font-size="14" fill="#555">Input: RefA/Msa read-based panSV PAV; presence counted from PAV value = 1</text>')
parts.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="#222" stroke-width="1.2"/>')
parts.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="#222" stroke-width="1.2"/>')
for tick in range(0, ymax + 1, step):
    y = y_pos(tick)
    parts.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left+plot_w}" y2="{y:.1f}" stroke="#e6e6e6" stroke-width="1"/>')
    parts.append(f'<text x="{left-10}" y="{y+4:.1f}" text-anchor="end" font-family="Arial" font-size="12" fill="#444">{fmt_k(tick)}</text>')
parts.append(f'<text x="{left-80}" y="{top+plot_h/2}" transform="rotate(-90 {left-80},{top+plot_h/2})" text-anchor="middle" font-family="Arial" font-size="15" fill="#333">Number of present read-mapping SVs</text>')
for i, sample in enumerate(samples):
    x = left + i * (bar_w + bar_gap)
    cumulative = 0
    for svtype in svtypes:
        val = counts[sample][svtype]
        if val <= 0:
            continue
        y1 = y_pos(cumulative + val)
        y0 = y_pos(cumulative)
        h = max(0, y0 - y1)
        parts.append(f'<rect x="{x:.1f}" y="{y1:.1f}" width="{bar_w:.1f}" height="{h:.1f}" fill="{COLORS.get(svtype, "#999999")}"/>')
        cumulative += val
    parts.append(f'<text x="{x+bar_w/2:.1f}" y="{y_pos(totals[sample])-6:.1f}" text-anchor="middle" font-family="Arial" font-size="10" fill="#222">{fmt_k(totals[sample])}</text>')
    label_x = x + bar_w/2
    label_y = top + plot_h + 18
    parts.append(f'<text x="{label_x:.1f}" y="{label_y:.1f}" transform="rotate(45 {label_x:.1f},{label_y:.1f})" text-anchor="start" font-family="Arial" font-size="12" fill="#333">{esc(sample)}</text>')
legend_x = left + plot_w + 35
legend_y = top + 20
parts.append(f'<text x="{legend_x}" y="{legend_y-10}" font-family="Arial" font-size="16" font-weight="700">SVTYPE</text>')
for j, svtype in enumerate(svtypes):
    y = legend_y + j * 28
    parts.append(f'<rect x="{legend_x}" y="{y}" width="18" height="18" fill="{COLORS.get(svtype, "#999999")}"/>')
    parts.append(f'<text x="{legend_x+28}" y="{y+14}" font-family="Arial" font-size="14" fill="#333">{esc(svtype)}</text>')
parts.append(f'<text x="{left}" y="{height-22}" font-family="Arial" font-size="12" fill="#666">Input: {esc(str(PAV))}</text>')
parts.append('</svg>')
SVG.write_text('\n'.join(parts) + '\n', encoding='utf-8')

README.write_text(
    'Read-mapping-only SVTYPE counts by species for genome_Msa single-reference panSV\n'
    f'Input PAV: {PAV}\n'
    'Counting rule: for each species/sample column, count records with PAV value exactly 1; values 0 and NA are not counted.\n'
    'This excludes all SVGAP-only events and uses only the read-based RefA/Msa panSV matrix.\n'
    f'Total read-mapping SV records scanned: {rows}\n'
    f'Wide table: {WIDE_TSV}\n'
    f'Long table: {LONG_TSV}\n'
    f'Total table: {TOTAL_TSV}\n'
    f'Stacked bar SVG: {SVG}\n',
    encoding='utf-8'
)
print('rows_scanned', rows)
print('samples', len(samples), ','.join(samples))
print('svtypes', ','.join(svtypes))
print('wide_tsv', WIDE_TSV)
print('svg', SVG)
