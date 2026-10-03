#!/usr/bin/env python3
from pathlib import Path
import csv
from collections import defaultdict

base = Path('path/to/project/N_3.call_SV/07.read_mapping_flankQC_panSV_Msa')
pav = base / 'results/Msa/panSV/panSV.Msa.flankQC.read_primary.PAV.matrix.tsv'
out_dir = base / 'results/Msa/tables'
out_dir.mkdir(parents=True, exist_ok=True)
exclude = 'genome_Msa'
svtypes = ['DEL', 'INS']

bins = []
# labels 1..10 mean <=1kb, 1-2kb, ..., 9-10kb; final label >10.
for i in range(1, 11):
    lo = 0 if i == 1 else (i - 1) * 1000
    hi = i * 1000
    label = str(i)
    display = f'<= {i} kb' if i == 1 else f'{i-1}-{i} kb'
    bins.append((label, display, lo, hi))
bins.append(('>10', '>10 kb', 10000, None))

def parse_len(row):
    svlen = row.get('SVLEN', '')
    try:
        return abs(int(float(svlen)))
    except Exception:
        pass
    try:
        pos = int(row.get('POS', 0))
        end = int(row.get('END', 0))
        return abs(end - pos)
    except Exception:
        return None

def bin_len(length):
    for label, display, lo, hi in bins:
        if hi is None:
            if length > lo:
                return label
        else:
            if length <= hi and length > lo:
                return label
    return None

counts = {sv: defaultdict(int) for sv in svtypes}
length_sum = defaultdict(int)
len_values = {sv: [] for sv in svtypes}
rows_total = rows_kept = 0
with pav.open(newline='') as fh:
    reader = csv.DictReader(fh, delimiter='\t')
    genome_cols = [c for c in reader.fieldnames if c.startswith('genome_')]
    keep_cols = [c for c in genome_cols if c != exclude]
    for row in reader:
        rows_total += 1
        sv = row.get('SVTYPE')
        if sv not in svtypes:
            continue
        # Exclude genome_Msa self; keep non-redundant event if any remaining 17 species carries it.
        present17 = any(row.get(c, '0') not in ('0', '.', '', 'NA') for c in keep_cols)
        if not present17:
            continue
        length = parse_len(row)
        if length is None or length < 50:
            continue
        b = bin_len(length)
        if b is None:
            continue
        counts[sv][b] += 1
        length_sum[sv] += length
        len_values[sv].append(length)
        rows_kept += 1

out_tsv = out_dir / 'panSV_DEL_INS_length_distribution.exclude_genome_Msa.tsv'
with out_tsv.open('w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['SVTYPE', 'bin_label', 'bin_display', 'length_min_bp_exclusive', 'length_max_bp_inclusive', 'event_count'])
    for sv in svtypes:
        for label, display, lo, hi in bins:
            w.writerow([sv, label, display, lo, 'Inf' if hi is None else hi, counts[sv][label]])

out_summary = out_dir / 'panSV_DEL_INS_length_distribution.exclude_genome_Msa.summary.tsv'
with out_summary.open('w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['SVTYPE','total_events','mean_bp','median_bp','min_bp','max_bp','le_1kb','gt_10kb'])
    for sv in svtypes:
        vals = sorted(len_values[sv])
        n = len(vals)
        if n:
            med = vals[n//2] if n % 2 else (vals[n//2 - 1] + vals[n//2]) / 2
            row = [sv, n, round(sum(vals)/n, 2), med, vals[0], vals[-1], counts[sv]['1'], counts[sv]['>10']]
        else:
            row = [sv, 0, 'NA', 'NA', 'NA', 'NA', 0, 0]
        w.writerow(row)

print('WROTE', out_tsv)
print('WROTE', out_summary)
for sv in svtypes:
    print('COUNTS', sv, '\t'.join(f'{label}={counts[sv][label]}' for label, _, _, _ in bins))
print('ROWS_SCANNED', rows_total)
print('DEL_INS_EVENTS_COUNTED', rows_kept)