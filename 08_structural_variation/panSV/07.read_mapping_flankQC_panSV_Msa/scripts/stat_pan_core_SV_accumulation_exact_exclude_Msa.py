#!/usr/bin/env python3
from pathlib import Path
import csv, argparse
from collections import Counter

ap = argparse.ArgumentParser()
ap.add_argument('--pav', required=True, type=Path)
ap.add_argument('--exclude', default='genome_Msa')
ap.add_argument('--out-prefix', required=True, type=Path)
args = ap.parse_args()

def comb(n, k):
    if k < 0 or k > n:
        return 0
    k = min(k, n-k)
    num = 1
    den = 1
    for i in range(1, k+1):
        num *= n - k + i
        den *= i
    return num // den

with args.pav.open(newline='') as fh:
    reader = csv.DictReader(fh, delimiter='\t')
    genome_cols = [c for c in reader.fieldnames if c.startswith('genome_')]
    keep_cols = [c for c in genome_cols if c != args.exclude]
    if args.exclude not in genome_cols:
        raise SystemExit('exclude sample not found: ' + args.exclude)
    n = len(keep_cols)
    freq_counter = Counter()
    total_rows = kept_rows = 0
    for row in reader:
        total_rows += 1
        s = 0
        for c in keep_cols:
            if row.get(c, '0') not in ('0', '.', '', 'NA'):
                s += 1
        if s == 0:
            continue
        kept_rows += 1
        freq_counter[s] += 1

rows = []
for k in range(1, n + 1):
    denom = float(comb(n, k))
    pan = 0.0
    core = 0.0
    for s, count in freq_counter.items():
        p_pan = 1.0
        if n - s >= k:
            p_pan -= comb(n - s, k) / denom
        p_core = comb(s, k) / denom if s >= k else 0.0
        pan += count * p_pan
        core += count * p_core
    rows.append({
        'sample_number': k,
        'pan_SVs_expected': pan,
        'core_SVs_expected': core,
        'species_count': n,
        'event_count': kept_rows,
    })

args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
out_tsv = args.out_prefix.with_suffix('.tsv')
with out_tsv.open('w', newline='') as fh:
    fields = ['sample_number','pan_SVs_expected','core_SVs_expected','species_count','event_count']
    w = csv.DictWriter(fh, delimiter='\t', fieldnames=fields)
    w.writeheader()
    for r in rows:
        rr = r.copy()
        rr['pan_SVs_expected'] = '%.3f' % r['pan_SVs_expected']
        rr['core_SVs_expected'] = '%.3f' % r['core_SVs_expected']
        w.writerow(rr)

out_freq = args.out_prefix.with_suffix('.frequency.tsv')
with out_freq.open('w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['presence_frequency','event_count'])
    for s in range(1, n + 1):
        w.writerow([s, freq_counter[s]])

out_samples = args.out_prefix.with_suffix('.samples.tsv')
with out_samples.open('w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['index','sample'])
    for i, c in enumerate(keep_cols, start=1):
        w.writerow([i,c])

print('WROTE', out_tsv)
print('WROTE', out_freq)
print('WROTE', out_samples)
print('SPECIES', n, ','.join(keep_cols))
print('EVENTS_TOTAL', total_rows)
print('EVENTS_17SPECIES', kept_rows)
print('FINAL_K', rows[-1])