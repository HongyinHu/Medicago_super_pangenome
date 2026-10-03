#!/usr/bin/env python3
from pathlib import Path
import csv, random, argparse
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument('--pav', required=True, type=Path)
ap.add_argument('--exclude', default='genome_Msa')
ap.add_argument('--repeats', type=int, default=1000)
ap.add_argument('--seed', type=int, default=20260702)
ap.add_argument('--out-prefix', required=True, type=Path)
args = ap.parse_args()

svtypes = ['DEL','INS','TRA','DUP','INV']
# Load PAV as compact bitmask over non-excluded genomes. Keep only events present in >=1 remaining genome.
with args.pav.open(newline='') as fh:
    reader = csv.DictReader(fh, delimiter='\t')
    genome_cols = [c for c in reader.fieldnames if c.startswith('genome_')]
    keep_cols = [c for c in genome_cols if c != args.exclude]
    if args.exclude not in genome_cols:
        raise SystemExit(f'exclude sample not found: {args.exclude}')
    n = len(keep_cols)
    masks = []
    type_masks = {sv: [] for sv in svtypes}
    for row in reader:
        mask = 0
        for i, c in enumerate(keep_cols):
            if row.get(c, '0') not in ('0', '.', '', 'NA'):
                mask |= (1 << i)
        if mask == 0:
            continue
        masks.append(mask)
        sv = row.get('SVTYPE')
        if sv in type_masks:
            type_masks[sv].append(mask)

rng = random.Random(args.seed)
all_bits = (1 << n) - 1
rows = []
for k in range(1, n + 1):
    pan_vals = []
    core_vals = []
    for _ in range(args.repeats):
        idx = rng.sample(range(n), k)
        subset = 0
        for i in idx:
            subset |= (1 << i)
        pan = 0
        core = 0
        for m in masks:
            if m & subset:
                pan += 1
            if (m & subset) == subset:
                core += 1
        pan_vals.append(pan)
        core_vals.append(core)
    rows.append({
        'sample_number': k,
        'pan_SVs_mean': sum(pan_vals) / args.repeats,
        'core_SVs_mean': sum(core_vals) / args.repeats,
        'pan_SVs_min': min(pan_vals),
        'pan_SVs_max': max(pan_vals),
        'core_SVs_min': min(core_vals),
        'core_SVs_max': max(core_vals),
        'repeats': args.repeats,
        'species_count': n,
        'event_count': len(masks),
    })

args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
out_tsv = args.out_prefix.with_suffix('.tsv')
with out_tsv.open('w', newline='') as fh:
    fields = ['sample_number','pan_SVs_mean','core_SVs_mean','pan_SVs_min','pan_SVs_max','core_SVs_min','core_SVs_max','repeats','species_count','event_count']
    w = csv.DictWriter(fh, delimiter='\t', fieldnames=fields)
    w.writeheader()
    for r in rows:
        rr = r.copy()
        rr['pan_SVs_mean'] = f"{r['pan_SVs_mean']:.3f}"
        rr['core_SVs_mean'] = f"{r['core_SVs_mean']:.3f}"
        w.writerow(rr)

out_samples = args.out_prefix.with_suffix('.samples.tsv')
with out_samples.open('w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['index','sample'])
    for i, c in enumerate(keep_cols, start=1):
        w.writerow([i,c])

print('WROTE', out_tsv)
print('WROTE', out_samples)
print('SPECIES', n, ','.join(keep_cols))
print('EVENTS', len(masks))
print('FINAL_K', rows[-1])