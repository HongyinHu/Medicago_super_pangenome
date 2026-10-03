#!/usr/bin/env python3
import sys, os, math, glob
from pathlib import Path


def read_hist(path):
    depths = []
    counts = []
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            try:
                d = int(float(parts[0]))
                c = int(float(parts[1]))
            except ValueError:
                continue
            if d >= 0 and c >= 0:
                depths.append(d)
                counts.append(c)
    return depths, counts


def smooth(vals, radius=2):
    out = []
    n = len(vals)
    for i in range(n):
        lo = max(0, i - radius)
        hi = min(n, i + radius + 1)
        out.append(sum(vals[lo:hi]) / (hi - lo))
    return out


def choose_peak(depths, counts, min_depth=5):
    pairs = [(d, c) for d, c in zip(depths, counts) if d >= min_depth and c > 0]
    if not pairs:
        return None, [], 'no_peak_after_min_depth'
    ds = [p[0] for p in pairs]
    cs = [p[1] for p in pairs]
    sm = smooth(cs, radius=2)
    candidates = []
    for i, (d, c) in enumerate(zip(ds, sm)):
        left = sm[i - 1] if i else -1
        right = sm[i + 1] if i + 1 < len(sm) else -1
        if c >= left and c >= right:
            candidates.append((d, c, cs[i]))
    if not candidates:
        idx = max(range(len(cs)), key=lambda i: cs[i])
        candidates = [(ds[idx], sm[idx], cs[idx])]
    candidates = sorted(candidates, key=lambda x: x[1], reverse=True)
    top = candidates[0]
    peak = top[0]
    note = 'highest_smoothed_peak'

    # In heterozygous genomes the half-depth heterozygous peak can be taller.
    # If a strong peak exists near 2x the tallest peak, use it as the homozygous coverage candidate.
    for d, smc, rawc in candidates[1:12]:
        if 1.75 * top[0] <= d <= 2.25 * top[0] and smc >= 0.08 * top[1]:
            peak = d
            note = f'higher_depth_peak_near_2x_top; top_peak={top[0]}'
            break
    compact = ';'.join([f'{d}:{int(rawc)}' for d, smc, rawc in candidates[:8]])
    return peak, compact, note


def estimate_one(hist_path):
    depths, counts = read_hist(hist_path)
    if not depths:
        return None
    total_all = sum(d * c for d, c in zip(depths, counts))
    total_ge2 = sum(d * c for d, c in zip(depths, counts) if d >= 2)
    total_ge3 = sum(d * c for d, c in zip(depths, counts) if d >= 3)
    distinct_all = sum(counts)
    distinct_ge2 = sum(c for d, c in zip(depths, counts) if d >= 2)
    peak, peaks, note = choose_peak(depths, counts)
    if not peak:
        peak = 0
    def div(x):
        return '' if peak == 0 else f'{x / peak:.0f}'
    p = Path(hist_path)
    name = p.name
    # Expected: sample.k21.histo
    sample = name
    k = ''
    parts = name.split('.')
    if len(parts) >= 3 and parts[-1] == 'histo':
        sample = '.'.join(parts[:-2])
        k = parts[-2].lstrip('k')
    return {
        'sample': sample,
        'k': k,
        'histogram': str(p),
        'peak_depth': str(peak),
        'genome_size_all_bp': div(total_all),
        'genome_size_ge2_bp': div(total_ge2),
        'genome_size_ge3_bp': div(total_ge3),
        'total_kmers_all': str(total_all),
        'total_kmers_ge2': str(total_ge2),
        'total_kmers_ge3': str(total_ge3),
        'distinct_kmers_all': str(distinct_all),
        'distinct_kmers_ge2': str(distinct_ge2),
        'candidate_peaks_depth:count': peaks if isinstance(peaks, str) else '',
        'peak_rule': note,
    }


def main():
    if len(sys.argv) < 3:
        print('Usage: 02_estimate_genome_size.py OUT.tsv HISTO...', file=sys.stderr)
        sys.exit(2)
    out = sys.argv[1]
    histos = []
    for arg in sys.argv[2:]:
        histos.extend(glob.glob(arg))
    rows = []
    for h in sorted(set(histos)):
        row = estimate_one(h)
        if row:
            rows.append(row)
    fields = ['sample','k','peak_depth','genome_size_all_bp','genome_size_ge2_bp','genome_size_ge3_bp','total_kmers_all','total_kmers_ge2','total_kmers_ge3','distinct_kmers_all','distinct_kmers_ge2','candidate_peaks_depth:count','peak_rule','histogram']
    with open(out, 'w', encoding='utf-8') as fh:
        fh.write('\t'.join(fields) + '\n')
        for r in rows:
            fh.write('\t'.join(r.get(f, '') for f in fields) + '\n')

if __name__ == '__main__':
    main()
