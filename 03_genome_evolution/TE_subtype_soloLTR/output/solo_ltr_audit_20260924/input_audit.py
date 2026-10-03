#!/usr/bin/env python3
"""Audit the inputs behind the 39.TE_type_soloLTR summary."""

import csv
from collections import Counter
from pathlib import Path

BASE = Path('path/to/project/39.TE_type_soloLTR')
SUMMARY = BASE / 'output/solo_ltr/solo_full_ltr_summary.tsv'
OUT = BASE / 'output/solo_ltr_audit_20260924/input_audit.tsv'


def features(path):
    counts = Counter()
    if path is None or not path.is_file():
        return counts
    with path.open(encoding='utf-8', errors='replace') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) >= 9:
                counts[fields[2]] += 1
    return counts


def find_raw_ltr(species_dir):
    paths = sorted(species_dir.glob('*.EDTA.raw/LTR/*.LTR.intact.raw.gff3'))
    if not paths:
        paths = sorted(species_dir.glob('*.EDTA.raw/LTR/*.LTR.intact.gff3'))
    return paths[0] if paths else None


rows = []
with SUMMARY.open(newline='') as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        species = row['species']
        species_dir = BASE / 'data' / (species + '_EDTA')
        final = Path(row['intact_gff3'])
        raw = find_raw_ltr(species_dir)
        final_counts = features(final)
        raw_counts = features(raw)
        final_ltr = sum(value for key, value in final_counts.items() if 'LTR_retrotransposon' in key)
        raw_ltr = sum(value for key, value in raw_counts.items() if 'LTR_retrotransposon' in key)
        rows.append({
            'species': species,
            'reported_full': row['full_length_ltr_rt_count'],
            'reported_solo': row['solo_ltr_count'],
            'final_ltr_features': final_ltr,
            'raw_ltr_features': raw_ltr,
            'final_raw_ratio': f'{final_ltr/raw_ltr:.4f}' if raw_ltr else 'NA',
            'final_path': str(final),
            'raw_path': str(raw) if raw else 'MISSING',
            'raw_ltr_group': raw_counts['LTR_retrotransposon'],
            'raw_copia': raw_counts['Copia_LTR_retrotransposon'],
            'raw_gypsy': raw_counts['Gypsy_LTR_retrotransposon'],
            'merged_ltr_hits': row['merged_ltr_hits'],
            'excluded_near_int': row['excluded_near_ltr_internal'],
        })

OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter='\t')
    writer.writeheader()
    writer.writerows(rows)
print(OUT)
