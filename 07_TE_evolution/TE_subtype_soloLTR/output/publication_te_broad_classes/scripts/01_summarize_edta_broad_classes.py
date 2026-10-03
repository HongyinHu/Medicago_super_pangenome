#!/usr/bin/env python3
import csv
import re
from pathlib import Path
from urllib.parse import unquote
from collections import defaultdict, Counter

BASE = Path('path/to/project/39.TE_type_soloLTR')
OLD = BASE / 'output' / 'species_totals.tsv'
OUT = BASE / 'output' / 'publication_te_broad_classes'
TABLE = OUT / 'tables' / 'edta_broad_class_percent.tsv'
WIDE = OUT / 'tables' / 'edta_broad_class_percent.wide.tsv'
CLASS_RAW = OUT / 'tables' / 'edta_raw_classification_counts.tsv'

CATEGORY_ORDER = [
    'LTR/Copia',
    'LTR/Gypsy',
    'LTR/unknown',
    'DNA transposons',
    'Helitron',
    'Other / unclassified',
]

# Stable order from the existing analysis; this keeps the figure comparable to earlier outputs.
PREFERRED_ORDER = [
    'genome_395','genome_436','genome_457','genome_468','genome_474','genome_482',
    'genome_M22','genome_Mar','genome_Mru','genome_Msa','genome_ZM4',
    'genome_410','genome_454','genome_461','genome_472','genome_474_T2T',
    'genome_A17','genome_M46','genome_MPO','genome_Msa_T2T','genome_R108'
]

LABELS = {
    'genome_395': 'genome_395',
    'genome_436': 'genome_436',
    'genome_457': 'genome_457',
    'genome_468': 'genome_468',
    'genome_474': 'genome_474',
    'genome_482': 'genome_482',
    'genome_M22': 'genome_M22',
    'genome_Mar': 'M. archiducis-nicolai',
    'genome_Mru': 'M. ruthenica',
    'genome_Msa': 'M. sativa',
    'genome_ZM4': 'M. sativa cv. Zhongmu-4',
    'genome_410': 'genome_410',
    'genome_454': 'genome_454',
    'genome_461': 'genome_461',
    'genome_472': 'genome_472',
    'genome_474_T2T': 'genome_474_T2T',
    'genome_A17': 'M. truncatula A17',
    'genome_M46': 'genome_M46',
    'genome_MPO': 'M. polymorpha',
    'genome_Msa_T2T': 'M. sativa T2T',
    'genome_R108': 'M. truncatula R108',
}

def parse_attrs(text):
    attrs = {}
    for part in text.strip().split(';'):
        if not part:
            continue
        if '=' in part:
            k, v = part.split('=', 1)
            attrs[k] = unquote(v)
    return attrs

def broad_category(feature, classification):
    c = (classification or '').strip()
    f = (feature or '').strip()
    low = f'{c} {f}'.lower()
    if c == 'LTR/Copia' or 'copia_ltr_retrotransposon' in low:
        return 'LTR/Copia'
    if c == 'LTR/Gypsy' or 'gypsy_ltr_retrotransposon' in low:
        return 'LTR/Gypsy'
    if c == 'LTR/unknown' or (c.startswith('LTR/') and c.lower().endswith('/unknown')) or f == 'LTR_retrotransposon':
        return 'LTR/unknown'
    if 'helitron' in low:
        return 'Helitron'
    if c.startswith('DNA/') or c.startswith('MITE/') or 'tir_transposon' in low or 'terminal_inverted_repeat' in low:
        return 'DNA transposons'
    return 'Other / unclassified'

def add_interval(store, category, chrom, start, end):
    if end < start:
        start, end = end, start
    store[category][chrom].append((start, end))

def merged_bp(chrom_intervals):
    total = 0
    for intervals in chrom_intervals.values():
        if not intervals:
            continue
        intervals.sort()
        s0, e0 = intervals[0]
        for s, e in intervals[1:]:
            if s <= e0 + 1:
                if e > e0:
                    e0 = e
            else:
                total += e0 - s0 + 1
                s0, e0 = s, e
        total += e0 - s0 + 1
    return total

def pct(n, d):
    return 0.0 if not d else n * 100.0 / d

species_rows = []
with OLD.open(newline='', encoding='utf-8') as fh:
    reader = csv.DictReader(fh, delimiter='\t')
    for row in reader:
        species_rows.append(row)

order_index = {s:i for i,s in enumerate(PREFERRED_ORDER)}
species_rows.sort(key=lambda r: (order_index.get(r['species'], 999), r['species']))

rows = []
raw_rows = []
for srow in species_rows:
    species = srow['species']
    gff = Path(srow['teanno_gff3'])
    genome_bp = int(float(srow['genome_bp']))
    intervals = defaultdict(lambda: defaultdict(list))
    raw_counter = Counter()
    feature_counter = Counter()
    if not gff.exists():
        raise FileNotFoundError(f'Missing GFF3 for {species}: {gff}')
    with gff.open(encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attrs_text = parts[:9]
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue
            attrs = parse_attrs(attrs_text)
            classification = attrs.get('Classification', '')
            raw_counter[classification or 'NA'] += 1
            feature_counter[feature or 'NA'] += 1
            category = broad_category(feature, classification)
            add_interval(intervals, category, chrom, start_i, end_i)
    for category in CATEGORY_ORDER:
        bp = merged_bp(intervals.get(category, {}))
        rows.append({
            'species': species,
            'label': LABELS.get(species, species),
            'category': category,
            'category_order': CATEGORY_ORDER.index(category) + 1,
            'species_order': order_index.get(species, 999),
            'bp_merged': bp,
            'genome_bp': genome_bp,
            'percent_of_genome': f'{pct(bp, genome_bp):.6f}',
        })
    for cls, count in raw_counter.most_common():
        raw_rows.append({'species': species, 'classification': cls, 'line_count': count})

OUT.joinpath('tables').mkdir(parents=True, exist_ok=True)
with TABLE.open('w', newline='', encoding='utf-8') as fh:
    fields = ['species','label','category','category_order','species_order','bp_merged','genome_bp','percent_of_genome']
    writer = csv.DictWriter(fh, delimiter='\t', fieldnames=fields)
    writer.writeheader(); writer.writerows(rows)

lookup = {(r['species'], r['category']): r['percent_of_genome'] for r in rows}
species_names = [r['species'] for r in species_rows]
with WIDE.open('w', newline='', encoding='utf-8') as fh:
    writer = csv.writer(fh, delimiter='\t')
    writer.writerow(['species','label'] + CATEGORY_ORDER)
    for sp in species_names:
        writer.writerow([sp, LABELS.get(sp, sp)] + [lookup.get((sp, c), '0.000000') for c in CATEGORY_ORDER])

with CLASS_RAW.open('w', newline='', encoding='utf-8') as fh:
    fields = ['species','classification','line_count']
    writer = csv.DictWriter(fh, delimiter='\t', fieldnames=fields)
    writer.writeheader(); writer.writerows(raw_rows)

print(f'Wrote {TABLE}')
print(f'Wrote {WIDE}')
