#!/usr/bin/env python3
from pathlib import Path
from collections import Counter
import csv

base = Path('path/to/project/N_3.call_SV/07.read_mapping_flankQC_panSV_Msa')
summary_dir = base / 'results/Msa/species_level_mergedSV/summary'
pan_dir = base / 'results/Msa/panSV'
fig_dir = base / 'results/Msa/figures'
table_dir = base / 'results/Msa/tables'
fig_dir.mkdir(parents=True, exist_ok=True)
table_dir.mkdir(parents=True, exist_ok=True)

exclude = 'genome_Msa'
overview_in = summary_dir / 'species_level_mergedSV.overview.tsv'
overview_out = fig_dir / 'species_level_mergedSV.overview.exclude_genome_Msa.tsv'
pav_in = pan_dir / 'panSV.Msa.flankQC.read_primary.PAV.matrix.tsv'
stats_out = pan_dir / 'panSV.Msa.flankQC.read_primary.exclude_genome_Msa.stats.tsv'
pav17_out = pan_dir / 'panSV.Msa.flankQC.read_primary.exclude_genome_Msa.PAV17.summary.tsv'

# Bars: remove the self-reference row.
with overview_in.open(newline='') as fin, overview_out.open('w', newline='') as fout:
    reader = csv.DictReader(fin, delimiter='\t')
    writer = csv.DictWriter(fout, delimiter='\t', fieldnames=reader.fieldnames)
    writer.writeheader()
    n_in = n_out = 0
    for row in reader:
        n_in += 1
        if row['species'] == exclude:
            continue
        n_out += 1
        writer.writerow(row)

# Pie: count non-redundant panSV events still present in at least one of the 17 non-Msa species.
svtypes = ['DEL', 'INS', 'TRA', 'DUP', 'INV']
counts = Counter()
source_counts = Counter()
rows_total = rows_kept = rows_msa_only = 0
with pav_in.open(newline='') as fin:
    reader = csv.DictReader(fin, delimiter='\t')
    genome_cols = [c for c in reader.fieldnames if c.startswith('genome_')]
    keep_cols = [c for c in genome_cols if c != exclude]
    if exclude not in genome_cols:
        raise SystemExit(f'{exclude} not found in PAV genome columns')
    for row in reader:
        rows_total += 1
        keep_present = any(row.get(c, '0') not in ('0', '.', '', 'NA') for c in keep_cols)
        msa_present = row.get(exclude, '0') not in ('0', '.', '', 'NA')
        if keep_present:
            rows_kept += 1
            counts[row['SVTYPE']] += 1
            source_counts[row.get('source_methods', 'NA')] += 1
        elif msa_present:
            rows_msa_only += 1

with stats_out.open('w', newline='') as fout:
    w = csv.writer(fout, delimiter='\t')
    w.writerow(['section', 'key', 'count'])
    w.writerow(['parameters', 'reference', 'genome_Msa'])
    w.writerow(['parameters', 'excluded_species', exclude])
    w.writerow(['parameters', 'event_policy', 'nonredundant_panSV_events_with_presence_in_remaining_17_species'])
    for sv in svtypes:
        w.writerow(['SVTYPE', sv, counts[sv]])
    w.writerow(['rows', 'original_final_events_18species', rows_total])
    w.writerow(['rows', 'final_events_exclude_genome_Msa', rows_kept])
    w.writerow(['rows', 'genome_Msa_only_events_removed', rows_msa_only])
    for k in sorted(source_counts):
        w.writerow(['source_methods', k, source_counts[k]])

with pav17_out.open('w', newline='') as fout:
    w = csv.writer(fout, delimiter='\t')
    w.writerow(['metric', 'value'])
    w.writerow(['overview_rows_in', n_in])
    w.writerow(['overview_rows_out', n_out])
    w.writerow(['pav_rows_in', rows_total])
    w.writerow(['pav_rows_kept_exclude_genome_Msa', rows_kept])
    w.writerow(['pav_rows_removed_genome_Msa_only', rows_msa_only])
    for sv in svtypes:
        w.writerow([f'SVTYPE_{sv}', counts[sv]])

print('WROTE', overview_out)
print('WROTE', stats_out)
print('WROTE', pav17_out)
print('BAR_SPECIES', n_out)
print('PIE_COUNTS', '\t'.join(f'{sv}={counts[sv]}' for sv in svtypes))
print('ROWS', f'original={rows_total}', f'kept17={rows_kept}', f'msa_only_removed={rows_msa_only}')