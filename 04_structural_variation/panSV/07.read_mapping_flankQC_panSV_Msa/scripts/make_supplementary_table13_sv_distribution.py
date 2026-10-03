#!/usr/bin/env python3
from pathlib import Path
import csv

base = Path('path/to/project/N_3.call_SV/07.read_mapping_flankQC_panSV_Msa')
overview_path = base / 'results/Msa/species_level_mergedSV/summary/species_level_mergedSV.overview.tsv'
out_dir = base / 'results/Msa/tables'
out_dir.mkdir(parents=True, exist_ok=True)

# Order and names follow the user's Medicago material table; genome_A17 is excluded.
rows_meta = [
    ('genome_474', 'Mcar', 'Medicago carstiensis'),
    ('genome_472', 'Msuf', 'Medicago suffruticosa'),
    ('genome_468', 'Mcre', 'Medicago cretacea'),
    ('genome_482', 'Medg', 'Medicago edgeworthii'),
    ('genome_M46', 'Mfis', 'Medicago fischeriana'),
    ('genome_457', 'Mmar', 'Medicago marina'),
    ('genome_454', 'Mlan', 'Medicago lanigera'),
    ('genome_M22', 'Morb', 'Medicago orbicularis'),
    ('genome_410', 'Mpra', 'Medicago praecox'),
    ('genome_436', 'Mrad', 'Medicago radiata'),
    ('genome_461', 'Msec', 'Medicago secundiflora'),
    ('genome_395', 'Mlup', 'Medicago lupulina'),
    ('genome_Msa', 'Msat_T2T', 'Medicago sativa'),
    ('genome_ZM4', 'Msat_ZM4', 'Medicago sativa cv. Zhongmu-4'),
    ('genome_Mpo', 'Mpol', 'Medicago polymorpha'),
    ('genome_R108', 'Mtru_R108', 'Medicago truncatula R108'),
    ('genome_Mar', 'Marc', 'Medicago archiducis-nicolai'),
    ('genome_Mru', 'Mrut', 'Medicago ruthenica'),
]

with overview_path.open(newline='') as fh:
    reader = csv.DictReader(fh, delimiter='\t')
    counts = {row['species']: row for row in reader}

missing = [g for g, _, _ in rows_meta if g not in counts]
if missing:
    raise SystemExit('Missing genome ids in overview: ' + ','.join(missing))

field_exact = ['Sample id', 'Species name', 'Deletion', 'Insertion', 'Translocation', 'Duplication', 'Inversion']
field_trace = ['Genome id'] + field_exact + ['Total', 'read_based_only', 'read_based_plus_SVGAP']

out_exact = out_dir / 'Supplementary_Table_13_SV_distribution_each_genome.tsv'
out_trace = out_dir / 'Supplementary_Table_13_SV_distribution_each_genome.with_genome_id.tsv'
out_total = out_dir / 'Supplementary_Table_13_SV_distribution_total_summary.tsv'

records_exact = []
records_trace = []
summary = {k: 0 for k in ['Deletion', 'Insertion', 'Translocation', 'Duplication', 'Inversion', 'Total']}

for genome, sample_id, species_name in rows_meta:
    c = counts[genome]
    rec = {
        'Sample id': sample_id,
        'Species name': species_name,
        'Deletion': int(c['DEL']),
        'Insertion': int(c['INS']),
        'Translocation': int(c['TRA']),
        'Duplication': int(c['DUP']),
        'Inversion': int(c['INV']),
    }
    total = rec['Deletion'] + rec['Insertion'] + rec['Translocation'] + rec['Duplication'] + rec['Inversion']
    for k in ['Deletion', 'Insertion', 'Translocation', 'Duplication', 'Inversion']:
        summary[k] += rec[k]
    summary['Total'] += total
    records_exact.append(rec)
    trace = {'Genome id': genome, **rec, 'Total': total,
             'read_based_only': int(c['read_based_only']),
             'read_based_plus_SVGAP': int(c['read_based_plus_SVGAP'])}
    records_trace.append(trace)

for path, fields, records in [(out_exact, field_exact, records_exact), (out_trace, field_trace, records_trace)]:
    with path.open('w', newline='') as fh:
        writer = csv.DictWriter(fh, delimiter='\t', fieldnames=fields)
        writer.writeheader()
        for rec in records:
            # Keep TSV numeric cells unformatted for downstream spreadsheet use.
            writer.writerow(rec)

with out_total.open('w', newline='') as fh:
    writer = csv.writer(fh, delimiter='\t')
    writer.writerow(['SVTYPE', 'count'])
    for k in ['Deletion', 'Insertion', 'Translocation', 'Duplication', 'Inversion', 'Total']:
        writer.writerow([k, summary[k]])

# Also write a markdown preview with thousands separators for quick reading.
out_md = out_dir / 'Supplementary_Table_13_SV_distribution_each_genome.preview.md'
with out_md.open('w') as fh:
    fh.write('| Sample id | Species name | Deletion | Insertion | Translocation | Duplication | Inversion |\n')
    fh.write('|---|---:|---:|---:|---:|---:|---:|\n')
    for rec in records_exact:
        fh.write('| ' + ' | '.join([
            rec['Sample id'],
            rec['Species name'],
            f"{rec['Deletion']:,}",
            f"{rec['Insertion']:,}",
            f"{rec['Translocation']:,}",
            f"{rec['Duplication']:,}",
            f"{rec['Inversion']:,}",
        ]) + ' |\n')

print('WROTE', out_exact)
print('WROTE', out_trace)
print('WROTE', out_total)
print('WROTE', out_md)
print('TOTALS', '\t'.join(f'{k}={summary[k]}' for k in ['Deletion','Insertion','Translocation','Duplication','Inversion','Total']))