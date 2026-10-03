#!/usr/bin/env python3
import csv
import re
from pathlib import Path
from urllib.parse import unquote

BASE = Path('path/to/project/39.TE_type_soloLTR')
OUT = BASE / 'output' / 'publication_ltr_insertion_time'
SPECIES_TOTALS = BASE / 'output' / 'species_totals.tsv'
MU = 1.3e-8

# Keep all current EDTA objects; labels can be revised without rerunning stats.
LABELS = {
    'genome_Mar': 'Marc',
    'genome_Mru': 'Mrut',
    'genome_MPO': 'Mpol',
    'genome_ZM4': 'Msat_ZM4',
    'genome_R108': 'Mtru_R108',
    'genome_Msa': 'Msat_Cae',
    'genome_Msa_T2T': 'Msat_T2T',
    'genome_A17': 'Mtru_A17',
}
PREFERRED_ORDER = [
    'genome_395','genome_436','genome_457','genome_468','genome_474','genome_482',
    'genome_M22','genome_Mar','genome_Mru','genome_Msa','genome_ZM4',
    'genome_410','genome_454','genome_461','genome_472','genome_474_T2T',
    'genome_A17','genome_M46','genome_MPO','genome_Msa_T2T','genome_R108'
]

def parse_attrs(text):
    attrs = {}
    for part in text.strip().split(';'):
        if not part:
            continue
        if '=' in part:
            k, v = part.split('=', 1)
            attrs[k] = unquote(v)
    return attrs

def classify_ltr(classification, feature):
    c = (classification or '').strip()
    f = (feature or '').strip().lower()
    low = f'{c} {feature}'.lower()
    if c == 'LTR/Gypsy' or 'gypsy_ltr_retrotransposon' in low:
        return 'Gypsy'
    if c == 'LTR/Copia' or 'copia_ltr_retrotransposon' in low:
        return 'Copia'
    if c.startswith('LTR/') or 'ltr_retrotransposon' in f or f == 'repeat_region':
        return 'Unknown'
    return None

def find_ltr_intact_gff(edta_dir):
    patterns = [
        '*.EDTA.raw/LTR/*.LTR.intact.gff3',
        '*.EDTA.raw/*.LTR.intact.gff3',
        '*.EDTA.anno/*.EDTA.intact.gff3',
        '*.EDTA.final/*.EDTA.intact.gff3',
        '*.EDTA.intact.gff3',
    ]
    hits = []
    for pat in patterns:
        hits.extend([p for p in edta_dir.glob(pat) if p.is_file() and p.stat().st_size > 0])
    if hits:
        # Prefer raw/LTR intact file because it contains one repeat_region per intact LTR.
        hits = sorted(set(hits), key=lambda p: (0 if '/LTR/' in str(p) else 1, len(str(p)), str(p)))
        return hits[0]
    all_hits = sorted([p for p in edta_dir.rglob('*LTR.intact.gff3') if p.is_file() and p.stat().st_size > 0], key=lambda p: str(p))
    return all_hits[0] if all_hits else None

species_rows = []
with SPECIES_TOTALS.open(encoding='utf-8', newline='') as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        species_rows.append(row)
order_index = {s:i for i,s in enumerate(PREFERRED_ORDER)}
species_rows.sort(key=lambda r: (order_index.get(r['species'], 999), r['species']))

records = []
manifest = []
for sr in species_rows:
    sp = sr['species']
    # Infer EDTA directory from the TEanno path stored in existing species_totals.tsv.
    teanno = Path(sr['teanno_gff3'])
    edta_dir = None
    for parent in [teanno.parent, teanno.parent.parent, teanno.parent.parent.parent]:
        if parent.name.endswith('_EDTA') or parent.name == 'EDTA' or 'EDTA_monoploid' in parent.name:
            edta_dir = parent
            break
    if edta_dir is None:
        # fallback to symlink under data
        d = BASE / 'data' / f'{sp}_EDTA'
        edta_dir = d if d.exists() else teanno.parent
    gff = find_ltr_intact_gff(edta_dir)
    label = LABELS.get(sp, sp)
    count = 0
    used = set()
    if not gff:
        manifest.append({'species': sp, 'label': label, 'ltr_intact_gff3': '', 'status': 'missing_ltr_intact_gff3', 'intact_ltr_count': 0})
        continue
    with gff.open(encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attr_text = parts[:9]
            attrs = parse_attrs(attr_text)
            # Prefer one row per intact LTR. repeat_region is ideal; fallback to LTRRT feature if no repeat_region duplicate exists.
            if feature != 'repeat_region':
                continue
            ltr_id = attrs.get('ID') or attrs.get('Name') or f'{chrom}:{start}-{end}'
            if ltr_id in used:
                continue
            used.add(ltr_id)
            classification = attrs.get('Classification', '')
            ltr_class = classify_ltr(classification, feature)
            if not ltr_class:
                continue
            try:
                identity = float(attrs.get('ltr_identity', 'nan'))
            except ValueError:
                continue
            if not (0 < identity <= 1.0):
                continue
            time_mya = (1.0 - identity) / (2.0 * MU) / 1e6
            records.append({
                'species': sp,
                'label': label,
                'species_order': order_index.get(sp, 999),
                'ltr_id': ltr_id,
                'chrom': chrom,
                'start': start,
                'end': end,
                'classification': classification,
                'ltr_class': ltr_class,
                'ltr_identity': f'{identity:.6f}',
                'insertion_time_mya': f'{time_mya:.6f}',
                'source_gff3': str(gff),
            })
            count += 1
    manifest.append({'species': sp, 'label': label, 'ltr_intact_gff3': str(gff), 'status': 'ok', 'intact_ltr_count': count})

OUT.joinpath('tables').mkdir(parents=True, exist_ok=True)
fields = ['species','label','species_order','ltr_id','chrom','start','end','classification','ltr_class','ltr_identity','insertion_time_mya','source_gff3']
with (OUT/'tables/ltr_insertion_time_records.tsv').open('w', encoding='utf-8', newline='') as fh:
    w=csv.DictWriter(fh, delimiter='\t', fieldnames=fields); w.writeheader(); w.writerows(records)
with (OUT/'tables/ltr_insertion_time_manifest.tsv').open('w', encoding='utf-8', newline='') as fh:
    fields2=['species','label','ltr_intact_gff3','status','intact_ltr_count']
    w=csv.DictWriter(fh, delimiter='\t', fieldnames=fields2); w.writeheader(); w.writerows(manifest)
print(f'wrote {len(records)} records')
