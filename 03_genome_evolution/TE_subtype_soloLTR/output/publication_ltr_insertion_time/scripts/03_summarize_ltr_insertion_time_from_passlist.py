#!/usr/bin/env python3
import csv
from pathlib import Path

BASE = Path('path/to/project/39.TE_type_soloLTR')
OUT = BASE / 'output' / 'publication_ltr_insertion_time'
SPECIES_TOTALS = BASE / 'output' / 'species_totals.tsv'
LABELS = {
    'genome_Mar': 'Marc', 'genome_Mru': 'Mrut', 'genome_MPO': 'Mpol',
    'genome_ZM4': 'Msat_ZM4', 'genome_R108': 'Mtru_R108', 'genome_Msa': 'Msat_Cae',
    'genome_Msa_T2T': 'Msat_T2T', 'genome_A17': 'Mtru_A17',
}
PREFERRED_ORDER = [
    'genome_395','genome_436','genome_457','genome_468','genome_474','genome_482',
    'genome_M22','genome_Mar','genome_Mru','genome_Msa','genome_ZM4',
    'genome_410','genome_454','genome_461','genome_472','genome_474_T2T',
    'genome_A17','genome_M46','genome_MPO','genome_Msa_T2T','genome_R108'
]

def find_pass_list(species):
    d = BASE / 'data' / f'{species}_EDTA'
    patterns = [
        'cal_LAI/*.pass.list',
        '*.pass.list',
        '*.EDTA.raw/LTR/*.pass.list',
        '*.EDTA.raw/*.pass.list',
        '*.EDTA.anno/*.pass.list',
    ]
    hits = []
    for pat in patterns:
        hits.extend([p for p in d.glob(pat) if p.is_file() and p.stat().st_size > 0 and 'nmtf' not in p.name])
        if hits:
            break
    if hits:
        return sorted(hits, key=lambda p: (0 if 'cal_LAI' in str(p) else 1, len(str(p)), str(p)))[0]
    return None

def norm_class(sf):
    sf = (sf or '').strip().lower()
    if sf == 'gypsy': return 'Gypsy'
    if sf == 'copia': return 'Copia'
    return 'Unknown'

species=[]
with SPECIES_TOTALS.open(encoding='utf-8', newline='') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        species.append(r['species'])
order_index={s:i for i,s in enumerate(PREFERRED_ORDER)}
species=sorted(species, key=lambda s:(order_index.get(s,999), s))
records=[]; manifest=[]
for sp in species:
    label=LABELS.get(sp, sp)
    p=find_pass_list(sp)
    if not p:
        manifest.append({'species':sp,'label':label,'pass_list':'','status':'missing','records':0})
        continue
    n=0
    with p.open(encoding='utf-8', errors='replace') as fh:
        header=None
        for line in fh:
            line=line.rstrip('\n')
            if not line:
                continue
            if line.startswith('#'):
                header=line.lstrip('#').split()
                continue
            parts=line.split()
            if header and len(parts) >= len(header):
                row=dict(zip(header, parts))
                loc=row.get('LTR_loc','')
                sf=row.get('SuperFamily','')
                identity=row.get('Identity','')
                it=row.get('Insertion_Time','')
            else:
                # EDTA pass.list standard columns: loc category motif tsd 5tsd 3tsd internal identity strand superfamily te_type insertion_time
                if len(parts) < 12:
                    continue
                loc, identity, sf, it = parts[0], parts[7], parts[9], parts[11]
            try:
                time_mya=float(it)/1e6
            except ValueError:
                continue
            try:
                ident=float(identity)
            except ValueError:
                ident=float('nan')
            ltr_class=norm_class(sf)
            records.append({'species':sp,'label':label,'species_order':order_index.get(sp,999),'ltr_id':loc,'superfamily':sf,'ltr_class':ltr_class,'identity':identity,'insertion_time_mya':f'{time_mya:.6f}','source_pass_list':str(p)})
            n += 1
    manifest.append({'species':sp,'label':label,'pass_list':str(p),'status':'ok','records':n})
OUT.joinpath('tables').mkdir(parents=True, exist_ok=True)
with (OUT/'tables/ltr_insertion_time_passlist_records.tsv').open('w', encoding='utf-8', newline='') as fh:
    fields=['species','label','species_order','ltr_id','superfamily','ltr_class','identity','insertion_time_mya','source_pass_list']
    w=csv.DictWriter(fh, delimiter='\t', fieldnames=fields); w.writeheader(); w.writerows(records)
with (OUT/'tables/ltr_insertion_time_passlist_manifest.tsv').open('w', encoding='utf-8', newline='') as fh:
    fields=['species','label','pass_list','status','records']
    w=csv.DictWriter(fh, delimiter='\t', fieldnames=fields); w.writeheader(); w.writerows(manifest)
print(f'wrote {len(records)} pass.list records')
