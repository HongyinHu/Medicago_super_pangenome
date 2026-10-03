#!/usr/bin/env python3
from pathlib import Path
import csv
OUT=Path('path/to/project/N_1.EDTA_single/02_EDTA_single/LAI_recalc_single_EDTA_IDfixed_20260625')
manifest=OUT/'prepare_manifest_single_idfixed.tsv'
summary=OUT/'LAI_single_EDTA_IDfixed_summary.tsv'
rows=[]
with manifest.open() as f:
    for r in csv.DictReader(f, delimiter='\t'):
        sample=r['sample']; work=Path(r['workdir']); lai=work/f'{sample}.singleEDTA.IDfixed.out.LAI'
        status='done' if (OUT/'logs/done'/f'{sample}.done').exists() else ('not_applicable' if (OUT/'logs/not_applicable'/f'{sample}.done').exists() else ('failed' if (OUT/'logs/failed'/f'{sample}.failed').exists() else ('running' if (OUT/'logs/locks'/f'{sample}.lockdir').exists() else ('pending' if r['status'].startswith('ok') else r['status']))))
        rec={k:'' for k in ['Chr','From','To','Intact','Total','raw_LAI','LAI']}
        if lai.exists() and lai.stat().st_size:
            with lai.open(errors='replace') as fh:
                header=None
                for line in fh:
                    parts=line.rstrip('\n').split('\t')
                    if not parts or not parts[0]: continue
                    if parts[0]=='Chr': header=parts
                    elif parts[0]=='whole_genome' and header:
                        vals=dict(zip(header,parts))
                        for k in rec: rec[k]=vals.get(k,'')
                        break
        rows.append(dict(sample=sample,status=status,**rec,lai_file=str(lai) if lai.exists() else '',original_genome=r['original_genome'],mod_genome=r['mod_genome'],pass_list=r['pass_list'],rmout=r['rmout'],idfixed_out=r['idfixed_out'],mapped_records=r['mapped_records'],unchanged_records=r['unchanged_records'],unmapped_records=r['unmapped_records'],prepare_status=r['status']))
fields=['sample','status','Chr','From','To','Intact','Total','raw_LAI','LAI','lai_file','original_genome','mod_genome','pass_list','rmout','idfixed_out','mapped_records','unchanged_records','unmapped_records','prepare_status']
with summary.open('w', newline='') as f:
    w=csv.DictWriter(f, fieldnames=fields, delimiter='\t')
    w.writeheader(); w.writerows(rows)
print(summary)
counts={}
for r in rows: counts[r['status']]=counts.get(r['status'],0)+1
print('counts', counts)
