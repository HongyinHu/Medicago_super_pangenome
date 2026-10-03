#!/usr/bin/env python3
from pathlib import Path
import csv, re, sys
out = Path('path/to/project/N_1.EDTA_single/02_EDTA_single/LAI_recalc_single_EDTA_20260625')
manifest = out/'manifest_single_EDTA_LAI.tsv'
summary = out/'LAI_single_EDTA_recalc_summary.tsv'
rows=[]
with manifest.open() as f:
    for r in csv.DictReader(f, delimiter='\t'):
        sample = r['sample']
        wd = out/'work'/sample
        lai = wd/'all.mod.out.LAI'
        status = 'done' if (out/'status'/f'{sample}.done').exists() else ('failed' if (out/'status'/f'{sample}.failed').exists() else ('running' if (out/'status'/f'{sample}.lock').exists() else ('pending' if r['status']=='ready' else r['status'])))
        rec = dict(sample=sample, status=status, Chr='', From='', To='', Intact='', Total='', raw_LAI='', LAI='', lai_file=str(lai) if lai.exists() else '', fa=r['fa'], pass_list=r['pass_list'], rmout=r['rmout'])
        if lai.exists() and lai.stat().st_size:
            with lai.open(errors='replace') as lf:
                header=None
                for line in lf:
                    line=line.rstrip('\n')
                    if not line: continue
                    parts=line.split('\t')
                    if parts[0]=='Chr':
                        header=parts
                    elif parts[0]=='whole_genome' and header:
                        vals=dict(zip(header, parts))
                        for k in ('Chr','From','To','Intact','Total','raw_LAI','LAI'):
                            rec[k]=vals.get(k,'')
                        break
        rows.append(rec)
fields=['sample','status','Chr','From','To','Intact','Total','raw_LAI','LAI','lai_file','fa','pass_list','rmout']
with summary.open('w', newline='') as f:
    w=csv.DictWriter(f, fieldnames=fields, delimiter='\t')
    w.writeheader(); w.writerows(rows)
print(summary)
counts={}
for r in rows: counts[r['status']]=counts.get(r['status'],0)+1
print('counts', counts)
for r in rows:
    if r['status']!='done': print('NOT_DONE', r['sample'], r['status'])
