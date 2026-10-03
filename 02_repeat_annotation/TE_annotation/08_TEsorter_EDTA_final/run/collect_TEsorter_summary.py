#!/usr/bin/env python3
from pathlib import Path
import csv, re
BASE=Path('path/to/project/N_1.EDTA_single')
OUT=BASE/'08_TEsorter_EDTA_final'
manifest=OUT/'manifest.tsv'
summary=OUT/'TEsorter_summary.tsv'
rows=[]
with manifest.open() as f:
    for r in csv.DictReader(f, delimiter='\t'):
        sample=r['sample']; outdir=Path(r['out_dir'])
        done=outdir/'.done'; failed=outdir/'.failed'
        status='pending'
        if done.exists(): status='ok'
        if failed.exists(): status='failed'
        cls_files=sorted(outdir.glob('*.cls.tsv'))
        n_total=n_classified=n_unknown=0
        if cls_files:
            cls=cls_files[0]
            with cls.open(errors='replace') as cf:
                header=cf.readline().rstrip('\n').split('\t')
                for line in cf:
                    if not line.strip(): continue
                    n_total += 1
                    parts=line.rstrip('\n').split('\t')
                    text='\t'.join(parts).lower()
                    if 'unknown' in text or 'unclassified' in text:
                        n_unknown += 1
                    else:
                        n_classified += 1
        rows.append({**r,'status':status,'cls_tsv':str(cls_files[0]) if cls_files else '', 'n_records_in_cls':n_total, 'n_classified_rough':n_classified, 'n_unknown_rough':n_unknown})
with summary.open('w', newline='') as f:
    fields=['sample','type','input_path','input_size','input_mtime','out_dir','status','cls_tsv','n_records_in_cls','n_classified_rough','n_unknown_rough']
    w=csv.DictWriter(f, fieldnames=fields, delimiter='\t')
    w.writeheader(); w.writerows(rows)
print(summary)
