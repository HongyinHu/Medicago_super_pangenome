#!/usr/bin/env python3
import csv
from pathlib import Path
REC=Path('path/to/project/N_4.pod_spiny/05_direction_free_sv_cosegregation_20260701/summary/Chr23997_caller_records_20260701.tsv')
OUT=Path('path/to/project/N_4.pod_spiny/19_Chr23997_callset_truth_calibration_20260708')
TCHROM='Chr4'; TSTART=88980831; TEND=88981052; TLEN=221
samples={'genome_474','genome_R108','genome_436','genome_457','genome_454','genome_A17','genome_395','genome_461','genome_468','genome_472','genome_M22','genome_Mar','genome_Mru'}
rows=[]
with REC.open() as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if r['sample'] not in samples or r['chrom'] != TCHROM:
            continue
        try:
            pos=int(r['pos']); end=int(r['end']); svlen=abs(int(float(r['svlen'])))
        except Exception:
            continue
        # keep records within target +/- 2 kb or overlapping target
        if end < TSTART-2000 or pos > TEND+2000:
            continue
        ov=max(0, min(end,TEND)-max(pos,TSTART)+1)
        rec_len=max(1,end-pos+1)
        target_len=TEND-TSTART+1
        ro=min(ov/rec_len, ov/target_len) if ov else 0.0
        bpdist=max(abs(pos-TSTART), abs(end-TEND))
        len_ratio=(svlen/TLEN) if TLEN else 0
        exact_like = (r['svtype']=='DEL' and ov>0 and ro>=0.5 and 0.5<=len_ratio<=2.0 and bpdist<=500)
        rows.append({
            'sample':r['sample'], 'phenotype':r['phenotype'], 'caller':r['caller'], 'id':r['id'],
            'chrom':r['chrom'], 'pos':pos, 'end':end, 'svtype':r['svtype'], 'svlen_abs':svlen,
            'overlap_bp':ov, 'reciprocal_overlap_min':round(ro,3), 'breakpoint_max_dist':bpdist, 'len_ratio_to_target':round(len_ratio,3),
            'exact_target_like_DEL': 'yes' if exact_like else 'no',
            'filter':r.get('filter',''), 'near_manual_del':r.get('near_manual_del',''), 'overlap_manual_window':r.get('overlap_manual_window','')
        })
rows.sort(key=lambda x:(x['sample'], x['caller'], x['pos'], x['end']))
fields=list(rows[0].keys()) if rows else []
out=OUT/'results/Chr23997.target_DEL_exactness_audit_from_raw_callers.tsv'
with out.open('w') as fo:
    w=csv.DictWriter(fo, delimiter='\t', fieldnames=fields)
    w.writeheader(); w.writerows(rows)
# summarize exact-like by sample/caller
summary={}
for r in rows:
    key=(r['sample'],r['caller'])
    s=summary.setdefault(key, {'sample':r['sample'],'caller':r['caller'],'near_records':0,'exact_target_like_DEL_records':0,'ids':[]})
    s['near_records']+=1
    if r['exact_target_like_DEL']=='yes':
        s['exact_target_like_DEL_records']+=1; s['ids'].append(r['id'])
srows=[]
for s in summary.values():
    s['exact_target_like_ids']=','.join(s['ids']) if s['ids'] else 'NA'
    del s['ids']
    srows.append(s)
srows.sort(key=lambda x:(x['sample'],x['caller']))
sumout=OUT/'summary/Chr23997.target_DEL_exactness_by_sample_caller.tsv'
with sumout.open('w') as fo:
    w=csv.DictWriter(fo, delimiter='\t', fieldnames=['sample','caller','near_records','exact_target_like_DEL_records','exact_target_like_ids'])
    w.writeheader(); w.writerows(srows)
print(out)
print(sumout)
