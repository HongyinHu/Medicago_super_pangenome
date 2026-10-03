#!/usr/bin/env python3
import os, glob, csv, re, math
BASE=os.environ.get('BASE','.')
target_chrom='Chr4'
target_start=88980831
target_end=88981052
target_len=target_end-target_start+1
max_bp_dist=1000
max_len_for_strict=2000
min_len_for_sv=50

def parse_info(s):
    d={}
    for part in s.split(';'):
        if not part: continue
        if '=' in part:
            k,v=part.split('=',1); d[k]=v
        else:
            d[part]=True
    return d

def parse_sample(fmt, sample):
    keys=fmt.split(':')
    vals=sample.split(':')
    return dict(zip(keys, vals))

def is_alt_gt(gt):
    if gt in ('0/0','0|0','./.','.'):
        return False
    return '1' in re.split(r'[\/|]', gt)

def overlap(a1,a2,b1,b2):
    return max(0, min(a2,b2)-max(a1,b1)+1)

# phenotype map
pheno={}
with open(os.path.join(BASE,'input','valid_phenotyped_bams.for_delly.tsv')) as f:
    for line in f:
        if not line.strip(): continue
        sid, phe, bam = line.rstrip('\n').split('\t')[:3]
        pheno[sid]=phe

out_rows=[]
strict_rows=[]
for path in sorted(glob.glob(os.path.join(BASE,'results','delly_DEL_per_sample','*','*.delly.DEL.Chr23997_window.vcf'))):
    sid=os.path.basename(path).split('.delly.DEL.Chr23997_window.vcf')[0]
    best_any=None
    best_bp=None
    strict=[]
    records=0
    alt_records=0
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            records+=1
            fields=line.rstrip('\n').split('\t')
            if len(fields)<10: continue
            chrom,pos,vid,ref,alt,qual,flt,info_s,fmt,samp=fields[:10]
            if chrom!=target_chrom or alt!='<DEL>': continue
            pos=int(pos)
            info=parse_info(info_s)
            end=int(info.get('END',pos))
            svlen=abs(int(info.get('SVLEN', end-pos+1))) if str(info.get('SVLEN','')).lstrip('-').isdigit() else abs(end-pos+1)
            sm=parse_sample(fmt,samp)
            gt=sm.get('GT','.')
            alt_gt=is_alt_gt(gt)
            if alt_gt: alt_records+=1
            ov=overlap(pos,end,target_start,target_end)
            ro_target=ov/target_len if target_len else 0
            ro_record=ov/svlen if svlen else 0
            bp_dist=max(abs(pos-target_start), abs(end-target_end))
            center_dist=abs(((pos+end)/2)-((target_start+target_end)/2))
            row={
                'sample_id':sid,'phenotype':pheno.get(sid,'NA'),'chrom':chrom,'pos':pos,'end':end,'id':vid,
                'svlen':svlen,'filter':flt,'qual':qual,'gt':gt,'sample_filter':sm.get('FT',''),
                'PE':info.get('PE','0'),'SR':info.get('SR','0'),'MAPQ':info.get('MAPQ',''),
                'DR':sm.get('DR',''),'DV':sm.get('DV',''),'RR':sm.get('RR',''),'RV':sm.get('RV',''),
                'overlap_bp':ov,'reciprocal_target':round(ro_target,4),'reciprocal_record':round(ro_record,4),
                'bp_dist':bp_dist,'center_dist':round(center_dist,1),'alt_gt':int(alt_gt),
                'raw':line.rstrip('\n')
            }
            # best by breakpoint distance among plausible SV-size calls
            if svlen>=min_len_for_sv and svlen<=max_len_for_strict and alt_gt:
                if best_bp is None or (bp_dist, -ov, svlen) < (best_bp['bp_dist'], -best_bp['overlap_bp'], best_bp['svlen']):
                    best_bp=row
            # strict match: alt genotype, reasonable length, both breakpoints close, reciprocal overlap with candidate not tiny
            if alt_gt and svlen>=min_len_for_sv and svlen<=max_len_for_strict and bp_dist<=max_bp_dist and ro_target>=0.5 and ro_record>=0.2:
                strict.append(row)
            # best any overlapping candidate but not huge: useful QC
            if alt_gt and svlen>=min_len_for_sv and svlen<=10000 and ov>0:
                if best_any is None or (-(min(ro_target,ro_record)), bp_dist) < (-(min(best_any['reciprocal_target'],best_any['reciprocal_record'])), best_any['bp_dist']):
                    best_any=row
    chosen = sorted(strict, key=lambda r:(r['bp_dist'], -r['overlap_bp'], r['svlen']))[0] if strict else None
    status = 'strict_match' if chosen else ('nearby_alt_no_strict' if best_bp else 'no_nearby_alt')
    out_rows.append({
        'sample_id':sid,'phenotype':pheno.get(sid,'NA'),'records_in_window':records,'alt_gt_records_in_window':alt_records,
        'status':status,
        'strict_match_count':len(strict),
        'best_pos': chosen['pos'] if chosen else (best_bp['pos'] if best_bp else ''),
        'best_end': chosen['end'] if chosen else (best_bp['end'] if best_bp else ''),
        'best_svlen': chosen['svlen'] if chosen else (best_bp['svlen'] if best_bp else ''),
        'best_filter': chosen['filter'] if chosen else (best_bp['filter'] if best_bp else ''),
        'best_gt': chosen['gt'] if chosen else (best_bp['gt'] if best_bp else ''),
        'best_PE': chosen['PE'] if chosen else (best_bp['PE'] if best_bp else ''),
        'best_SR': chosen['SR'] if chosen else (best_bp['SR'] if best_bp else ''),
        'best_overlap_bp': chosen['overlap_bp'] if chosen else (best_bp['overlap_bp'] if best_bp else ''),
        'best_reciprocal_target': chosen['reciprocal_target'] if chosen else (best_bp['reciprocal_target'] if best_bp else ''),
        'best_reciprocal_record': chosen['reciprocal_record'] if chosen else (best_bp['reciprocal_record'] if best_bp else ''),
        'best_bp_dist': chosen['bp_dist'] if chosen else (best_bp['bp_dist'] if best_bp else ''),
    })
    for r in strict:
        strict_rows.append(r)

outdir=os.path.join(BASE,'results','Chr23997_short_read_DEL_summary')
os.makedirs(outdir, exist_ok=True)
summary_path=os.path.join(outdir,'Chr23997_DEL_221.delly_per_sample_match.tsv')
with open(summary_path,'w',newline='') as f:
    cols=list(out_rows[0].keys()) if out_rows else []
    w=csv.DictWriter(f, fieldnames=cols, delimiter='\t')
    w.writeheader(); w.writerows(out_rows)
strict_path=os.path.join(outdir,'Chr23997_DEL_221.delly_strict_records.tsv')
with open(strict_path,'w',newline='') as f:
    cols=['sample_id','phenotype','chrom','pos','end','id','svlen','filter','qual','gt','sample_filter','PE','SR','MAPQ','DR','DV','RR','RV','overlap_bp','reciprocal_target','reciprocal_record','bp_dist','center_dist','alt_gt','raw']
    w=csv.DictWriter(f, fieldnames=cols, delimiter='\t', extrasaction='ignore')
    w.writeheader(); w.writerows(strict_rows)
# summary counts
from collections import Counter,defaultdict
cnt=Counter(r['status'] for r in out_rows)
by=defaultdict(Counter)
for r in out_rows:
    by[r['phenotype']][r['status']]+=1
text=[]
text.append(f'target\t{target_chrom}:{target_start}-{target_end}\n')
text.append(f'target_len\t{target_len}\n')
text.append(f'total_samples\t{len(out_rows)}\n')
for k in sorted(cnt): text.append(f'status\t{k}\t{cnt[k]}\n')
for phe in sorted(by):
    total=sum(by[phe].values())
    for st in sorted(by[phe]):
        text.append(f'phenotype_status\t{phe}\t{st}\t{by[phe][st]}\t{by[phe][st]/total:.4f}\n')
with open(os.path.join(outdir,'Chr23997_DEL_221.delly_match.summary.txt'),'w') as f:
    f.writelines(text)
print(summary_path)
print(strict_path)
print(os.path.join(outdir,'Chr23997_DEL_221.delly_match.summary.txt'))
