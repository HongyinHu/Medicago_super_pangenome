#!/usr/bin/env python3
from pathlib import Path
import csv, subprocess, math, statistics, sys
samtools=sys.argv[1]
bamdir=Path(sys.argv[2])
pheno=Path(sys.argv[3])
qc=Path(sys.argv[4])
out=Path(sys.argv[5])
chrom='Chr4'; start=88980831; end=88981052; flank=500
# read pheno, strip CR
ph={}
with pheno.open(newline='') as f:
    for line in f:
        line=line.strip()
        if not line or line.startswith('#'): continue
        parts=line.split()
        if len(parts)>=3 and parts[2] != '-9':
            ph[parts[1]]=parts[2]
# pass quickcheck
valid=set()
with qc.open(newline='') as f:
    rd=csv.DictReader(f, delimiter='\t')
    for r in rd:
        sid=(r.get('sample_id') or '').strip()
        status=(r.get('quickcheck_status') or '').strip()
        if sid and status=='PASS': valid.add(sid)
def mean_depth(bam,a,b):
    cmd=[samtools,'depth','-aa','-r',f'{chrom}:{a}-{b}',str(bam)]
    p=subprocess.run(cmd, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if p.returncode!=0:
        return math.nan,0,p.stderr.replace('\n','; ')
    vals=[]
    for line in p.stdout.splitlines():
        sp=line.split('\t')
        if len(sp)>=3:
            try: vals.append(int(sp[2]))
            except: pass
    return (sum(vals)/len(vals) if vals else 0.0), len(vals), ''
rows=[]
for sid in sorted(ph):
    if sid not in valid: continue
    code=ph[sid]
    phen='spiny' if code=='2' else 'spineless' if code=='1' else 'unknown'
    bam=bamdir/f'{sid}.dedup.bam'
    if not bam.exists(): continue
    lm,lb,le=mean_depth(bam,max(1,start-flank),start-1)
    sm,sb,se=mean_depth(bam,start,end)
    rm,rb,re=mean_depth(bam,end+1,end+flank)
    flank_mean=(lm+rm)/2 if lb and rb else math.nan
    ratio=sm/flank_mean if flank_mean and not math.isnan(flank_mean) else math.nan
    if math.isnan(ratio): gt='NA'
    elif ratio<=0.35: gt='DEL_hom_or_strong'
    elif ratio<=0.70: gt='DEL_het_or_partial'
    else: gt='REF_or_no_depth_loss'
    rows.append(dict(sample_id=sid, phenotype_code=code, phenotype=phen, bam=str(bam), chrom=chrom, start=start, end=end, flank=flank, left_mean=lm, sv_mean=sm, right_mean=rm, flank_mean=flank_mean, sv_to_flank_ratio=ratio, left_bp=lb, sv_bp=sb, right_bp=rb, depth_gt=gt, depth_error=le+se+re))
fields=['sample_id','phenotype_code','phenotype','bam','chrom','start','end','flank','left_mean','sv_mean','right_mean','flank_mean','sv_to_flank_ratio','left_bp','sv_bp','right_bp','depth_gt','depth_error']
out.parent.mkdir(parents=True, exist_ok=True)
with out.open('w', newline='') as f:
    w=csv.DictWriter(f, delimiter='\t', fieldnames=fields)
    w.writeheader(); w.writerows(rows)
summary=out.with_suffix('.summary.txt')
with summary.open('w') as f:
    f.write(f'valid_rows\t{len(rows)}\n')
    for phen in ['spiny','spineless']:
        sub=[r for r in rows if r['phenotype']==phen]
        ratios=[r['sv_to_flank_ratio'] for r in sub if not math.isnan(r['sv_to_flank_ratio'])]
        counts={}
        for r in sub: counts[r['depth_gt']]=counts.get(r['depth_gt'],0)+1
        f.write(f'{phen}\tn={len(sub)}\tmean_ratio={(statistics.mean(ratios) if ratios else float("nan")):.4f}\tmedian_ratio={(statistics.median(ratios) if ratios else float("nan")):.4f}\tcounts={counts}\n')
print(summary.read_text())
