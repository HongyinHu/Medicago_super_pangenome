#!/usr/bin/env python3
import sys, subprocess, math
samtools, bam, chrom, start, end = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
flank=int(sys.argv[6]) if len(sys.argv)>6 else 500
regions={
 'left':(max(1,start-flank), start-1),
 'sv':(start,end),
 'right':(end+1,end+flank),
}
def mean_depth(a,b):
    if b<a: return 0.0,0
    cmd=[samtools,'depth','-aa','-r','{}:{}-{}'.format(chrom,a,b),bam]
    vals=[]
    p=subprocess.run(cmd, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if p.returncode != 0:
        sys.stderr.write(p.stderr)
        return float('nan'),0
    for line in p.stdout.splitlines():
        parts=line.split('\t')
        if len(parts)>=3:
            try: vals.append(int(parts[2]))
            except Exception: pass
    return (sum(vals)/len(vals) if vals else 0.0), len(vals)
res={}
for k,(a,b) in regions.items(): res[k]=mean_depth(a,b)
fl=(res['left'][0]+res['right'][0])/2.0 if res['left'][1] and res['right'][1] else float('nan')
ratio=res['sv'][0]/fl if fl and not math.isnan(fl) else float('nan')
sid=bam.split('/')[-1].replace('.dedup.bam','')
print('\t'.join(map(str,[sid,bam,chrom,start,end,flank,res['left'][0],res['sv'][0],res['right'][0],fl,ratio,res['left'][1],res['sv'][1],res['right'][1]])))
