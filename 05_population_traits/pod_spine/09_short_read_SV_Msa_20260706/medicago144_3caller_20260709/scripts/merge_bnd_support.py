#!/usr/bin/env python3
import sys,gzip,re,math
from pathlib import Path

def op(path): return gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)
def info_dict(s):
    d={}
    for p in s.split(';') if s and s!='.' else []:
        if '=' in p:
            k,v=p.split('=',1); d[k]=v
        else: d[p]=True
    return d
def mate_from_alt(alt):
    m=re.search(r'[][]([^:\[\]]+):(\d+)[][]', alt)
    if m: return m.group(1), int(m.group(2))
    return None,None
def recs(path, caller):
    for line in op(path):
        if not line or line.startswith('#'): continue
        c=line.rstrip('\n').split('\t')
        if len(c)<8: continue
        chrom,pos,rid,ref,alt,qual,filt,info=c[:8]
        d=info_dict(info); svt=d.get('SVTYPE','').upper()
        if svt not in ('BND','TRA'): continue
        chr2=d.get('CHR2') or d.get('CHROM2')
        end=d.get('END')
        a_chr,a_pos=mate_from_alt(alt)
        if not chr2 and a_chr: chr2=a_chr
        if not end and a_pos: end=str(a_pos)
        if not chr2: chr2=chrom
        if not end: end=pos
        try: p1=int(pos); p2=int(end)
        except Exception: continue
        key=(chrom,p1,chr2,p2)
        # order pair so reciprocal breakends can cluster
        if (key[2],key[3]) < (key[0],key[1]):
            key=(key[2],key[3],key[0],key[1])
        yield {'key':key,'caller':caller,'id':rid,'chrom':key[0],'pos':key[1],'chr2':key[2],'end':key[3]}

def close(a,b,dist):
    return a['chrom']==b['chrom'] and a['chr2']==b['chr2'] and abs(a['pos']-b['pos'])<=dist and abs(a['end']-b['end'])<=dist

if len(sys.argv)<5:
    sys.exit('usage: merge_bnd_support.py sample out.vcf max_dist caller=vcf ...')
sample,out,maxdist=sys.argv[1],sys.argv[2],int(sys.argv[3])
items=[]
for arg in sys.argv[4:]:
    caller,path=arg.split('=',1)
    if Path(path).exists(): items += list(recs(path, caller))
clusters=[]
for r in items:
    placed=False
    for cl in clusters:
        if close(cl[0],r,maxdist):
            cl.append(r); placed=True; break
    if not placed: clusters.append([r])
with open(out,'w') as fo:
    fo.write('##fileformat=VCFv4.2\n')
    for contig in sorted({x['chrom'] for x in items} | {x['chr2'] for x in items}):
        fo.write(f'##contig=<ID={contig}>\n')
    fo.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n')
    fo.write('##INFO=<ID=CHR2,Number=1,Type=String,Description="Mate chromosome">\n')
    fo.write('##INFO=<ID=POS2,Number=1,Type=Integer,Description="Mate breakpoint position">\n')
    fo.write('##INFO=<ID=SUPP,Number=1,Type=Integer,Description="Number of supporting callers">\n')
    fo.write('##INFO=<ID=CALLERS,Number=.,Type=String,Description="Supporting callers">\n')
    fo.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    fo.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample+'\n')
    n=0
    for cl in clusters:
        callers=sorted(set(x['caller'] for x in cl))
        if len(callers)<2: continue
        # use median breakpoints
        chrom=cl[0]['chrom']; chr2=cl[0]['chr2']
        pos=sorted(x['pos'] for x in cl)[len(cl)//2]
        end=sorted(x['end'] for x in cl)[len(cl)//2]
        n+=1
        rid=f'{sample}_BND_2caller_{n}'
        info=f'SVTYPE=BND;CHR2={chr2};POS2={end};SUPP={len(callers)};CALLERS={",".join(callers)}'
        fo.write(f'{chrom}\t{pos}\t{rid}\tN\t<BND>\t.\tPASS\t{info}\tGT\t0/1\n')
