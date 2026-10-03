#!/usr/bin/env python3
import sys, gzip, re
from pathlib import Path

def open_text(path):
    path=str(path)
    return gzip.open(path,'rt') if path.endswith('.gz') else open(path)

def info_get(info, key):
    for part in info.split(';'):
        if part.startswith(key+'='):
            return part.split('=',1)[1]
    return None

def set_info(info, updates):
    parts=[]
    seen=set()
    if info and info!='.':
        for p in info.split(';'):
            if not p: continue
            k=p.split('=',1)[0]
            if k in updates:
                if updates[k] is not None: parts.append(f'{k}={updates[k]}')
                seen.add(k)
            else:
                parts.append(p)
    for k,v in updates.items():
        if k not in seen and v is not None:
            parts.append(f'{k}={v}')
    return ';'.join(parts) if parts else '.'

def infer_svtype(alt, info):
    svt=info_get(info,'SVTYPE')
    if svt: return svt.upper()
    m=re.search(r'<([^>]+)>', alt)
    if m: return m.group(1).upper()
    if '[' in alt or ']' in alt: return 'BND'
    return None

def parse_svlen(info, pos, svt):
    v=info_get(info,'SVLEN')
    if v:
        try: return max(abs(int(float(x))) for x in v.split(','))
        except Exception: pass
    end=info_get(info,'END')
    if end:
        try: return abs(int(end)-int(pos))+1
        except Exception: pass
    return None

if len(sys.argv) != 5:
    sys.exit('usage: normalize_vcf.py in.vcf[.gz] caller sample out.vcf')
infile, caller, sample, outfile = sys.argv[1:]
keep_types={'DEL','INS','DUP','INV','BND','TRA'}
header_extra=[
'##INFO=<ID=CALLER,Number=1,Type=String,Description="Source short-read SV caller">',
'##INFO=<ID=ORIG_ID,Number=1,Type=String,Description="Original caller record ID">',
'##INFO=<ID=ORIG_SVTYPE,Number=1,Type=String,Description="Original caller SVTYPE">'
]
records=0
with open_text(infile) as fi, open(outfile,'w') as fo:
    inserted=False
    for line in fi:
        if line.startswith('##'):
            fo.write(line)
            continue
        if line.startswith('#CHROM'):
            for h in header_extra: fo.write(h+'\n')
            cols=line.rstrip('\n').split('\t')
            if len(cols) >= 10:
                cols[9]=sample
            fo.write('\t'.join(cols)+'\n')
            inserted=True
            continue
        if not line.strip():
            continue
        cols=line.rstrip('\n').split('\t')
        if len(cols) < 8: continue
        chrom,pos,rid,ref,alt,qual,filt,info=cols[:8]
        if filt not in ('PASS','.'):
            continue
        svt=infer_svtype(alt, info)
        if not svt: continue
        orig_svt=svt
        if svt in ('TRA','CTX'): svt='BND'
        if svt not in keep_types: continue
        svlen=parse_svlen(info, pos, svt)
        if svt not in ('BND','TRA') and (svlen is None or svlen < 50):
            continue
        # Keep BND as BND. Downstream stats can label it as Translocation.
        updates={'CALLER':caller,'ORIG_ID':rid.replace(';',','),'ORIG_SVTYPE':orig_svt,'SVTYPE':svt}
        if svt in ('BND','TRA') and info_get(info,'END') is None:
            updates['END']=pos
        cols[2]=f'{sample}_{caller}_{records+1}'
        cols[6]='PASS'
        cols[7]=set_info(info, updates)
        if len(cols) >= 10:
            cols[9]='0/1'
        fo.write('\t'.join(cols)+'\n')
        records += 1
if records == 0:
    Path(outfile).with_suffix(Path(outfile).suffix+'.empty').write_text('no_records_after_filter\n')
