#!/usr/bin/env python3
import sys, math, csv, os

def read_gwas(path):
    if not os.path.exists(path) or os.path.getsize(path)==0:
        return []
    rows=[]
    with open(path) as f:
        header=f.readline().strip().split()
        idx={h:i for i,h in enumerate(header)}
        for line in f:
            if not line.strip(): continue
            parts=line.strip().split()
            try:
                p=parts[idx['P']]
                if p=='NA': continue
                p=float(p)
                pos=int(float(parts[idx['POS']]))
            except Exception:
                continue
            rows.append({
                'CHROM':parts[idx.get('#CHROM', idx.get('CHROM',0))],
                'POS':pos,
                'ID':parts[idx.get('ID',2)] if idx.get('ID',2)<len(parts) else '.',
                'REF':parts[idx.get('REF',3)] if idx.get('REF',3)<len(parts) else '.',
                'ALT':parts[idx.get('ALT',4)] if idx.get('ALT',4)<len(parts) else '.',
                'OBS_CT':parts[idx.get('OBS_CT',7)] if 'OBS_CT' in idx else '.',
                'OR':parts[idx.get('OR',8)] if 'OR' in idx else '.',
                'P':p,
                'minus_log10P': -math.log10(p) if p>0 else float('inf')
            })
    return rows

def best(rows, start, end):
    sub=[r for r in rows if start <= r['POS'] <= end]
    if not sub: return None,0
    return min(sub, key=lambda r: r['P']), len(sub)

def emit(out, dataset, model, rows):
    regions=[
        ('Chr4_88_90Mb', 88000000, 90000000),
        ('Chr23997_pm100kb', 88979818-100000, 88982331+100000),
        ('Chr23997_pm10kb', 88979818-10000, 88982331+10000),
        ('Chr23997_gene_body', 88979818, 88982331),
    ]
    for name,s,e in regions:
        r,n=best(rows,s,e)
        if r is None:
            out.writerow([dataset,model,name,s,e,0,'.','.','.','.','.','.','.','.','.'])
        else:
            out.writerow([dataset,model,name,s,e,n,r['CHROM'],r['POS'],r['ID'],r['REF'],r['ALT'],r['OBS_CT'],r['OR'],f"{r['P']:.6g}",f"{r['minus_log10P']:.6g}"])

if __name__ == '__main__':
    out=csv.writer(sys.stdout, delimiter='\t', lineterminator='\n')
    out.writerow(['dataset','model','region','start','end','n_tested','top_CHROM','top_POS','top_ID','REF','ALT','OBS_CT','OR','P','minus_log10P'])
    for spec in sys.argv[1:]:
        dataset,model,path=spec.split('=',2)
        emit(out,dataset,model,read_gwas(path))
