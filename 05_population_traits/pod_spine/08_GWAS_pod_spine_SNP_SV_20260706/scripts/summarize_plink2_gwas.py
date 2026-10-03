#!/usr/bin/env python3
import sys, math, csv
inp, chrom, start, end, out_tsv, out_sum = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6]
rows=[]
with open(inp, newline='') as f:
    reader=csv.DictReader(f, delimiter=' ', skipinitialspace=True)
    # DictReader with delimiter space may produce empty keys if multiple spaces; fallback manual split
    if reader.fieldnames and None not in reader.fieldnames and len(reader.fieldnames)>5:
        fields=reader.fieldnames
        for r in reader:
            if not r: continue
            chrval=r.get('#CHROM') or r.get('CHROM')
            if chrval == chrom and start <= int(r['POS']) <= end:
                rows.append(r)
    else:
        f.seek(0)
        header=f.readline().strip().split()
        for line in f:
            if not line.strip(): continue
            vals=line.strip().split()
            if len(vals)<len(header): continue
            r=dict(zip(header, vals))
            chrval=r.get('#CHROM') or r.get('CHROM')
            if chrval == chrom and start <= int(r['POS']) <= end:
                rows.append(r)
if rows:
    for r in rows:
        try: p=float(r.get('P','nan'))
        except Exception: p=float('nan')
        r['minus_log10_p'] = 'inf' if p==0 else (f'{-math.log10(p):.8g}' if p>0 else 'nan')
    rows.sort(key=lambda r: (float(r.get('P','nan')) if r.get('P','nan')!='NA' else float('inf'), int(r['POS'])))
header=list(rows[0].keys()) if rows else []
with open(out_tsv,'w', newline='') as f:
    if rows:
        w=csv.DictWriter(f, delimiter='\t', fieldnames=header)
        w.writeheader(); w.writerows(rows)
with open(out_sum,'w') as f:
    f.write(f'input\t{inp}\n')
    f.write(f'window\t{chrom}:{start}-{end}\n')
    f.write(f'n_variants_in_window\t{len(rows)}\n')
    if rows:
        top=rows[0]
        f.write('top_variant\n')
        for k in ['#CHROM','CHROM','POS','ID','REF','ALT','A1','TEST','OBS_CT','OR','Z_STAT','P','minus_log10_p']:
            if k in top: f.write(f'{k}\t{top[k]}\n')
        gstart,gend=88979818,88982331
        sub=[r for r in rows if gstart-10000 <= int(r['POS']) <= gend+10000]
        f.write(f'n_variants_Chr23997_plus10kb\t{len(sub)}\n')
        if sub:
            st=sub[0]
            f.write('top_Chr23997_plus10kb\n')
            for k in ['#CHROM','CHROM','POS','ID','REF','ALT','A1','TEST','OBS_CT','OR','Z_STAT','P','minus_log10_p']:
                if k in st: f.write(f'{k}\t{st[k]}\n')
