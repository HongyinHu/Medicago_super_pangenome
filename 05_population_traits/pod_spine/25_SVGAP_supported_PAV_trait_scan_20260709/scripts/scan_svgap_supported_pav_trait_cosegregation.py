#!/usr/bin/env python3
import csv, math, os
from collections import Counter

SRC = "path/to/project/N_3.call_SV/06.integrated_panSV_Msa_paper_style_unified/results/paper_style_unified_panSV.Msa_ref.PAV.matrix.tsv"
OUT = "path/to/project/N_4.pod_spiny/25_SVGAP_supported_PAV_trait_scan_20260709"
CHROMS = [f"Chr{i}" for i in range(1, 9)]
# 06.integrated_panSV_Msa_paper_style_unified intentionally excludes genome_A17; genome_R108 represents M. truncatula.
INTACT_GROUP = ["genome_410", "genome_436", "genome_457", "genome_454", "genome_474", "genome_Mpo", "genome_R108"]
DEL_GROUP = ["genome_395", "genome_461", "genome_468", "genome_472", "genome_M22", "genome_Mar", "genome_Mru"]
TARGET_CHROM = "Chr4"
TARGET_POS = 88_980_831
TARGET_END = 88_981_052
TARGET_LEN = 221
TARGET_ID = "Chr23997_second_intron_DEL_manual_IGV_corrected"

def fisher_two_sided(a,b,c,d):
    n1=a+b; m1=a+c; n=a+b+c+d
    def lc(nv,kv):
        if kv<0 or kv>nv: return float('-inf')
        return math.lgamma(nv+1)-math.lgamma(kv+1)-math.lgamma(nv-kv+1)
    def prob(x):
        if x<0 or x>n1 or x>m1 or (n1-x)>(n-m1): return 0.0
        return math.exp(lc(m1,x)+lc(n-m1,n1-x)-lc(n,n1))
    obs=prob(a); lo=max(0,n1-(n-m1)); hi=min(n1,m1)
    return min(1.0, sum(prob(x) for x in range(lo,hi+1) if prob(x)<=obs+1e-15))

def gt_present(v):
    return str(v).strip() == "1"

def target_like(row):
    if row['CHROM'] != TARGET_CHROM or row['SVTYPE'] != 'DEL':
        return False
    try:
        pos=int(row['POS']); end=int(row['END']); svlen=abs(int(float(row['SVLEN'])))
    except Exception:
        return False
    start=min(pos,end); stop=max(pos,end)
    ratio=min(max(svlen,1),TARGET_LEN)/max(max(svlen,1),TARGET_LEN)
    if ratio < 0.5:
        return False
    bp_ok = abs(start - TARGET_POS) <= 500 and abs(stop - TARGET_END) <= 500
    ov=max(0,min(stop,TARGET_END)-max(start,TARGET_POS)+1)
    ro=min(ov/max(1,stop-start+1), ov/max(1,TARGET_END-TARGET_POS+1))
    return bp_ok or ro >= 0.5

os.makedirs(os.path.join(OUT,'results'), exist_ok=True)
os.makedirs(os.path.join(OUT,'summary'), exist_ok=True)
rows=[]; n_total=0; n_hc=0; svtype_counts=Counter(); target_rows=[]
with open(SRC) as f:
    r=csv.DictReader(f, delimiter='\t')
    missing=[s for s in INTACT_GROUP+DEL_GROUP if s not in r.fieldnames]
    if missing:
        raise SystemExit(f"Missing samples in matrix: {missing}")
    for row in r:
        n_total += 1
        if row.get('CHROM') not in CHROMS:
            continue
        if row.get('final_confidence') != 'read_svgap_supported':
            continue
        n_hc += 1
        svtype_counts[row.get('SVTYPE','NA')] += 1
        ip=sum(gt_present(row[s]) for s in INTACT_GROUP)
        dp=sum(gt_present(row[s]) for s in DEL_GROUP)
        ia=len(INTACT_GROUP)-ip
        da=len(DEL_GROUP)-dp
        ifreq=ip/len(INTACT_GROUP)
        dfreq=dp/len(DEL_GROUP)
        if dfreq>ifreq:
            direction='DEL_group_enriched'
        elif ifreq>dfreq:
            direction='INTACT_group_enriched'
        else:
            direction='tie'
        p=fisher_two_sided(dp,da,ip,ia)
        outrow={
            'FinalSV_ID':row['FinalSV_ID'],'CHROM':row['CHROM'],'POS':row['POS'],'END':row['END'],'SVTYPE':row['SVTYPE'],'SVLEN':row['SVLEN'],
            'read_SV_ID':row.get('read_SV_ID','.'),'svgap_ids':row.get('svgap_ids','.'),'final_confidence':row.get('final_confidence','.'),
            'intact_present':ip,'intact_absent':ia,'del_present':dp,'del_absent':da,
            'intact_freq':ifreq,'del_freq':dfreq,'direction':direction,'fisher_p':p,'neglog10p':-math.log10(max(p,1e-300)),
            'complete_segregation': ((dp==len(DEL_GROUP) and ip==0) or (ip==len(INTACT_GROUP) and dp==0)),
            'group_specific80': ((dfreq>=0.8 and ifreq<=0.2) or (ifreq>=0.8 and dfreq<=0.2)),
            'target_like_Chr23997': target_like(row),
        }
        rows.append(outrow)
        if outrow['target_like_Chr23997']:
            target_rows.append(outrow)
fields=['FinalSV_ID','CHROM','POS','END','SVTYPE','SVLEN','read_SV_ID','svgap_ids','final_confidence','intact_present','intact_absent','del_present','del_absent','intact_freq','del_freq','direction','fisher_p','neglog10p','complete_segregation','group_specific80','target_like_Chr23997']
res=os.path.join(OUT,'results/svgap_supported_PAV_trait_cosegregation.tsv')
with open(res,'w',newline='') as f:
    w=csv.DictWriter(f, delimiter='\t', fieldnames=fields)
    w.writeheader(); w.writerows(rows)
with open(os.path.join(OUT,'results/Chr23997_target_like_svgap_supported_PAV.tsv'),'w',newline='') as f:
    w=csv.DictWriter(f, delimiter='\t', fieldnames=fields)
    w.writeheader(); w.writerows(target_rows)
p_manual=fisher_two_sided(len(DEL_GROUP),0,0,len(INTACT_GROUP))
with open(os.path.join(OUT,'results/Chr23997_manual_IGV_corrected_target_point.tsv'),'w',newline='') as f:
    w=csv.writer(f, delimiter='\t')
    w.writerow(['id','chrom','pos','end','svtype','svlen','intact_present','intact_absent','del_present','del_absent','fisher_p_core','neglog10p_core','note'])
    w.writerow([TARGET_ID,TARGET_CHROM,TARGET_POS,TARGET_END,'DEL',TARGET_LEN,0,len(INTACT_GROUP),len(DEL_GROUP),0,p_manual,-math.log10(max(p_manual,1e-300)),'Manual IGV-corrected Chr23997 intron-2 DEL; high-confidence PAV scan excludes genome_A17 because 06 final matrix excludes A17.'])
with open(os.path.join(OUT,'summary/run_summary.tsv'),'w',newline='') as f:
    w=csv.writer(f, delimiter='\t')
    w.writerow(['item','value'])
    w.writerow(['analysis','Msa_ref_read_svgap_supported_PAV_trait_cosegregation'])
    w.writerow(['source_matrix',SRC])
    w.writerow(['total_input_PAV',n_total])
    w.writerow(['read_svgap_supported_PAV_chr1_8',n_hc])
    w.writerow(['intact_group',','.join(INTACT_GROUP)])
    w.writerow(['del_group',','.join(DEL_GROUP)])
    w.writerow(['excluded_by_matrix_design','genome_A17 absent from 06 integrated final matrix; genome_Msa self; genome_M46/genome_482/genome_ZM4 not used for phenotype contrast'])
    w.writerow(['group_specific80_sites',sum(1 for x in rows if x['group_specific80'])])
    w.writerow(['complete_segregation_sites',sum(1 for x in rows if x['complete_segregation'])])
    w.writerow(['p_lt_0.05_sites',sum(1 for x in rows if x['fisher_p']<0.05)])
    w.writerow(['target_like_Chr23997_sites',len(target_rows)])
    w.writerow(['Chr23997_manual_core_p',p_manual])
    w.writerow(['Chr23997_manual_core_neglog10p',-math.log10(max(p_manual,1e-300))])
    for k,v in svtype_counts.most_common():
        w.writerow([f'svtype_{k}',v])
print(f"high_confidence_pav={n_hc} rows={len(rows)} target_like={len(target_rows)} out={res}")
