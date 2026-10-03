#!/usr/bin/env python3
import os, gzip, math, re, csv
from collections import defaultdict, Counter
from statistics import median

RUN = "path/to/project/N_3.call_SV/01.read_based_dualref_hifi"
OUT = "path/to/project/N_4.pod_spiny/22_readmapping_rawSV_trait_cosegregation_scan_20260708"
FAI = os.path.join(RUN, "02_ref_prepare/Msa/ref.fa.fai")
CHROMS = [f"Chr{i}" for i in range(1, 9)]
MIN_SVLEN = 50
BIN = 500
LEN_RATIO_MIN = 0.50
INTACT_GROUP = ["genome_410", "genome_436", "genome_457", "genome_454", "genome_474", "genome_Mpo", "genome_R108", "genome_A17"]
DEL_GROUP = ["genome_395", "genome_461", "genome_468", "genome_472", "genome_M22", "genome_Mar", "genome_Mru"]
SAMPLES = INTACT_GROUP + DEL_GROUP
CALLERS = {
    "pbsv": os.path.join(RUN, "03_per_sample/Msa/pbsv/{sample}.pbsv.vcf.gz"),
    "sniffles2": os.path.join(RUN, "03_per_sample/Msa/sniffles2/{sample}.sniffles2.vcf.gz"),
    "cuteSV": os.path.join(RUN, "03_per_sample/Msa/cutesv/{sample}.cutesv.vcf.gz"),
}
TARGET_CHROM = "Chr4"
TARGET_POS = 88_980_831
TARGET_END = 88_981_052
TARGET_LEN = 221
TARGET_ID = "Chr23997_second_intron_DEL_manual_IGV_corrected"

def open_vcf(path):
    return gzip.open(path, 'rt', errors='replace') if path.endswith('.gz') else open(path, 'rt', errors='replace')

def parse_info(info):
    d = {}
    for item in info.split(';'):
        if not item: continue
        if '=' in item:
            k, v = item.split('=', 1); d[k] = v
        else:
            d[item] = True
    return d

def first_int(v):
    if v is None: return None
    m = re.search(r'-?\d+', str(v).split(',')[0])
    return int(m.group(0)) if m else None

def parse_svtype(info, alt):
    svt = info.get('SVTYPE')
    if not svt:
        m = re.search(r'<([^>]+)>', alt)
        svt = m.group(1) if m else None
    if not svt and ('[' in alt or ']' in alt): svt = 'BND'
    if not svt: return 'UNK'
    svt = svt.upper().replace('SNV', 'UNK')
    if svt in ('TRA','TRANSLOCATION'): return 'TRA'
    return svt

def parse_svlen(info, pos, end, svtype):
    vals = []
    if 'SVLEN' in info:
        for x in str(info['SVLEN']).split(','):
            try: vals.append(abs(int(float(x))))
            except Exception: pass
    if vals: return max(vals)
    if end and end >= pos and svtype not in ('INS','BND','TRA'):
        return abs(end - pos + 1)
    return 0

def parse_rec(line, sample, caller):
    a = line.rstrip('\n').split('\t')
    if len(a) < 8: return None
    chrom, pos_s, rid, ref, alt, qual, flt, info_s = a[:8]
    if chrom not in CHROMS or flt not in ('.','PASS'): return None
    try: pos = int(pos_s)
    except Exception: return None
    info = parse_info(info_s)
    svtype = parse_svtype(info, alt)
    if svtype in ('UNK','SNV'): return None
    end = first_int(info.get('END'))
    svlen = parse_svlen(info, pos, end, svtype)
    if svtype == 'INS':
        if end is None: end = pos
    elif svtype in ('BND','TRA'):
        end = pos
        svlen = max(svlen, 50)
    elif end is None:
        if not svlen: return None
        end = pos + svlen - 1
    if svtype not in ('BND','TRA') and svlen < MIN_SVLEN: return None
    start, stop = min(pos,end), max(pos,end)
    return {'sample':sample,'caller':caller,'chrom':chrom,'pos':pos,'end':end,'start':start,'stop':stop,'svtype':svtype,'svlen':abs(svlen)}

def key_for(r):
    if r['svtype'] in ('INS','BND','TRA'):
        return (r['chrom'], r['svtype'], int(round(r['pos']/BIN)))
    # Fast genome-wide clustering for plotting: same type, same approximate breakpoints, and broad length class.
    # The final local candidate interpretation still uses the stricter pairwise clustering and manual IGV correction.
    len_class = int(round(math.log10(max(r['svlen'], 50)) * 10))
    return (r['chrom'], r['svtype'], int(round(r['start']/BIN)), int(round(r['stop']/BIN)), len_class)

def target_like(r):
    if r['chrom'] != TARGET_CHROM or r['svtype'] != 'DEL': return False
    ratio = min(max(r['svlen'],1), TARGET_LEN) / max(max(r['svlen'],1), TARGET_LEN)
    if ratio < LEN_RATIO_MIN: return False
    bp_ok = abs(r['start'] - TARGET_POS) <= BIN and abs(r['stop'] - TARGET_END) <= BIN
    ov = max(0, min(r['stop'], TARGET_END) - max(r['start'], TARGET_POS) + 1)
    ro = min(ov / max(1, r['stop']-r['start']+1), ov / max(1, TARGET_END-TARGET_POS+1))
    return bp_ok or ro >= 0.5

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

def main():
    os.makedirs(os.path.join(OUT,'results'), exist_ok=True)
    os.makedirs(os.path.join(OUT,'summary'), exist_ok=True)
    groups = {}
    status=[]; raw_n=0
    for sample in SAMPLES:
        for caller, tmpl in CALLERS.items():
            path=tmpl.format(sample=sample); n=0
            if not os.path.exists(path):
                status.append([sample,caller,path,'missing',0]); continue
            with open_vcf(path) as fh:
                for line in fh:
                    if not line or line.startswith('#'): continue
                    r=parse_rec(line,sample,caller)
                    if not r: continue
                    raw_n += 1; n += 1
                    k=key_for(r)
                    g=groups.get(k)
                    if g is None:
                        g={'chrom':r['chrom'],'svtype':r['svtype'],'starts':[],'stops':[],'poss':[],'ends':[],'lens':[],'records':0,'sample_callers':defaultdict(set),'target_like':False}
                        groups[k]=g
                    g['starts'].append(r['start']); g['stops'].append(r['stop']); g['poss'].append(r['pos']); g['ends'].append(r['end']); g['lens'].append(r['svlen']); g['records'] += 1
                    g['sample_callers'][sample].add(caller)
                    if target_like(r): g['target_like'] = True
            status.append([sample,caller,path,'ok',n])
    with open(os.path.join(OUT,'summary/raw_vcf_file_status_genomewide_fast.tsv'),'w',newline='') as f:
        w=csv.writer(f, delimiter='\t'); w.writerow(['sample','caller','path','status','records_chr1_8_pass_svlen50']); w.writerows(status)
    rows=[]
    intact_set=set(INTACT_GROUP); del_set=set(DEL_GROUP)
    for g in groups.values():
        samples=set(g['sample_callers'])
        ip=len(samples & intact_set); dp=len(samples & del_set)
        ia=len(INTACT_GROUP)-ip; da=len(DEL_GROUP)-dp
        ifreq=ip/len(INTACT_GROUP); dfreq=dp/len(DEL_GROUP)
        if dfreq>ifreq: direction='DEL_group_enriched'
        elif ifreq>dfreq: direction='INTACT_group_enriched'
        else: direction='tie'
        p=fisher_two_sided(dp,da,ip,ia)
        callers_by_sample=';'.join(f"{s}:{','.join(sorted(g['sample_callers'][s]))}" for s in sorted(samples))
        rows.append({'cluster_id':'','chrom':g['chrom'],'pos':int(median(g['poss'])),'end':int(median(g['ends'])),'start':min(g['starts']),'stop':max(g['stops']),'svtype':g['svtype'],'svlen_median':int(median(g['lens'])),'n_records':g['records'],'n_samples_present':len(samples),'samples_present':','.join(sorted(samples)),'callers_by_sample':callers_by_sample,'intact_present':ip,'intact_absent':ia,'del_present':dp,'del_absent':da,'intact_freq':ifreq,'del_freq':dfreq,'direction':direction,'fisher_p':p,'neglog10p':-math.log10(max(p,1e-300)),'complete_segregation':((dp==len(DEL_GROUP) and ip==0) or (ip==len(INTACT_GROUP) and dp==0)),'group_specific80':((dfreq>=0.8 and ifreq<=0.2) or (ifreq>=0.8 and dfreq<=0.2)),'target_like_raw_cluster':g['target_like']})
    order={c:i for i,c in enumerate(CHROMS)}
    rows.sort(key=lambda r:(order.get(r['chrom'],999), r['pos'], r['svtype']))
    for i,r in enumerate(rows,1): r['cluster_id']=f"RawSVfast_{i:07d}"
    fields=['cluster_id','chrom','pos','end','start','stop','svtype','svlen_median','n_records','n_samples_present','samples_present','callers_by_sample','intact_present','intact_absent','del_present','del_absent','intact_freq','del_freq','direction','fisher_p','neglog10p','complete_segregation','group_specific80','target_like_raw_cluster']
    out_table=os.path.join(OUT,'results/readmapping_rawSV_genomewide_fast.breakpoint_bin_clusters_with_fisher.tsv')
    with open(out_table,'w',newline='') as f:
        w=csv.DictWriter(f, delimiter='\t', fieldnames=fields); w.writeheader(); w.writerows(rows)
    p=fisher_two_sided(len(DEL_GROUP),0,0,len(INTACT_GROUP))
    with open(os.path.join(OUT,'results/Chr23997_manual_IGV_corrected_target_point.tsv'),'w',newline='') as f:
        w=csv.writer(f, delimiter='\t')
        w.writerow(['id','chrom','pos','end','svtype','svlen','intact_present','intact_absent','del_present','del_absent','fisher_p_core','neglog10p_core','note'])
        w.writerow([TARGET_ID,TARGET_CHROM,TARGET_POS,TARGET_END,'DEL',TARGET_LEN,0,len(INTACT_GROUP),len(DEL_GROUP),0,p,-math.log10(p),'Manual IGV-corrected Chr23997 intron-2 DEL'])
    with open(os.path.join(OUT,'summary/run_summary_genomewide_fast.tsv'),'w',newline='') as f:
        w=csv.writer(f, delimiter='\t')
        w.writerow(['item','value'])
        for k,v in [('analysis','raw_readmapping_per_species_no_panSV_genomewide_fast_breakpoint_bin'),('chromosomes',','.join(CHROMS)),('raw_records',raw_n),('fast_clusters',len(rows)),('group_specific80_clusters',sum(1 for r in rows if r['group_specific80'])),('complete_segregation_clusters',sum(1 for r in rows if r['complete_segregation'])),('p_lt_0.05_clusters',sum(1 for r in rows if r['fisher_p']<0.05)),('target_like_raw_clusters',sum(1 for r in rows if r['target_like_raw_cluster'])),('Chr23997_manual_core_p',p),('Chr23997_manual_core_neglog10p',-math.log10(p)),('cluster_table',out_table),('note','Genomewide plot uses fast 500-bp breakpoint-bin clustering for display; local Chr4 candidate plot uses strict pairwise matching plus manual IGV correction.')]:
            w.writerow([k,v])
    print(f"raw_records={raw_n} fast_clusters={len(rows)} out={out_table}")
if __name__ == '__main__': main()
