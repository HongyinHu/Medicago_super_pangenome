#!/usr/bin/env python3
import csv, math, os, sys
from pathlib import Path
from collections import defaultdict

RUN=Path('path/to/project/N_4.pod_spiny/18_Chr23997_genomewide_cosegregation_rank_20260708')
MATRIX=Path('path/to/project/N_3.call_SV/06.integrated_panSV_Msa_paper_style_unified/results/paper_style_unified_panSV.Msa_ref.PAV.matrix.tsv')
EVENTS=Path('path/to/project/N_3.call_SV/06.integrated_panSV_Msa_paper_style_unified/results/paper_style_unified_panSV.Msa_ref.events.tsv')
CHR23997_RECODE=Path('path/to/project/N_4.pod_spiny/03_candidate_gene_sv_cosegregation_Chr23997/results/Chr23997.intron23_DEL.read_based_discovery_vs_genotype.manual_recode.tsv')
CHR23997_REGION=("Chr4",88979818,88982331)
CHR23997_INTRON23_EVENT_ID='0_0_pbsv.DEL.130205'

MODELS={
    # Conservative model: unambiguous spiny vs unambiguous spineless, excludes hairy/weak/low-quality/uncertain/reference.
    'core_clear': {
        'spiny':['genome_410','genome_Mpo','genome_R108'],
        'spineless':['genome_395','genome_461','genome_468','genome_472','genome_M22','genome_Mar','genome_Mru','genome_ZM4'],
        'exclude_note':'exclude genome_436, genome_457, genome_454, genome_474, genome_482, genome_M46, genome_Msa; A17 absent from final Msa panSV matrix'
    },
    # User biological model: weak/short-spine/hairy/small-indel accessions treated as spiny-like/intact after manual IGV review.
    'spiny_like_inclusive_rawPAV': {
        'spiny':['genome_410','genome_474','genome_436','genome_457','genome_454','genome_Mpo','genome_R108'],
        'spineless':['genome_395','genome_461','genome_468','genome_472','genome_M22','genome_Mar','genome_Mru','genome_ZM4'],
        'exclude_note':'exclude genome_482 uncertain, genome_M46 low quality, genome_Msa reference; A17 absent from final Msa panSV matrix'
    }
}

def present_value(x):
    x=(x or '').strip()
    if x in ('NA','.','./.',''):
        return None
    if x in ('1','1.0'):
        return 1
    if x in ('0','0.0'):
        return 0
    # fall back for GT-like values if encountered
    if x in ('0/1','1/0','1/1','0|1','1|0','1|1'):
        return 1
    if x in ('0/0','0|0'):
        return 0
    return None

def fisher_two_sided(a,b,c,d):
    # table [[a,b],[c,d]] with fixed margins; two-sided by probability <= observed.
    n=a+b+c+d
    if n==0: return 1.0
    row1=a+b; row2=c+d; col1=a+c
    def prob(x):
        return math.comb(col1,x)*math.comb(n-col1,row1-x)/math.comb(n,row1)
    lo=max(0,row1-(n-col1)); hi=min(row1,col1)
    pobs=prob(a)
    return min(1.0, sum(prob(x) for x in range(lo,hi+1) if prob(x) <= pobs + 1e-15))

def bh_qvalues(rows):
    indexed=[(i,r['fisher_p']) for i,r in enumerate(rows) if r['fisher_p'] is not None]
    indexed.sort(key=lambda x:x[1])
    m=len(indexed)
    q=[None]*len(rows)
    prev=1.0
    for rank,(i,p) in enumerate(reversed(indexed), start=1):
        k=m-rank+1
        val=min(prev, p*m/k)
        prev=val
        q[i]=min(1.0,val)
    for i,v in enumerate(q): rows[i]['fdr_bh_q']=v

def score_row(row, spiny, spineless):
    sp=[]; sl=[]
    for s in spiny:
        if s in row: sp.append(present_value(row[s]))
    for s in spineless:
        if s in row: sl.append(present_value(row[s]))
    sp_p=sum(1 for x in sp if x==1); sp_a=sum(1 for x in sp if x==0); sp_m=sum(1 for x in sp if x is None)
    sl_p=sum(1 for x in sl if x==1); sl_a=sum(1 for x in sl if x==0); sl_m=sum(1 for x in sl if x is None)
    sp_n=sp_p+sp_a; sl_n=sl_p+sl_a
    sp_rate=sp_p/sp_n if sp_n else None
    sl_rate=sl_p/sl_n if sl_n else None
    if sp_rate is None or sl_rate is None:
        delta=None; p=1.0; direction='NA'
    else:
        delta=sp_rate-sl_rate
        p=fisher_two_sided(sp_p,sp_a,sl_p,sl_a)
        if delta>0: direction='spiny_present_enriched'
        elif delta<0: direction='spineless_present_enriched'
        else: direction='no_delta'
    nonmissing=sp_n+sl_n
    abs_delta=abs(delta) if delta is not None else 0.0
    # Tier-like score: emphasize effect size, completeness, then p-value, then SVGAP evidence.
    pscore=-math.log10(max(p,1e-300))
    support_bonus=0.15 if row.get('final_confidence','') != 'read_based' else 0.0
    score=abs_delta*100 + min(nonmissing,20)*1.5 + min(pscore,50)*2 + support_bonus
    return dict(spiny_present=sp_p,spiny_absent=sp_a,spiny_missing=sp_m,spineless_present=sl_p,spineless_absent=sl_a,spineless_missing=sl_m,spiny_rate=sp_rate,spineless_rate=sl_rate,delta_spiny_minus_spineless=delta,abs_delta=abs_delta,fisher_p=p,nonmissing=nonmissing,direction=direction,score=score)

def read_events_by_read_id():
    m={}
    with EVENTS.open() as f:
        r=csv.DictReader(f,delimiter='\t')
        for row in r:
            m[row['read_SV_ID']]=row
    return m

def rank_model(model_name, cfg):
    rows=[]; chr_hits=[]
    with MATRIX.open() as f:
        reader=csv.DictReader(f, delimiter='\t')
        for row in reader:
            sc=score_row(row,cfg['spiny'],cfg['spineless'])
            out={k:row.get(k,'') for k in ['FinalSV_ID','source_methods','reference','CHROM','POS','END','SVTYPE','SVLEN','read_SV_ID','svgap_ids','final_confidence']}
            out.update(sc)
            # Flags for Chr23997 region and target intron23 event.
            try:
                pos=int(row['POS']); end=int(row['END'])
            except Exception:
                pos=end=-1
            overlap_chr23997=(row.get('CHROM')==CHR23997_REGION[0] and end>=CHR23997_REGION[1] and pos<=CHR23997_REGION[2])
            out['overlap_Chr23997_gene']='yes' if overlap_chr23997 else 'no'
            out['is_Chr23997_intron23_large_DEL']='yes' if row.get('read_SV_ID')==CHR23997_INTRON23_EVENT_ID else 'no'
            rows.append(out)
            if overlap_chr23997 or row.get('read_SV_ID')==CHR23997_INTRON23_EVENT_ID:
                chr_hits.append(out)
    bh_qvalues(rows)
    rows.sort(key=lambda r:(r['fisher_p'], -r['abs_delta'], -r['nonmissing'], -r['score']))
    for i,r in enumerate(rows, start=1): r['rank_by_p']=i
    rows_by_score=sorted(rows, key=lambda r:(-r['score'], r['fisher_p']))
    score_rank={id(r):i for i,r in enumerate(rows_by_score, start=1)}
    for r in rows: r['rank_by_score']=score_rank[id(r)]
    out=RUN/'results'/f'{model_name}.genomewide_sv_trait_rank.tsv'
    fields=['rank_by_p','rank_by_score','FinalSV_ID','read_SV_ID','CHROM','POS','END','SVTYPE','SVLEN','source_methods','final_confidence','spiny_present','spiny_absent','spiny_missing','spineless_present','spineless_absent','spineless_missing','spiny_rate','spineless_rate','delta_spiny_minus_spineless','abs_delta','direction','fisher_p','fdr_bh_q','nonmissing','score','overlap_Chr23997_gene','is_Chr23997_intron23_large_DEL','svgap_ids']
    with out.open('w') as fo:
        w=csv.DictWriter(fo,delimiter='\t',fieldnames=fields,extrasaction='ignore')
        w.writeheader(); w.writerows(rows)
    top=RUN/'results'/f'{model_name}.top200_by_p.tsv'
    with top.open('w') as fo:
        w=csv.DictWriter(fo,delimiter='\t',fieldnames=fields,extrasaction='ignore')
        w.writeheader(); w.writerows(rows[:200])
    chrout=RUN/'results'/f'{model_name}.Chr23997_gene_region_rows.tsv'
    chr_hits_sorted=sorted(chr_hits, key=lambda r:(r['fisher_p'], -r['abs_delta'], -r['nonmissing']))
    with chrout.open('w') as fo:
        w=csv.DictWriter(fo,delimiter='\t',fieldnames=fields,extrasaction='ignore')
        w.writeheader(); w.writerows(chr_hits_sorted)
    target=[r for r in rows if r.get('is_Chr23997_intron23_large_DEL')=='yes']
    return rows, target

def manual_chr23997_score():
    # Manual IGV/read-level recoding for the intron2-3 large DEL: DEL=1, INTACT=0, uncertain/excluded/reference not used.
    calls={}
    with CHR23997_RECODE.open() as f:
        r=csv.DictReader(f,delimiter='\t')
        for row in r:
            call=row['recommended_event_call']
            if call=='DEL': v=1
            elif call.startswith('INTACT'): v=0
            else: v=None
            calls[row['sample']]=v
    models={
        'manual_core_clear': MODELS['core_clear'],
        'manual_spiny_like_inclusive': MODELS['spiny_like_inclusive_rawPAV'],
    }
    rows=[]
    for name,cfg in models.items():
        row={s:('NA' if calls.get(s) is None else str(calls.get(s))) for s in set(cfg['spiny']+cfg['spineless'])}
        sc=score_row(row,cfg['spiny'],cfg['spineless'])
        out={'model':name,'event':'Chr23997_intron2_3_large_DEL_manual_IGV_recode','CHROM':'Chr4','POS':'88980831','END':'88981052','SVTYPE':'DEL','SVLEN':'-221','note':'manual read-level recode from IGV/depth; 474/436/457 treated intact despite raw genotype overcall; 482 missing; M46/Msa excluded'}
        out.update(sc)
        rows.append(out)
    out=RUN/'results'/'Chr23997.manual_IGV_recode_score.tsv'
    fields=['model','event','CHROM','POS','END','SVTYPE','SVLEN','spiny_present','spiny_absent','spiny_missing','spineless_present','spineless_absent','spineless_missing','spiny_rate','spineless_rate','delta_spiny_minus_spineless','abs_delta','direction','fisher_p','nonmissing','score','note']
    with out.open('w') as fo:
        w=csv.DictWriter(fo,delimiter='\t',fieldnames=fields,extrasaction='ignore')
        w.writeheader(); w.writerows(rows)
    return rows

def main():
    RUN.joinpath('results').mkdir(exist_ok=True)
    RUN.joinpath('summary').mkdir(exist_ok=True)
    summary=[]
    for model,cfg in MODELS.items():
        rows,target=rank_model(model,cfg)
        summary.append((model,cfg,target))
    manual=manual_chr23997_score()
    with (RUN/'summary'/'analysis_summary.tsv').open('w') as fo:
        fo.write('section\tkey\tvalue\n')
        fo.write(f'input\tmatrix\t{MATRIX}\n')
        fo.write(f'input\tchr23997_manual_recode\t{CHR23997_RECODE}\n')
        for model,cfg,target in summary:
            fo.write(f'model\t{model}.spiny\t{",".join(cfg["spiny"])}\n')
            fo.write(f'model\t{model}.spineless\t{",".join(cfg["spineless"])}\n')
            fo.write(f'model\t{model}.exclude_note\t{cfg["exclude_note"]}\n')
            if target:
                t=target[0]
                fo.write(f'chr23997_raw\t{model}.rank_by_p\t{t["rank_by_p"]}\n')
                fo.write(f'chr23997_raw\t{model}.rank_by_score\t{t["rank_by_score"]}\n')
                fo.write(f'chr23997_raw\t{model}.spiny_rate\t{t["spiny_rate"]}\n')
                fo.write(f'chr23997_raw\t{model}.spineless_rate\t{t["spineless_rate"]}\n')
                fo.write(f'chr23997_raw\t{model}.fisher_p\t{t["fisher_p"]}\n')
                fo.write(f'chr23997_raw\t{model}.fdr_bh_q\t{t["fdr_bh_q"]}\n')
            else:
                fo.write(f'chr23997_raw\t{model}.target_found\tno\n')
        for m in manual:
            fo.write(f'chr23997_manual\t{m["model"]}.spiny_rate\t{m["spiny_rate"]}\n')
            fo.write(f'chr23997_manual\t{m["model"]}.spineless_rate\t{m["spineless_rate"]}\n')
            fo.write(f'chr23997_manual\t{m["model"]}.fisher_p\t{m["fisher_p"]}\n')
            fo.write(f'chr23997_manual\t{m["model"]}.score\t{m["score"]}\n')
    with (RUN/'summary'/'README.md').open('w') as fo:
        fo.write('''# Chr23997 genome-wide SV-trait cosegregation ranking\n\nPurpose: test whether the Chr23997 intron2-3 large deletion stands out under a reproducible SV-phenotype cosegregation screen.\n\nInputs:\n- Msa single-reference integrated panSV matrix: read-mapping-primary events with SVGAP used only as supporting evidence; SVGAP-only events removed.\n- Phenotype models:\n  - `core_clear`: only unambiguous spiny and spineless species.\n  - `spiny_like_inclusive_rawPAV`: includes genome_436, genome_457, genome_454 and genome_474 as spiny-like/intact group, but uses raw PAV matrix.\n- Manual Chr23997 recode: IGV/read-level correction of the intron2-3 large deletion.\n\nStatistics:\nFor each SV, presence rates were compared between spiny and spineless groups. Fisher exact test was applied to the 2x2 table: spiny present/absent vs spineless present/absent. Both directions were allowed. BH-FDR was calculated within the genome-wide matrix.\n\nImportant interpretation:\nThe raw genotyped PAV matrix is intentionally reported separately from the manual read-level recode. Chr23997 is expected to be penalized in raw PAV because genome_436, genome_457 and genome_474 were overcalled as DEL in the genotyped matrix, while IGV/depth inspection supports intact or only small-indel status.\n''')

if __name__=='__main__': main()
