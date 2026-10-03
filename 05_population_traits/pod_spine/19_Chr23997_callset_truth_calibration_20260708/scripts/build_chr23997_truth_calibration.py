#!/usr/bin/env python3
import csv, math
from pathlib import Path

POD=Path('path/to/project/N_4.pod_spiny')
RUN=POD/'19_Chr23997_callset_truth_calibration_20260708'
SUPPORT=POD/'05_direction_free_sv_cosegregation_20260701/summary/Chr23997_per_sample_caller_support_20260701.tsv'
MANUAL=POD/'03_candidate_gene_sv_cosegregation_Chr23997/results/Chr23997.intron23_DEL.read_based_discovery_vs_genotype.manual_recode.tsv'
REGION=POD/'03_candidate_gene_sv_cosegregation_Chr23997/results/Chr23997.region_sv_gt_strings.tsv'
TARGET='0_0_pbsv.DEL.130205'
LOCAL_IGV='local_windows_IGV_manual_review_folder__SV_trait_cosegregation_IGV_check'

RUN.joinpath('results').mkdir(exist_ok=True)
RUN.joinpath('summary').mkdir(exist_ok=True)

# Curated interpretation used for the validation-first analysis.
CURATED_EXTRA={
    'genome_A17': {
        'phenotype_group_curated':'spiny_extra_A17_intact_manual',
        'read_based_discovery_SUPP_VEC':'NA',
        'genotyped_GT':'NA','DR':'NA','DV':'NA',
        'recommended_event_call':'INTACT_manual',
        'note':'extra A17 accession; local IGV snapshot in spiny group shows no large intron DEL; keep outside final 18-species Msa-reference statistics unless explicitly needed'
    }
}

support={}
with SUPPORT.open() as f:
    for r in csv.DictReader(f, delimiter='\t'):
        support[r['sample']]=r

manual={}
with MANUAL.open() as f:
    for r in csv.DictReader(f, delimiter='\t'):
        manual[r['sample']]=r
manual.update(CURATED_EXTRA)

target_gt={}
if REGION.exists():
    with REGION.open() as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r.get('sv_id')==TARGET:
                target_gt={k:v for k,v in r.items() if k.startswith('genome_')}
                break

order=['genome_410','genome_474','genome_Mpo','genome_R108','genome_436','genome_457','genome_454','genome_A17',
       'genome_395','genome_461','genome_468','genome_472','genome_M22','genome_Mar','genome_Mru','genome_ZM4',
       'genome_482','genome_M46','genome_Msa']

def final_class(call):
    if call.startswith('DEL'):
        return 'DEL_present'
    if call.startswith('INTACT'):
        return 'INTACT_absent'
    if call == 'UNCERTAIN':
        return 'UNCERTAIN'
    if call == 'EXCLUDE':
        return 'EXCLUDE'
    return 'NA'

def eval_group(sample, group, call):
    if sample == 'genome_Msa': return 'reference_self_exclude'
    if sample == 'genome_M46': return 'low_quality_exclude'
    if sample == 'genome_482': return 'uncertain_exclude'
    if sample == 'genome_ZM4': return 'spineless_exception_intact'
    if sample in {'genome_410','genome_Mpo','genome_R108','genome_436','genome_457','genome_454','genome_474','genome_A17'}:
        return 'validated_spiny_or_spiny_like_intact_set'
    if sample in {'genome_395','genome_461','genome_468','genome_472','genome_M22','genome_Mar','genome_Mru'}:
        return 'validated_spineless_DEL_set'
    return 'other'

def is_present_gt(gt):
    return gt in {'0/1','1/0','1/1','0|1','1|0','1|1'}

rows=[]
for sample in order:
    m=manual.get(sample, {})
    s=support.get(sample, {})
    call=m.get('recommended_event_call','NA')
    raw_gt=m.get('genotyped_GT') or target_gt.get(sample,'NA')
    raw_present=is_present_gt(raw_gt)
    truth=final_class(call)
    truth_present = truth == 'DEL_present'
    fp='no'
    if raw_present and truth == 'INTACT_absent':
        fp='yes_genotype_or_nearby_caller_overcall'
    if int(s.get('n_callers_support','0') or 0)>0 and truth == 'INTACT_absent':
        if fp == 'no': fp='yes_nearby_caller_support_not_target_large_DEL'
        else: fp += ';nearby_caller_support_not_target_large_DEL'
    rows.append({
        'sample':sample,
        'phenotype_group_curated':m.get('phenotype_group_curated', s.get('phenotype','NA')),
        'eval_group':eval_group(sample,m.get('phenotype_group_curated',''),call),
        'raw_near_region_pbsv_count':s.get('pbsv','NA'),
        'raw_near_region_sniffles2_count':s.get('sniffles2','NA'),
        'raw_near_region_cutesv_count':s.get('cutesv','NA'),
        'raw_near_region_n_callers_support':s.get('n_callers_support','NA'),
        'raw_near_region_support_callers':s.get('support_callers','NA'),
        'target_genotyped_GT':raw_gt,
        'target_genotyped_DR':m.get('DR','NA'),
        'target_genotyped_DV':m.get('DV','NA'),
        'read_based_discovery_SUPP_VEC':m.get('read_based_discovery_SUPP_VEC','NA'),
        'validated_truth_call':call,
        'validated_truth_class':truth,
        'truth_DEL_binary': '1' if truth_present else ('0' if truth=='INTACT_absent' else 'NA'),
        'raw_false_positive_flag':fp,
        'manual_igv_source':LOCAL_IGV,
        'note':m.get('note','NA')
    })

fields=list(rows[0].keys())
out=RUN/'results/Chr23997.target_large_intron_DEL.truth_calibrated_per_species.tsv'
with out.open('w') as fo:
    w=csv.DictWriter(fo, delimiter='\t', fieldnames=fields)
    w.writeheader(); w.writerows(rows)

def nCk(n,k):
    if k<0 or k>n: return 0
    fn=getattr(math,'comb',None)
    if fn: return fn(n,k)
    k=min(k,n-k); out=1
    for i in range(1,k+1): out=out*(n-k+i)//i
    return out

def fisher(a,b,c,d):
    n=a+b+c+d
    if n==0: return 1.0
    row1=a+b; col1=a+c
    lo=max(0,row1-(n-col1)); hi=min(row1,col1)
    def prob(x): return nCk(col1,x)*nCk(n-col1,row1-x)/float(nCk(n,row1))
    pobs=prob(a)
    return min(1.0, sum(prob(x) for x in range(lo,hi+1) if prob(x)<=pobs+1e-15))

def summarize_model(name, spiny, spineless):
    by={r['sample']:r for r in rows}
    def count(samples):
        pres=absn=miss=0
        detail=[]
        for sm in samples:
            cl=by[sm]['validated_truth_class']
            if cl=='DEL_present': pres+=1
            elif cl=='INTACT_absent': absn+=1
            else: miss+=1
            detail.append(sm+':'+by[sm]['validated_truth_call'])
        return pres,absn,miss,';'.join(detail)
    sp_p,sp_a,sp_m,sp_detail=count(spiny)
    sl_p,sl_a,sl_m,sl_detail=count(spineless)
    return {
        'model':name,
        'spiny_set':','.join(spiny),
        'spineless_set':','.join(spineless),
        'spiny_DEL_present':sp_p,'spiny_intact_absent':sp_a,'spiny_missing_or_excluded':sp_m,
        'spineless_DEL_present':sl_p,'spineless_intact_absent':sl_a,'spineless_missing_or_excluded':sl_m,
        'spiny_DEL_rate': sp_p/(sp_p+sp_a) if sp_p+sp_a else 'NA',
        'spineless_DEL_rate': sl_p/(sl_p+sl_a) if sl_p+sl_a else 'NA',
        'fisher_p_DEL_presence': fisher(sp_p,sp_a,sl_p,sl_a),
        'spiny_detail':sp_detail,
        'spineless_detail':sl_detail
    }

models=[
    summarize_model('core_clear_no_ambiguous_no_A17', ['genome_410','genome_Mpo','genome_R108'], ['genome_395','genome_461','genome_468','genome_472','genome_M22','genome_Mar','genome_Mru']),
    summarize_model('validated_spiny_like_vs_core_spineless_with_A17', ['genome_410','genome_474','genome_Mpo','genome_R108','genome_436','genome_457','genome_454','genome_A17'], ['genome_395','genome_461','genome_468','genome_472','genome_M22','genome_Mar','genome_Mru']),
    summarize_model('validated_spiny_like_vs_spineless_plus_ZM4_exception', ['genome_410','genome_474','genome_Mpo','genome_R108','genome_436','genome_457','genome_454','genome_A17'], ['genome_395','genome_461','genome_468','genome_472','genome_M22','genome_Mar','genome_Mru','genome_ZM4']),
]
summary=RUN/'summary/Chr23997.truth_calibration_summary.tsv'
with summary.open('w') as fo:
    w=csv.DictWriter(fo, delimiter='\t', fieldnames=list(models[0].keys()))
    w.writeheader(); w.writerows(models)

# False-positive audit for raw calls near/at the target event.
fp=[r for r in rows if r['raw_false_positive_flag']!='no']
with (RUN/'summary/Chr23997.raw_call_false_positive_audit.tsv').open('w') as fo:
    w=csv.DictWriter(fo, delimiter='\t', fieldnames=fields)
    w.writeheader(); w.writerows(fp)

md=RUN/'summary/Chr23997.validation_first_candidate_standard.md'
md.write_text(f'''# Chr23997 validation-first candidate standard\n\nTarget event: `{TARGET}`, Msa reference `Chr4:88,980,831-88,981,052`, large DEL in the Chr23997 intron 2-3 region.\n\nRationale: genome-wide genotype/PAV matrices alone over-call this locus in several spiny-like samples. Therefore, the candidate is evaluated by a validation-first rule: a sample is counted as DEL-present only when the large intronic deletion is supported by read-level IGV/depth evidence; nearby small indels or shifted breakpoints are not counted as the target large DEL.\n\nKey outputs:\n\n- `results/Chr23997.target_large_intron_DEL.truth_calibrated_per_species.tsv`: per-species raw caller/genotype/manual truth table.\n- `summary/Chr23997.truth_calibration_summary.tsv`: Fisher exact tests for curated phenotype models.\n- `summary/Chr23997.raw_call_false_positive_audit.tsv`: samples where raw genotype or nearby caller support conflicts with IGV truth.\n\nRecommended manuscript wording: Chr23997 was recovered by an SV-phenotype co-segregation screen only after a read-level validation step. The raw caller/genotype signal alone is insufficient because repeat/small-indel signals near the intron produce false positives in genome_436, genome_457 and genome_474.\n''', encoding='utf-8')

print(out)
print(summary)
