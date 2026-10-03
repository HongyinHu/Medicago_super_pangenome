#!/usr/bin/env python3
import os, re, sys, math
import pandas as pd
import numpy as np

WORK = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
BASE = 'path/to/project'
N3 = os.path.join(BASE, 'N_3.call_SV')
SPECIES_MAP = os.path.join(WORK, 'config', 'species_map.tsv')
SAMPLES = os.path.join(WORK, 'config', 'samples.tsv')
ORTHO = os.path.join(BASE, '1.orthology_family_2/output/1.orthofinder/orthofinder_for_pangenome/Results_Dec07/Orthogroups/Orthogroups.tsv')
COUNTS = os.path.join(BASE, '1.orthology_family_2/output/1.orthofinder/orthofinder_for_pangenome/Results_Dec07/Orthogroups/Orthogroups.GeneCount.tsv')
PAV_MATRIX = os.path.join(N3, '07.read_mapping_flankQC_panSV_Msa/results/Msa/panSV/panSV.Msa.flankQC.read_primary.PAV.matrix.tsv')
PAV_CONTEXT = os.path.join(N3, '07.read_mapping_flankQC_panSV_Msa/results/Msa/tables/panSV_PAV_gene_context.exclude_genome_Msa.events.tsv')
OLD08 = os.path.join(N3, '08_mixed_illumina_expression_SV_distance_20260705/tables/genome_Msa_vs_targets_expression_divergence.expressed_genes.tsv')
TABLES = os.path.join(WORK, 'tables')
os.makedirs(TABLES, exist_ok=True)

def clean_msa_gene(full_id):
    x = str(full_id).strip()
    if '|' in x:
        x = x.split('|', 1)[1]
    # Msa_T2T orthofinder IDs are mostly Chr00001.1; PAV context uses Chr00001.
    x = re.sub(r'\.\d+$', '', x)
    return x

def split_gene_cell(x):
    if pd.isna(x):
        return []
    s = str(x).strip()
    if not s:
        return []
    return [g.strip() for g in s.split(',') if g.strip()]

def yes_series(s):
    return s.astype(str).str.strip().isin(['1','1.0','TRUE','True','true','yes','YES'])

species = pd.read_csv(SPECIES_MAP, sep='\t')
samples = pd.read_csv(SAMPLES, sep='\t')
# Validate quant outputs.
quant_records = []
for _, row in samples.iterrows():
    q = os.path.join(WORK, 'quants', row['sample_id'], 'quant.sf')
    if not os.path.exists(q) or os.path.getsize(q) == 0:
        raise SystemExit('missing quant.sf: %s' % q)
    df = pd.read_csv(q, sep='\t', usecols=['Name','TPM'])
    df['sample_id'] = row['sample_id']
    df['group'] = row['group']
    quant_records.append(df)
expr = pd.concat(quant_records, ignore_index=True)
expr.to_csv(os.path.join(TABLES, 'strict_all_sample_gene_TPM.long.tsv'), sep='\t', index=False)
mean_expr = expr.groupby(['group','Name'], as_index=False)['TPM'].mean().rename(columns={'TPM':'mean_TPM'})
mean_expr.to_csv(os.path.join(TABLES, 'strict_gene_expression_group_mean_TPM.tsv'), sep='\t', index=False)
expr_lookup = dict(zip(zip(mean_expr['group'], mean_expr['Name']), mean_expr['mean_TPM']))

og = pd.read_csv(ORTHO, sep='\t', dtype=str)
gc = pd.read_csv(COUNTS, sep='\t')
for col in gc.columns:
    if col not in ['Orthogroup']:
        gc[col] = pd.to_numeric(gc[col], errors='coerce').fillna(0).astype(int)

msa_col = 'Msa_T2T'
all_orth = []
count_rows = []
for _, sp in species.iterrows():
    group = sp['group']
    target_col = sp['orthofinder_col']
    if group == 'genome_Msa':
        continue
    if target_col not in og.columns or target_col not in gc.columns:
        raise SystemExit('orthofinder column missing: %s' % target_col)
    sub_counts = gc[(gc[msa_col] == 1) & (gc[target_col] == 1)][['Orthogroup', msa_col, target_col]]
    sub = sub_counts.merge(og[['Orthogroup', msa_col, target_col]], on='Orthogroup', how='inner', suffixes=('_n',''))
    n_total = len(sub)
    good = []
    for _, r in sub.iterrows():
        msa_genes = split_gene_cell(r[msa_col])
        tgt_genes = split_gene_cell(r[target_col])
        if len(msa_genes) == 1 and len(tgt_genes) == 1:
            good.append((group, target_col, r['Orthogroup'], msa_genes[0], tgt_genes[0], clean_msa_gene(msa_genes[0])))
    tmp = pd.DataFrame(good, columns=['target_group','target_orthofinder_col','Orthogroup','msa_full_id','target_full_id','gene_id'])
    all_orth.append(tmp)
    count_rows.append({'target_group':group, 'target_orthofinder_col':target_col, 'one_to_one_orthogroups':n_total, 'usable_single_gene_pairs':len(tmp)})
orth = pd.concat(all_orth, ignore_index=True)
orth.to_csv(os.path.join(TABLES, 'strict_Msa_to_target_one_to_one_orthologs.tsv'), sep='\t', index=False)
pd.DataFrame(count_rows).to_csv(os.path.join(TABLES, 'strict_ortholog_pair_counts.tsv'), sep='\t', index=False)

ctx = pd.read_csv(PAV_CONTEXT, sep='\t', dtype={'SV_ID':str, 'nearest_gene_id':str})
ctx['distance_to_gene_bp'] = pd.to_numeric(ctx['distance_to_gene_bp'], errors='coerce')
ctx_near = ctx[(ctx['nearest_gene_id'].notna()) & (ctx['nearest_gene_id'] != '') & (ctx['distance_to_gene_bp'] <= 2000)].copy()
ctx_near = ctx_near[['SV_ID','nearest_gene_id','category4','gene_body_detail','distance_to_gene_bp']].drop_duplicates()
need_cols = ['SV_ID'] + [x for x in species['pav_col'].tolist() if x != 'genome_Msa']
pav = pd.read_csv(PAV_MATRIX, sep='\t', usecols=need_cols)
near_pav = ctx_near.merge(pav, on='SV_ID', how='inner')

sv_status = {}
cat_records = []
priority = ['Exon', 'Intron', 'Upstream', 'Downstream', 'Other']
for _, sp in species.iterrows():
    group = sp['group']
    pav_col = sp['pav_col']
    if group == 'genome_Msa':
        continue
    present = near_pav[yes_series(near_pav[pav_col])].copy()
    genes = set(present['nearest_gene_id'].astype(str))
    sv_status[group] = genes
    if len(present) == 0:
        continue
    def event_cat(row):
        if str(row.get('gene_body_detail','')) == 'exon_overlap': return 'Exon'
        if str(row.get('gene_body_detail','')) == 'intron_only': return 'Intron'
        if str(row.get('category4','')) == 'Upstream': return 'Upstream'
        if str(row.get('category4','')) == 'Downstream': return 'Downstream'
        return 'Other'
    present['PAV_gene_structure_category'] = present.apply(event_cat, axis=1)
    rank = {c:i for i,c in enumerate(priority)}
    present['rank'] = present['PAV_gene_structure_category'].map(rank).fillna(99)
    best = present.sort_values(['nearest_gene_id','rank','distance_to_gene_bp']).drop_duplicates(['nearest_gene_id'])
    for _, rr in best.iterrows():
        cat_records.append({'target_group':group, 'gene_id':str(rr['nearest_gene_id']), 'PAV_gene_structure_category':rr['PAV_gene_structure_category']})
cat_df = pd.DataFrame(cat_records)
cat_df.to_csv(os.path.join(TABLES, 'strict_gene_species_best_PAV_gene_structure_category.tsv'), sep='\t', index=False)

rows = []
cat_rows = []
for group in [g for g in species['group'].tolist() if g != 'genome_Msa']:
    m = orth[orth['target_group'] == group]
    cat_map = {}
    if len(cat_df):
        subcat = cat_df[cat_df['target_group'] == group]
        cat_map = dict(zip(subcat['gene_id'], subcat['PAV_gene_structure_category']))
    with_genes = sv_status.get(group, set())
    for _, r in m.iterrows():
        ref_tpm = float(expr_lookup.get(('genome_Msa', r['msa_full_id']), 0.0))
        tgt_tpm = float(expr_lookup.get((group, r['target_full_id']), 0.0))
        log2fc = math.log(tgt_tpm + 1.0, 2) - math.log(ref_tpm + 1.0, 2)
        div = abs(log2fc)
        expressed_any = (ref_tpm >= 1.0) or (tgt_tpm >= 1.0)
        status = 'With SV' if r['gene_id'] in with_genes else 'Without SV'
        base = {
            'gene_id': r['gene_id'],
            'comparison': 'genome_Msa vs %s' % group,
            'target_group': group,
            'Orthogroup': r['Orthogroup'],
            'msa_full_id': r['msa_full_id'],
            'target_full_id': r['target_full_id'],
            'mean_TPM_ref': ref_tpm,
            'mean_TPM_target': tgt_tpm,
            'log2FC_target_vs_ref': log2fc,
            'expression_divergence': div,
            'expressed_any': expressed_any,
            'SV_within_2kb': status == 'With SV',
            'sv_status': status,
        }
        rows.append(base)
        cat = cat_map.get(r['gene_id'], 'Without SV')
        b2 = dict(base)
        b2['PAV_gene_structure_category'] = cat
        cat_rows.append(b2)
all_df = pd.DataFrame(rows)
all_df.to_csv(os.path.join(TABLES, 'strict_ortholog_expression_divergence.all_one_to_one_genes.tsv'), sep='\t', index=False)
exp_df = all_df[all_df['expressed_any']].copy()
exp_df.to_csv(os.path.join(TABLES, 'strict_ortholog_expression_divergence.expressed_genes.tsv'), sep='\t', index=False)
cat_all = pd.DataFrame(cat_rows)
cat_exp = cat_all[cat_all['expressed_any']].copy()
cat_exp.to_csv(os.path.join(TABLES, 'strict_PAV_gene_structure_expression_divergence.expressed_genes.tsv'), sep='\t', index=False)

summary = exp_df.groupby(['comparison','target_group','sv_status']).agg(
    n_genes=('gene_id','count'),
    median_expression_divergence=('expression_divergence','median'),
    mean_expression_divergence=('expression_divergence','mean'),
    q25=('expression_divergence', lambda x: np.quantile(x, 0.25)),
    q75=('expression_divergence', lambda x: np.quantile(x, 0.75))
).reset_index()
summary.to_csv(os.path.join(TABLES, 'strict_ortholog_expression_divergence.summary.tsv'), sep='\t', index=False)
cat_summary = cat_exp.groupby(['comparison','target_group','PAV_gene_structure_category']).agg(
    n_genes=('gene_id','count'),
    median_expression_divergence=('expression_divergence','median'),
    mean_expression_divergence=('expression_divergence','mean')
).reset_index()
cat_summary.to_csv(os.path.join(TABLES, 'strict_PAV_gene_structure_expression_divergence.summary.tsv'), sep='\t', index=False)

if os.path.exists(OLD08):
    old = pd.read_csv(OLD08, sep='\t')
    old = old[old['target_group'].isin([g for g in species['group'].tolist() if g != 'genome_Msa'])].copy()
    old['method'] = '08_reference_projected'
    strict = exp_df.copy()
    strict['method'] = '09_strict_one_to_one_ortholog'
    common_cols = ['gene_id','comparison','target_group','mean_TPM_ref','mean_TPM_target','log2FC_target_vs_ref','expression_divergence','expressed_any','SV_within_2kb','sv_status','method']
    comp = pd.concat([old[common_cols], strict[common_cols]], ignore_index=True)
    comp.to_csv(os.path.join(TABLES, 'strict_vs_reference_projected.expression_divergence.long.tsv'), sep='\t', index=False)
    comp_sum = comp.groupby(['method','comparison','target_group','sv_status']).agg(
        n_genes=('gene_id','count'),
        median_expression_divergence=('expression_divergence','median'),
        mean_expression_divergence=('expression_divergence','mean')
    ).reset_index()
    comp_sum.to_csv(os.path.join(TABLES, 'strict_vs_reference_projected.expression_divergence.summary.tsv'), sep='\t', index=False)
    shared = old.merge(strict, on=['gene_id','target_group'], suffixes=('_08','_09'))
    corr_rows = []
    for group, sub in shared.groupby('target_group'):
        corr = sub[['expression_divergence_08','expression_divergence_09']].corr(method='spearman').iloc[0,1] if len(sub) > 2 else np.nan
        corr_rows.append({'target_group':group, 'shared_gene_count':len(sub), 'spearman_expression_divergence_08_vs_09':corr})
    pd.DataFrame(corr_rows).to_csv(os.path.join(TABLES, 'strict_vs_reference_projected.shared_gene_correlation.tsv'), sep='\t', index=False)
print('strict expression analysis complete')
print('expressed strict rows:', len(exp_df))
