#!/usr/bin/env python3
import csv
from pathlib import Path
from urllib.parse import unquote
from collections import defaultdict, Counter

BASE = Path('path/to/project/39.TE_type_soloLTR')
OLD = BASE / 'output' / 'species_totals.tsv'
OUT = BASE / 'output' / 'publication_te_broad_classes'
TABLE = OUT / 'tables' / 'edta_broad_class_percent.nonoverlap_priority.tsv'
WIDE = OUT / 'tables' / 'edta_broad_class_percent.nonoverlap_priority.wide.tsv'
CHECK = OUT / 'tables' / 'nonoverlap_sum_vs_all_te_union.tsv'

LEFT_TO_RIGHT = ['Other / unclassified','Helitron','DNA transposons','LTR/unknown','LTR/Gypsy','LTR/Copia']
# In overlaps, prefer biologically specific LTR classifications, then Helitron/DNA, then residual classes.
PRIORITY = ['LTR/Gypsy','LTR/Copia','LTR/unknown','Helitron','DNA transposons','Other / unclassified']
PRANK = {c:i for i,c in enumerate(PRIORITY)}
PREFERRED_ORDER = [
    'genome_395','genome_436','genome_457','genome_468','genome_474','genome_482',
    'genome_M22','genome_Mar','genome_Mru','genome_Msa','genome_ZM4',
    'genome_410','genome_454','genome_461','genome_472','genome_474_T2T',
    'genome_A17','genome_M46','genome_MPO','genome_Msa_T2T','genome_R108'
]
LABELS = {
    'genome_Mar': 'M. archiducis-nicolai',
    'genome_Mru': 'M. ruthenica',
    'genome_Msa': 'M. sativa',
    'genome_ZM4': 'M. sativa cv. Zhongmu-4',
    'genome_A17': 'M. truncatula A17',
    'genome_MPO': 'M. polymorpha',
    'genome_Msa_T2T': 'M. sativa T2T',
    'genome_R108': 'M. truncatula R108',
}

def parse_attrs(text):
    attrs = {}
    for part in text.strip().split(';'):
        if not part: continue
        if '=' in part:
            k,v = part.split('=',1); attrs[k]=unquote(v)
    return attrs

def broad_category(feature, classification):
    c=(classification or '').strip(); f=(feature or '').strip(); low=f'{c} {f}'.lower()
    if c == 'LTR/Copia' or 'copia_ltr_retrotransposon' in low: return 'LTR/Copia'
    if c == 'LTR/Gypsy' or 'gypsy_ltr_retrotransposon' in low: return 'LTR/Gypsy'
    if c == 'LTR/unknown' or (c.startswith('LTR/') and c.lower().endswith('/unknown')) or f == 'LTR_retrotransposon': return 'LTR/unknown'
    if 'helitron' in low: return 'Helitron'
    if c.startswith('DNA/') or c.startswith('MITE/') or 'tir_transposon' in low or 'terminal_inverted_repeat' in low: return 'DNA transposons'
    return 'Other / unclassified'

def choose(active):
    if not active: return None
    return min(active.keys(), key=lambda c: PRANK.get(c, 999))

def sweep_bp(events_by_chrom):
    bp = Counter()
    for chrom, events in events_by_chrom.items():
        events.sort(key=lambda x: x[0])
        active = Counter()
        prev = None
        i = 0
        n = len(events)
        while i < n:
            pos = events[i][0]
            if prev is not None and pos > prev and active:
                bp[choose(active)] += pos - prev
            while i < n and events[i][0] == pos:
                _, cat, delta = events[i]
                active[cat] += delta
                if active[cat] <= 0:
                    active.pop(cat, None)
                i += 1
            prev = pos
    return bp

def pct(n,d): return 0.0 if not d else n*100.0/d

species_rows=[]
with OLD.open(newline='', encoding='utf-8') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        species_rows.append(r)
order_index={s:i for i,s in enumerate(PREFERRED_ORDER)}
species_rows.sort(key=lambda r: (order_index.get(r['species'],999), r['species']))
rows=[]; checks=[]
for sr in species_rows:
    sp=sr['species']; label=LABELS.get(sp, sp); genome_bp=int(float(sr['genome_bp']))
    events=defaultdict(list)
    gff=Path(sr['teanno_gff3'])
    with gff.open(encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if not line or line.startswith('#'): continue
            p=line.rstrip('\n').split('\t')
            if len(p)<9: continue
            chrom, feature, start, end, attrs = p[0], p[2], p[3], p[4], p[8]
            try:
                s=int(start); e=int(end)
            except ValueError:
                continue
            if e < s: s,e=e,s
            cat=broad_category(feature, parse_attrs(attrs).get('Classification',''))
            events[chrom].append((s, cat, 1))
            events[chrom].append((e+1, cat, -1))
    bp=sweep_bp(events)
    total=sum(bp.values())
    for cat in LEFT_TO_RIGHT:
        val=bp.get(cat,0)
        rows.append({'species':sp,'label':label,'category':cat,'category_order':LEFT_TO_RIGHT.index(cat)+1,'species_order':order_index.get(sp,999),'bp_nonoverlap':val,'genome_bp':genome_bp,'percent_of_genome':f'{pct(val,genome_bp):.6f}'})
    old_union=float(sr['total_te_merged_bp'])
    checks.append({'species':sp,'label':label,'nonoverlap_bp':total,'nonoverlap_percent':f'{pct(total,genome_bp):.6f}','old_all_te_union_bp':int(old_union),'old_all_te_union_percent':sr['te_percent_of_genome'],'delta_percent_points':f'{pct(total-old_union,genome_bp):.6f}'})

with TABLE.open('w', newline='', encoding='utf-8') as fh:
    fields=['species','label','category','category_order','species_order','bp_nonoverlap','genome_bp','percent_of_genome']
    w=csv.DictWriter(fh, delimiter='\t', fieldnames=fields); w.writeheader(); w.writerows(rows)
lookup={(r['species'],r['category']):r['percent_of_genome'] for r in rows}
with WIDE.open('w', newline='', encoding='utf-8') as fh:
    w=csv.writer(fh, delimiter='\t'); w.writerow(['species','label']+LEFT_TO_RIGHT)
    for sr in species_rows:
        sp=sr['species']; w.writerow([sp,LABELS.get(sp,sp)]+[lookup.get((sp,c),'0.000000') for c in LEFT_TO_RIGHT])
with CHECK.open('w', newline='', encoding='utf-8') as fh:
    fields=['species','label','nonoverlap_bp','nonoverlap_percent','old_all_te_union_bp','old_all_te_union_percent','delta_percent_points']
    w=csv.DictWriter(fh, delimiter='\t', fieldnames=fields); w.writeheader(); w.writerows(checks)
print(TABLE)
print(CHECK)
