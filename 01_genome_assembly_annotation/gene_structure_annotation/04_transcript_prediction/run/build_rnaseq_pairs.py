#!/usr/bin/env python3
from pathlib import Path
import re
base=Path('path/to/project/N_1.coding_gene_anno')
rna=base/'00_inputs/my_species_RNA'
out=base/'04_transcript_prediction/rnaseq_pairs.tsv'
rows=[]
def strip_ext(name):
    return re.sub(r'(?:\.f(?:ast)?q\.gz|\.fq\.gz|\.fastq\.gz|\.fq|\.fastq)$','',name,flags=re.I)
def readno(name):
    stem=strip_ext(name).lower()
    if re.search(r'(^|[_\.\-])(?:r1|read1|f1)$', stem): return '1'
    if re.search(r'(^|[_\.\-])(?:r2|read2)$', stem): return '2'
    if re.search(r'(^|[_\.\-])1$', stem): return '1'
    if re.search(r'(^|[_\.\-])2$', stem): return '2'
    return ''
def pair_key(name):
    stem=strip_ext(name)
    stem=re.sub(r'([_\.\-])(?:R1|R2|r1|r2|read1|read2|f1)$','',stem)
    stem=re.sub(r'([_\.\-])(?:1|2)$','',stem)
    stem=re.sub(r'_clean$','',stem,flags=re.I)
    return stem
for sample_dir in sorted(p for p in rna.glob('genome_*') if p.is_dir()):
    illum=sample_dir/'illumina-seq'
    if not illum.is_dir():
        continue
    groups={}
    for p in sorted(illum.iterdir()):
        if not p.is_file() and not p.is_symlink():
            continue
        if not re.search(r'\.(fq|fastq)(\.gz)?$', p.name, re.I):
            continue
        rn=readno(p.name)
        if rn not in {'1','2'}:
            continue
        key=pair_key(p.name)
        groups.setdefault(key,{})[rn]=p.resolve()
    for key,d in sorted(groups.items()):
        if '1' in d and '2' in d:
            rows.append((sample_dir.name,key,str(d['1']),str(d['2'])))
with out.open('w') as fh:
    fh.write('sample\tpair_id\tread1\tread2\n')
    for r in rows:
        fh.write('\t'.join(r)+'\n')
print(f'wrote {len(rows)} paired RNA-seq libraries to {out}')
