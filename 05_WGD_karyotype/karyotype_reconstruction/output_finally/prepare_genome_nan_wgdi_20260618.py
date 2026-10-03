#!/usr/bin/env python3
import re
import sys
from collections import defaultdict
from pathlib import Path

genome_fa = Path('data/genome_nan/genome_Mpo.T2T.chr.final.fa')
gff_in = Path('data/genome_nan/Mpo_genome.anno.gff')
cds_in = Path('data/genome_nan/Mpo_genome.anno.cds')
pep_in = Path('data/genome_nan/Mpo_genome.anno.pep')
prefix = 'genome_nan'
outdir = Path('data/genome_nan')

def chr_num(name):
    token = name.strip().split()[0]
    m = re.search(r'(?:^|[^0-9])Chr0*([0-9]+)\b|^0*([0-9]+)$', token, re.I)
    if not m:
        return None
    n = int(m.group(1) or m.group(2))
    if 1 <= n <= 8:
        return str(n)
    return None

def parse_attrs(s):
    d = {}
    for part in s.strip().split(';'):
        if not part:
            continue
        if '=' in part:
            k, v = part.split('=', 1)
            d[k] = v
    return d

def fasta_lengths(path):
    lengths = {}
    current = None
    total = 0
    with path.open() as f:
        for line in f:
            if line.startswith('>'):
                if current is not None:
                    c = chr_num(current)
                    if c is not None:
                        lengths[c] = total
                current = line[1:].strip().split()[0]
                total = 0
            else:
                total += len(line.strip())
        if current is not None:
            c = chr_num(current)
            if c is not None:
                lengths[c] = total
    return lengths

def rewrite_fasta(infile, outfile, idmap, is_pep=False):
    written = 0
    keep = False
    with infile.open() as inf, outfile.open('w') as out:
        for line in inf:
            if line.startswith('>'):
                old = line[1:].strip().split()[0]
                new = idmap.get(old)
                if new:
                    keep = True
                    written += 1
                    out.write(f'>{new}\n')
                else:
                    keep = False
            elif keep:
                seq = line.strip()
                if is_pep:
                    seq = seq.replace('.', '')
                if seq:
                    out.write(seq + '\n')
    return written

chr_lengths = fasta_lengths(genome_fa)
if not chr_lengths:
    raise SystemExit('No chromosome lengths parsed from FASTA')

idmap = {}
gene_count = defaultdict(int)
gff_rows = []
with gff_in.open() as inf:
    for line in inf:
        if not line.strip() or line.startswith('#'):
            continue
        cols = line.rstrip('\n').split('\t')
        if len(cols) < 9 or cols[2] != 'mRNA':
            continue
        c = chr_num(cols[0])
        if c is None:
            continue
        attrs = parse_attrs(cols[8])
        old_id = attrs.get('ID')
        if not old_id:
            continue
        gene_count[c] += 1
        order = gene_count[c]
        new_id = f'{prefix}_{int(c)}g{order:05d}'
        idmap[old_id] = new_id
        gene_name = attrs.get('Parent') or attrs.get('Name') or old_id
        gff_rows.append([c, new_id, cols[3], cols[4], cols[6], str(order), gene_name])

with (outdir/f'{prefix}.gff').open('w') as out:
    for row in gff_rows:
        out.write('\t'.join(row) + '\n')
with (outdir/f'{prefix}.lens').open('w') as out:
    for c in sorted(chr_lengths, key=int):
        out.write(f'{c}\t{chr_lengths[c]}\t{gene_count.get(c,0)}\n')
cds_n = rewrite_fasta(cds_in, outdir/f'{prefix}.cds', idmap, is_pep=False)
pep_n = rewrite_fasta(pep_in, outdir/f'{prefix}.pep', idmap, is_pep=True)

with (Path('output_finally')/'prepare_genome_nan_wgdi_20260618.summary.txt').open('w') as out:
    out.write(f'chromosomes={len(chr_lengths)}\n')
    out.write(f'mRNA_in_gff={len(idmap)}\n')
    out.write(f'cds_written={cds_n}\n')
    out.write(f'pep_written={pep_n}\n')
    for c in sorted(chr_lengths, key=int):
        out.write(f'chr{c}\tlength={chr_lengths[c]}\tgenes={gene_count.get(c,0)}\n')
print('DONE', len(chr_lengths), len(idmap), cds_n, pep_n)
