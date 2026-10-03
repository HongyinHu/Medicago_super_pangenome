#!/usr/bin/env python3
from pathlib import Path
import gzip
vcf = Path('path/to/project/38.medicago_resequence/test_population/01.population_analysis/04_tree/sativa182.pruned.vcf')
out_phy = Path('path/to/project/38.medicago_resequence/test_population/01.population_analysis/04_tree/sativa182.pruned.phy')
out_fa = Path('path/to/project/38.medicago_resequence/test_population/01.population_analysis/04_tree/sativa182.pruned.fasta')
iupac = {frozenset('AG'):'R', frozenset('CT'):'Y', frozenset('CG'):'S', frozenset('AT'):'W', frozenset('GT'):'K', frozenset('AC'):'M'}
openfn = gzip.open if vcf.suffix == '.gz' else open
samples = []
seqs = []
nvar = 0
with openfn(vcf, 'rt', errors='replace') as fh:
    for line in fh:
        if line.startswith('##'):
            continue
        if line.startswith('#CHROM'):
            samples = line.rstrip('\n').split('\t')[9:]
            seqs = [[] for _ in samples]
            continue
        if not samples:
            continue
        cols = line.rstrip('\n').split('\t')
        ref, alt = cols[3].upper(), cols[4].upper()
        if len(ref) != 1 or len(alt) != 1 or ref not in 'ACGT' or alt not in 'ACGT':
            continue
        alleles = [ref, alt]
        for i, gtfield in enumerate(cols[9:]):
            gt = gtfield.split(':', 1)[0].replace('|', '/')
            if gt in ('0/0', '0'):
                base = ref
            elif gt in ('1/1', '1'):
                base = alt
            elif gt in ('0/1', '1/0'):
                base = iupac.get(frozenset((ref, alt)), 'N')
            else:
                base = 'N'
            seqs[i].append(base)
        nvar += 1
if not samples or nvar == 0:
    raise SystemExit('No usable SNPs found for tree alignment')
with out_phy.open('w') as fo:
    fo.write(f'{len(samples)} {nvar}\n')
    for name, seq in zip(samples, seqs):
        safe = ''.join(c if c.isalnum() or c in '_.-' else '_' for c in name)[:30]
        fo.write(f'{safe:<32} {"".join(seq)}\n')
with out_fa.open('w') as fo:
    for name, seq in zip(samples, seqs):
        safe = ''.join(c if c.isalnum() or c in '_.-' else '_' for c in name)
        fo.write(f'>{safe}\n')
        s = ''.join(seq)
        for j in range(0, len(s), 80):
            fo.write(s[j:j+80] + '\n')
print(f'Wrote {len(samples)} samples x {nvar} SNPs')
