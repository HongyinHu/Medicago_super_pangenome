#!/usr/bin/env python3
from pathlib import Path
inp=Path('path/to/project/N_4.pod_spiny/24_Chr23997_strict_RBH_gene_tree_with_outgroup_20260709/alignments/Chr23997_strict_RBH_plus_outgroup.cds.codon.aln.fa')
out=Path('path/to/project/N_4.pod_spiny/24_Chr23997_strict_RBH_gene_tree_with_outgroup_20260709/alignments/Chr23997_strict_RBH_plus_outgroup.codon12.aln.fa')
def read_fasta(p):
    arr=[]; h=None; seq=[]
    for line in open(p):
        line=line.strip()
        if not line: continue
        if line.startswith('>'):
            if h is not None: arr.append((h,''.join(seq)))
            h=line[1:]; seq=[]
        else: seq.append(line)
    if h is not None: arr.append((h,''.join(seq)))
    return arr
with open(out,'w') as fo:
    for h,s in read_fasta(inp):
        if len(s)%3 != 0: raise SystemExit('codon alignment length not divisible by 3: '+h)
        ns=''.join(s[i:i+2] for i in range(0,len(s),3))
        fo.write('>'+h+'\n'+'\n'.join(ns[i:i+70] for i in range(0,len(ns),70))+'\n')
