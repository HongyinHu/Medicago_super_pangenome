#!/usr/bin/env python3
from pathlib import Path
import csv
out = Path('path/to/project/N_4.pod_spiny/24_Chr23997_strict_RBH_gene_tree_with_outgroup_20260709')
src = Path('path/to/project/N_4.pod_spiny/14_Chr23997_ortholog_sequences_20260708_with_A17/results')
msa_cds = Path('path/to/projects_all/pan_genome/P_genome_Medicago_sativa_T2T/9.genome_annotation_gene_predict/output/5.evidencemodeler_combind/finally.genome.anno.cds')
msa_pep = Path('path/to/projects_all/pan_genome/P_genome_Medicago_sativa_T2T/9.genome_annotation_gene_predict/output/5.evidencemodeler_combind/finally.genome.anno.pep')

def read_fasta(p):
    d={}; order=[]; h=None; seq=[]
    for line in open(p):
        line=line.strip()
        if not line: continue
        if line.startswith('>'):
            if h is not None:
                d[h]=''.join(seq); order.append(h)
            h=line[1:]; seq=[]
        else:
            seq.append(line)
    if h is not None:
        d[h]=''.join(seq); order.append(h)
    return d, order

def get_by_prefix(path, prefix):
    d, order = read_fasta(path)
    for h in order:
        if h.split()[0] == prefix or h.startswith(prefix+' ') or h == prefix:
            return h, d[h]
    raise SystemExit('not found: '+prefix+' in '+str(path))

def short_id(h): return h.split('|')[0]
cds, order = read_fasta(src/'Chr23997_strict_RBH_orthologs.cds.fa')
pep, orderp = read_fasta(src/'Chr23997_strict_RBH_orthologs.pep.fa')
out_cds_h, out_cds_seq = get_by_prefix(msa_cds, 'Chr23998.1')
out_pep_h, out_pep_seq = get_by_prefix(msa_pep, 'Chr23998.1')
rows=[]
with open(out/'data/Chr23997_strict_RBH_plus_outgroup.id_map.tsv','w', newline='') as fh:
    w=csv.DictWriter(fh, fieldnames=['tree_id','source','transcript','gene','cds_len','pep_len','cds_mod3','original_header'], delimiter='\t')
    w.writeheader()
    with open(out/'data/Chr23997_strict_RBH_plus_outgroup.cds.simple.fa','w') as fc, open(out/'data/Chr23997_strict_RBH_plus_outgroup.pep.simple.fa','w') as fp:
        for h in order:
            sid=short_id(h)
            ph=next(x for x in orderp if short_id(x)==sid)
            cs=cds[h].upper().replace('U','T').replace(' ','')
            ps=pep[ph].upper().replace(' ','')
            fc.write('>'+sid+'\n'+'\n'.join(cs[i:i+70] for i in range(0,len(cs),70))+'\n')
            fp.write('>'+sid+'\n'+'\n'.join(ps[i:i+70] for i in range(0,len(ps),70))+'\n')
            transcript=''; gene=''
            for part in h.split('|'):
                if part.startswith('transcript='): transcript=part.split('=',1)[1]
                if part.startswith('gene='): gene=part.split('=',1)[1]
            w.writerow({'tree_id':sid,'source':'Chr23997_strict_RBH','transcript':transcript,'gene':gene,'cds_len':len(cs),'pep_len':len(ps),'cds_mod3':len(cs)%3,'original_header':h})
        ogid='outgroup_Msa_Chr23998'
        cs=out_cds_seq.upper().replace('U','T').replace(' ','')
        ps=out_pep_seq.upper().replace(' ','')
        fc.write('>'+ogid+'\n'+'\n'.join(cs[i:i+70] for i in range(0,len(cs),70))+'\n')
        fp.write('>'+ogid+'\n'+'\n'.join(ps[i:i+70] for i in range(0,len(ps),70))+'\n')
        w.writerow({'tree_id':ogid,'source':'Msa_paralog_outgroup','transcript':'Chr23998.1','gene':'Chr23998','cds_len':len(cs),'pep_len':len(ps),'cds_mod3':len(cs)%3,'original_header':out_cds_h+' / '+out_pep_h})
