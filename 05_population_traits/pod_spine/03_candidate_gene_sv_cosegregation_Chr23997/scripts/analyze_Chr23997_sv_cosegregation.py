#!/usr/bin/env python3
import gzip, os, re, sys, math, json
from collections import Counter, defaultdict

root = "path/to/project"
base = root + "/N_4.pod_spiny"
run = base + "/03_candidate_gene_sv_cosegregation_Chr23997"
vcf = root + "/N_3.call_SV/01.read_based_dualref_hifi/05_dualref_panSV_20260626/results/RefA/panSV/panSV.RefA.genotyped.18species.vcf.gz"
gff = base + "/00_data/4.reference_anno/genome_Msa.gff"
ann = root + "/16.T2T_ref_function_anno/output"
chrom = "Chr4"
gene_id = "Chr23997"
transcript_id = "Chr23997.1"
gene_start = 88979818
gene_end = 88982331
strand = "-"
windows = [0, 10_000, 50_000, 100_000]

phenotype = {
    "genome_410": "spiny_main",
    "genome_474": "spiny_main",
    "genome_Mpo": "spiny_main",
    "genome_R108": "spiny_main",
    "genome_395": "spineless_main",
    "genome_454": "spineless_main_hairy",
    "genome_461": "spineless_main",
    "genome_468": "spineless_main",
    "genome_472": "spineless_main",
    "genome_482": "spineless_main",
    "genome_M22": "spineless_main",
    "genome_Mar": "spineless_main",
    "genome_Mru": "spineless_main",
    "genome_Msa": "spineless_main_ref",
    "genome_ZM4": "spineless_main",
    "genome_436": "weak_spiny_ambiguous",
    "genome_457": "spiny_low_assembly_quality",
    "genome_M46": "spineless_low_assembly_quality",
}
main_spiny = [s for s,p in phenotype.items() if p == "spiny_main"]
main_spineless = [s for s,p in phenotype.items() if p.startswith("spineless_main")]
all_spiny_sensitive = main_spiny + ["genome_436", "genome_457"]
all_spineless_sensitive = main_spineless + ["genome_M46"]

os.makedirs(run + "/results", exist_ok=True)
os.makedirs(run + "/summary", exist_ok=True)
os.makedirs(run + "/logs", exist_ok=True)

def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def parse_info(info):
    d = {}
    for item in info.split(';'):
        if not item: continue
        if '=' in item:
            k,v = item.split('=',1)
            d[k]=v
        else:
            d[item]=True
    return d

def get_svtype(row, info):
    alt = row[4]
    if 'SVTYPE' in info: return info['SVTYPE']
    m = re.search(r'<([^>]+)>', alt)
    return m.group(1) if m else alt

def get_end(row, info, svtype):
    pos = int(row[1])
    if 'END' in info:
        try: return int(float(info['END']))
        except Exception: pass
    if 'SVLEN' in info:
        vals = info['SVLEN'].split(',')
        try: return pos + abs(int(float(vals[0])))
        except Exception: pass
    return pos

def get_svlen(row, info, start, end):
    if 'SVLEN' in info:
        try:
            return int(float(info['SVLEN'].split(',')[0]))
        except Exception:
            pass
    return end-start+1

def gt_presence(sample_field, fmt_keys):
    vals = sample_field.split(':')
    fm = dict(zip(fmt_keys, vals))
    gt = fm.get('GT', './.')
    if gt in ('./.', '.', '.|.'):
        return None, gt
    alleles = re.split(r'[|/]', gt)
    if any(a not in ('.','0') for a in alleles):
        return 1, gt
    if all(a == '0' for a in alleles if a != '.'):
        return 0, gt
    return None, gt

def overlap(a,b,c,d):
    return max(0, min(b,d)-max(a,c)+1)

def context_of(start, end, features):
    exon_ol = 0; cds_ol = 0; gene_ol = overlap(start,end,gene_start,gene_end)
    for typ,s,e in features:
        ol = overlap(start,end,s,e)
        if typ == 'exon': exon_ol += ol
        if typ == 'CDS': cds_ol += ol
    if cds_ol > 0 and gene_ol > cds_ol:
        return 'CDS+intron_boundary'
    if cds_ol > 0:
        return 'CDS'
    if exon_ol > 0:
        return 'UTR_or_exon_nonCDS'
    if gene_ol > 0:
        return 'intron'
    if end < gene_start:
        return 'downstream' if strand == '-' else 'upstream'
    if start > gene_end:
        return 'upstream' if strand == '-' else 'downstream'
    return 'intergenic'

features=[]
with open(gff) as fh:
    for line in fh:
        if line.startswith('#'): continue
        parts=line.rstrip('\n').split('\t')
        if len(parts)<9: continue
        seq,src,typ,s,e,score,st,phase,attrs=parts
        if seq != chrom: continue
        if f"Parent={transcript_id}" in attrs or f"ID={transcript_id}" in attrs or f"ID={gene_id}" in attrs:
            if typ in ('exon','CDS'):
                features.append((typ,int(s),int(e)))

samples=[]
records=[]
region_min = gene_start - max(windows)
region_max = gene_end + max(windows)
with open_text(vcf) as fh:
    for line in fh:
        if line.startswith('##'):
            continue
        if line.startswith('#CHROM'):
            samples = line.rstrip('\n').split('\t')[9:]
            continue
        row=line.rstrip('\n').split('\t')
        if row[0] != chrom: continue
        pos=int(row[1])
        if pos > region_max + 1000000:  # VCF sorted; enough slack for long events starting nearby
            break
        info=parse_info(row[7])
        svtype=get_svtype(row, info)
        end=get_end(row, info, svtype)
        svlen=get_svlen(row, info, pos, end)
        if end < region_min or pos > region_max:
            continue
        fmt=row[8].split(':')
        calls={}
        gts={}
        for sm, sf in zip(samples, row[9:]):
            pres, gt = gt_presence(sf, fmt)
            calls[sm]=pres
            gts[sm]=gt
        rec_id=row[2] if row[2] != '.' else f"{row[0]}:{pos}-{end}:{svtype}"
        # window category by shortest window containing the SV interval start/end overlap
        min_window=None
        for w in windows:
            if end >= gene_start-w and pos <= gene_end+w:
                min_window=w
                break
        if min_window is None:
            continue
        sp_present=sum(1 for s in main_spiny if calls.get(s)==1)
        sp_missing=sum(1 for s in main_spiny if calls.get(s) is None)
        sp_absent=sum(1 for s in main_spiny if calls.get(s)==0)
        sl_present=sum(1 for s in main_spineless if calls.get(s)==1)
        sl_missing=sum(1 for s in main_spineless if calls.get(s) is None)
        sl_absent=sum(1 for s in main_spineless if calls.get(s)==0)
        sens_sp_present=sum(1 for s in all_spiny_sensitive if calls.get(s)==1)
        sens_sl_present=sum(1 for s in all_spineless_sensitive if calls.get(s)==1)
        # orientation: spiny-enriched or spineless-enriched
        sp_rate=sp_present/max(1, (len(main_spiny)-sp_missing))
        sl_rate=sl_present/max(1, (len(main_spineless)-sl_missing))
        if sp_present==len(main_spiny) and sl_present==0:
            status='perfect_main_spiny_presence'
        elif sp_present>=3 and sl_present<=1:
            status='near_spiny_enriched'
        elif sl_present==len(main_spineless) and sp_present==0:
            status='perfect_main_spineless_presence'
        elif sl_present>=8 and sp_present<=1:
            status='spineless_enriched'
        else:
            status='not_cosegregating'
        records.append({
            'id':rec_id,'chrom':row[0],'pos':pos,'end':end,'svtype':svtype,'svlen':svlen,
            'min_window_bp':min_window,'context':context_of(pos,end,features),
            'overlap_gene_bp':overlap(pos,end,gene_start,gene_end),
            'spiny_present':sp_present,'spiny_absent':sp_absent,'spiny_missing':sp_missing,
            'spineless_present':sl_present,'spineless_absent':sl_absent,'spineless_missing':sl_missing,
            'spiny_rate':sp_rate,'spineless_rate':sl_rate,'delta_spiny_minus_spineless':sp_rate-sl_rate,
            'sensitive_spiny_present':sens_sp_present,'sensitive_spineless_present':sens_sl_present,
            'status':status,
            'calls':calls,'gts':gts,
            'raw_info':row[7]
        })

records.sort(key=lambda r: (r['min_window_bp'], -abs(r['delta_spiny_minus_spineless']), r['pos']))

sample_order = samples
with open(run + '/results/Chr23997.region_sv_genotype_matrix.tsv','w') as out:
    out.write('\t'.join(['sv_id','chrom','pos','end','svtype','svlen','min_window_bp','context','overlap_gene_bp','status','spiny_present','spiny_absent','spiny_missing','spineless_present','spineless_absent','spineless_missing','spiny_rate','spineless_rate','delta_spiny_minus_spineless','sensitive_spiny_present','sensitive_spineless_present'] + sample_order) + '\n')
    for r in records:
        vals=[r['id'],r['chrom'],str(r['pos']),str(r['end']),r['svtype'],str(r['svlen']),str(r['min_window_bp']),r['context'],str(r['overlap_gene_bp']),r['status'],str(r['spiny_present']),str(r['spiny_absent']),str(r['spiny_missing']),str(r['spineless_present']),str(r['spineless_absent']),str(r['spineless_missing']),f"{r['spiny_rate']:.4f}",f"{r['spineless_rate']:.4f}",f"{r['delta_spiny_minus_spineless']:.4f}",str(r['sensitive_spiny_present']),str(r['sensitive_spineless_present'])]
        vals += [('1' if r['calls'].get(s)==1 else '0' if r['calls'].get(s)==0 else 'NA') for s in sample_order]
        out.write('\t'.join(vals)+'\n')

with open(run + '/results/Chr23997.region_sv_gt_strings.tsv','w') as out:
    out.write('\t'.join(['sv_id','chrom','pos','end','svtype','svlen'] + sample_order) + '\n')
    for r in records:
        out.write('\t'.join([r['id'],r['chrom'],str(r['pos']),str(r['end']),r['svtype'],str(r['svlen'])] + [r['gts'].get(s,'NA') for s in sample_order])+'\n')

# summarize by window and status
with open(run + '/summary/Chr23997.summary.tsv','w') as out:
    out.write('item\tvalue\n')
    out.write(f'gene_id\t{gene_id}\n')
    out.write(f'transcript_id\t{transcript_id}\n')
    out.write(f'coordinate\t{chrom}:{gene_start}-{gene_end}({strand})\n')
    out.write('candidate_sequence_files\t00_data/6.candicated_gene/Msa_T2T_Chr23997.1_CDS.fa;00_data/6.candicated_gene/Msa_T2T_Chr23997.1_PEP.fa\n')
    out.write('main_spiny\t' + ','.join(main_spiny) + '\n')
    out.write('main_spineless\t' + ','.join(main_spineless) + '\n')
    out.write('ambiguous_or_low_quality\tgenome_436,genome_457,genome_M46\n')
    out.write(f'total_sv_in_100kb_window\t{len(records)}\n')
    for w in windows:
        sub=[r for r in records if r['min_window_bp']<=w]
        out.write(f'sv_count_gene_plus_{w}bp\t{len(sub)}\n')
    for status,count in Counter(r['status'] for r in records).most_common():
        out.write(f'status_{status}\t{count}\n')
    perfect=[r for r in records if r['status']=='perfect_main_spiny_presence']
    near=[r for r in records if r['status']=='near_spiny_enriched']
    out.write('perfect_main_spiny_presence_ids\t' + (';'.join(r['id'] for r in perfect) if perfect else 'NA') + '\n')
    out.write('near_spiny_enriched_ids\t' + (';'.join(r['id'] for r in near) if near else 'NA') + '\n')

with open(run + '/results/Chr23997.cosegregating_candidates.tsv','w') as out:
    out.write('\t'.join(['sv_id','chrom','pos','end','svtype','svlen','min_window_bp','context','overlap_gene_bp','status','spiny_present','spineless_present','spiny_rate','spineless_rate','delta','sensitive_spiny_present','sensitive_spineless_present'])+'\n')
    for r in records:
        if r['status'] != 'not_cosegregating':
            out.write('\t'.join([r['id'],r['chrom'],str(r['pos']),str(r['end']),r['svtype'],str(r['svlen']),str(r['min_window_bp']),r['context'],str(r['overlap_gene_bp']),r['status'],str(r['spiny_present']),str(r['spineless_present']),f"{r['spiny_rate']:.4f}",f"{r['spineless_rate']:.4f}",f"{r['delta_spiny_minus_spineless']:.4f}",str(r['sensitive_spiny_present']),str(r['sensitive_spineless_present'])])+'\n')

# collect annotation snippets
ann_lines=[]
def grep_file(path, pattern, maxn=20):
    res=[]
    try:
        with open(path, errors='ignore') as fh:
            for line in fh:
                if pattern in line:
                    res.append(line.rstrip('\n'))
                    if len(res)>=maxn: break
    except FileNotFoundError:
        pass
    return res

# write annotation summary manually from known files if present
with open(run + '/summary/Chr23997.functional_annotation.txt','w') as out:
    out.write('Gene: Chr23997 / transcript Chr23997.1 / protein Chr239971\n')
    out.write(f'Coordinate: {chrom}:{gene_start}-{gene_end}({strand})\n')
    out.write('Candidate protein length: 435 aa; CDS length: 1308 bp (from candidate FASTA headers)\n\n')
    for label,path in [
        ('Pfam', ann + '/pfam/ref_Msa.T2T_ctg.pep.pfam.out'),
        ('InterPro domains', ann + '/interproscan/ref_Msa.T2T_ctg.pep.domains.txt'),
        ('InterPro full', ann + '/interproscan/ref_Msa.T2T_ctg.pep.interproscan.out'),
        ('TAIR blast annotation', ann + '/blast_tair/T2T_tait/ref_Msa.T2T_ctg.pep.tab.anno'),
        ('SwissProt annotation', ann + '/uniprot_swissprot/ref_Msa.T2T_ctg.pep_blastp.tab.anno.txt'),
        ('TrEMBL plants annotation', ann + '/uniprot_trembl_plants/ref_Msa.T2T_ctg.pep_blastp.tab.anno.txt'),
    ]:
        out.write(f'## {label}: {path}\n')
        for line in grep_file(path, 'Chr239971', 15) + grep_file(path, 'Chr23997', 5):
            out.write(line + '\n')
        out.write('\n')

# README
with open(run + '/README.md','w') as out:
    out.write('# Chr23997 SV-phenotype cosegregation analysis\n\n')
    out.write('This directory was generated under N_4.pod_spiny for candidate gene Chr23997.\n\n')
    out.write('Inputs:\n')
    out.write(f'- Candidate gene FASTA: {base}/00_data/6.candicated_gene/Msa_T2T_Chr23997.1_CDS.fa and PEP.fa\n')
    out.write(f'- Msa GFF: {gff}\n')
    out.write(f'- RefA 18-sample genotyped pan-SV VCF: {vcf}\n\n')
    out.write('Phenotype grouping:\n')
    out.write('- Main spiny: genome_410, genome_474, genome_Mpo, genome_R108\n')
    out.write('- Main spineless: genome_395, genome_454, genome_461, genome_468, genome_472, genome_482, genome_M22, genome_Mar, genome_Mru, genome_Msa, genome_ZM4\n')
    out.write('- Ambiguous/low-quality sensitivity: genome_436, genome_457, genome_M46\n\n')
    out.write('Outputs:\n')
    out.write('- results/Chr23997.region_sv_genotype_matrix.tsv: all SVs in gene +/-100 kb with presence/absence matrix\n')
    out.write('- results/Chr23997.cosegregating_candidates.tsv: SVs enriched or perfectly co-segregating with phenotype under the main grouping\n')
    out.write('- results/Chr23997.region_sv_gt_strings.tsv: raw GT strings for each SV/sample\n')
    out.write('- summary/Chr23997.summary.tsv: compact counts and candidate IDs\n')
    out.write('- summary/Chr23997.functional_annotation.txt: annotation evidence snippets\n')
