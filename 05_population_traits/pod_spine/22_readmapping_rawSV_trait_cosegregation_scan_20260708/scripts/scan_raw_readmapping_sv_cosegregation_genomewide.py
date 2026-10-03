#!/usr/bin/env python3
import os, gzip, math, re, csv
from collections import defaultdict, Counter
from statistics import median

RUN = "path/to/project/N_3.call_SV/01.read_based_dualref_hifi"
OUT = "path/to/project/N_4.pod_spiny/22_readmapping_rawSV_trait_cosegregation_scan_20260708"
FAI = os.path.join(RUN, "02_ref_prepare/Msa/ref.fa.fai")
CHROMS = [f"Chr{i}" for i in range(1, 9)]
MIN_SVLEN = 50
BP_DIST = 500
INS_DIST = 500
MIN_RO = 0.50
LEN_RATIO_MIN = 0.50

INTACT_GROUP = ["genome_410", "genome_436", "genome_457", "genome_454", "genome_474", "genome_Mpo", "genome_R108", "genome_A17"]
DEL_GROUP = ["genome_395", "genome_461", "genome_468", "genome_472", "genome_M22", "genome_Mar", "genome_Mru"]
EXCLUDED = {"genome_Msa", "genome_M46", "genome_482", "genome_ZM4"}
SAMPLES = INTACT_GROUP + DEL_GROUP
CALLERS = {
    "pbsv": os.path.join(RUN, "03_per_sample/Msa/pbsv/{sample}.pbsv.vcf.gz"),
    "sniffles2": os.path.join(RUN, "03_per_sample/Msa/sniffles2/{sample}.sniffles2.vcf.gz"),
    "cuteSV": os.path.join(RUN, "03_per_sample/Msa/cutesv/{sample}.cutesv.vcf.gz"),
}
TARGET_CHROM = "Chr4"
TARGET_POS = 88_980_831
TARGET_END = 88_981_052
TARGET_LEN = 221
TARGET_ID = "Chr23997_second_intron_DEL_manual_IGV_corrected"

def read_chrom_lengths():
    d = {}
    with open(FAI) as f:
        for line in f:
            if not line.strip():
                continue
            a = line.rstrip('\n').split('\t')
            if a[0] in CHROMS:
                d[a[0]] = int(a[1])
    missing = [c for c in CHROMS if c not in d]
    if missing:
        raise SystemExit(f"Missing chromosomes in FAI: {missing}")
    return d

def open_vcf(path):
    return gzip.open(path, 'rt', errors='replace') if path.endswith('.gz') else open(path, 'rt', errors='replace')

def parse_info(info):
    d = {}
    for item in info.split(';'):
        if not item:
            continue
        if '=' in item:
            k, v = item.split('=', 1)
            d[k] = v
        else:
            d[item] = True
    return d

def first_int(v):
    if v is None:
        return None
    m = re.search(r'-?\d+', str(v).split(',')[0])
    return int(m.group(0)) if m else None

def parse_svtype(info, alt):
    svt = info.get('SVTYPE')
    if not svt:
        m = re.search(r'<([^>]+)>', alt)
        svt = m.group(1) if m else None
    if not svt and ('[' in alt or ']' in alt):
        svt = 'BND'
    if not svt:
        return 'UNK'
    svt = svt.upper().replace('SNV', 'UNK')
    if svt in ('TRA', 'TRANSLOCATION'):
        return 'TRA'
    if svt == 'BND':
        return 'BND'
    return svt

def parse_svlen(info, pos, end, svtype):
    vals = []
    if 'SVLEN' in info:
        for x in str(info['SVLEN']).split(','):
            try:
                vals.append(abs(int(float(x))))
            except Exception:
                pass
    if vals:
        return max(vals)
    if end and end >= pos and svtype not in ('INS', 'BND', 'TRA'):
        return abs(end - pos + 1)
    return 0

def record_from_vcf_line(line, sample, caller):
    fields = line.rstrip('\n').split('\t')
    if len(fields) < 8:
        return None
    chrom, pos_s, rid, ref, alt, qual, flt, info_s = fields[:8]
    if chrom not in CHROMS:
        return None
    try:
        pos = int(pos_s)
    except Exception:
        return None
    if flt not in ('.', 'PASS'):
        return None
    info = parse_info(info_s)
    svtype = parse_svtype(info, alt)
    if svtype in ('UNK', 'SNV'):
        return None
    end = first_int(info.get('END'))
    svlen = parse_svlen(info, pos, end, svtype)
    if svtype == 'INS':
        if end is None:
            end = pos
    elif svtype in ('BND', 'TRA'):
        end = pos
    elif end is None:
        if svlen:
            end = pos + svlen - 1
        else:
            return None
    rec_start, rec_end = min(pos, end), max(pos, end)
    if svtype not in ('BND', 'TRA') and svlen < MIN_SVLEN:
        return None
    if svtype in ('BND', 'TRA'):
        svlen = max(svlen, 50)
    return {
        'sample': sample, 'caller': caller, 'chrom': chrom, 'pos': pos, 'end': end,
        'start': rec_start, 'stop': rec_end, 'id': rid if rid else '.', 'alt': alt,
        'svtype': svtype, 'svlen': svlen, 'abs_svlen': abs(svlen), 'filter': flt,
    }

def reciprocal_overlap(a, b):
    if a['svtype'] in ('INS', 'BND', 'TRA') or b['svtype'] in ('INS', 'BND', 'TRA'):
        return 0.0
    s = max(a['start'], b['start'])
    e = min(a['stop'], b['stop'])
    ov = max(0, e - s + 1)
    la = max(1, a['stop'] - a['start'] + 1)
    lb = max(1, b['stop'] - b['start'] + 1)
    return min(ov / la, ov / lb)

def same_event(a, b):
    if a['chrom'] != b['chrom'] or a['svtype'] != b['svtype']:
        return False
    la, lb = max(1, a['abs_svlen']), max(1, b['abs_svlen'])
    ratio = min(la, lb) / max(la, lb)
    if ratio < LEN_RATIO_MIN:
        return False
    if a['svtype'] in ('INS', 'BND', 'TRA'):
        return abs(a['pos'] - b['pos']) <= INS_DIST
    bp_ok = abs(a['start'] - b['start']) <= BP_DIST and abs(a['stop'] - b['stop']) <= BP_DIST
    ro_ok = reciprocal_overlap(a, b) >= MIN_RO
    return bp_ok or ro_ok

def fisher_two_sided(a, b, c, d):
    n1 = a + b; m1 = a + c; n = a + b + c + d
    def logchoose(nv, kv):
        if kv < 0 or kv > nv:
            return float('-inf')
        return math.lgamma(nv + 1) - math.lgamma(kv + 1) - math.lgamma(nv - kv + 1)
    def prob(x):
        if x < 0 or x > n1 or x > m1 or (n1 - x) > (n - m1):
            return 0.0
        return math.exp(logchoose(m1, x) + logchoose(n - m1, n1 - x) - logchoose(n, n1))
    obs = prob(a)
    lo = max(0, n1 - (n - m1)); hi = min(n1, m1)
    return min(1.0, sum(prob(x) for x in range(lo, hi + 1) if prob(x) <= obs + 1e-15))

def cluster_records(records):
    groups = defaultdict(list)
    for idx, r in enumerate(records):
        groups[(r['chrom'], r['svtype'])].append(idx)
    parent = list(range(len(records)))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
    for key, idxs in groups.items():
        idxs.sort(key=lambda i: records[i]['start'])
        active = []
        for i in idxs:
            r = records[i]
            active = [j for j in active if records[j]['stop'] >= r['start'] - BP_DIST * 2]
            for j in active:
                if same_event(r, records[j]):
                    union(i, j)
            active.append(i)
    clusters = defaultdict(list)
    for i in range(len(records)):
        clusters[find(i)].append(i)
    return list(clusters.values())

def target_like(r):
    if r['chrom'] != TARGET_CHROM or r['svtype'] != 'DEL':
        return False
    if r['abs_svlen'] < 50:
        return False
    ratio = min(r['abs_svlen'], TARGET_LEN) / max(r['abs_svlen'], TARGET_LEN)
    if ratio < 0.50:
        return False
    fake = {'chrom': TARGET_CHROM, 'svtype': 'DEL', 'start': TARGET_POS, 'stop': TARGET_END, 'pos': TARGET_POS, 'end': TARGET_END, 'abs_svlen': TARGET_LEN}
    return abs(r['start'] - TARGET_POS) <= BP_DIST and abs(r['stop'] - TARGET_END) <= BP_DIST or reciprocal_overlap(r, fake) >= MIN_RO

def main():
    os.makedirs(os.path.join(OUT, 'results'), exist_ok=True)
    os.makedirs(os.path.join(OUT, 'summary'), exist_ok=True)
    chrom_lengths = read_chrom_lengths()
    records = []
    status_rows = []
    for sample in SAMPLES:
        for caller, template in CALLERS.items():
            path = template.format(sample=sample)
            n = 0
            if not os.path.exists(path):
                status_rows.append([sample, caller, path, 'missing', 0])
                continue
            try:
                with open_vcf(path) as fh:
                    for line in fh:
                        if not line or line.startswith('#'):
                            continue
                        rec = record_from_vcf_line(line, sample, caller)
                        if rec:
                            records.append(rec)
                            n += 1
                status_rows.append([sample, caller, path, 'ok', n])
            except Exception as e:
                status_rows.append([sample, caller, path, f'error:{e}', n])
    status_path = os.path.join(OUT, 'summary/raw_vcf_file_status_genomewide.tsv')
    with open(status_path, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['sample','caller','path','status','records_chr1_8_pass_svlen50'])
        w.writerows(status_rows)
    clusters = cluster_records(records)
    rows = []
    for members in clusters:
        recs = [records[i] for i in members]
        samples_present = sorted({r['sample'] for r in recs})
        callers_by_sample = []
        for s in samples_present:
            callers = sorted({r['caller'] for r in recs if r['sample'] == s})
            callers_by_sample.append(f"{s}:{','.join(callers)}")
        intact_present = len(set(samples_present) & set(INTACT_GROUP))
        del_present = len(set(samples_present) & set(DEL_GROUP))
        intact_total = len(INTACT_GROUP); del_total = len(DEL_GROUP)
        intact_absent = intact_total - intact_present
        del_absent = del_total - del_present
        intact_freq = intact_present / intact_total
        del_freq = del_present / del_total
        if del_freq > intact_freq:
            direction = 'DEL_group_enriched'
        elif intact_freq > del_freq:
            direction = 'INTACT_group_enriched'
        else:
            direction = 'tie'
        p = fisher_two_sided(del_present, del_absent, intact_present, intact_absent)
        svtype = Counter(r['svtype'] for r in recs).most_common(1)[0][0]
        chrom = Counter(r['chrom'] for r in recs).most_common(1)[0][0]
        starts = [r['start'] for r in recs]
        stops = [r['stop'] for r in recs]
        poss = [r['pos'] for r in recs]
        ends = [r['end'] for r in recs]
        lens = [r['abs_svlen'] for r in recs]
        complete = (del_present == del_total and intact_present == 0) or (intact_present == intact_total and del_present == 0)
        group80 = (del_freq >= 0.8 and intact_freq <= 0.2) or (intact_freq >= 0.8 and del_freq <= 0.2)
        tgt_like = any(target_like(r) for r in recs)
        rows.append({
            'cluster_id': '', 'chrom': chrom, 'pos': int(median(poss)), 'end': int(median(ends)),
            'start': min(starts), 'stop': max(stops), 'svtype': svtype, 'svlen_median': int(median(lens)),
            'n_records': len(recs), 'n_samples_present': len(samples_present), 'samples_present': ','.join(samples_present),
            'callers_by_sample': ';'.join(callers_by_sample), 'intact_present': intact_present, 'intact_absent': intact_absent,
            'del_present': del_present, 'del_absent': del_absent, 'intact_freq': intact_freq, 'del_freq': del_freq,
            'direction': direction, 'fisher_p': p, 'neglog10p': -math.log10(max(p, 1e-300)),
            'complete_segregation': complete, 'group_specific80': group80, 'target_like_raw_cluster': tgt_like,
        })
    chrom_order = {c:i for i,c in enumerate(CHROMS)}
    rows.sort(key=lambda r: (chrom_order.get(r['chrom'], 999), r['pos'], r['svtype']))
    for i, r in enumerate(rows, 1):
        r['cluster_id'] = f"RawSVg_{i:07d}"
    out_table = os.path.join(OUT, 'results/readmapping_rawSV_genomewide.strict_clusters_with_fisher.tsv')
    fields = ['cluster_id','chrom','pos','end','start','stop','svtype','svlen_median','n_records','n_samples_present','samples_present','callers_by_sample','intact_present','intact_absent','del_present','del_absent','intact_freq','del_freq','direction','fisher_p','neglog10p','complete_segregation','group_specific80','target_like_raw_cluster']
    with open(out_table, 'w', newline='') as f:
        w = csv.DictWriter(f, delimiter='\t', fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    manual_path = os.path.join(OUT, 'results/Chr23997_manual_IGV_corrected_target_point.tsv')
    with open(manual_path, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['id','chrom','pos','end','svtype','svlen','intact_present','intact_absent','del_present','del_absent','fisher_p_core','neglog10p_core','note'])
        p = fisher_two_sided(len(DEL_GROUP), 0, 0, len(INTACT_GROUP))
        w.writerow([TARGET_ID, TARGET_CHROM, TARGET_POS, TARGET_END, 'DEL', TARGET_LEN, 0, len(INTACT_GROUP), len(DEL_GROUP), 0, p, -math.log10(max(p,1e-300)), 'Manual IGV-corrected Chr23997 intron-2 DEL: DEL only in deletion/spineless group, intact in spiny-like group'])
    with open(os.path.join(OUT, 'summary/run_summary_genomewide.tsv'), 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['item','value'])
        w.writerow(['analysis','raw_readmapping_per_species_no_panSV_genomewide'])
        w.writerow(['chromosomes', ','.join(CHROMS)])
        w.writerow(['raw_records', len(records)])
        w.writerow(['strict_clusters', len(rows)])
        w.writerow(['group_specific80_clusters', sum(1 for r in rows if r['group_specific80'])])
        w.writerow(['complete_segregation_clusters', sum(1 for r in rows if r['complete_segregation'])])
        w.writerow(['p_lt_0.05_clusters', sum(1 for r in rows if r['fisher_p'] < 0.05)])
        w.writerow(['target_like_raw_clusters', sum(1 for r in rows if r['target_like_raw_cluster'])])
        w.writerow(['Chr23997_manual_core_p', p])
        w.writerow(['Chr23997_manual_core_neglog10p', -math.log10(max(p,1e-300))])
        w.writerow(['matching_rule_DEL_DUP_INV', 'same SVTYPE+chrom; length ratio >= 0.5; reciprocal overlap >= 0.5 OR both breakpoints <= 500 bp'])
        w.writerow(['matching_rule_INS_BND_TRA', 'same SVTYPE+chrom; length ratio >= 0.5; breakpoint distance <= 500 bp'])
        w.writerow(['excluded_samples', ','.join(sorted(EXCLUDED))])
        w.writerow(['cluster_table', out_table])
    print(f"records={len(records)} clusters={len(rows)} out={out_table}")

if __name__ == '__main__':
    main()
