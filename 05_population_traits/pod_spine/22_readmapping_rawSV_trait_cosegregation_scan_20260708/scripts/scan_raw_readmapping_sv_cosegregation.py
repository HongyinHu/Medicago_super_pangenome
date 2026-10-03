#!/usr/bin/env python3
import os, sys, gzip, math, re, csv
from collections import defaultdict, Counter
from itertools import combinations
from statistics import median

RUN = "path/to/project/N_3.call_SV/01.read_based_dualref_hifi"
OUT = "path/to/project/N_4.pod_spiny/22_readmapping_rawSV_trait_cosegregation_scan_20260708"
REGION_CHROM = "Chr4"
REGION_START = 86_000_000
REGION_END = 92_000_000
MIN_SVLEN = 50
BP_DIST = 500
INS_DIST = 500
MIN_RO = 0.50
LEN_RATIO_MIN = 0.50
LEN_RATIO_MAX = 2.00

# Manual phenotype/truth grouping for the Chr23997 validation context.
# These are intentionally not taken from panSV.
INTACT_GROUP = ["genome_410", "genome_436", "genome_457", "genome_454", "genome_474", "genome_Mpo", "genome_R108", "genome_A17"]
DEL_GROUP = ["genome_395", "genome_461", "genome_468", "genome_472", "genome_M22", "genome_Mar", "genome_Mru"]
EXCLUDED = {"genome_Msa", "genome_M46", "genome_482", "genome_ZM4"}
SAMPLES = INTACT_GROUP + DEL_GROUP
CALLERS = {
    "pbsv": os.path.join(RUN, "03_per_sample/Msa/pbsv/{sample}.pbsv.vcf.gz"),
    "sniffles2": os.path.join(RUN, "03_per_sample/Msa/sniffles2/{sample}.sniffles2.vcf.gz"),
    "cuteSV": os.path.join(RUN, "03_per_sample/Msa/cutesv/{sample}.cutesv.vcf.gz"),
}

# Manual IGV-corrected truth for the candidate second-intron DEL in Chr23997.
TARGET_CHROM = "Chr4"
TARGET_POS = 88_980_831
TARGET_END = 88_981_052
TARGET_LEN = 221
TARGET_ID = "Chr23997_second_intron_DEL_manual_IGV_corrected"


def open_vcf(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt', errors='replace')
    return open(path, 'rt', errors='replace')

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
    if chrom != REGION_CHROM:
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
    # overlap target plotting region
    rec_start, rec_end = min(pos, end), max(pos, end)
    if rec_end < REGION_START or rec_start > REGION_END:
        return None
    if svtype not in ('BND', 'TRA') and svlen < MIN_SVLEN:
        return None
    if svtype in ('BND', 'TRA'):
        svlen = max(svlen, 50)
    return {
        'sample': sample,
        'caller': caller,
        'chrom': chrom,
        'pos': pos,
        'end': end,
        'start': rec_start,
        'stop': rec_end,
        'id': rid if rid else '.',
        'alt': alt,
        'svtype': svtype,
        'svlen': svlen,
        'abs_svlen': abs(svlen),
        'filter': flt,
        'raw_info': info_s,
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
    if a['sample'] == b['sample'] and a['caller'] == b['caller'] and a['id'] == b['id'] and a['pos'] == b['pos']:
        return True
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
    # table [[a,b],[c,d]]; row/column sums fixed.
    # Uses lgamma so it works on older Python versions without math.comb.
    n1 = a + b; n2 = c + d; m1 = a + c; n = n1 + n2
    def logchoose(nv, kv):
        if kv < 0 or kv > nv:
            return float('-inf')
        return math.lgamma(nv + 1) - math.lgamma(kv + 1) - math.lgamma(nv - kv + 1)
    def prob(x):
        if x < 0 or x > n1 or x > m1 or (n1 - x) > (n - m1):
            return 0.0
        lp = logchoose(m1, x) + logchoose(n - m1, n1 - x) - logchoose(n, n1)
        return math.exp(lp)
    obs = prob(a)
    lo = max(0, n1 - (n - m1))
    hi = min(n1, m1)
    p = sum(prob(x) for x in range(lo, hi + 1) if prob(x) <= obs + 1e-15)
    return min(1.0, p)

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
    # Sweep by position to avoid all-vs-all across chromosome.
    for key, idxs in groups.items():
        idxs.sort(key=lambda i: records[i]['start'])
        active = []
        for i in idxs:
            r = records[i]
            # keep possible overlapping/bp-near records
            active = [j for j in active if records[j]['stop'] >= r['start'] - BP_DIST * 2]
            for j in active:
                if same_event(r, records[j]):
                    union(i, j)
            active.append(i)
    clusters = defaultdict(list)
    for i in range(len(records)):
        clusters[find(i)].append(i)
    return list(clusters.values())

def summarize_cluster(cid, idxs, records):
    recs = [records[i] for i in idxs]
    samples_present = sorted(set(r['sample'] for r in recs))
    callers_by_sample = {s: sorted(set(r['caller'] for r in recs if r['sample'] == s)) for s in samples_present}
    poss = [r['pos'] for r in recs]
    ends = [r['end'] for r in recs]
    svlens = [r['abs_svlen'] for r in recs if r['abs_svlen']]
    svtype = Counter(r['svtype'] for r in recs).most_common(1)[0][0]
    chrom = recs[0]['chrom']
    pos = int(round(median(poss)))
    end = int(round(median(ends)))
    start, stop = min(pos, end), max(pos, end)
    svlen = int(round(median(svlens))) if svlens else max(1, stop - start + 1)
    intact_present = sum(1 for s in INTACT_GROUP if s in samples_present)
    intact_absent = len(INTACT_GROUP) - intact_present
    del_present = sum(1 for s in DEL_GROUP if s in samples_present)
    del_absent = len(DEL_GROUP) - del_present
    p = fisher_two_sided(intact_present, intact_absent, del_present, del_absent)
    f_intact = intact_present / len(INTACT_GROUP)
    f_del = del_present / len(DEL_GROUP)
    complete = ((intact_present == 0 and del_present == len(DEL_GROUP)) or
                (intact_present == len(INTACT_GROUP) and del_present == 0))
    group_specific80 = ((f_intact <= 0.2 and f_del >= 0.8) or (f_intact >= 0.8 and f_del <= 0.2))
    direction = 'DEL_group_enriched' if f_del > f_intact else ('INTACT_group_enriched' if f_intact > f_del else 'tie')
    # Does this cluster represent the target neighborhood? strict match against target interval.
    target_rec = {'chrom': TARGET_CHROM, 'svtype': 'DEL', 'pos': TARGET_POS, 'end': TARGET_END,
                  'start': TARGET_POS, 'stop': TARGET_END, 'abs_svlen': TARGET_LEN, 'sample': 'target', 'caller': 'manual', 'id': 'target'}
    target_like = any(same_event(r, target_rec) for r in recs)
    return {
        'cluster_id': f'RawSV_{cid:06d}',
        'chrom': chrom,
        'pos': pos,
        'end': end,
        'start': start,
        'stop': stop,
        'svtype': svtype,
        'svlen_median': svlen,
        'n_records': len(recs),
        'n_samples_present': len(samples_present),
        'samples_present': ','.join(samples_present),
        'callers_by_sample': ';'.join(f'{s}:{"|".join(callers_by_sample[s])}' for s in samples_present),
        'intact_present': intact_present,
        'intact_absent': intact_absent,
        'del_present': del_present,
        'del_absent': del_absent,
        'intact_freq': f_intact,
        'del_freq': f_del,
        'direction': direction,
        'fisher_p': p,
        'neglog10p': -math.log10(max(p, 1e-300)),
        'complete_segregation': complete,
        'group_specific80': group_specific80,
        'target_like_raw_cluster': target_like,
    }

def write_tsv(path, rows, fields):
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow(r)

def main():
    records = []
    file_status = []
    for sample in SAMPLES:
        for caller, pattern in CALLERS.items():
            path = pattern.format(sample=sample)
            status = {'sample': sample, 'caller': caller, 'path': path, 'exists': os.path.exists(path), 'records_in_region': 0}
            if not os.path.exists(path):
                file_status.append(status)
                continue
            with open_vcf(path) as fh:
                for line in fh:
                    if not line or line.startswith('#'):
                        continue
                    rec = record_from_vcf_line(line, sample, caller)
                    if rec:
                        records.append(rec)
                        status['records_in_region'] += 1
            file_status.append(status)
    rec_fields = ['sample','caller','chrom','pos','end','start','stop','id','svtype','svlen','abs_svlen','filter','alt']
    write_tsv(os.path.join(OUT, 'results/raw_readmapping_SV_records_Chr4_86_92Mb.tsv'), records, rec_fields)
    write_tsv(os.path.join(OUT, 'summary/raw_vcf_file_status.tsv'), file_status, ['sample','caller','path','exists','records_in_region'])
    clusters = cluster_records(records)
    summaries = [summarize_cluster(i+1, idxs, records) for i, idxs in enumerate(clusters)]
    summaries.sort(key=lambda r: (r['chrom'], r['start'], r['svtype'], r['svlen_median']))
    fields = ['cluster_id','chrom','pos','end','start','stop','svtype','svlen_median','n_records','n_samples_present','samples_present','callers_by_sample','intact_present','intact_absent','del_present','del_absent','intact_freq','del_freq','direction','fisher_p','neglog10p','complete_segregation','group_specific80','target_like_raw_cluster']
    write_tsv(os.path.join(OUT, 'results/readmapping_rawSV_Chr4_86_92Mb.strict_clusters_with_fisher.tsv'), summaries, fields)
    # Manual corrected target point and support table from raw calls near the target.
    support_rows = []
    for sample in SAMPLES:
        for caller, pattern in CALLERS.items():
            path = pattern.format(sample=sample)
            hits = []
            if os.path.exists(path):
                with open_vcf(path) as fh:
                    for line in fh:
                        if line.startswith('#'):
                            continue
                        rec = record_from_vcf_line(line, sample, caller)
                        if rec and rec['chrom'] == TARGET_CHROM and rec['svtype'] == 'DEL':
                            if rec['stop'] >= TARGET_POS - 2000 and rec['start'] <= TARGET_END + 2000:
                                target_rec = {'chrom': TARGET_CHROM, 'svtype': 'DEL', 'pos': TARGET_POS, 'end': TARGET_END,
                                              'start': TARGET_POS, 'stop': TARGET_END, 'abs_svlen': TARGET_LEN,
                                              'sample': 'target', 'caller': 'manual', 'id': 'target'}
                                ro = reciprocal_overlap(rec, target_rec)
                                bp = max(abs(rec['start'] - TARGET_POS), abs(rec['stop'] - TARGET_END))
                                ratio = min(max(1, rec['abs_svlen']), TARGET_LEN) / max(max(1, rec['abs_svlen']), TARGET_LEN)
                                hits.append(f"{rec['id']}:{rec['start']}-{rec['stop']}:len{rec['abs_svlen']}:bp{bp}:ro{ro:.2f}:ratio{ratio:.2f}")
            support_rows.append({'sample': sample, 'caller': caller, 'nearby_DEL_hits': len(hits), 'hit_details': '|'.join(hits)})
    write_tsv(os.path.join(OUT, 'results/Chr23997_target_region_raw_DEL_support_by_species_caller.tsv'), support_rows, ['sample','caller','nearby_DEL_hits','hit_details'])
    manual_p_core = fisher_two_sided(0, len(INTACT_GROUP), len(DEL_GROUP), 0)
    manual_p_with_zm4 = fisher_two_sided(0, len(INTACT_GROUP), len(DEL_GROUP), 1)  # spineless includes one intact exception
    manual = [{
        'target_id': TARGET_ID,
        'chrom': TARGET_CHROM,
        'pos': TARGET_POS,
        'end': TARGET_END,
        'svtype': 'DEL',
        'svlen': TARGET_LEN,
        'manual_intact_group_n': len(INTACT_GROUP),
        'manual_intact_group_DEL_present': 0,
        'manual_del_group_n': len(DEL_GROUP),
        'manual_del_group_DEL_present': len(DEL_GROUP),
        'fisher_p_core': manual_p_core,
        'neglog10p_core': -math.log10(manual_p_core),
        'fisher_p_with_ZM4_exception': manual_p_with_zm4,
        'neglog10p_with_ZM4_exception': -math.log10(manual_p_with_zm4),
        'note': 'Manual IGV truth: intact group has no target second-intron DEL; core spineless group has target DEL. genome_ZM4 is a spineless exception with intact locus and is not used in the core scan.'
    }]
    write_tsv(os.path.join(OUT, 'results/Chr23997_manual_IGV_corrected_target_point.tsv'), manual, list(manual[0].keys()))
    total = len(summaries)
    n80 = sum(1 for r in summaries if r['group_specific80'])
    ncomp = sum(1 for r in summaries if r['complete_segregation'])
    np05 = sum(1 for r in summaries if r['fisher_p'] < 0.05)
    target_clusters = [r for r in summaries if r['target_like_raw_cluster']]
    with open(os.path.join(OUT, 'summary/run_summary.tsv'), 'w') as fh:
        fh.write('item\tvalue\n')
        fh.write(f'analysis\traw_readmapping_per_species_no_panSV\n')
        fh.write(f'region\t{REGION_CHROM}:{REGION_START}-{REGION_END}\n')
        fh.write(f'raw_records\t{len(records)}\n')
        fh.write(f'strict_clusters\t{total}\n')
        fh.write(f'group_specific80_clusters\t{n80}\n')
        fh.write(f'complete_segregation_clusters\t{ncomp}\n')
        fh.write(f'p_lt_0.05_clusters\t{np05}\n')
        fh.write(f'target_like_raw_clusters\t{len(target_clusters)}\n')
        fh.write(f'Chr23997_manual_core_p\t{manual_p_core}\n')
        fh.write(f'Chr23997_manual_core_neglog10p\t{-math.log10(manual_p_core)}\n')
        fh.write(f'matching_rule_DEL_DUP_INV\tsame SVTYPE+chrom; length ratio >= {LEN_RATIO_MIN}; reciprocal overlap >= {MIN_RO} OR both breakpoints <= {BP_DIST} bp\n')
        fh.write(f'matching_rule_INS_BND_TRA\tsame SVTYPE+chrom; length ratio >= {LEN_RATIO_MIN}; breakpoint distance <= {INS_DIST} bp\n')
        fh.write(f'excluded_samples\t{",".join(sorted(EXCLUDED))}\n')
    with open(os.path.join(OUT, 'summary/README.md'), 'w') as fh:
        fh.write('# Raw read-mapping SV co-segregation scan, no panSV matrix\n\n')
        fh.write('This run re-scans the Chr4 86-92 Mb window using raw per-species read-mapping SV VCFs under the genome_Msa reference. It does not use the panSV PAV matrix.\n\n')
        fh.write('Groups used for the primary validation scan:\n\n')
        fh.write('- Intact/spiny-like group: ' + ', '.join(INTACT_GROUP) + '\n')
        fh.write('- Core DEL/spineless group: ' + ', '.join(DEL_GROUP) + '\n')
        fh.write('- Excluded from primary Fisher test: genome_Msa self, genome_M46, genome_482, genome_ZM4 exception.\n\n')
        fh.write('Strict event matching rules are recorded in summary/run_summary.tsv.\n')

if __name__ == '__main__':
    main()
