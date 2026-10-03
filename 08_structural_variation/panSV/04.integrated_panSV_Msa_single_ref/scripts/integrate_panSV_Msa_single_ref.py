#!/usr/bin/env python3
import csv
import shutil
import time
from pathlib import Path
from collections import defaultdict, Counter

BASE = Path('path/to/project/N_3.call_SV')
READ_BASE = BASE / '01.read_based_dualref_hifi' / '05_dualref_panSV_20260626'
SVG_BASE = BASE / '02.assembly_svgap_dualref'
OUT = BASE / '04.integrated_panSV_Msa_single_ref'
RESULTS = OUT / 'results'
SUMMARY = OUT / 'summary'

READ_PAV = READ_BASE / 'results/RefA/panSV/panSV.RefA.PAV.matrix.tsv'
SVG_DIR = SVG_BASE / 'Msa_ref/CombinedSV'

BIN_SIZE = 100000
DIST = 500
MAX_MATCHES_PER_SVGAP = 5
MAX_IDS_PER_READ = 20

EVENTS_OUT = RESULTS / 'Msa_single_ref.final_panSV.events.tsv'
PAV_OUT = RESULTS / 'Msa_single_ref.final_panSV.PAV.matrix.tsv'
EVIDENCE_OUT = RESULTS / 'Msa_single_ref.svgap_readbased.overlap_evidence.tsv'
STATS_OUT = RESULTS / 'Msa_single_ref.final_panSV.stats.tsv'
TYPE_OUT = RESULTS / 'Msa_single_ref.final_panSV.by_type_source.tsv'
SMALL_OUT = RESULTS / 'Msa_single_ref.excluded_svgap_small_indels.tsv'
TMP_SVGAP_ONLY = RESULTS / '.Msa_single_ref.svgap_only.tmp.tsv'
DONE = SUMMARY / 'Msa_single_ref.final_panSV.done'


def log(msg):
    print(f"[{time.strftime('%F %T')}] {msg}", flush=True)


def open_tsv(path):
    return open(path, newline='', encoding='utf-8', errors='replace')


def safe_int(x, default=None):
    try:
        if x is None or x == '' or x == '.':
            return default
        return int(float(x))
    except Exception:
        return default


def clean_str(x):
    return '.' if x is None or x == '' else str(x)


def norm_interval(start, end):
    if start is None or end is None:
        return None, None
    if start > end:
        start, end = end, start
    return start, end


def length_from(start, end, svlen=None):
    if svlen is not None:
        try:
            return abs(int(float(svlen)))
        except Exception:
            pass
    if start is None or end is None:
        return 0
    return abs(end - start) + 1


def len_ratio_ok(a, b, max_ratio=1.5):
    a = abs(a or 0)
    b = abs(b or 0)
    if a == 0 or b == 0:
        return True
    return max(a, b) / max(1, min(a, b)) <= max_ratio


def interval_overlap(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def norm_chrom(chrom):
    if not chrom or chrom == '.':
        return chrom
    if chrom.startswith('Msa.'):
        return chrom[4:]
    return chrom


def bins_for(start, end):
    if start is None or end is None:
        return []
    lo = max(0, start - DIST)
    hi = end + DIST
    return range(lo // BIN_SIZE, hi // BIN_SIZE + 1)


def sample_from_query_token(token):
    if not token or token == '.':
        return None
    token = token.split(',')[0].split('_')[0]
    prefix = token.split('.')[0]
    if not prefix:
        return None
    if prefix.startswith('G') and prefix[1:].isdigit():
        return 'genome_' + prefix[1:]
    return 'genome_' + prefix


def samples_from_merge_samples(text):
    samples = set()
    if not text or text == '.':
        return samples
    for tok in text.split(','):
        sample = sample_from_query_token(tok)
        if sample:
            samples.add(sample)
    return samples


def validate_inputs():
    missing = []
    for path in [READ_PAV, SVG_DIR]:
        if not path.exists():
            missing.append(str(path))
    if missing:
        raise SystemExit('Missing required inputs:\n' + '\n'.join(missing))


def backup_existing_outputs():
    existing = [p for p in [EVENTS_OUT, PAV_OUT, EVIDENCE_OUT, STATS_OUT, TYPE_OUT, SMALL_OUT, DONE] if p.exists()]
    if not existing:
        if TMP_SVGAP_ONLY.exists():
            TMP_SVGAP_ONLY.unlink()
        return None
    bdir = OUT / '99_rerun_backup' / time.strftime('%Y%m%d_%H%M%S')
    bdir.mkdir(parents=True, exist_ok=True)
    for path in existing:
        shutil.move(str(path), str(bdir / path.name))
    if TMP_SVGAP_ONLY.exists():
        TMP_SVGAP_ONLY.unlink()
    return bdir


def read_readbased_events():
    log('loading read-based Msa single-ref PAV/events')
    events = []
    by_id = {}
    index = defaultdict(list)
    with open_tsv(READ_PAV) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        sample_cols = reader.fieldnames[6:]
        for row in reader:
            rid = row['SV_ID']
            start, end = norm_interval(safe_int(row.get('POS')), safe_int(row.get('END')))
            svtype = row.get('SVTYPE', '.')
            svlen = safe_int(row.get('SVLEN'), length_from(start, end))
            rec = {
                'row': row,
                'pav': [row.get(s, 'NA') for s in sample_cols],
                'svgap_count': 0,
                'svgap_ids': [],
                'svgap_samples': set(),
            }
            events.append(rec)
            by_id[rid] = rec
            if start is None or end is None:
                continue
            types_to_index = [svtype]
            if svtype == 'DUP':
                types_to_index.append('CNV')
            item = (rid, start, end, abs(svlen or length_from(start, end)), svtype)
            chrom = norm_chrom(row.get('CHROM'))
            for indexed_type in types_to_index:
                for b in bins_for(start, end):
                    index[(chrom, indexed_type, b)].append(item)
    log(f'loaded read events={len(events)} samples={len(sample_cols)} index_bins={len(index)}')
    return sample_cols, events, by_id, index


def match_read_events(index, chrom, start, end, svtype, svlen):
    if svtype == 'TRL':
        return []
    chrom = norm_chrom(chrom)
    candidates = {}
    for b in bins_for(start, end):
        for item in index.get((chrom, svtype, b), []):
            candidates[item[0]] = item
    matches = []
    qlen = max(1, abs(svlen or length_from(start, end)))
    for rid, rstart, rend, rlen, _rtype in candidates.values():
        if svtype == 'INS':
            dist = abs(start - rstart)
            if dist <= DIST and len_ratio_ok(qlen, rlen):
                matches.append((dist, rid, 'ins_breakpoint_len'))
        else:
            ov = interval_overlap(start, end, rstart, rend)
            alen = max(1, end - start + 1)
            blen = max(1, rend - rstart + 1)
            ro_ok = ov > 0 and ov / alen >= 0.5 and ov / blen >= 0.5
            end_ok = abs(start - rstart) <= DIST and abs(end - rend) <= DIST and len_ratio_ok(qlen, rlen)
            if ro_ok or end_ok:
                score = -ov if ov > 0 else abs(start - rstart) + abs(end - rend)
                reason = 'reciprocal_overlap' if ov > 0 else 'endpoint_distance_len'
                matches.append((score, rid, reason))
    matches.sort(key=lambda x: x[0])
    return matches[:MAX_MATCHES_PER_SVGAP]


def iter_svgap_large_events():
    files = [
        (SVG_DIR / 'All.DELs.50bplarge.bed.combined.sorted.txt', 'combined'),
        (SVG_DIR / 'All.INTs.50bplarge.bed.combined.sorted.txt', 'combined'),
        (SVG_DIR / 'All.CNVs.bed', 'raw_cnv'),
        (SVG_DIR / 'All.INVs.bed', 'raw_inv'),
        (SVG_DIR / 'All.TRLs.bed', 'raw_trl'),
    ]
    seq = 0
    for path, kind in files:
        log(f'processing SVGAP Msa {path.name}')
        seen_raw = set()
        with open_tsv(path) as fh:
            for line in fh:
                line = line.rstrip('\n')
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                source_file = 'Msa_ref/CombinedSV/' + path.name
                if kind == 'combined':
                    if len(parts) < 15:
                        continue
                    chrom = parts[0]
                    start = safe_int(parts[1])
                    end = safe_int(parts[2])
                    source_id = parts[3]
                    svtype = parts[4]
                    svlen = abs(safe_int(parts[5], 0) or 0)
                    merge_n = safe_int(parts[6], 0) or 0
                    support_samples = samples_from_merge_samples(parts[7])
                elif kind == 'raw_cnv':
                    if len(parts) < 8:
                        continue
                    key = tuple(parts[:8])
                    if key in seen_raw:
                        continue
                    seen_raw.add(key)
                    query, q_start, q_end = parts[0], parts[1], parts[2]
                    chrom = parts[4]
                    start = safe_int(parts[5])
                    end = safe_int(parts[6])
                    svtype = 'CNV'
                    svlen = abs(safe_int(parts[7], 0) or length_from(start, end))
                    source_id = f'{query}_{q_start}_{q_end}_{svlen}'
                    support_samples = {sample_from_query_token(query)} - {None}
                    merge_n = len(support_samples)
                elif kind == 'raw_inv':
                    if len(parts) < 8:
                        continue
                    key = tuple(parts[:8])
                    if key in seen_raw:
                        continue
                    seen_raw.add(key)
                    chrom = parts[0]
                    start = safe_int(parts[1])
                    end = safe_int(parts[2])
                    query, q_start, q_end = parts[4], parts[5], parts[6]
                    svtype = 'INV'
                    svlen = abs(safe_int(parts[7], 0) or length_from(start, end))
                    source_id = f'{query}_{q_start}_{q_end}_{svlen}'
                    support_samples = {sample_from_query_token(query)} - {None}
                    merge_n = len(support_samples)
                elif kind == 'raw_trl':
                    if len(parts) < 9:
                        continue
                    key = tuple(parts[:9])
                    if key in seen_raw:
                        continue
                    seen_raw.add(key)
                    query, q_start, q_end = parts[0], parts[1], parts[2]
                    chrom = parts[3]
                    start = safe_int(parts[4])
                    end = safe_int(parts[5])
                    svtype = 'TRL'
                    svlen = max(abs(safe_int(parts[7], 0) or 0), abs(safe_int(parts[8], 0) or 0), length_from(start, end))
                    source_id = f'{query}_{q_start}_{q_end}_{chrom}_{start}_{end}'
                    support_samples = {sample_from_query_token(query)} - {None}
                    merge_n = len(support_samples)
                else:
                    continue
                start, end = norm_interval(start, end)
                if start is None or end is None or not chrom:
                    continue
                if svlen < 50 and svtype in {'DEL', 'INS'}:
                    continue
                seq += 1
                yield {
                    'svgap_id': f'SVGAP_Msa_{svtype}_{seq:09d}',
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'svtype': svtype,
                    'svlen': svlen,
                    'source_id': source_id,
                    'source_file': source_file,
                    'support_samples': support_samples,
                    'merge_n': merge_n,
                }


def final_confidence_for_read(svgap_count):
    return 'read_svgap_supported' if svgap_count > 0 else 'read_based'


def final_confidence_for_svgap(svtype, support_n):
    if svtype in {'DEL', 'INS'} and support_n >= 2:
        return 'medium_svgap_combined'
    if svtype in {'CNV', 'INV', 'TRL'}:
        return 'medium_svgap_only'
    return 'low_svgap_only'


def process_svgap(events, by_id, index, sample_cols):
    log('streaming SVGAP Msa large SV events and matching read-based catalog')
    stats = Counter()
    evidence_n = 0
    svgap_only_n = 0
    with open(TMP_SVGAP_ONLY, 'w', newline='', encoding='utf-8') as tmp_fh, open(EVIDENCE_OUT, 'w', newline='', encoding='utf-8') as ev_fh:
        tmp_w = csv.writer(tmp_fh, delimiter='\t')
        ev_w = csv.writer(ev_fh, delimiter='\t')
        tmp_w.writerow(['chrom', 'start', 'end', 'svtype', 'svlen', 'svgap_id', 'svgap_source_id', 'svgap_source_file', 'support_samples', 'pav_values', 'final_confidence'])
        ev_w.writerow(['svgap_id', 'read_SV_ID', 'chrom', 'svgap_start', 'svgap_end', 'read_start', 'read_end', 'svtype', 'match_reason', 'svgap_source_file', 'svgap_source_id'])
        for i, ev in enumerate(iter_svgap_large_events(), 1):
            stats['svgap_large_processed'] += 1
            stats[f'svgap_large_{ev["svtype"]}'] += 1
            matches = match_read_events(index, ev['chrom'], ev['start'], ev['end'], ev['svtype'], ev['svlen'])
            if matches:
                for _score, rid, reason in matches:
                    rec = by_id.get(rid)
                    if rec is None:
                        continue
                    rec['svgap_count'] += 1
                    if len(rec['svgap_ids']) < MAX_IDS_PER_READ:
                        rec['svgap_ids'].append(ev['svgap_id'])
                    rec['svgap_samples'].update(ev['support_samples'])
                    row = rec['row']
                    ev_w.writerow([ev['svgap_id'], rid, ev['chrom'], ev['start'], ev['end'], row.get('POS', '.'), row.get('END', '.'), ev['svtype'], reason, ev['source_file'], ev['source_id']])
                    evidence_n += 1
            else:
                present = ev['support_samples']
                pav_values = ['1' if s in present else '0' for s in sample_cols]
                tmp_w.writerow([
                    ev['chrom'], ev['start'], ev['end'], ev['svtype'], ev['svlen'], ev['svgap_id'], ev['source_id'],
                    ev['source_file'], ','.join(sorted(present)) or '.', ','.join(pav_values), final_confidence_for_svgap(ev['svtype'], ev['merge_n'])
                ])
                svgap_only_n += 1
            if i % 100000 == 0:
                log(f'processed SVGAP events={i} svgap_only={svgap_only_n} evidence={evidence_n}')
    stats['svgap_only_main'] = svgap_only_n
    stats['svgap_read_overlap_evidence'] = evidence_n
    stats['read_events_with_svgap_support'] = sum(1 for event in events if event['svgap_count'] > 0)
    log(f'finished SVGAP stream processed={stats["svgap_large_processed"]} svgap_only={svgap_only_n} read_supported={stats["read_events_with_svgap_support"]} evidence={evidence_n}')
    return stats


def write_final_outputs(events, sample_cols, stats):
    log('writing Msa single-ref final catalog and PAV matrix')
    event_header = [
        'FinalSV_ID', 'source_methods', 'reference', 'CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN',
        'read_SV_ID', 'svgap_support_count', 'svgap_ids', 'svgap_support_samples', 'svgap_source_id',
        'svgap_source_file', 'final_confidence', 'notes'
    ]
    pav_header = ['FinalSV_ID', 'source_methods', 'reference', 'CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN', 'read_SV_ID', 'svgap_ids', 'final_confidence'] + sample_cols
    final_n = 0
    source_type_counter = Counter()
    with open(EVENTS_OUT, 'w', newline='', encoding='utf-8') as ev_fh, open(PAV_OUT, 'w', newline='', encoding='utf-8') as pav_fh:
        ev_w = csv.writer(ev_fh, delimiter='\t')
        pav_w = csv.writer(pav_fh, delimiter='\t')
        ev_w.writerow(event_header)
        pav_w.writerow(pav_header)
        for rec in events:
            final_n += 1
            fid = f'Msa_panSV_{final_n:09d}'
            row = rec['row']
            source = 'read_based+SVGAP' if rec['svgap_count'] > 0 else 'read_based'
            svgap_ids = ','.join(rec['svgap_ids']) if rec['svgap_ids'] else '.'
            svgap_samples = ','.join(sorted(rec['svgap_samples'])) if rec['svgap_samples'] else '.'
            conf = final_confidence_for_read(rec['svgap_count'])
            ev_w.writerow([
                fid, source, 'genome_Msa', clean_str(row.get('CHROM')), clean_str(row.get('POS')), clean_str(row.get('END')),
                clean_str(row.get('SVTYPE')), clean_str(row.get('SVLEN')), clean_str(row.get('SV_ID')),
                rec['svgap_count'], svgap_ids, svgap_samples, '.', '.', conf,
                'read_based_Msa_single_ref_catalog;svgap_overlap_support_added' if rec['svgap_count'] > 0 else 'read_based_Msa_single_ref_catalog'
            ])
            pav_w.writerow([
                fid, source, 'genome_Msa', clean_str(row.get('CHROM')), clean_str(row.get('POS')), clean_str(row.get('END')),
                clean_str(row.get('SVTYPE')), clean_str(row.get('SVLEN')), clean_str(row.get('SV_ID')), svgap_ids, conf
            ] + rec['pav'])
            source_type_counter[(source, clean_str(row.get('SVTYPE')))] += 1
        with open(TMP_SVGAP_ONLY, newline='', encoding='utf-8', errors='replace') as tmp_fh:
            reader = csv.DictReader(tmp_fh, delimiter='\t')
            for row in reader:
                final_n += 1
                fid = f'Msa_panSV_{final_n:09d}'
                pav_values = row['pav_values'].split(',') if row['pav_values'] else ['0'] * len(sample_cols)
                ev_w.writerow([
                    fid, 'SVGAP', 'genome_Msa', row['chrom'], row['start'], row['end'], row['svtype'], row['svlen'],
                    '.', 1, row['svgap_id'], row['support_samples'], row['svgap_source_id'], row['svgap_source_file'],
                    row['final_confidence'], 'assembly_svgap_Msa_single_ref_only_no_read_based_overlap;main_panSV_SVLEN_ge_50'
                ])
                pav_w.writerow([
                    fid, 'SVGAP', 'genome_Msa', row['chrom'], row['start'], row['end'], row['svtype'], row['svlen'],
                    '.', row['svgap_id'], row['final_confidence']
                ] + pav_values)
                source_type_counter[('SVGAP', row['svtype'])] += 1
    stats['final_total'] = final_n
    stats['final_read_based_rows'] = len(events)
    stats['final_svgap_only_rows'] = max(0, final_n - len(events))
    with open(TYPE_OUT, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['source_methods', 'SVTYPE', 'count'])
        for (source, svtype), n in sorted(source_type_counter.items()):
            writer.writerow([source, svtype, n])
    log(f'wrote final rows={final_n}')


def count_small_excluded():
    small_files = [
        ('DEL', SVG_DIR / 'All.DELs.50bpsmall.bed.combined.sorted.txt'),
        ('INS', SVG_DIR / 'All.INTs.50bpsmall.bed.combined.sorted.txt'),
    ]
    with open(SMALL_OUT, 'w', newline='', encoding='utf-8') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['reference', 'svtype', 'source_file', 'excluded_reason', 'records_excluding_header'])
        for svtype, path in small_files:
            n = 0
            with open(path, encoding='utf-8', errors='replace') as fh:
                for line in fh:
                    if line.strip() and not line.startswith('#'):
                        n += 1
            writer.writerow(['genome_Msa', svtype, str(path), 'SVLEN_lt_50_small_indel_not_in_main_panSV', n])
            log(f'excluded small Msa {svtype} records={n}')


def write_stats(stats, sample_cols):
    with open(STATS_OUT, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['section', 'key', 'count'])
        writer.writerow(['input', 'read_Msa_single_ref_events', len(open(READ_PAV, encoding='utf-8', errors='replace').read().splitlines()) - 1])
        writer.writerow(['input', 'samples', len(sample_cols)])
        for key, value in sorted(stats.items()):
            writer.writerow(['integrated', key, value])
        writer.writerow(['parameters', 'reference', 'genome_Msa'])
        writer.writerow(['parameters', 'main_panSV_min_svlen', 50])
        writer.writerow(['parameters', 'read_svgap_overlap_distance_bp', DIST])
        writer.writerow(['parameters', 'read_svgap_reciprocal_overlap', 0.5])
        writer.writerow(['parameters', 'max_matches_per_svgap_event', MAX_MATCHES_PER_SVGAP])
        writer.writerow(['paths', 'read_pav', str(READ_PAV)])
        writer.writerow(['paths', 'svgap_msa_combinedsv', str(SVG_DIR)])
    DONE.write_text(time.strftime('%F %T') + '\n', encoding='utf-8')


def write_readme():
    readme = SUMMARY / 'Msa_single_ref.final_panSV.README.txt'
    readme.write_text(
        'Msa single-reference final panSV integration\n'
        f'Date: {time.strftime("%F %T")}\n\n'
        'Inputs:\n'
        f'- Read-based Msa single-ref PAV/events: {READ_PAV}\n'
        f'- SVGAP Msa_ref CombinedSV: {SVG_DIR}\n\n'
        'Rules:\n'
        '- Main final panSV uses SVLEN >= 50 bp.\n'
        '- Read-based RefA/Msa panSV records are retained as the primary read-backed single-reference catalog.\n'
        '- SVGAP Msa_ref large DEL/INS combined.sorted.txt plus CNV/INV/TRL records are integrated.\n'
        '- SVGAP 50bpsmall DEL/INS records are excluded from the main panSV and counted separately.\n'
        '- SVGAP Msa. chromosome prefixes are normalized only for matching; original coordinates are retained in outputs.\n'
        '- Read-based and SVGAP records are linked by same chromosome/SVTYPE and reciprocal overlap >= 0.5 or endpoint/breakpoint distance <= 500 bp with compatible length.\n\n'
        'Main outputs:\n'
        f'- {EVENTS_OUT}\n'
        f'- {PAV_OUT}\n'
        f'- {EVIDENCE_OUT}\n'
        f'- {STATS_OUT}\n',
        encoding='utf-8'
    )


def main():
    t0 = time.time()
    RESULTS.mkdir(parents=True, exist_ok=True)
    SUMMARY.mkdir(parents=True, exist_ok=True)
    validate_inputs()
    backup_dir = backup_existing_outputs()
    if backup_dir:
        log(f'backed up existing outputs to {backup_dir}')
    sample_cols, events, by_id, index = read_readbased_events()
    stats = process_svgap(events, by_id, index, sample_cols)
    write_final_outputs(events, sample_cols, stats)
    count_small_excluded()
    write_stats(stats, sample_cols)
    write_readme()
    if TMP_SVGAP_ONLY.exists():
        TMP_SVGAP_ONLY.unlink()
    log(f'done elapsed_sec={int(time.time() - t0)}')


if __name__ == '__main__':
    main()
