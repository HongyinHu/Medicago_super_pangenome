#!/usr/bin/env python3
import os, sys, re, argparse, collections, bisect

BIN = 100000

def parse_attrs(s):
    d = {}
    for part in s.strip().split(';'):
        if not part:
            continue
        if '=' in part:
            k,v = part.split('=',1)
            d[k] = v
    return d

def add_interval(idx, seq, start, end):
    if start > end:
        start, end = end, start
    b1 = start // BIN
    b2 = end // BIN
    for b in range(b1, b2+1):
        idx[(seq,b)].append((start,end))

def has_overlap(idx, seq, start, end, min_bp=1):
    b1 = start // BIN
    b2 = end // BIN
    for b in range(b1, b2+1):
        for a,c in idx.get((seq,b), ()):
            ov = min(end,c) - max(start,a) + 1
            if ov >= min_bp:
                return True
    return False

def build_support_indices(gene_predictions, transcript_alignments):
    br = collections.defaultdict(list)
    gm = collections.defaultdict(list)
    pa = collections.defaultdict(list)
    with open(gene_predictions) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            seq, src, typ = cols[0], cols[1], cols[2]
            if typ != 'gene':
                continue
            try:
                st, en = int(cols[3]), int(cols[4])
            except ValueError:
                continue
            if src == 'BRAKER3':
                add_interval(br, seq, st, en)
            elif src == 'GeMoMa':
                add_interval(gm, seq, st, en)
    with open(transcript_alignments) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            seq, src, typ = cols[0], cols[1], cols[2]
            if src == 'PASA' and typ == 'cDNA_match':
                try:
                    st, en = int(cols[3]), int(cols[4])
                except ValueError:
                    continue
                add_interval(pa, seq, st, en)
    return br, gm, pa

def read_fasta_keep(in_fa, out_fa, keep_ids):
    n = 0
    kept = 0
    write = False
    with open(in_fa) as inp, open(out_fa, 'w') as out:
        for line in inp:
            if line.startswith('>'):
                n += 1
                ident = line[1:].strip().split()[0]
                write = ident in keep_ids or ident.replace('cds.','',1) in keep_ids
                if write:
                    kept += 1
                    out.write(line)
            else:
                if write:
                    out.write(line)
    return n, kept

def filter_sample(base, sample, out_root):
    evm = os.path.join(base, '05_EVM_integration')
    final = os.path.join(evm, 'work', sample, '02_final')
    inputs = os.path.join(evm, 'work', sample, '00_inputs')
    gff = os.path.join(final, sample + '.evm.gff3')
    pep = os.path.join(final, sample + '.evm.pep.fa')
    cdsfa = os.path.join(final, sample + '.evm.cds.fa')
    gp = os.path.join(inputs, 'gene_predictions.gff3')
    ta = os.path.join(inputs, 'transcript_alignments.gff3')
    outdir = os.path.join(out_root, sample)
    os.makedirs(outdir, exist_ok=True)
    br_idx, gm_idx, pa_idx = build_support_indices(gp, ta)

    gene_info = {}
    mrna_to_gene = {}
    gene_to_mrna = collections.defaultdict(list)
    gene_cds_bp = collections.Counter()
    gene_cds_count = collections.Counter()
    feature_gene = []
    header_lines = []

    with open(gff) as fh:
        for line in fh:
            if not line.strip():
                feature_gene.append((line, None))
                continue
            if line.startswith('#'):
                header_lines.append(line)
                feature_gene.append((line, None))
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                feature_gene.append((line, None))
                continue
            attrs = parse_attrs(cols[8])
            typ = cols[2]
            gid = None
            if typ == 'gene':
                gid = attrs.get('ID')
                if gid:
                    st, en = int(cols[3]), int(cols[4])
                    gene_info[gid] = {'seq': cols[0], 'start': st, 'end': en, 'strand': cols[6]}
            elif typ == 'mRNA':
                mid = attrs.get('ID')
                parents = attrs.get('Parent','').split(',') if attrs.get('Parent') else []
                gid = parents[0] if parents else None
                if mid and gid:
                    mrna_to_gene[mid] = gid
                    gene_to_mrna[gid].append(mid)
            else:
                parents = attrs.get('Parent','').split(',') if attrs.get('Parent') else []
                if parents:
                    parent = parents[0]
                    gid = mrna_to_gene.get(parent)
                    if gid is None and parent.startswith('evm.model.'):
                        gid = parent.replace('evm.model.', 'evm.TU.', 1)
                if typ == 'CDS' and gid:
                    gene_cds_bp[gid] += int(cols[4]) - int(cols[3]) + 1
                    gene_cds_count[gid] += 1
            feature_gene.append((line, gid))

    keep = set()
    removed_rows = []
    class_counts = collections.Counter()
    for gid, inf in gene_info.items():
        seq, st, en = inf['seq'], inf['start'], inf['end']
        aa_len = gene_cds_bp[gid] // 3
        br = has_overlap(br_idx, seq, st, en, 1)
        gm = has_overlap(gm_idx, seq, st, en, 1)
        pa = has_overlap(pa_idx, seq, st, en, 1)
        support = ('B' if br else '') + ('G' if gm else '') + ('P' if pa else '')
        if not support:
            support = 'none'
        # v2 stricter model filter:
        # - strong path: >=100 aa with at least two evidence classes among BRAKER3, GeMoMa, PASA
        # - single BRAKER3/PASA rescue: >=200 aa
        # - GeMoMa-only rescue: >=500 aa
        # - small ORF rescue: 50-99 aa only when both BRAKER3 and PASA overlap
        decision = False
        reason = 'removed'
        support_count = (1 if br else 0) + (1 if gm else 0) + (1 if pa else 0)
        if aa_len >= 100 and support_count >= 2:
            decision = True; reason = 'keep_multisupport_ge100'
        elif aa_len >= 200 and (br or pa):
            decision = True; reason = 'keep_single_braker_or_pasa_ge200'
        elif aa_len >= 500 and gm:
            decision = True; reason = 'keep_gemoma_only_ge500'
        elif aa_len >= 50 and br and pa:
            decision = True; reason = 'keep_small_braker_pasa_ge50'
        else:
            reason = 'drop_single_source_or_short'
        class_counts[(support, reason)] += 1
        if decision:
            keep.add(gid)
        else:
            removed_rows.append((gid, seq, st, en, aa_len, gene_cds_count[gid], support, reason))

    keep_mrna = set()
    for gid in keep:
        keep_mrna.update(gene_to_mrna.get(gid, []))

    out_gff = os.path.join(outdir, sample + '.filtered_v2.evm.gff3')
    with open(out_gff, 'w') as out:
        for line, gid in feature_gene:
            if line.startswith('#') or not line.strip():
                out.write(line)
            elif gid in keep:
                out.write(line)
            else:
                # gene and mRNA features have gid resolved later? For mRNA lines seen before mapping is set during same pass;
                # if gid was missing on an mRNA due to order, infer from attributes here.
                cols = line.rstrip('\n').split('\t')
                attrs = parse_attrs(cols[8]) if len(cols) >= 9 else {}
                typ = cols[2] if len(cols) >= 3 else ''
                outgid = gid
                if typ == 'gene': outgid = attrs.get('ID')
                elif typ == 'mRNA': outgid = (attrs.get('Parent','').split(',')[0] if attrs.get('Parent') else None)
                elif attrs.get('Parent'):
                    outgid = mrna_to_gene.get(attrs.get('Parent','').split(',')[0], outgid)
                if outgid in keep:
                    out.write(line)

    out_pep = os.path.join(outdir, sample + '.filtered_v2.evm.pep.fa')
    out_cds = os.path.join(outdir, sample + '.filtered_v2.evm.cds.fa')
    _, kept_pep = read_fasta_keep(pep, out_pep, keep_mrna)
    _, kept_cds = read_fasta_keep(cdsfa, out_cds, keep_mrna)

    out_removed = os.path.join(outdir, sample + '.filtered_v2.removed.tsv')
    with open(out_removed, 'w') as out:
        out.write('gene_id\tseqid\tstart\tend\taa_len\tcds_features\tsupport\treason\n')
        for row in removed_rows:
            out.write('\t'.join(map(str,row)) + '\n')

    out_stats = os.path.join(outdir, sample + '.filtered_v2.stats.tsv')
    original = len(gene_info)
    kept = len(keep)
    with open(out_stats, 'w') as out:
        out.write('sample\toriginal_genes\tfiltered_genes\tremoved_genes\tpep\tcds\tremoved_pct\n')
        out.write('%s\t%d\t%d\t%d\t%d\t%d\t%.2f\n' % (sample, original, kept, original-kept, kept_pep, kept_cds, (100.0*(original-kept)/original if original else 0)))
        out.write('\n# support_reason_counts\n')
        out.write('support\treason\tcount\n')
        for (support, reason), count in sorted(class_counts.items()):
            out.write('%s\t%s\t%d\n' % (support, reason, count))
    return sample, original, kept, original-kept, kept_pep, kept_cds

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--base', required=True)
    ap.add_argument('--out-root', required=True)
    ap.add_argument('samples', nargs='*')
    args = ap.parse_args()
    samples = args.samples
    if not samples:
        status = os.path.join(args.base, '05_EVM_integration', 'status')
        samples = sorted(p[:-5] for p in os.listdir(status) if p.endswith('.done'))
    os.makedirs(args.out_root, exist_ok=True)
    summary = os.path.join(args.out_root, 'summary.tsv')
    with open(summary, 'w') as out:
        out.write('sample\toriginal_genes\tfiltered_genes\tremoved_genes\tpep\tcds\n')
        for s in samples:
            sys.stderr.write('filtering %s\n' % s); sys.stderr.flush()
            row = filter_sample(args.base, s, args.out_root)
            out.write('%s\t%d\t%d\t%d\t%d\t%d\n' % row)
            out.flush()
    print(summary)
if __name__ == '__main__':
    main()
