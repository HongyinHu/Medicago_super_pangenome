#!/usr/bin/env python3
import argparse
import csv
import gzip
import os
from collections import Counter, defaultdict


INTERVAL_TYPES = {"DEL", "DUP", "INV"}
MATCHABLE_TYPES = {"DEL", "DUP", "INV", "INS"}
RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def parse_info(info):
    out = {}
    for item in info.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value
    return out


def normalize_ins_seq(ref, alt):
    if not alt or alt.startswith("<"):
        return None
    seq = alt.upper()
    ref = (ref or "").upper()
    if ref and seq.startswith(ref) and len(seq) > len(ref):
        seq = seq[len(ref) :]
    elif len(seq) > 1:
        seq = seq[1:]
    seq = "".join(base for base in seq if base in "ACGTN")
    return seq if len(seq) >= 30 else None


def read_vcf(path):
    records = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            info = parse_info(fields[7])
            svtype = info.get("SVTYPE", "NA").split(",", 1)[0]
            try:
                pos = int(fields[1])
                end = int(info.get("END", fields[1]))
            except ValueError:
                continue
            if end < pos:
                pos, end = end, pos
            svlen_s = str(info.get("SVLEN", end - pos)).split(",", 1)[0]
            try:
                svlen = abs(int(float(svlen_s)))
            except ValueError:
                svlen = max(1, abs(end - pos))
            seq = normalize_ins_seq(fields[3], fields[4]) if svtype == "INS" else None
            records.append(
                {
                    "id": fields[2],
                    "chrom": fields[0],
                    "pos": pos,
                    "end": end,
                    "type": svtype,
                    "len": max(1, svlen),
                    "seq": seq,
                }
            )
    return records


def parse_paf_tags(fields):
    tags = {}
    for tag in fields[12:]:
        parts = tag.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def add_block(bins, chrom, start, end, block, bin_size):
    first = max(0, start // bin_size)
    last = max(0, end // bin_size)
    for bin_id in range(first, last + 1):
        bins[chrom][bin_id].append(block)


def read_paf(path, args):
    query_bins = defaultdict(lambda: defaultdict(list))
    target_bins = defaultdict(lambda: defaultdict(list))
    kept = 0
    total = 0
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue
            total += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            tags = parse_paf_tags(fields)
            if args.require_primary and tags.get("tp") not in {None, "P"}:
                continue
            qchrom = fields[0]
            qstart, qend = int(fields[2]), int(fields[3])
            strand = fields[4]
            tchrom = fields[5]
            tstart, tend = int(fields[7]), int(fields[8])
            matches = int(fields[9])
            block_len = max(1, int(fields[10]))
            mapq = int(fields[11])
            dv = float(tags.get("dv", "1.0"))
            if mapq < args.min_mapq or block_len < args.min_synteny_block or dv > args.max_dv:
                continue
            block = {
                "id": f"{qchrom}:{qstart}-{qend}:{strand}:{tchrom}:{tstart}-{tend}",
                "qchrom": qchrom,
                "qstart": qstart,
                "qend": qend,
                "strand": strand,
                "tchrom": tchrom,
                "tstart": tstart,
                "tend": tend,
                "matches": matches,
                "block_len": block_len,
                "mapq": mapq,
                "dv": dv,
            }
            kept += 1
            add_block(query_bins, qchrom, qstart, qend, block, args.paf_bin)
            add_block(target_bins, tchrom, tstart, tend, block, args.paf_bin)
    return query_bins, target_bins, {"total": total, "kept": kept}


def block_quality(block):
    return (block["mapq"], -block["dv"], block["block_len"])


def blocks_for_interval(bins, chrom, start0, end0, bin_size):
    out = []
    seen = set()
    for bin_id in range(max(0, start0 // bin_size), max(0, end0 // bin_size) + 1):
        for block in bins.get(chrom, {}).get(bin_id, []):
            key = block["id"]
            if key not in seen:
                out.append(block)
                seen.add(key)
    return out


def choose_unique_block(candidates, args):
    if not candidates:
        return None, "no_high_quality_synteny", 0
    candidates.sort(key=block_quality, reverse=True)
    if len(candidates) > args.max_synteny_blocks:
        return None, "ambiguous_multi_synteny_blocks", len(candidates)
    return candidates[0], "ok", len(candidates)


def lift_query_to_target(query_bins, chrom, pos, end, args):
    qpos = max(0, pos - 1)
    qend = max(0, end - 1)
    containing = []
    for block in blocks_for_interval(query_bins, chrom, qpos, qend, args.paf_bin):
        if block["qchrom"] != chrom:
            continue
        if block["qstart"] <= qpos <= block["qend"] and block["qstart"] <= qend <= block["qend"]:
            containing.append(block)
    block, status, nblocks = choose_unique_block(containing, args)
    if not block:
        return None, status, nblocks
    if block["strand"] == "+":
        a = block["tstart"] + (qpos - block["qstart"]) + 1
        b = block["tstart"] + (qend - block["qstart"]) + 1
    else:
        a = block["tend"] - (qpos - block["qstart"])
        b = block["tend"] - (qend - block["qstart"])
    return (block["tchrom"], min(a, b), max(a, b), block["strand"], block), "ok", nblocks


def lift_target_to_query(target_bins, chrom, pos, end, args):
    tpos = max(0, pos - 1)
    tend0 = max(0, end - 1)
    containing = []
    for block in blocks_for_interval(target_bins, chrom, tpos, tend0, args.paf_bin):
        if block["tchrom"] != chrom:
            continue
        if block["tstart"] <= tpos <= block["tend"] and block["tstart"] <= tend0 <= block["tend"]:
            containing.append(block)
    block, status, nblocks = choose_unique_block(containing, args)
    if not block:
        return None, status, nblocks
    if block["strand"] == "+":
        a = block["qstart"] + (tpos - block["tstart"]) + 1
        b = block["qstart"] + (tend0 - block["tstart"]) + 1
    else:
        a = block["qstart"] + (block["tend"] - (tpos + 1)) + 1
        b = block["qstart"] + (block["tend"] - (tend0 + 1)) + 1
    return (block["qchrom"], min(a, b), max(a, b), block["strand"], block), "ok", nblocks


def reciprocal_overlap(a1, a2, b1, b2):
    overlap = max(0, min(a2, b2) - max(a1, b1) + 1)
    la = max(1, a2 - a1 + 1)
    lb = max(1, b2 - b1 + 1)
    return min(float(overlap) / la, float(overlap) / lb)


def build_record_index(records, bin_size):
    index = defaultdict(list)
    by_id = {}
    for record in records:
        by_id[record["id"]] = record
        key = (record["chrom"], record["type"])
        first = record["pos"] // bin_size - 1
        last = record["pos"] // bin_size + 1
        for bin_id in range(first, last + 1):
            index[(key, bin_id)].append(record)
    return index, by_id


def load_pav(path):
    data = {}
    with open(path, encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        samples = header[6:]
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            values = []
            for value in fields[6 : 6 + len(samples)]:
                if value == "1":
                    values.append("1")
                elif value == "0":
                    values.append("0")
                else:
                    values.append("N")
            while len(values) < len(samples):
                values.append("N")
            data[fields[0]] = "".join(values)
    return samples, data


def pav_jaccard(refa_id, refb_id, refa_pav, refb_pav, min_presence):
    a = refa_pav.get(refa_id)
    b = refb_pav.get(refb_id)
    if a is None or b is None:
        return None
    aset = {i for i, value in enumerate(a) if value == "1"}
    bset = {i for i, value in enumerate(b) if value == "1"}
    if len(aset) < min_presence or len(bset) < min_presence:
        return None
    union = len(aset | bset)
    if union == 0:
        return None
    return float(len(aset & bset)) / union


def kmer_set(seq, k):
    if not seq:
        return set()
    if len(seq) < k:
        return {seq}
    return {seq[i : i + k] for i in range(0, len(seq) - k + 1)}


def reverse_complement(seq):
    return seq.translate(RC_TABLE)[::-1].upper()


def sequence_similarity(seq_a, seq_b, k):
    if not seq_a or not seq_b:
        return None
    k = max(5, min(k, max(5, min(len(seq_a), len(seq_b)) // 3)))
    aset = kmer_set(seq_a.upper(), k)
    bset = kmer_set(seq_b.upper(), k)
    rcset = kmer_set(reverse_complement(seq_b), k)

    def score(query, target):
        inter = len(query & target)
        union = len(query | target)
        min_size = max(1, min(len(query), len(target)))
        return (float(inter) / union if union else 0.0, float(inter) / min_size)

    j1, c1 = score(aset, bset)
    j2, c2 = score(aset, rcset)
    if (j2, c2) > (j1, c1):
        return {"jaccard": j2, "containment": c2, "orientation": "reverse_complement", "k": k}
    return {"jaccard": j1, "containment": c1, "orientation": "same", "k": k}


def combine_pav(refa_id, refb_id, refa_pav, refb_pav, nsamples):
    a = refa_pav.get(refa_id) if refa_id and refa_id != "." else None
    b = refb_pav.get(refb_id) if refb_id and refb_id != "." else None
    values = []
    for i in range(nsamples):
        av = a[i] if a and i < len(a) else "N"
        bv = b[i] if b and i < len(b) else "N"
        if av == "1" or bv == "1":
            values.append("1")
        elif av == "0" or bv == "0":
            values.append("0")
        else:
            values.append("NA")
    return values


def evidence_block(block):
    return f"synteny_block={block['id']};mapq={block['mapq']};dv={block['dv']:.4f};block_len={block['block_len']}"


def score_interval(record, other, tpos, tend, ro, pav_score):
    dist = abs(other["pos"] - tpos) + abs(other["end"] - tend)
    ratio = float(max(record["len"], other["len"])) / max(1, min(record["len"], other["len"]))
    pav_bonus = 0.0 if pav_score is None else pav_score
    return (ro * 1000.0) + (pav_bonus * 50.0) - (dist / 1000.0) - (abs(ratio - 1.0) * 20.0)


def score_ins(record, other, tpos, sim, pav_score):
    dist = abs(other["pos"] - tpos)
    ratio = float(max(record["len"], other["len"])) / max(1, min(record["len"], other["len"]))
    pav_bonus = 0.0 if pav_score is None else pav_score
    return (sim["jaccard"] * 1000.0) + (sim["containment"] * 250.0) + (pav_bonus * 50.0) - (dist / 100.0) - (abs(ratio - 1.0) * 30.0)


def find_candidates(records, other_index, lift_func, args, refa_pav, refb_pav, direction):
    best = {}
    all_candidates = {}
    lift_status = {}
    for rec in records:
        if rec["type"] not in MATCHABLE_TYPES:
            lift_status[rec["id"]] = "unsupported_svtype"
            all_candidates[rec["id"]] = []
            continue
        lifted, status, nblocks = lift_func(rec["chrom"], rec["pos"], rec["end"])
        lift_status[rec["id"]] = status
        if not lifted:
            all_candidates[rec["id"]] = []
            continue
        tchrom, tpos, tend, strand, block = lifted
        candidates = []
        for bin_id in range(tpos // args.record_bin - 1, tpos // args.record_bin + 2):
            for other in other_index.get(((tchrom, rec["type"]), bin_id), []):
                if rec["type"] == "INS":
                    if not rec["seq"] or not other["seq"]:
                        continue
                    if abs(other["pos"] - tpos) > args.ins_distance:
                        continue
                    ratio = float(max(rec["len"], other["len"])) / max(1, min(rec["len"], other["len"]))
                    if ratio > args.ins_max_len_ratio:
                        continue
                    sim = sequence_similarity(rec["seq"], other["seq"], args.ins_kmer)
                    if not sim:
                        continue
                    if sim["jaccard"] < args.ins_min_jaccard or sim["containment"] < args.ins_min_containment:
                        continue
                    pav_score = pav_jaccard(rec["id"], other["id"], refa_pav, refb_pav, args.pav_min_presence) if direction == "AtoB" else pav_jaccard(other["id"], rec["id"], refa_pav, refb_pav, args.pav_min_presence)
                    if pav_score is not None and pav_score < args.min_pav_jaccard:
                        continue
                    candidates.append(
                        {
                            "other_id": other["id"],
                            "score": score_ins(rec, other, tpos, sim, pav_score),
                            "lift_chrom": tchrom,
                            "lift_pos": tpos,
                            "lift_end": tend,
                            "strand": strand,
                            "block": block,
                            "nblocks": nblocks,
                            "len_ratio": ratio,
                            "pav_jaccard": pav_score,
                            "seq_jaccard": sim["jaccard"],
                            "seq_containment": sim["containment"],
                            "seq_orientation": sim["orientation"],
                            "seq_k": sim["k"],
                            "distance": abs(other["pos"] - tpos),
                            "match_class": "INS",
                        }
                    )
                else:
                    if abs(other["pos"] - tpos) > args.distance or abs(other["end"] - tend) > args.distance:
                        continue
                    ratio = float(max(rec["len"], other["len"])) / max(1, min(rec["len"], other["len"]))
                    if ratio > args.max_len_ratio:
                        continue
                    ro = reciprocal_overlap(tpos, tend, other["pos"], other["end"])
                    if ro < args.min_ro:
                        continue
                    pav_score = pav_jaccard(rec["id"], other["id"], refa_pav, refb_pav, args.pav_min_presence) if direction == "AtoB" else pav_jaccard(other["id"], rec["id"], refa_pav, refb_pav, args.pav_min_presence)
                    if pav_score is not None and pav_score < args.min_pav_jaccard:
                        continue
                    candidates.append(
                        {
                            "other_id": other["id"],
                            "score": score_interval(rec, other, tpos, tend, ro, pav_score),
                            "lift_chrom": tchrom,
                            "lift_pos": tpos,
                            "lift_end": tend,
                            "strand": strand,
                            "block": block,
                            "nblocks": nblocks,
                            "len_ratio": ratio,
                            "pav_jaccard": pav_score,
                            "ro": ro,
                            "distance": abs(other["pos"] - tpos) + abs(other["end"] - tend),
                            "match_class": "INTERVAL",
                        }
                    )
        candidates.sort(key=lambda item: item["score"], reverse=True)
        all_candidates[rec["id"]] = candidates
        if candidates:
            best[rec["id"]] = candidates[0]
    return best, all_candidates, lift_status


def event_row(meta_id, refa, refb, svtype, svlen, status, method, confidence):
    if refa:
        a = [refa["id"], refa["chrom"], str(refa["pos"]), str(refa["end"])]
    else:
        a = [".", ".", ".", "."]
    if refb:
        b = [refb["id"], refb["chrom"], str(refb["pos"]), str(refb["end"])]
    else:
        b = [".", ".", ".", "."]
    return [meta_id] + a + b + [svtype, str(svlen), status, method, confidence]


def write_stats(path, samples, rows):
    counters = {sample: Counter() for sample in samples}
    type_counts = {sample: Counter() for sample in samples}
    for row in rows:
        present = []
        for sample, value in zip(samples, row["pav"]):
            if value == "1":
                counters[sample]["present"] += 1
                type_counts[sample][row["svtype"]] += 1
                present.append(sample)
            elif value == "0":
                counters[sample]["absent"] += 1
            else:
                counters[sample]["missing"] += 1
        if len(present) == 1:
            counters[present[0]]["singleton_present"] += 1
    all_types = sorted({svtype for values in type_counts.values() for svtype in values})
    with open(path, "w", encoding="utf-8") as out:
        out.write("sample\tpresent\tabsent\tmissing\tmissing_rate\tsingleton_present")
        for svtype in all_types:
            out.write("\tSVTYPE_" + svtype)
        out.write("\n")
        denom = float(len(rows)) if rows else 1.0
        for sample in samples:
            counts = counters[sample]
            out.write(
                f"{sample}\t{counts['present']}\t{counts['absent']}\t{counts['missing']}\t"
                f"{counts['missing'] / denom:.6f}\t{counts['singleton_present']}"
            )
            for svtype in all_types:
                out.write("\t" + str(type_counts[sample][svtype]))
            out.write("\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--refa-vcf", required=True)
    parser.add_argument("--refb-vcf", required=True)
    parser.add_argument("--refa-pav", required=True)
    parser.add_argument("--refb-pav", required=True)
    parser.add_argument("--paf", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--distance", type=int, default=500)
    parser.add_argument("--min-ro", type=float, default=0.7)
    parser.add_argument("--max-len-ratio", type=float, default=1.5)
    parser.add_argument("--ins-distance", type=int, default=500)
    parser.add_argument("--ins-max-len-ratio", type=float, default=1.5)
    parser.add_argument("--ins-min-jaccard", type=float, default=0.35)
    parser.add_argument("--ins-min-containment", type=float, default=0.70)
    parser.add_argument("--ins-kmer", type=int, default=15)
    parser.add_argument("--min-pav-jaccard", type=float, default=0.50)
    parser.add_argument("--pav-min-presence", type=int, default=2)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--max-dv", type=float, default=0.12)
    parser.add_argument("--min-synteny-block", type=int, default=10000)
    parser.add_argument("--require-primary", action="store_true", default=True)
    parser.add_argument("--max-synteny-blocks", type=int, default=1)
    parser.add_argument("--paf-bin", type=int, default=1000000)
    parser.add_argument("--record-bin", type=int, default=100000)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    refa_records = read_vcf(args.refa_vcf)
    refb_records = read_vcf(args.refb_vcf)
    refa_samples, refa_pav = load_pav(args.refa_pav)
    refb_samples, refb_pav = load_pav(args.refb_pav)
    if refa_samples != refb_samples:
        raise SystemExit("RefA and RefB PAV matrices have different sample order")
    query_bins, target_bins, paf_stats = read_paf(args.paf, args)
    refa_index, refa_by_id = build_record_index(refa_records, args.record_bin)
    refb_index, refb_by_id = build_record_index(refb_records, args.record_bin)

    best_b_for_a, candidates_b_for_a, lift_a = find_candidates(
        refa_records,
        refb_index,
        lambda chrom, pos, end: lift_query_to_target(query_bins, chrom, pos, end, args),
        args,
        refa_pav,
        refb_pav,
        "AtoB",
    )
    best_a_for_b, candidates_a_for_b, lift_b = find_candidates(
        refb_records,
        refa_index,
        lambda chrom, pos, end: lift_target_to_query(target_bins, chrom, pos, end, args),
        args,
        refa_pav,
        refb_pav,
        "BtoA",
    )

    best_b_counts = Counter(item["other_id"] for item in best_b_for_a.values())
    best_a_counts = Counter(item["other_id"] for item in best_a_for_b.values())

    publication_pairs = {}
    matched_a = set()
    matched_b = set()
    for a_id, cand in best_b_for_a.items():
        b_id = cand["other_id"]
        reciprocal = best_a_for_b.get(b_id)
        if reciprocal and reciprocal["other_id"] == a_id and best_b_counts[b_id] == 1 and best_a_counts[a_id] == 1:
            publication_pairs[a_id] = b_id
            matched_a.add(a_id)
            matched_b.add(b_id)

    event_header = [
        "MetaSV_ID",
        "RefA_SV_ID",
        "RefA_CHROM",
        "RefA_POS",
        "RefA_END",
        "RefB_SV_ID",
        "RefB_CHROM",
        "RefB_POS",
        "RefB_END",
        "SVTYPE",
        "SVLEN",
        "match_status",
        "match_method",
        "confidence",
    ]
    rows = []
    evidence_rows = []
    ambiguous_rows = []
    unresolved_ins_rows = []
    meta_i = 0

    def next_id():
        nonlocal meta_i
        meta_i += 1
        return "PubMetaSV_%08d" % meta_i

    def add_publication(refa, refb, status, method, confidence, evidence):
        meta_id = next_id()
        svtype = refa["type"] if refa else refb["type"]
        svlen = refa["len"] if refa else refb["len"]
        row = event_row(meta_id, refa, refb, svtype, svlen, status, method, confidence)
        pav = combine_pav(refa["id"] if refa else ".", refb["id"] if refb else ".", refa_pav, refb_pav, len(refa_samples))
        rows.append({"row": row, "pav": pav, "svtype": svtype})
        evidence_rows.append([meta_id, method, evidence])

    def add_ambiguous(refa, refb, status, method, confidence, reason):
        meta_id = next_id()
        svtype = refa["type"] if refa else refb["type"]
        svlen = refa["len"] if refa else refb["len"]
        ambiguous_rows.append(event_row(meta_id, refa, refb, svtype, svlen, status, method, confidence) + [reason])

    def add_unresolved_ins(refa, refb, reason):
        meta_id = next_id()
        svtype = refa["type"] if refa else refb["type"]
        svlen = refa["len"] if refa else refb["len"]
        unresolved_ins_rows.append(event_row(meta_id, refa, refb, svtype, svlen, "ins_unresolved", "flank_or_sequence_evidence_insufficient", "low") + [reason])

    for a in refa_records:
        if a["type"] not in MATCHABLE_TYPES:
            add_ambiguous(a, None, "unresolved_unsupported_svtype", "unsupported_svtype", "low", "not_in_DEL_DUP_INV_INS")
            continue
        if a["id"] in publication_pairs:
            b = refb_by_id[publication_pairs[a["id"]]]
            cand = best_b_for_a[a["id"]]
            if a["type"] == "INS":
                status = "shared_strict_sequence"
                method = "reciprocal_best_liftover_ins_kmer"
                confidence = "high" if cand["seq_jaccard"] >= 0.50 and cand["seq_containment"] >= 0.80 else "medium"
                evidence = (
                    f"lifted={cand['lift_chrom']}:{cand['lift_pos']};strand={cand['strand']};"
                    f"breakpoint_distance={cand['distance']};len_ratio={cand['len_ratio']:.4f};"
                    f"seq_jaccard={cand['seq_jaccard']:.4f};seq_containment={cand['seq_containment']:.4f};"
                    f"seq_orientation={cand['seq_orientation']};k={cand['seq_k']};"
                    f"pav_jaccard={'NA' if cand['pav_jaccard'] is None else '%.4f' % cand['pav_jaccard']};"
                    + evidence_block(cand["block"])
                )
            else:
                status = "shared_strict_liftover"
                method = "reciprocal_best_synteny_liftover_ro"
                confidence = "high" if cand["ro"] >= 0.80 and (cand["pav_jaccard"] is None or cand["pav_jaccard"] >= 0.70) else "medium"
                evidence = (
                    f"lifted={cand['lift_chrom']}:{cand['lift_pos']}-{cand['lift_end']};strand={cand['strand']};"
                    f"breakpoint_distance_sum={cand['distance']};reciprocal_overlap={cand['ro']:.4f};"
                    f"len_ratio={cand['len_ratio']:.4f};pav_jaccard={'NA' if cand['pav_jaccard'] is None else '%.4f' % cand['pav_jaccard']};"
                    + evidence_block(cand["block"])
                )
            add_publication(a, b, status, method, confidence, evidence)
            continue
        if a["type"] == "INS":
            reason = lift_a.get(a["id"], "unknown")
            if reason == "ok" and candidates_b_for_a.get(a["id"]):
                reason = "INS_candidate_not_reciprocal_or_not_one_to_one"
            elif reason == "ok":
                reason = "INS_no_sequence_match_in_syntenic_target_window"
            add_unresolved_ins(a, None, reason)
            continue
        if a["type"] not in INTERVAL_TYPES:
            add_ambiguous(a, None, "unresolved_unsupported_svtype", "unsupported_svtype", "low", "not_in_DEL_DUP_INV_INS")
            continue
        if lift_a.get(a["id"]) != "ok":
            add_ambiguous(a, None, "ambiguous_synteny", "synteny_filter_failed", "low", lift_a.get(a["id"], "unknown"))
        elif not candidates_b_for_a.get(a["id"]):
            add_publication(a, None, "RefA_only_publication_lifted_no_match", "unique_synteny_no_counterpart", "medium", "unique_high_quality_synteny_liftover_but_no_matching_RefB_SV")
        else:
            best = best_b_for_a[a["id"]]
            b = refb_by_id[best["other_id"]]
            if best_b_counts[best["other_id"]] > 1:
                status, reason = "ambiguous_many_to_one", "multiple_RefA_records_choose_same_RefB_best"
            elif best_a_counts[a["id"]] > 1:
                status, reason = "ambiguous_one_to_many", "multiple_RefB_records_choose_same_RefA_best"
            else:
                status, reason = "ambiguous_not_reciprocal", "best_hit_not_reciprocal"
            add_ambiguous(a, b, status, "publication_candidate_not_strict", "low", reason)

    for b in refb_records:
        if b["id"] in matched_b:
            continue
        if b["type"] == "INS":
            reason = lift_b.get(b["id"], "unknown")
            if reason == "ok" and candidates_a_for_b.get(b["id"]):
                reason = "INS_candidate_not_reciprocal_or_not_one_to_one"
            elif reason == "ok":
                reason = "INS_no_sequence_match_in_syntenic_target_window"
            add_unresolved_ins(None, b, reason)
            continue
        if b["type"] not in INTERVAL_TYPES:
            add_ambiguous(None, b, "unresolved_unsupported_svtype", "unsupported_svtype", "low", "not_in_DEL_DUP_INV_INS")
            continue
        if lift_b.get(b["id"]) != "ok":
            add_ambiguous(None, b, "ambiguous_synteny", "synteny_filter_failed", "low", lift_b.get(b["id"], "unknown"))
        elif not candidates_a_for_b.get(b["id"]):
            add_publication(None, b, "RefB_only_publication_lifted_no_match", "unique_synteny_no_counterpart", "medium", "unique_high_quality_synteny_liftover_but_no_matching_RefA_SV")
        else:
            best = best_a_for_b[b["id"]]
            a = refa_by_id[best["other_id"]]
            if best_a_counts[best["other_id"]] > 1:
                status, reason = "ambiguous_one_to_many", "multiple_RefB_records_choose_same_RefA_best"
            elif best_b_counts[b["id"]] > 1:
                status, reason = "ambiguous_many_to_one", "multiple_RefA_records_choose_same_RefB_best"
            else:
                status, reason = "ambiguous_not_reciprocal", "best_hit_not_reciprocal"
            add_ambiguous(a, b, status, "publication_candidate_not_strict", "low", reason)

    paths = {
        "events": os.path.join(args.out_dir, "meta-panSV.publication.events.tsv"),
        "pav": os.path.join(args.out_dir, "meta-panSV.publication.PAV.matrix.tsv"),
        "stats": os.path.join(args.out_dir, "meta-panSV.publication.stats.tsv"),
        "evidence": os.path.join(args.out_dir, "meta-panSV.publication.matching.evidence.tsv"),
        "ambiguous": os.path.join(args.out_dir, "meta-panSV.publication.ambiguous.events.tsv"),
        "unresolved_ins": os.path.join(args.out_dir, "meta-panSV.publication.unresolved_INS.events.tsv"),
        "summary": os.path.join(args.out_dir, "meta-panSV.publication.summary.tsv"),
    }
    with open(paths["events"], "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(event_header)
        for item in rows:
            writer.writerow(item["row"])
    with open(paths["pav"], "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["MetaSV_ID", "RefA_SV_ID", "RefB_SV_ID", "SVTYPE", "SVLEN", "match_status", "match_method", "confidence"] + refa_samples)
        for item in rows:
            row = item["row"]
            writer.writerow([row[0], row[1], row[5], row[9], row[10], row[11], row[12], row[13]] + item["pav"])
    with open(paths["evidence"], "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["MetaSV_ID", "method", "details"])
        writer.writerows(evidence_rows)
    with open(paths["ambiguous"], "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(event_header + ["reason"])
        writer.writerows(ambiguous_rows)
    with open(paths["unresolved_ins"], "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(event_header + ["reason"])
        writer.writerows(unresolved_ins_rows)
    write_stats(paths["stats"], refa_samples, rows)

    with open(paths["summary"], "w", encoding="utf-8") as out:
        out.write("section\tkey\tcount\n")
        out.write(f"input\tRefA_records\t{len(refa_records)}\n")
        out.write(f"input\tRefB_records\t{len(refb_records)}\n")
        out.write(f"paf\ttotal_blocks\t{paf_stats['total']}\n")
        out.write(f"paf\tkept_high_quality_blocks\t{paf_stats['kept']}\n")
        out.write(f"publication\ttotal\t{len(rows)}\n")
        for key, value in sorted(Counter(item["row"][11] for item in rows).items()):
            out.write(f"publication_status\t{key}\t{value}\n")
        for key, value in sorted(Counter(item["svtype"] for item in rows).items()):
            out.write(f"publication_svtype\t{key}\t{value}\n")
        out.write(f"ambiguous\ttotal\t{len(ambiguous_rows)}\n")
        for key, value in sorted(Counter(row[11] for row in ambiguous_rows).items()):
            out.write(f"ambiguous_status\t{key}\t{value}\n")
        out.write(f"unresolved_INS\ttotal\t{len(unresolved_ins_rows)}\n")
        for key, value in sorted(Counter(row[11] for row in unresolved_ins_rows).items()):
            out.write(f"unresolved_status\t{key}\t{value}\n")
        for name in [
            "distance",
            "min_ro",
            "max_len_ratio",
            "ins_distance",
            "ins_max_len_ratio",
            "ins_min_jaccard",
            "ins_min_containment",
            "min_pav_jaccard",
            "min_mapq",
            "max_dv",
            "min_synteny_block",
            "max_synteny_blocks",
        ]:
            out.write(f"parameters\t{name}\t{getattr(args, name)}\n")


if __name__ == "__main__":
    main()
