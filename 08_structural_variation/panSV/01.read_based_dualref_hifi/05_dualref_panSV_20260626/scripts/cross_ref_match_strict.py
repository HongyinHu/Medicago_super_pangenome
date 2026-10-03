#!/usr/bin/env python3
import argparse
import csv
import gzip
import os
from bisect import bisect_right
from collections import Counter, defaultdict


STRICT_TYPES = {"DEL", "DUP", "INV"}


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


def read_vcf(path, prefix):
    records = []
    seen = Counter()
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
            raw_id = fields[2] if fields[2] and fields[2] != "." else f"{prefix}_{fields[0]}_{pos}_{end}_{svtype}"
            seen[raw_id] += 1
            sv_id = raw_id if seen[raw_id] == 1 else f"{raw_id}__dup{seen[raw_id]}"
            records.append(
                {
                    "id": sv_id,
                    "raw_id": raw_id,
                    "chrom": fields[0],
                    "pos": pos,
                    "end": end,
                    "type": svtype,
                    "len": max(1, svlen),
                }
            )
    return records


def add_block(bins, starts, chrom, start, end, block, bin_size):
    bins_for_chrom = bins[chrom]
    first = max(0, start // bin_size)
    last = max(0, end // bin_size)
    for bin_id in range(first, last + 1):
        bins_for_chrom[bin_id].append(block)
    starts[chrom].append((start, block))


def read_paf(path, bin_size):
    query_bins = defaultdict(lambda: defaultdict(list))
    target_bins = defaultdict(lambda: defaultdict(list))
    query_starts = defaultdict(list)
    target_starts = defaultdict(list)
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            qchrom = fields[0]
            qstart, qend = int(fields[2]), int(fields[3])
            strand = fields[4]
            tchrom = fields[5]
            tstart, tend = int(fields[7]), int(fields[8])
            block = (qchrom, qstart, qend, strand, tchrom, tstart, tend)
            add_block(query_bins, query_starts, qchrom, qstart, qend, block, bin_size)
            add_block(target_bins, target_starts, tchrom, tstart, tend, block, bin_size)
    for chrom in query_starts:
        query_starts[chrom].sort(key=lambda x: x[0])
    for chrom in target_starts:
        target_starts[chrom].sort(key=lambda x: x[0])
    return query_bins, target_bins


def blocks_for_interval(bins, chrom, pos0, end0, bin_size):
    candidates = []
    seen = set()
    for bin_id in {max(0, pos0 // bin_size), max(0, end0 // bin_size)}:
        for block in bins.get(chrom, {}).get(bin_id, []):
            key = id(block)
            if key not in seen:
                candidates.append(block)
                seen.add(key)
    return candidates


def lift_query_to_target(query_bins, chrom, pos, end, bin_size):
    qpos = max(0, pos - 1)
    qend = max(0, end - 1)
    for qchrom, qs, qe, strand, tchrom, ts, te in blocks_for_interval(query_bins, chrom, qpos, qend, bin_size):
        if qchrom != chrom:
            continue
        if qs <= qpos <= qe and qs <= qend <= qe:
            if strand == "+":
                a = ts + (qpos - qs) + 1
                b = ts + (qend - qs) + 1
            else:
                a = te - (qpos - qs)
                b = te - (qend - qs)
            return tchrom, min(a, b), max(a, b), strand
    return None


def lift_target_to_query(target_bins, chrom, pos, end, bin_size):
    tpos = max(0, pos - 1)
    tend0 = max(0, end - 1)
    for qchrom, qs, qe, strand, tchrom, ts, te in blocks_for_interval(target_bins, chrom, tpos, tend0, bin_size):
        if tchrom != chrom:
            continue
        if ts <= tpos <= te and ts <= tend0 <= te:
            if strand == "+":
                a = qs + (tpos - ts) + 1
                b = qs + (tend0 - ts) + 1
            else:
                a = qs + (te - (tpos + 1)) + 1
                b = qs + (te - (tend0 + 1)) + 1
            return qchrom, min(a, b), max(a, b), strand
    return None


def reciprocal_overlap(a1, a2, b1, b2):
    overlap = max(0, min(a2, b2) - max(a1, b1) + 1)
    la = max(1, a2 - a1 + 1)
    lb = max(1, b2 - b1 + 1)
    return min(float(overlap) / la, float(overlap) / lb)


def build_record_index(records, bin_size):
    index = defaultdict(list)
    by_id = {}
    for rec in records:
        by_id[rec["id"]] = rec
        key = (rec["chrom"], rec["type"])
        first = rec["pos"] // bin_size - 1
        last = rec["pos"] // bin_size + 1
        for bin_id in range(first, last + 1):
            index[(key, bin_id)].append(rec)
    return index, by_id


def load_pav(path):
    data = {}
    with open(path, encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if len(header) <= 6 or header[0] != "SV_ID":
            raise SystemExit(f"Unexpected PAV header in {path}")
        samples = header[6:]
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            values = fields[6:]
            encoded = []
            for value in values[: len(samples)]:
                if value == "1":
                    encoded.append("1")
                elif value == "0":
                    encoded.append("0")
                else:
                    encoded.append("N")
            if len(encoded) < len(samples):
                encoded.extend("N" for _ in range(len(samples) - len(encoded)))
            data[fields[0]] = "".join(encoded)
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


def combine_pav(refa_id, refb_id, refa_pav, refb_pav, sample_count):
    a = refa_pav.get(refa_id) if refa_id and refa_id != "." else None
    b = refb_pav.get(refb_id) if refb_id and refb_id != "." else None
    values = []
    for i in range(sample_count):
        av = a[i] if a and i < len(a) else "N"
        bv = b[i] if b and i < len(b) else "N"
        if av == "1" or bv == "1":
            values.append("1")
        elif av == "0" or bv == "0":
            values.append("0")
        else:
            values.append("NA")
    return values


def candidate_score(record, candidate, lifted_pos, lifted_end, ro, pav_score):
    dist = abs(candidate["pos"] - lifted_pos) + abs(candidate["end"] - lifted_end)
    ratio = float(max(record["len"], candidate["len"])) / max(1, min(record["len"], candidate["len"]))
    pav_bonus = 0.0 if pav_score is None else pav_score
    return (ro * 1000.0) + (pav_bonus * 50.0) - (dist / 1000.0) - (abs(ratio - 1.0) * 20.0)


def find_candidates(records, other_index, lift_func, args, refa_pav, refb_pav, direction):
    best = {}
    all_candidates = {}
    lifted_intervals = {}
    lifted_ok = set()
    for rec in records:
        if rec["type"] not in STRICT_TYPES:
            all_candidates[rec["id"]] = []
            continue
        lifted = lift_func(rec["chrom"], rec["pos"], rec["end"])
        if not lifted:
            all_candidates[rec["id"]] = []
            continue
        tchrom, tpos, tend, strand = lifted
        lifted_intervals[rec["id"]] = lifted
        lifted_ok.add(rec["id"])
        candidates = []
        for bin_id in range(tpos // args.record_bin - 1, tpos // args.record_bin + 2):
            for other in other_index.get(((tchrom, rec["type"]), bin_id), []):
                if abs(other["pos"] - tpos) > args.distance or abs(other["end"] - tend) > args.distance:
                    continue
                ratio = float(max(rec["len"], other["len"])) / max(1, min(rec["len"], other["len"]))
                if ratio > args.max_len_ratio:
                    continue
                ro = reciprocal_overlap(tpos, tend, other["pos"], other["end"])
                if ro < args.min_ro:
                    continue
                if direction == "AtoB":
                    pav_score = pav_jaccard(rec["id"], other["id"], refa_pav, refb_pav, args.pav_min_presence)
                else:
                    pav_score = pav_jaccard(other["id"], rec["id"], refa_pav, refb_pav, args.pav_min_presence)
                if pav_score is not None and pav_score < args.min_pav_jaccard:
                    continue
                score = candidate_score(rec, other, tpos, tend, ro, pav_score)
                candidates.append(
                    {
                        "other_id": other["id"],
                        "score": score,
                        "ro": ro,
                        "pav_jaccard": pav_score,
                        "strand": strand,
                        "lift_chrom": tchrom,
                        "lift_pos": tpos,
                        "lift_end": tend,
                        "len_ratio": ratio,
                    }
                )
        candidates.sort(key=lambda x: x["score"], reverse=True)
        all_candidates[rec["id"]] = candidates
        if candidates:
            best[rec["id"]] = candidates[0]
    return best, all_candidates, lifted_intervals, lifted_ok


def event_row(meta_id, refa, refb, svtype, svlen, status, method, confidence):
    if refa:
        a_id, a_chrom, a_pos, a_end = refa["id"], refa["chrom"], str(refa["pos"]), str(refa["end"])
    else:
        a_id, a_chrom, a_pos, a_end = ".", ".", ".", "."
    if refb:
        b_id, b_chrom, b_pos, b_end = refb["id"], refb["chrom"], str(refb["pos"]), str(refb["end"])
    else:
        b_id, b_chrom, b_pos, b_end = ".", ".", ".", "."
    return [
        meta_id,
        a_id,
        a_chrom,
        a_pos,
        a_end,
        b_id,
        b_chrom,
        b_pos,
        b_end,
        svtype,
        str(svlen),
        status,
        method,
        confidence,
    ]


def write_stats(path, samples, rows):
    counters = {sample: Counter() for sample in samples}
    type_counts = {sample: Counter() for sample in samples}
    for row in rows:
        svtype = row["svtype"]
        present = [sample for sample, value in zip(samples, row["pav"]) if value == "1"]
        for sample, value in zip(samples, row["pav"]):
            if value == "1":
                counters[sample]["present"] += 1
                type_counts[sample][svtype] += 1
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
                f"{sample}\t{counts['present']}\t{counts['absent']}\t"
                f"{counts['missing']}\t{counts['missing'] / denom:.6f}\t"
                f"{counts['singleton_present']}"
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
    parser.add_argument("--distance", type=int, default=1000)
    parser.add_argument("--min-ro", type=float, default=0.5)
    parser.add_argument("--max-len-ratio", type=float, default=2.0)
    parser.add_argument("--min-pav-jaccard", type=float, default=0.30)
    parser.add_argument("--pav-min-presence", type=int, default=2)
    parser.add_argument("--paf-bin", type=int, default=1000000)
    parser.add_argument("--record-bin", type=int, default=100000)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    refa_records = read_vcf(args.refa_vcf, "RefA")
    refb_records = read_vcf(args.refb_vcf, "RefB")
    refa_samples, refa_pav = load_pav(args.refa_pav)
    refb_samples, refb_pav = load_pav(args.refb_pav)
    if refa_samples != refb_samples:
        raise SystemExit("RefA and RefB PAV matrices have different sample order")

    query_bins, target_bins = read_paf(args.paf, args.paf_bin)
    refa_index, refa_by_id = build_record_index(refa_records, args.record_bin)
    refb_index, refb_by_id = build_record_index(refb_records, args.record_bin)

    best_b_for_a, candidates_b_for_a, lifted_a, lifted_a_ok = find_candidates(
        refa_records,
        refb_index,
        lambda chrom, pos, end: lift_query_to_target(query_bins, chrom, pos, end, args.paf_bin),
        args,
        refa_pav,
        refb_pav,
        "AtoB",
    )
    best_a_for_b, candidates_a_for_b, lifted_b, lifted_b_ok = find_candidates(
        refb_records,
        refa_index,
        lambda chrom, pos, end: lift_target_to_query(target_bins, chrom, pos, end, args.paf_bin),
        args,
        refa_pav,
        refb_pav,
        "BtoA",
    )

    best_b_counts = Counter(item["other_id"] for item in best_b_for_a.values())
    best_a_counts = Counter(item["other_id"] for item in best_a_for_b.values())

    strict_pairs = {}
    matched_a = set()
    matched_b = set()
    for a_id, cand in best_b_for_a.items():
        b_id = cand["other_id"]
        reciprocal = best_a_for_b.get(b_id)
        if (
            reciprocal
            and reciprocal["other_id"] == a_id
            and best_b_counts[b_id] == 1
            and best_a_counts[a_id] == 1
        ):
            strict_pairs[a_id] = b_id
            matched_a.add(a_id)
            matched_b.add(b_id)

    events_path = os.path.join(args.out_dir, "meta-panSV.strict.events.tsv")
    pav_path = os.path.join(args.out_dir, "meta-panSV.strict.PAV.matrix.tsv")
    stats_path = os.path.join(args.out_dir, "meta-panSV.strict.stats.tsv")
    evidence_path = os.path.join(args.out_dir, "meta-panSV.strict.matching.evidence.tsv")
    ambiguous_path = os.path.join(args.out_dir, "meta-panSV.ambiguous.events.tsv")
    unresolved_ins_path = os.path.join(args.out_dir, "meta-panSV.unresolved_INS.events.tsv")
    summary_path = os.path.join(args.out_dir, "meta-panSV.strict.summary.tsv")

    strict_rows = []
    ambiguous_rows = []
    unresolved_ins_rows = []
    evidence_rows = []
    meta_i = 0

    def next_id():
        nonlocal meta_i
        meta_i += 1
        return "StrictMetaSV_%08d" % meta_i

    def add_strict(refa, refb, status, method, confidence, evidence):
        meta_id = next_id()
        svtype = refa["type"] if refa else refb["type"]
        svlen = refa["len"] if refa else refb["len"]
        row = event_row(meta_id, refa, refb, svtype, svlen, status, method, confidence)
        pav = combine_pav(refa["id"] if refa else ".", refb["id"] if refb else ".", refa_pav, refb_pav, len(refa_samples))
        strict_rows.append({"row": row, "pav": pav, "svtype": svtype})
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
        unresolved_ins_rows.append(event_row(meta_id, refa, refb, svtype, svlen, "ins_unresolved", "sequence_matching_pending", "low") + [reason])

    for a in refa_records:
        if a["type"] == "INS":
            add_unresolved_ins(a, None, "RefA_INS_requires_flank_and_insertion_sequence_matching")
            continue
        if a["type"] not in STRICT_TYPES:
            add_ambiguous(a, None, "unresolved_unsupported_svtype", "unsupported_svtype", "low", "not_in_DEL_DUP_INV_INS")
            continue
        if a["id"] in strict_pairs:
            b = refb_by_id[strict_pairs[a["id"]]]
            cand = best_b_for_a[a["id"]]
            pav_score = cand["pav_jaccard"]
            confidence = "high" if cand["ro"] >= 0.70 and (pav_score is None or pav_score >= 0.50) else "medium"
            add_strict(
                a,
                b,
                "shared_strict_liftover",
                "reciprocal_best_liftover_ro",
                confidence,
                "lifted=%s:%d-%d;strand=%s;reciprocal_overlap=%.4f;len_ratio=%.4f;pav_jaccard=%s;score=%.4f"
                % (
                    cand["lift_chrom"],
                    cand["lift_pos"],
                    cand["lift_end"],
                    cand["strand"],
                    cand["ro"],
                    cand["len_ratio"],
                    "NA" if pav_score is None else "%.4f" % pav_score,
                    cand["score"],
                ),
            )
            continue
        if a["id"] not in lifted_a_ok:
            add_ambiguous(a, None, "unresolved_no_liftover", "no_full_interval_liftover", "low", "RefA_interval_not_fully_lifted_to_RefB")
        elif not candidates_b_for_a.get(a["id"]):
            add_strict(a, None, "RefA_only_strict_lifted_no_match", "lifted_no_counterpart", "medium", "RefA_interval_lifted_to_RefB_but_no_matching_SV")
        else:
            best = best_b_for_a[a["id"]]
            b = refb_by_id[best["other_id"]]
            if best_b_counts[best["other_id"]] > 1:
                reason = "multiple_RefA_records_choose_same_RefB_best"
                status = "ambiguous_many_to_one"
            elif best_a_counts[a["id"]] > 1:
                reason = "multiple_RefB_records_choose_same_RefA_best"
                status = "ambiguous_one_to_many"
            else:
                reason = "best_hit_not_reciprocal"
                status = "ambiguous_not_reciprocal"
            add_ambiguous(a, b, status, "candidate_liftover_not_strict", "low", reason)

    for b in refb_records:
        if b["id"] in matched_b:
            continue
        if b["type"] == "INS":
            add_unresolved_ins(None, b, "RefB_INS_requires_flank_and_insertion_sequence_matching")
            continue
        if b["type"] not in STRICT_TYPES:
            add_ambiguous(None, b, "unresolved_unsupported_svtype", "unsupported_svtype", "low", "not_in_DEL_DUP_INV_INS")
            continue
        if b["id"] not in lifted_b_ok:
            add_ambiguous(None, b, "unresolved_no_liftover", "no_full_interval_liftover", "low", "RefB_interval_not_fully_lifted_to_RefA")
        elif not candidates_a_for_b.get(b["id"]):
            add_strict(None, b, "RefB_only_strict_lifted_no_match", "lifted_no_counterpart", "medium", "RefB_interval_lifted_to_RefA_but_no_matching_SV")
        else:
            best = best_a_for_b[b["id"]]
            a = refa_by_id[best["other_id"]]
            if best_a_counts[best["other_id"]] > 1:
                reason = "multiple_RefB_records_choose_same_RefA_best"
                status = "ambiguous_one_to_many"
            elif best_b_counts[b["id"]] > 1:
                reason = "multiple_RefA_records_choose_same_RefB_best"
                status = "ambiguous_many_to_one"
            else:
                reason = "best_hit_not_reciprocal"
                status = "ambiguous_not_reciprocal"
            add_ambiguous(a, b, status, "candidate_liftover_not_strict", "low", reason)

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
    with open(events_path, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(event_header)
        for item in strict_rows:
            writer.writerow(item["row"])
    with open(pav_path, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["MetaSV_ID", "RefA_SV_ID", "RefB_SV_ID", "SVTYPE", "SVLEN", "match_status", "match_method", "confidence"] + refa_samples)
        for item in strict_rows:
            row = item["row"]
            writer.writerow([row[0], row[1], row[5], row[9], row[10], row[11], row[12], row[13]] + item["pav"])
    with open(evidence_path, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["MetaSV_ID", "method", "details"])
        writer.writerows(evidence_rows)
    with open(ambiguous_path, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(event_header + ["reason"])
        writer.writerows(ambiguous_rows)
    with open(unresolved_ins_path, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(event_header + ["reason"])
        writer.writerows(unresolved_ins_rows)
    write_stats(stats_path, refa_samples, strict_rows)

    summary = Counter(item["row"][11] for item in strict_rows)
    ambiguous_summary = Counter(row[11] for row in ambiguous_rows)
    unresolved_summary = Counter(row[11] for row in unresolved_ins_rows)
    type_summary = Counter(item["svtype"] for item in strict_rows)
    with open(summary_path, "w", encoding="utf-8") as out:
        out.write("section\tkey\tcount\n")
        out.write(f"input\tRefA_records\t{len(refa_records)}\n")
        out.write(f"input\tRefB_records\t{len(refb_records)}\n")
        out.write(f"strict\ttotal\t{len(strict_rows)}\n")
        for key, value in sorted(summary.items()):
            out.write(f"strict_status\t{key}\t{value}\n")
        for key, value in sorted(type_summary.items()):
            out.write(f"strict_svtype\t{key}\t{value}\n")
        out.write(f"ambiguous\ttotal\t{len(ambiguous_rows)}\n")
        for key, value in sorted(ambiguous_summary.items()):
            out.write(f"ambiguous_status\t{key}\t{value}\n")
        out.write(f"unresolved_INS\ttotal\t{len(unresolved_ins_rows)}\n")
        for key, value in sorted(unresolved_summary.items()):
            out.write(f"unresolved_status\t{key}\t{value}\n")
        out.write(f"parameters\tdistance\t{args.distance}\n")
        out.write(f"parameters\tmin_ro\t{args.min_ro}\n")
        out.write(f"parameters\tmax_len_ratio\t{args.max_len_ratio}\n")
        out.write(f"parameters\tmin_pav_jaccard\t{args.min_pav_jaccard}\n")
        out.write(f"parameters\tpav_min_presence\t{args.pav_min_presence}\n")


if __name__ == "__main__":
    main()
