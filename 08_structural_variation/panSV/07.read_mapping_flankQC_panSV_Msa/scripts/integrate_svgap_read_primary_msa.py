#!/usr/bin/env python3
import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


BIN_SIZE = 100000
DIST = 500
INS_DIST = 1000
MIN_RO = 0.5
MAX_IDS_PER_READ = 50


FINAL_SAMPLES = [
    "genome_395", "genome_410", "genome_436", "genome_454", "genome_457", "genome_461",
    "genome_468", "genome_472", "genome_474", "genome_482", "genome_M22", "genome_M46",
    "genome_Mar", "genome_Mpo", "genome_Mru", "genome_Msa", "genome_R108", "genome_ZM4",
]
DROP_LABELS = set(["genome_A17", "A17", "genome_G474b", "G474b"])
MAP_LABELS = {"genome_G474a": "genome_474", "G474a": "genome_474"}


def split_labels(value):
    if not value or value == ".":
        return []
    out = []
    for token in str(value).replace(";", ",").split(","):
        token = token.strip()
        if token and token != ".":
            out.append(token)
    return out


def normalize_sample(label):
    if label == "genome_474a":
        label = "genome_474"
    if label == "genome_474b":
        return None
    label = MAP_LABELS.get(label, label)
    if label in DROP_LABELS:
        return None
    if label not in FINAL_SAMPLES:
        return None
    return label


def normalize_samples(value):
    out = []
    seen = set()
    for label in split_labels(value):
        sample = normalize_sample(label)
        if sample and sample not in seen:
            out.append(sample)
            seen.add(sample)
    return out


def sample_from_query_token(token):
    if not token:
        return None
    base = token.split(".", 1)[0]
    base = base.split("_Chr", 1)[0]
    aliases = {
        "395": "genome_395",
        "410": "genome_410",
        "436": "genome_436",
        "454": "genome_454",
        "457": "genome_457",
        "461": "genome_461",
        "468": "genome_468",
        "472": "genome_472",
        "482": "genome_482",
        "G395": "genome_395",
        "G410": "genome_410",
        "G436": "genome_436",
        "G454": "genome_454",
        "G457": "genome_457",
        "G461": "genome_461",
        "G468": "genome_468",
        "G472": "genome_472",
        "G474a": "genome_474",
        "G474b": None,
        "G482": "genome_482",
        "A17": None,
        "M22": "genome_M22",
        "M46": "genome_M46",
        "Mar": "genome_Mar",
        "Mpo": "genome_Mpo",
        "Mru": "genome_Mru",
        "R108": "genome_R108",
        "ZM4": "genome_ZM4",
        "Msa": "genome_Msa",
    }
    if base in aliases:
        return aliases[base]
    if base.startswith("genome_"):
        return normalize_sample(base)
    return normalize_sample("genome_" + base)


def samples_from_svgap_tokens(value):
    out = set()
    for token in split_labels(value):
        sample = sample_from_query_token(token)
        if sample:
            out.add(sample)
    return out


def safe_int(value, default=None):
    try:
        return int(float(str(value).split(",", 1)[0]))
    except Exception:
        return default


def norm_chrom(chrom):
    return chrom[4:] if chrom.startswith("Msa.") else chrom


def interval_len(start, end):
    if start is None or end is None:
        return 0
    return max(1, end - start + 1)


def reciprocal_overlap(a1, a2, b1, b2):
    ov = max(0, min(a2, b2) - max(a1, b1) + 1)
    if ov <= 0:
        return 0.0
    return min(float(ov) / interval_len(a1, a2), float(ov) / interval_len(b1, b2))


def length_ratio(a, b):
    a = abs(safe_int(a, 0) or 0)
    b = abs(safe_int(b, 0) or 0)
    if a == 0 or b == 0:
        return 1.0
    return min(float(a) / b, float(b) / a)


def bins_for(chrom, start, end):
    if start is None or end is None:
        return []
    a = int(start) // BIN_SIZE
    b = int(end) // BIN_SIZE
    return [(chrom, i) for i in range(a, b + 1)]


def add_index(index, event):
    for key in bins_for(event["CHROM"], event["POS"], event["END"]):
        index[key].append(event)


def read_pav(path):
    events = []
    index = defaultdict(list)
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            svtype = row["SVTYPE"]
            if svtype == "BND":
                svtype = "TRA"
            event = {
                "SV_ID": row["SV_ID"],
                "CHROM": norm_chrom(row["CHROM"]),
                "POS": safe_int(row["POS"], 0) or 0,
                "END": safe_int(row.get("END", row["POS"]), safe_int(row["POS"], 0)) or safe_int(row["POS"], 0) or 0,
                "SVTYPE": svtype,
                "SVLEN": safe_int(row.get("SVLEN", 0), 0) or 0,
                "CHR2": row.get("CHR2", "."),
                "MATEPOS": row.get("MATEPOS", "."),
                "pav": [row.get(s, "0") for s in FINAL_SAMPLES],
                "svgap_count": 0,
                "svgap_ids": [],
                "svgap_samples": set(),
                "svgap_source_ids": [],
            }
            events.append(event)
            add_index(index, event)
    return events, index


def samples_from_merge_samples(value):
    samples = set(normalize_samples(value))
    samples.update(samples_from_svgap_tokens(value))
    return samples


def iter_svgap_events(svg_dir):
    files = [
        (svg_dir / "All.DELs.50bplarge.bed.combined.sorted.txt", "combined"),
        (svg_dir / "All.INTs.50bplarge.bed.combined.sorted.txt", "combined"),
        (svg_dir / "All.CNVs.bed", "raw_cnv"),
        (svg_dir / "All.INVs.bed", "raw_inv"),
        (svg_dir / "All.TRLs.bed", "raw_trl"),
    ]
    seq = 0
    for source_file, kind in files:
        if not source_file.exists():
            continue
        with open(source_file, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t", 14)
                try:
                    if kind == "combined":
                        if len(parts) < 8:
                            continue
                        chrom = norm_chrom(parts[0])
                        start = safe_int(parts[1])
                        end = safe_int(parts[2])
                        svtype = parts[4]
                        if svtype == "INT":
                            svtype = "INS"
                        svlen = abs(safe_int(parts[5], 0) or interval_len(start, end))
                        source_id = "%s:%s-%s:%s" % (chrom, start, end, parts[3])
                        support_samples = samples_from_merge_samples(parts[7])
                    elif kind == "raw_cnv":
                        if len(parts) < 8:
                            continue
                        chrom = norm_chrom(parts[0])
                        start = safe_int(parts[1])
                        end = safe_int(parts[2])
                        query = parts[3]
                        svtype = "DUP"
                        svlen = abs(safe_int(parts[7], 0) or interval_len(start, end))
                        source_id = "%s_%s_%s_%s" % (query, parts[5], parts[6], svlen)
                        sample = sample_from_query_token(query)
                        support_samples = set([sample]) if sample else set()
                    elif kind == "raw_inv":
                        if len(parts) < 8:
                            continue
                        chrom = norm_chrom(parts[0])
                        start = safe_int(parts[1])
                        end = safe_int(parts[2])
                        query = parts[3]
                        svtype = "INV"
                        svlen = abs(safe_int(parts[7], 0) or interval_len(start, end))
                        source_id = "%s_%s_%s_%s" % (query, parts[5], parts[6], svlen)
                        sample = sample_from_query_token(query)
                        support_samples = set([sample]) if sample else set()
                    elif kind == "raw_trl":
                        if len(parts) < 9:
                            continue
                        chrom = norm_chrom(parts[0])
                        start = safe_int(parts[1])
                        end = safe_int(parts[2])
                        query = parts[3]
                        svtype = "TRA"
                        svlen = max(abs(safe_int(parts[7], 0) or 0), abs(safe_int(parts[8], 0) or 0), interval_len(start, end))
                        source_id = "%s_%s_%s_%s_%s_%s" % (query, parts[5], parts[6], chrom, start, end)
                        sample = sample_from_query_token(query)
                        support_samples = set([sample]) if sample else set()
                    else:
                        continue
                except Exception:
                    continue
                if start is None or end is None or not support_samples:
                    continue
                if start > end:
                    start, end = end, start
                seq += 1
                yield {
                    "svgap_id": "SVGAP_Msa_%s_%09d" % (svtype, seq),
                    "CHROM": chrom,
                    "POS": start,
                    "END": end,
                    "SVTYPE": svtype,
                    "SVLEN": svlen,
                    "source_file": str(source_file),
                    "source_id": source_id.replace("G474a", "genome_474"),
                    "support_samples": support_samples,
                }


def match_reason(read_ev, svg_ev):
    if read_ev["SVTYPE"] != svg_ev["SVTYPE"]:
        return None
    if read_ev["CHROM"] != svg_ev["CHROM"]:
        return None
    if read_ev["SVTYPE"] == "INS":
        if abs(read_ev["POS"] - svg_ev["POS"]) <= INS_DIST:
            return "ins_position_distance"
        return None
    if read_ev["SVTYPE"] == "TRA":
        if abs(read_ev["POS"] - svg_ev["POS"]) <= INS_DIST:
            return "tra_breakpoint_distance"
        return None
    ro = reciprocal_overlap(read_ev["POS"], read_ev["END"], svg_ev["POS"], svg_ev["END"])
    if ro >= MIN_RO and length_ratio(read_ev["SVLEN"], svg_ev["SVLEN"]) >= 0.5:
        return "reciprocal_overlap"
    if abs(read_ev["POS"] - svg_ev["POS"]) <= DIST or abs(read_ev["END"] - svg_ev["END"]) <= DIST:
        if length_ratio(read_ev["SVLEN"], svg_ev["SVLEN"]) >= 0.5:
            return "breakpoint_distance"
    return None


def match_svgap(index, svg_ev):
    seen = set()
    matches = []
    for key in bins_for(svg_ev["CHROM"], max(0, svg_ev["POS"] - INS_DIST), svg_ev["END"] + INS_DIST):
        for read_ev in index.get(key, []):
            if read_ev["SV_ID"] in seen:
                continue
            seen.add(read_ev["SV_ID"])
            reason = match_reason(read_ev, svg_ev)
            if reason:
                matches.append((read_ev, reason))
    return matches


def write_outputs(events, pav_path, svg_dir, out_events, out_pav, out_evidence, out_stats):
    index = defaultdict(list)
    for event in events:
        add_index(index, event)
    stats = Counter()
    with open(out_evidence, "w", newline="") as ev_fh:
        ev_w = csv.writer(ev_fh, delimiter="\t", lineterminator="\n")
        ev_w.writerow(["svgap_id", "read_SV_ID", "chrom", "svgap_start", "svgap_end", "read_start", "read_end", "svtype", "match_reason", "svgap_source_file", "svgap_source_id", "svgap_support_samples"])
        for svg_ev in iter_svgap_events(svg_dir):
            stats["svgap_events_seen"] += 1
            matches = match_svgap(index, svg_ev)
            if not matches:
                stats["svgap_only_not_retained"] += 1
                continue
            for read_ev, reason in matches:
                read_ev["svgap_count"] += 1
                if len(read_ev["svgap_ids"]) < MAX_IDS_PER_READ:
                    read_ev["svgap_ids"].append(svg_ev["svgap_id"])
                if len(read_ev["svgap_source_ids"]) < MAX_IDS_PER_READ:
                    read_ev["svgap_source_ids"].append(svg_ev["source_id"])
                read_ev["svgap_samples"].update(svg_ev["support_samples"])
                ev_w.writerow([
                    svg_ev["svgap_id"], read_ev["SV_ID"], svg_ev["CHROM"], svg_ev["POS"], svg_ev["END"],
                    read_ev["POS"], read_ev["END"], svg_ev["SVTYPE"], reason, svg_ev["source_file"],
                    svg_ev["source_id"], ",".join(sorted(svg_ev["support_samples"]))
                ])
                stats["svgap_overlap_evidence"] += 1

    with open(out_events, "w", newline="") as ev_out, open(out_pav, "w", newline="") as pav_out:
        ev_w = csv.writer(ev_out, delimiter="\t", lineterminator="\n")
        pav_w = csv.writer(pav_out, delimiter="\t", lineterminator="\n")
        ev_w.writerow(["FinalSV_ID", "source_methods", "reference", "CHROM", "POS", "END", "SVTYPE", "SVLEN", "read_SV_ID", "svgap_support_count", "svgap_ids", "svgap_support_samples", "svgap_source_id", "final_confidence", "notes"])
        pav_w.writerow(["SV_ID", "source_methods", "reference", "CHROM", "POS", "END", "SVTYPE", "SVLEN", "read_SV_ID", "svgap_ids", "final_confidence"] + FINAL_SAMPLES)
        for i, ev in enumerate(events, start=1):
            final_id = "Msa_flankQC_panSV_%09d" % i
            source = "read_based+SVGAP" if ev["svgap_count"] else "read_based"
            conf = "read_svgap_supported" if ev["svgap_count"] else "read_based"
            svg_ids = ",".join(ev["svgap_ids"]) if ev["svgap_ids"] else "."
            svg_samples = ",".join(sorted(ev["svgap_samples"])) if ev["svgap_samples"] else "."
            svg_src = ",".join(ev["svgap_source_ids"]) if ev["svgap_source_ids"] else "."
            ev_w.writerow([final_id, source, "genome_Msa", ev["CHROM"], ev["POS"], ev["END"], ev["SVTYPE"], ev["SVLEN"], ev["SV_ID"], ev["svgap_count"], svg_ids, svg_samples, svg_src, conf, "read_mapping_primary;svgap_evidence_only;flankQC_pass"])
            pav_w.writerow([final_id, source, "genome_Msa", ev["CHROM"], ev["POS"], ev["END"], ev["SVTYPE"], ev["SVLEN"], ev["SV_ID"], svg_ids, conf] + ev["pav"])
            stats["final_events"] += 1
            stats["source_methods:" + source] += 1
            stats["SVTYPE:" + ev["SVTYPE"]] += 1

    with open(out_stats, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["section", "key", "count"])
        w.writerow(["parameters", "reference", "genome_Msa"])
        w.writerow(["parameters", "event_policy", "read_mapping_primary_svgap_evidence_only_no_svgap_only"])
        w.writerow(["parameters", "flankQC", "applied_before_panSV_merge"])
        for key in sorted(stats):
            if ":" in key:
                section, sub = key.split(":", 1)
                w.writerow([section, sub, stats[key]])
            else:
                w.writerow(["rows", key, stats[key]])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--read-pav", required=True)
    ap.add_argument("--svgap-dir", required=True)
    ap.add_argument("--out-events", required=True)
    ap.add_argument("--out-pav", required=True)
    ap.add_argument("--out-evidence", required=True)
    ap.add_argument("--out-stats", required=True)
    args = ap.parse_args()
    events, _index = read_pav(args.read_pav)
    write_outputs(events, Path(args.read_pav), Path(args.svgap_dir), args.out_events, args.out_pav, args.out_evidence, args.out_stats)


if __name__ == "__main__":
    main()
