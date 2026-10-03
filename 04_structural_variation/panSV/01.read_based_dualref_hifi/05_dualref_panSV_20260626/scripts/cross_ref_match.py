#!/usr/bin/env python3
import argparse
import gzip
import os
from collections import defaultdict


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def parse_info(info):
    d = {}
    for item in info.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d


def read_vcf(path):
    recs = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            info = parse_info(c[7])
            svtype = info.get("SVTYPE", "NA").split(",", 1)[0]
            try:
                pos = int(c[1])
                end = int(info.get("END", c[1]))
            except ValueError:
                continue
            svlen_s = str(info.get("SVLEN", end - pos)).split(",", 1)[0]
            try:
                svlen = abs(int(float(svlen_s)))
            except ValueError:
                svlen = abs(end - pos)
            recs.append({"id": c[2], "chrom": c[0], "pos": pos, "end": end,
                         "type": svtype, "len": svlen})
    return recs


def read_paf(path):
    by_query = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            q, qstart, qend = p[0], int(p[2]), int(p[3])
            strand = p[4]
            t, tstart, tend = p[5], int(p[7]), int(p[8])
            by_query[q].append((qstart, qend, strand, t, tstart, tend))
    for q in by_query:
        by_query[q].sort()
    return by_query


def lift_interval(paf, chrom, pos, end):
    qpos = max(0, pos - 1)
    qend = max(0, end - 1)
    for qs, qe, strand, t, ts, te in paf.get(chrom, []):
        if qs <= qpos <= qe and qs <= qend <= qe:
            if strand == "+":
                a = ts + (qpos - qs) + 1
                b = ts + (qend - qs) + 1
            else:
                a = te - (qpos - qs)
                b = te - (qend - qs)
            return t, min(a, b), max(a, b), strand
    return None


def reciprocal_overlap(a1, a2, b1, b2):
    ov = max(0, min(a2, b2) - max(a1, b1) + 1)
    la = max(1, a2 - a1 + 1)
    lb = max(1, b2 - b1 + 1)
    return min(float(ov) / la, float(ov) / lb)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refa-vcf", required=True)
    ap.add_argument("--refb-vcf", required=True)
    ap.add_argument("--paf", required=True)
    ap.add_argument("--out-events", required=True)
    ap.add_argument("--out-evidence", required=True)
    ap.add_argument("--distance", type=int, default=1000)
    ap.add_argument("--min-ro", type=float, default=0.5)
    args = ap.parse_args()

    a_recs = read_vcf(args.refa_vcf)
    b_recs = read_vcf(args.refb_vcf)
    paf = read_paf(args.paf)
    b_index = defaultdict(list)
    for b in b_recs:
        key = (b["chrom"], b["type"])
        for bin_id in range(b["pos"] // 100000 - 1, b["pos"] // 100000 + 2):
            b_index[(key, bin_id)].append(b)

    os.makedirs(os.path.dirname(args.out_events), exist_ok=True)
    matched_b = set()
    meta_i = 0
    with open(args.out_events, "w") as ev, open(args.out_evidence, "w") as evid:
        ev.write("MetaSV_ID\tRefA_SV_ID\tRefA_CHROM\tRefA_POS\tRefA_END\tRefB_SV_ID\tRefB_CHROM\tRefB_POS\tRefB_END\tSVTYPE\tSVLEN\tmatch_status\tmatch_method\tconfidence\n")
        evid.write("MetaSV_ID\tmethod\tdetails\n")
        for a in a_recs:
            meta_i += 1
            meta_id = "MetaSV_%08d" % meta_i
            if a["type"] == "INS":
                ev.write("%s\t%s\t%s\t%d\t%d\t.\t.\t.\t.\t%s\t%d\tRefA_only\tins_sequence_matching_pending\tlow\n" %
                         (meta_id, a["id"], a["chrom"], a["pos"], a["end"], a["type"], a["len"]))
                continue
            lifted = lift_interval(paf, a["chrom"], a["pos"], a["end"])
            best = None
            best_score = -1
            if lifted:
                tchrom, tpos, tend, strand = lifted
                for bin_id in range(tpos // 100000 - 1, tpos // 100000 + 2):
                    for b in b_index.get(((tchrom, a["type"]), bin_id), []):
                        if abs(b["pos"] - tpos) > args.distance or abs(b["end"] - tend) > args.distance:
                            continue
                        ratio = float(max(a["len"], b["len"])) / max(1, min(a["len"], b["len"]))
                        if ratio > 2.0:
                            continue
                        ro = reciprocal_overlap(tpos, tend, b["pos"], b["end"])
                        if ro < args.min_ro:
                            continue
                        score = ro - (abs(b["pos"] - tpos) + abs(b["end"] - tend)) / 1000000.0
                        if score > best_score:
                            best_score = score
                            best = (b, ro, strand, tpos, tend)
            if best:
                b, ro, strand, tpos, tend = best
                matched_b.add(b["id"])
                ev.write("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%d\tshared_by_liftover\tbreakpoint_liftover_ro\tmedium\n" %
                         (meta_id, a["id"], a["chrom"], a["pos"], a["end"], b["id"], b["chrom"], b["pos"], b["end"], a["type"], a["len"]))
                evid.write("%s\tliftover\tlifted=%s:%d-%d;strand=%s;reciprocal_overlap=%.4f\n" %
                           (meta_id, best[0]["chrom"], tpos, tend, strand, ro))
            else:
                ev.write("%s\t%s\t%s\t%d\t%d\t.\t.\t.\t.\t%s\t%d\tRefA_only\tno_liftover_match\tlow\n" %
                         (meta_id, a["id"], a["chrom"], a["pos"], a["end"], a["type"], a["len"]))
        for b in b_recs:
            if b["id"] in matched_b:
                continue
            meta_i += 1
            meta_id = "MetaSV_%08d" % meta_i
            ev.write("%s\t.\t.\t.\t.\t%s\t%s\t%d\t%d\t%s\t%d\tRefB_only\tno_liftover_match\tlow\n" %
                     (meta_id, b["id"], b["chrom"], b["pos"], b["end"], b["type"], b["len"]))


if __name__ == "__main__":
    main()
