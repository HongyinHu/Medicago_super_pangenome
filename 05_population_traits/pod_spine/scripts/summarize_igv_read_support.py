#!/usr/bin/env python3
import csv
import sys
from collections import Counter, defaultdict


def ratio_group(n, d):
    return f"{n}/{d}"


def main():
    if len(sys.argv) != 4:
        raise SystemExit("usage: summarize_igv_read_support.py <candidate_tsv> <metrics_tsv> <out_tsv>")
    candidate_tsv, metrics_tsv, out_tsv = sys.argv[1:]
    candidates = {}
    with open(candidate_tsv, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            candidates[row["candidate_id"]] = row

    grouped = defaultdict(list)
    with open(metrics_tsv, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            grouped[row["candidate_id"]].append(row)

    out_rows = []
    for cid, rows in grouped.items():
        c = candidates.get(cid, {})
        calls = Counter(r["support_call"] for r in rows)
        phen_calls = defaultdict(Counter)
        del_reads = defaultdict(list)
        cov_ratio = defaultdict(list)
        for r in rows:
            phen_calls[r["phenotype"]][r["support_call"]] += 1
            del_reads[r["phenotype"]].append(int(r["del_cigar_support_reads"]))
            try:
                cov_ratio[r["phenotype"]].append(float(r["event_to_flank_cov_ratio"]))
            except ValueError:
                pass

        direction = rows[0]["expected_pattern"]
        if direction == "spiny_present":
            expected_present_pheno = "spiny"
            expected_absent_pheno = "spineless"
        else:
            expected_present_pheno = "spineless"
            expected_absent_pheno = "spiny"

        expected_present_rows = [r for r in rows if r["phenotype"] == expected_present_pheno]
        expected_absent_rows = [r for r in rows if r["phenotype"] == expected_absent_pheno]
        present_ok = sum(r["support_call"] == "supports_expected_presence" for r in expected_present_rows)
        absent_ok = sum(r["support_call"] == "supports_expected_absence" for r in expected_absent_rows)
        present_fail = sum(r["support_call"] == "lacks_expected_presence_support" for r in expected_present_rows)
        absent_fail = sum(r["support_call"] == "unexpected_presence_support" for r in expected_absent_rows)
        present_n = len(expected_present_rows)
        absent_n = len(expected_absent_rows)

        if present_ok >= max(3, present_n - 1) and absent_fail == 0:
            verdict = "IGV_priority_confirm_likely"
        elif present_ok >= 2 and absent_fail <= 1:
            verdict = "IGV_manual_review_possible"
        elif absent_fail >= 2:
            verdict = "IGV_suspect_unexpected_absent_group_support"
        else:
            verdict = "IGV_suspect_low_expected_support"

        def mean(vals):
            return round(sum(vals) / len(vals), 3) if vals else "NA"

        out_rows.append(
            {
                "candidate_id": cid,
                "MetaSV_ID": rows[0]["MetaSV_ID"],
                "ref": rows[0]["ref"],
                "locus": f"{rows[0]['chrom']}:{rows[0]['start']}-{rows[0]['end']}",
                "SVTYPE": rows[0]["SVTYPE"],
                "SVLEN": rows[0]["SVLEN"],
                "association_direction": direction,
                "expected_present_group": expected_present_pheno,
                "expected_presence_supported": ratio_group(present_ok, present_n),
                "expected_absence_supported": ratio_group(absent_ok, absent_n),
                "expected_presence_fail": present_fail,
                "unexpected_absent_group_presence": absent_fail,
                "spiny_mean_del_support_reads": mean(del_reads["spiny"]),
                "spineless_mean_del_support_reads": mean(del_reads["spineless"]),
                "spiny_mean_event_to_flank_cov_ratio": mean(cov_ratio["spiny"]),
                "spineless_mean_event_to_flank_cov_ratio": mean(cov_ratio["spineless"]),
                "support_verdict": verdict,
                "support_call_counts": ";".join(f"{k}={v}" for k, v in sorted(calls.items())),
            }
        )

    rank_order = {
        "IGV_priority_confirm_likely": 0,
        "IGV_manual_review_possible": 1,
        "IGV_suspect_low_expected_support": 2,
        "IGV_suspect_unexpected_absent_group_support": 3,
    }
    out_rows.sort(key=lambda r: (rank_order.get(r["support_verdict"], 9), int(r["candidate_id"].split("_", 1)[0])))
    with open(out_tsv, "w", newline="") as f:
        fields = list(out_rows[0].keys())
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out_rows)


if __name__ == "__main__":
    main()
