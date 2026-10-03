#!/usr/bin/env python3
import argparse
import csv
import os
from collections import defaultdict

from sv_utils import iter_vcf_records, overlap_len


def read_table(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_phenotype(path):
    data = {}
    for row in read_table(path):
        row["include_primary"] = str(row.get("include_primary", "1")) == "1"
        data[row["sample"]] = row
    return data


def expected_match(state, phenotype, expected_state, positive, negative):
    if expected_state == "deleted_in_spineless":
        if phenotype == negative:
            return state == "large_DEL"
        if phenotype == positive:
            return state in {"intact", "small_indel_only"}
    if expected_state == "deleted_in_spiny":
        if phenotype == positive:
            return state == "large_DEL"
        if phenotype == negative:
            return state in {"intact", "small_indel_only"}
    return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--targets", required=True)
    parser.add_argument("--vcf-manifest", required=True)
    parser.add_argument("--positive", required=True)
    parser.add_argument("--negative", required=True)
    parser.add_argument("--min-callers", type=int, required=True)
    parser.add_argument("--min-svlen", type=int, required=True)
    parser.add_argument("--small-indel-max-len", type=int, required=True)
    parser.add_argument("--min-overlap-fraction", type=float, required=True)
    parser.add_argument("--min-match", type=float, required=True)
    parser.add_argument("--states-out", required=True)
    parser.add_argument("--candidates-out", required=True)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.states_out), exist_ok=True)
    os.makedirs(os.path.dirname(args.candidates_out), exist_ok=True)

    phenotypes = read_phenotype(args.phenotype)
    targets = read_table(args.targets)
    manifest = read_table(args.vcf_manifest)

    vcf_by_sample_ref = defaultdict(list)
    for row in manifest:
        vcf_by_sample_ref[(row["sample"], row["ref"])].append(row)

    state_rows = []
    candidate_rows = []

    for target in targets:
        target_id = target["target_id"]
        ref = target["ref"]
        chrom = target["chrom"]
        start = int(target["start"])
        end = int(target["end"])
        target_len = end - start + 1
        expected_state = target.get("expected_state", "deleted_in_spineless")

        per_sample = []
        for sample, pheno in phenotypes.items():
            caller_hits = {}
            small_hits = {}
            detail = []
            for vcf_row in vcf_by_sample_ref.get((sample, ref), []):
                caller = vcf_row["caller"]
                vcf = vcf_row["vcf"]
                if not os.path.exists(vcf):
                    continue
                for rec in iter_vcf_records(vcf, chrom=chrom, start=start, end=end):
                    if rec["filter"] not in {"PASS", "."}:
                        continue
                    if rec["svtype"] != "DEL":
                        continue
                    svlen = abs(int(rec["svlen"]))
                    ov = overlap_len(start, end, rec["pos"], rec["end"])
                    ov_frac = ov / target_len if target_len else 0.0
                    if ov == 0:
                        continue
                    rec_desc = f"{caller}:{rec['chrom']}:{rec['pos']}-{rec['end']}:{svlen}:ov={ov_frac:.2f}"
                    if svlen >= args.min_svlen and ov_frac >= args.min_overlap_fraction:
                        caller_hits[caller] = max(caller_hits.get(caller, 0), svlen)
                        detail.append(rec_desc)
                    elif svlen <= args.small_indel_max_len:
                        small_hits[caller] = max(small_hits.get(caller, 0), svlen)

            n_callers = len(caller_hits)
            if n_callers >= args.min_callers:
                state = "large_DEL"
            elif small_hits:
                state = "small_indel_only"
            else:
                state = "intact"

            match = expected_match(state, pheno["phenotype"], expected_state, args.positive, args.negative)
            row = {
                "target_id": target_id,
                "sample": sample,
                "phenotype": pheno["phenotype"],
                "role": pheno.get("role", ""),
                "include_primary": int(pheno.get("include_primary", False)),
                "state": state,
                "n_large_del_callers": n_callers,
                "large_del_callers": ",".join(sorted(caller_hits)) if caller_hits else "NA",
                "small_indel_callers": ",".join(sorted(small_hits)) if small_hits else "NA",
                "match_expected": int(match),
                "records": ";".join(detail) if detail else "NA",
            }
            state_rows.append(row)
            per_sample.append(row)

        group_counts = {}
        for phenotype in (args.positive, args.negative):
            primary = [r for r in per_sample if r["include_primary"] and r["phenotype"] == phenotype]
            total = len(primary)
            match = sum(int(r["match_expected"]) for r in primary)
            group_counts[phenotype] = (match, total, match / total if total else 0.0)

        pos_match, pos_total, pos_rate = group_counts[args.positive]
        neg_match, neg_total, neg_rate = group_counts[args.negative]
        pass_flag = pos_rate >= args.min_match and neg_rate >= args.min_match
        candidate_rows.append({
            "candidate_id": target_id,
            "evidence_track": "gene_interval",
            "ref": ref,
            "chrom": chrom,
            "start": start,
            "end": end,
            "gene_id": target.get("gene_id", ""),
            "region_type": target.get("region_type", ""),
            "expected_state": expected_state,
            "positive_match": pos_match,
            "positive_total": pos_total,
            "positive_match_rate": f"{pos_rate:.3f}",
            "negative_match": neg_match,
            "negative_total": neg_total,
            "negative_match_rate": f"{neg_rate:.3f}",
            "pass_80pct": int(pass_flag),
            "priority_note": "gene_interval_pass" if pass_flag else "gene_interval_fail",
        })

    state_fields = [
        "target_id", "sample", "phenotype", "role", "include_primary", "state",
        "n_large_del_callers", "large_del_callers", "small_indel_callers",
        "match_expected", "records",
    ]
    with open(args.states_out, "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=state_fields)
        writer.writeheader()
        writer.writerows(state_rows)

    candidate_fields = [
        "candidate_id", "evidence_track", "ref", "chrom", "start", "end",
        "gene_id", "region_type", "expected_state",
        "positive_match", "positive_total", "positive_match_rate",
        "negative_match", "negative_total", "negative_match_rate",
        "pass_80pct", "priority_note",
    ]
    with open(args.candidates_out, "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=candidate_fields)
        writer.writeheader()
        writer.writerows(candidate_rows)


if __name__ == "__main__":
    main()
