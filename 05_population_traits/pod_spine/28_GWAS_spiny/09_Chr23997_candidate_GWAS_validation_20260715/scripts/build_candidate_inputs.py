#!/usr/bin/env python3
"""Build typed Chr23997 candidate-SV inputs without imputing ambiguous calls."""

import argparse
import csv
from collections import Counter
from pathlib import Path


TARGET_CHROM = "Chr4"
TARGET_POS = 88980831
TARGET_END = 88981052
TARGET_SVLEN = -221


def state_to_gt(state):
    """Map evidence state to the candidate VCF convention."""
    normalized = state.strip().upper()
    if normalized == "DEL":
        return "1/1"
    if normalized in {"REF", "PRESENT"}:
        return "0/0"
    if normalized == "MISSING":
        return "./."
    raise ValueError("unsupported evidence state: {}".format(state))


def assert_same_sample_order(left, right):
    if list(left) != list(right):
        raise ValueError("sample order mismatch between required inputs")


def read_tsv(path):
    with open(path, "r", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fieldnames, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sample_order_from_keep(path):
    samples = []
    with open(path, "r") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError("invalid keep-file row: {}".format(line.rstrip()))
            samples.append(fields[1])
    if len(samples) != len(set(samples)):
        raise ValueError("duplicate IID in sample-order file")
    return samples


def phenotype_map(path):
    records = read_tsv(path)
    required = {"Sample_ID", "Latin_name", "Section", "Pod_spine", "Pod_spine_binary_spiny1_spineless0"}
    if not records or not required.issubset(records[0]):
        raise ValueError("phenotype table lacks required fields")
    mapping = {}
    for row in records:
        sample = row["Sample_ID"]
        if sample in mapping:
            raise ValueError("duplicate sample in phenotype table: {}".format(sample))
        mapping[sample] = row
    return mapping


def canonical_state(value):
    value = value.strip().upper()
    if value == "DEL":
        return "DEL"
    if value in {"REF", "PRESENT"}:
        return "PRESENT"
    if value == "MISSING":
        return "MISSING"
    raise ValueError("unsupported state in target genotype table: {}".format(value))


def build_rows(target_records, order, phenotypes, state_column, reason_column):
    by_sample = {}
    for row in target_records:
        sample = row["sample"]
        if sample in by_sample:
            raise ValueError("duplicate sample in target genotype table: {}".format(sample))
        by_sample[sample] = row
    if set(order) != set(by_sample):
        missing = sorted(set(order) - set(by_sample))
        extra = sorted(set(by_sample) - set(order))
        raise ValueError("target genotype samples do not match analysis samples; missing={}, extra={}".format(missing[:5], extra[:5]))
    if not set(order).issubset(phenotypes):
        missing = sorted(set(order) - set(phenotypes))
        raise ValueError("phenotypes missing analysis samples: {}".format(missing[:5]))

    evidence_fields = [
        "svtyper_GT", "svtyper_GQ", "svtyper_PE", "svtyper_SR", "svtyper_SU",
        "delly_record_count", "delly_GT", "delly_GQ", "delly_DR", "delly_DV",
        "flank_mean_depth", "target_to_flank_depth_ratio", "target_fraction_depth_ge5",
        "target_del_cigar_reads", "target_split_alignment_reads", "target_direct_alt_reads",
        "target_discordant_pair_reads", "left_breakpoint_ref_reads", "right_breakpoint_ref_reads",
        "target_breakpoint_clip_reads", "evidence_call",
    ]
    rows = []
    for sample in order:
        source = by_sample[sample]
        meta = phenotypes[sample]
        state = canonical_state(source[state_column])
        gt = state_to_gt(state)
        row = {
            "sample": sample,
            "phenotype": meta["Pod_spine_binary_spiny1_spineless0"],
            "pod_spine": meta["Pod_spine"],
            "latin_name": meta["Latin_name"],
            "population_section": meta["Section"],
            "state": state,
            "genotype": "1" if state == "DEL" else ("0" if state == "PRESENT" else "NA"),
            "GT": gt,
            "call_reason": source.get(reason_column, ""),
        }
        for field in evidence_fields:
            row[field] = source.get(field, "")
        rows.append(row)
    return rows


def write_candidate_vcf(path, rows, callset):
    with open(path, "w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##source=Chr23997_candidate_targeted_evidence_{}\n".format(callset))
        handle.write("##contig=<ID=Chr4>\n")
        handle.write("##INFO=<ID=END,Number=1,Type=Integer,Description=End position of deletion>\n")
        handle.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=Structural variant type>\n")
        handle.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=Length of deletion>\n")
        handle.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=Candidate genotype: 0/0=presence,1/1=deletion>\n")
        samples = [row["sample"] for row in rows]
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n".format("\t".join(samples)))
        info = "END={};SVTYPE=DEL;SVLEN={}".format(TARGET_END, TARGET_SVLEN)
        gts = "\t".join(row["GT"] for row in rows)
        handle.write("{}\t{}\tChr23997_intron2_DEL221\tN\t<DEL>\t.\tPASS\t{}\tGT\t{}\n".format(TARGET_CHROM, TARGET_POS, info, gts))


def write_qc(outdir, callset, rows):
    counts = Counter(row["state"] for row in rows)
    by_pheno = {}
    for phenotype in ("0", "1"):
        subset = [row for row in rows if row["phenotype"] == phenotype]
        subset_counts = Counter(row["state"] for row in subset)
        by_pheno[phenotype] = subset_counts
    called = counts["DEL"] + counts["PRESENT"]
    qc_rows = [
        {"callset": callset, "metric": "total_samples", "value": len(rows)},
        {"callset": callset, "metric": "called_samples", "value": called},
        {"callset": callset, "metric": "missing_samples", "value": counts["MISSING"]},
        {"callset": callset, "metric": "call_rate", "value": "{:.8f}".format(called / float(len(rows)))},
        {"callset": callset, "metric": "DEL_samples", "value": counts["DEL"]},
        {"callset": callset, "metric": "PRESENT_samples", "value": counts["PRESENT"]},
        {"callset": callset, "metric": "DEL_frequency_among_called", "value": "{:.8f}".format(counts["DEL"] / float(called)) if called else "NA"},
    ]
    for phenotype, label in (("0", "spineless"), ("1", "spiny")):
        for state in ("DEL", "PRESENT", "MISSING"):
            qc_rows.append({"callset": callset, "metric": "{}_{}".format(label, state), "value": by_pheno[phenotype][state]})
    return qc_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-genotypes", required=True)
    parser.add_argument("--sample-order", required=True)
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    order = sample_order_from_keep(args.sample_order)
    phenotypes = phenotype_map(args.phenotype)
    target_records = read_tsv(args.target_genotypes)
    if not target_records or "sample" not in target_records[0]:
        raise ValueError("target genotype table is empty or lacks sample column")

    all_qc = []
    combined_rows = []
    for callset, state_column, reason_column in (
        ("strict", "strict_state", "strict_reason"),
        ("sensitivity", "sensitivity_state", "sensitivity_reason"),
    ):
        rows = build_rows(target_records, order, phenotypes, state_column, reason_column)
        fields = list(rows[0])
        write_tsv(outdir / "Chr23997_221bp_PAV.{}.tsv".format(callset), fields, rows)
        write_candidate_vcf(outdir / "Chr23997_221bp_PAV.{}.vcf".format(callset), rows, callset)
        all_qc.extend(write_qc(outdir, callset, rows))
        combined_rows.extend([dict(row, callset=callset) for row in rows])

    write_tsv(outdir / "genotype_qc.tsv", ["callset", "metric", "value"], all_qc)
    write_tsv(outdir / "genotype_by_sample.tsv", ["callset"] + list(combined_rows[0]), combined_rows)
    contingency = []
    for row in all_qc:
        if row["metric"].startswith(("spineless_", "spiny_")):
            group, state = row["metric"].split("_", 1)
            contingency.append({"callset": row["callset"], "phenotype_group": group, "state": state, "n": row["value"]})
    write_tsv(outdir / "genotype_case_control_distribution.tsv", ["callset", "phenotype_group", "state", "n"], contingency)


if __name__ == "__main__":
    main()
