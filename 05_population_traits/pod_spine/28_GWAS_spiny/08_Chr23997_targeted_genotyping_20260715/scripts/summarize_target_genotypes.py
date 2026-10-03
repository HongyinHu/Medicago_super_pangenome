#!/usr/bin/env python3
"""Combine exact-site callers and independent BAM evidence without imputing failures."""

import argparse
import csv
import math
import os
from collections import Counter

import pysam
from scipy.stats import fisher_exact


def numeric(value, default=0.0):
    if value is None:
        return default
    if isinstance(value, (tuple, list)):
        value = value[0] if value else None
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def gt_text(gt):
    if gt is None:
        return "./."
    return "/".join("." if allele is None else str(allele) for allele in gt)


def is_alt(gt):
    return bool(gt) and any(allele is not None and allele > 0 for allele in gt)


def is_ref(gt):
    return bool(gt) and all(allele == 0 for allele in gt if allele is not None) and all(
        allele is not None for allele in gt
    )


def read_variant(path):
    with pysam.VariantFile(path) as handle:
        records = list(handle)
        samples = list(handle.header.samples)
    if not records:
        return None, {}, 0
    if len(records) != 1:
        raise ValueError(f"Expected <=1 record in {path}, found {len(records)}")
    record = records[0]
    sample_data = dict(record.samples[samples[0]].items()) if samples else {}
    return record, sample_data, 1


def classify(sv, de, evidence):
    sv_gt = sv.get("GT")
    de_gt = de.get("GT")
    sv_su = numeric(sv.get("SU"))
    de_alt_reads = numeric(de.get("DV")) + numeric(de.get("RV"))
    direct = int(evidence["target_direct_alt_reads"])
    discordant = int(evidence["target_discordant_pair_reads"])
    clips = int(evidence["target_breakpoint_clip_reads"])
    depth_ratio = numeric(evidence["target_to_flank_depth_ratio"], float("nan"))
    evidence_call = evidence["evidence_call"]

    high_reasons = []
    if direct >= 2:
        high_reasons.append("direct_read>=2")
    if is_alt(sv_gt) and sv_su >= 3:
        high_reasons.append("svtyper_nonref_SU>=3")
    if is_alt(de_gt) and de_alt_reads >= 3:
        high_reasons.append("delly_nonref_alt_reads>=3")

    moderate_reasons = []
    if evidence_call in {
        "DEL_MULTIMODAL_SUPPORT",
        "DEL_PAIR_DEPTH_CLIP_SUPPORT",
    }:
        moderate_reasons.append(evidence_call)
    if direct >= 1 and (discordant >= 2 or clips >= 2):
        moderate_reasons.append("direct_plus_pair_or_clip")
    if discordant >= 3 and clips >= 2 and not math.isnan(depth_ratio) and depth_ratio <= 0.75:
        moderate_reasons.append("pair_clip_depth")

    flank_depth = numeric(evidence["flank_mean_depth"])
    left_ref = int(evidence["left_breakpoint_ref_reads"])
    right_ref = int(evidence["right_breakpoint_ref_reads"])
    sv_gq = numeric(sv.get("GQ"))
    strong_ref = (
        flank_depth >= 5
        and direct == 0
        and not high_reasons
        and not moderate_reasons
        and is_ref(sv_gt)
        and sv_gq >= 20
        and left_ref >= 3
        and right_ref >= 3
    )

    if high_reasons:
        strict = "DEL"
        strict_reason = ";".join(high_reasons)
    elif strong_ref:
        strict = "REF"
        strict_reason = "svtyper_0/0_GQ>=20_and_both_breakpoints_intact"
    else:
        strict = "MISSING"
        strict_reason = "insufficient_or_conflicting_evidence"

    if high_reasons or moderate_reasons:
        sensitivity = "DEL"
        sensitivity_reason = ";".join(high_reasons + moderate_reasons)
    elif strong_ref:
        sensitivity = "REF"
        sensitivity_reason = strict_reason
    else:
        sensitivity = "MISSING"
        sensitivity_reason = strict_reason
    return strict, strict_reason, sensitivity, sensitivity_reason


def write_proxy_vcf(path, rows, state_column):
    samples = [row["sample"] for row in rows]
    with open(path, "w", newline="") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##contig=<ID=4>\n")
        handle.write('##INFO=<ID=EVENT,Number=1,Type=String,Description="Validated biological event">\n')
        handle.write('##INFO=<ID=SOURCE_COORD,Number=1,Type=String,Description="Exact source coordinate on Msa reference">\n')
        handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Consensus presence genotype; A/T is a PLINK-compatible proxy encoding">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        handle.write("\t".join(samples) + "\n")
        calls = {"DEL": "1/1", "REF": "0/0", "MISSING": "./."}
        genotypes = [calls[row[state_column]] for row in rows]
        handle.write(
            "4\t88980831\tChr23997_intron2_DEL_221bp\tA\tT\t.\tPASS\t"
            "EVENT=Chr23997_intron2_DEL_221bp;SOURCE_COORD=Chr4:88980831-88981052\tGT\t"
            + "\t".join(genotypes)
            + "\n"
        )


def association(rows, state_column, label):
    counts = Counter((int(row["phenotype"]), row[state_column]) for row in rows)
    spineless_del = counts[(0, "DEL")]
    spineless_ref = counts[(0, "REF")]
    spiny_del = counts[(1, "DEL")]
    spiny_ref = counts[(1, "REF")]
    called = spineless_del + spineless_ref + spiny_del + spiny_ref
    missing = len(rows) - called
    spineless_total = sum(int(row["phenotype"]) == 0 for row in rows)
    spiny_total = sum(int(row["phenotype"]) == 1 for row in rows)
    spineless_called = spineless_del + spineless_ref
    spiny_called = spiny_del + spiny_ref
    spineless_missing = spineless_total - spineless_called
    spiny_missing = spiny_total - spiny_called
    carriers = spineless_del + spiny_del
    refs = spineless_ref + spiny_ref
    if carriers and refs and (spineless_del + spineless_ref) and (spiny_del + spiny_ref):
        odds_ratio, pvalue = fisher_exact(
            [[spineless_del, spineless_ref], [spiny_del, spiny_ref]],
            alternative="two-sided",
        )
    else:
        odds_ratio, pvalue = float("nan"), float("nan")
    callability_odds, callability_p = fisher_exact(
        [[spineless_called, spineless_missing], [spiny_called, spiny_missing]],
        alternative="two-sided",
    )
    return {
        "callset": label,
        "total_samples": len(rows),
        "called_samples": called,
        "missing_samples": missing,
        "call_rate": f"{called / len(rows):.8f}",
        "DEL_carriers": carriers,
        "REF_samples": refs,
        "spineless_DEL": spineless_del,
        "spineless_REF": spineless_ref,
        "spiny_DEL": spiny_del,
        "spiny_REF": spiny_ref,
        "spineless_missing": spineless_missing,
        "spiny_missing": spiny_missing,
        "fisher_odds_DEL_spineless_vs_spiny": f"{odds_ratio:.12g}",
        "fisher_two_sided_p": f"{pvalue:.12g}",
        "callability_fisher_odds_spineless_vs_spiny": f"{callability_odds:.12g}",
        "callability_fisher_two_sided_p": f"{callability_p:.12g}",
        "gemma_eligible": int(carriers > 0 and refs > 0 and called >= 20),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--per-sample-dir", required=True)
    parser.add_argument("--summary-dir", required=True)
    parser.add_argument("--vcf-dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.summary_dir, exist_ok=True)
    os.makedirs(args.vcf_dir, exist_ok=True)
    with open(args.manifest, newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    if len(manifest) != 143:
        raise ValueError(f"Expected 143 samples, found {len(manifest)}")

    rows = []
    for item in manifest:
        sample = item["sample"]
        sample_dir = os.path.join(args.per_sample_dir, sample)
        done = os.path.join(sample_dir, f"{sample}.done")
        if not os.path.getsize(done):
            raise ValueError(f"Missing completion marker: {done}")
        _, sv, sv_records = read_variant(
            os.path.join(sample_dir, f"{sample}-smoove.genotyped.vcf.gz")
        )
        _, de, de_records = read_variant(os.path.join(sample_dir, f"{sample}.delly.bcf"))
        if sv_records != 1:
            raise ValueError(f"svtyper did not emit exactly one record for {sample}")
        with open(os.path.join(sample_dir, f"{sample}.evidence.tsv"), newline="") as handle:
            evidence_rows = list(csv.DictReader(handle, delimiter="\t"))
        if len(evidence_rows) != 1:
            raise ValueError(f"Expected one evidence row for {sample}")
        evidence = evidence_rows[0]

        strict, strict_reason, sensitivity, sensitivity_reason = classify(sv, de, evidence)
        row = {
            "task_index": item["task_index"],
            "sample": sample,
            "phenotype": item["phenotype"],
            "strict_state": strict,
            "strict_reason": strict_reason,
            "sensitivity_state": sensitivity,
            "sensitivity_reason": sensitivity_reason,
            "svtyper_GT": gt_text(sv.get("GT")),
            "svtyper_GQ": numeric(sv.get("GQ"), float("nan")),
            "svtyper_PE": numeric(sv.get("PE")),
            "svtyper_SR": numeric(sv.get("SR")),
            "svtyper_SU": numeric(sv.get("SU")),
            "delly_record_count": de_records,
            "delly_GT": gt_text(de.get("GT")),
            "delly_GQ": numeric(de.get("GQ"), float("nan")),
            "delly_DR": numeric(de.get("DR")),
            "delly_DV": numeric(de.get("DV")),
            "delly_RR": numeric(de.get("RR")),
            "delly_RV": numeric(de.get("RV")),
        }
        row.update(evidence)
        rows.append(row)

    table_path = os.path.join(args.summary_dir, "Chr23997.targeted_genotypes.143.tsv")
    with open(table_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    strict_vcf = os.path.join(args.vcf_dir, "Chr23997.target.strict.proxy.vcf")
    sensitivity_vcf = os.path.join(args.vcf_dir, "Chr23997.target.sensitivity.proxy.vcf")
    write_proxy_vcf(strict_vcf, rows, "strict_state")
    write_proxy_vcf(sensitivity_vcf, rows, "sensitivity_state")

    stats = [
        association(rows, "strict_state", "strict"),
        association(rows, "sensitivity_state", "sensitivity"),
    ]
    stats_path = os.path.join(args.summary_dir, "Chr23997.targeted_association.fisher.tsv")
    with open(stats_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(stats[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(stats)

    reason_counts = Counter()
    for row in rows:
        reason_counts[(row["strict_state"], row["strict_reason"])] += 1
    with open(os.path.join(args.summary_dir, "strict_call_reason_counts.tsv"), "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["strict_state", "strict_reason", "sample_count"])
        for (state, reason), count in sorted(reason_counts.items()):
            writer.writerow([state, reason, count])

    for stat in stats:
        state_column = f"{stat['callset']}_state"
        phenotype_path = os.path.join(
            args.summary_dir, f"{stat['callset']}.phenotype.called_only.txt"
        )
        with open(phenotype_path, "w") as handle:
            for row in rows:
                value = row["phenotype"] if row[state_column] != "MISSING" else "NA"
                handle.write(value + "\n")
        carrier_path = os.path.join(
            args.summary_dir, f"{stat['callset']}.DEL_carriers.txt"
        )
        with open(carrier_path, "w") as handle:
            for row in rows:
                if row[state_column] == "DEL":
                    handle.write(row["sample"] + "\n")
        eligible = os.path.join(args.summary_dir, f"{stat['callset']}.gemma_eligible")
        if stat["gemma_eligible"]:
            with open(eligible, "w") as handle:
                handle.write("eligible\n")


if __name__ == "__main__":
    main()
