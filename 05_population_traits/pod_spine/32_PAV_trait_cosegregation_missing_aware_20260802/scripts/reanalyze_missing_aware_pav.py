#!/usr/bin/env python3
"""Missing-aware PAV/phenotype association reanalysis.

This preserves NA as missing, applies an explicit call-rate filter, replaces the
three documented false-positive Chr23997 calls only in the corrected analysis,
and recalculates Fisher exact P, BH-FDR, Bonferroni and exhaustive 7-vs-7 maxT.
"""

import argparse
import csv
import hashlib
import itertools
import json
import math
from collections import Counter
from datetime import datetime
from pathlib import Path

import numpy as np


CHROMS = ["Chr%d" % i for i in range(1, 9)]
DEL_GROUP = [
    "genome_395", "genome_461", "genome_468", "genome_472",
    "genome_M22", "genome_Mar", "genome_Mru",
]
INTACT_GROUP = [
    "genome_410", "genome_436", "genome_457", "genome_454",
    "genome_474", "genome_Mpo", "genome_R108",
]
SAMPLES = DEL_GROUP + INTACT_GROUP
TARGET_ID = "Msa_panSV_000339055"
TARGET_CHROM = "Chr4"
TARGET_POS = 88980831
TARGET_END = 88981052
MANUAL_ZERO_SAMPLES = {"genome_436", "genome_457", "genome_474"}
MODELS = [
    ("complete_case", 14, 7),
    ("primary_callrate_ge80", 12, 5),
    ("relaxed_callrate_ge70", 10, 4),
]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def comb(n, k):
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)


def fisher_two_sided(a, b, c, d):
    """Probability-ordering two-sided Fisher exact test for [[a,b],[c,d]]."""
    n = a + b + c + d
    if n == 0:
        return 1.0
    row1 = a + b
    col1 = a + c
    lo = max(0, row1 - (n - col1))
    hi = min(row1, col1)
    denom = float(comb(n, row1))

    def prob(x):
        return comb(col1, x) * comb(n - col1, row1 - x) / denom

    observed = prob(a)
    return min(1.0, sum(prob(x) for x in range(lo, hi + 1)
                        if prob(x) <= observed + 1e-15))


def bh_adjust(pvalues):
    m = len(pvalues)
    order = sorted(range(m), key=lambda i: pvalues[i])
    qvalues = [1.0] * m
    running = 1.0
    for pos in range(m - 1, -1, -1):
        idx = order[pos]
        rank = pos + 1
        running = min(running, pvalues[idx] * m / float(rank))
        qvalues[idx] = min(1.0, running)
    return qvalues


def parse_call(value):
    value = str(value).strip()
    if value in {"1", "1.0"}:
        return 1
    if value in {"0", "0.0"}:
        return 0
    return None


def calculate_counts(calls):
    del_calls = [calls[s] for s in DEL_GROUP if calls[s] is not None]
    intact_calls = [calls[s] for s in INTACT_GROUP if calls[s] is not None]
    dp = sum(del_calls)
    ip = sum(intact_calls)
    dn = len(del_calls)
    inn = len(intact_calls)
    da = dn - dp
    ia = inn - ip
    return dp, da, ip, ia, dn, inn


def test_calls(calls):
    dp, da, ip, ia, dn, inn = calculate_counts(calls)
    called_total = dn + inn
    missing_total = len(SAMPLES) - called_total
    del_freq = dp / float(dn) if dn else None
    intact_freq = ip / float(inn) if inn else None
    pvalue = fisher_two_sided(dp, da, ip, ia) if dn and inn else 1.0
    if del_freq is None or intact_freq is None:
        direction = "NA"
        abs_delta = 0.0
    elif del_freq > intact_freq:
        direction = "DEL_group_enriched"
        abs_delta = del_freq - intact_freq
    elif intact_freq > del_freq:
        direction = "INTACT_group_enriched"
        abs_delta = intact_freq - del_freq
    else:
        direction = "tie"
        abs_delta = 0.0
    complete = bool(dn and inn and ((dp == dn and ip == 0) or (ip == inn and dp == 0)))
    group80 = bool(
        del_freq is not None and intact_freq is not None and
        ((del_freq >= 0.8 and intact_freq <= 0.2) or
         (intact_freq >= 0.8 and del_freq <= 0.2))
    )
    return {
        "del_present": dp,
        "del_absent": da,
        "del_called": dn,
        "intact_present": ip,
        "intact_absent": ia,
        "intact_called": inn,
        "called_total": called_total,
        "missing_total": missing_total,
        "del_freq": del_freq,
        "intact_freq": intact_freq,
        "abs_frequency_difference": abs_delta,
        "direction": direction,
        "fisher_p": pvalue,
        "neglog10p": -math.log10(max(pvalue, 1e-300)),
        "complete_separation": complete,
        "group_specific80": group80,
    }


def passes_model(stats, min_total, min_group):
    return (
        stats["called_total"] >= min_total and
        stats["del_called"] >= min_group and
        stats["intact_called"] >= min_group
    )


def calls_to_masks(calls):
    called_mask = 0
    present_mask = 0
    for i, sample in enumerate(SAMPLES):
        value = calls[sample]
        if value is None:
            continue
        called_mask |= (1 << i)
        if value == 1:
            present_mask |= (1 << i)
    return called_mask, present_mask


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_permutation_max_t(patterns, target_pattern, nominal_p_values):
    all_mask = (1 << 14) - 1
    bitcount = np.array([bin(int(i)).count("1") for i in range(1 << 14)], dtype=np.int8)
    called_masks = np.array([p[0] for p in patterns], dtype=np.uint16)
    present_masks = np.array([p[1] for p in patterns], dtype=np.uint16)

    p_table = np.ones((8, 8, 8, 8), dtype=np.float64)
    for dp in range(8):
        for da in range(8):
            for ip in range(8):
                for ia in range(8):
                    if dp + da > 0 and ip + ia > 0:
                        p_table[dp, da, ip, ia] = fisher_two_sided(dp, da, ip, ia)

    label_masks = []
    for combo in itertools.combinations(range(14), 7):
        mask = 0
        for i in combo:
            mask |= (1 << i)
        label_masks.append(mask)

    min_p = np.empty(len(label_masks), dtype=np.float64)
    target_p = np.empty(len(label_masks), dtype=np.float64)
    target_called, target_present = target_pattern

    for i, label_mask in enumerate(label_masks):
        complement = all_mask ^ label_mask
        dc = bitcount[called_masks & label_mask]
        ic = bitcount[called_masks & complement]
        dp = bitcount[present_masks & label_mask]
        ip = bitcount[present_masks & complement]
        pvals = p_table[dp, dc - dp, ip, ic - ip]
        min_p[i] = float(np.min(pvals))

        tdc = int(bitcount[target_called & label_mask])
        tic = int(bitcount[target_called & complement])
        tdp = int(bitcount[target_present & label_mask])
        tip = int(bitcount[target_present & complement])
        target_p[i] = p_table[tdp, tdc - tdp, tip, tic - tip]

    observed_mask = (1 << 7) - 1
    observed_index = label_masks.index(observed_mask)
    observed_target_p = float(target_p[observed_index])
    point_empirical_p = float(np.mean(target_p <= observed_target_p + 1e-15))
    max_t_empirical_p = float(np.mean(min_p <= observed_target_p + 1e-15))

    unique_nominal = sorted(set(nominal_p_values))
    adjusted_lookup = {}
    for pvalue in unique_nominal:
        adjusted_lookup[pvalue] = float(np.mean(min_p <= pvalue + 1e-15))
    max_t_adjusted = [adjusted_lookup[p] for p in nominal_p_values]
    attainable = [p for p in unique_nominal if adjusted_lookup[p] <= 0.05 + 1e-15]

    min_p_counts = Counter("%.15g" % value for value in min_p.tolist())
    return {
        "permutations": len(label_masks),
        "unique_patterns": len(patterns),
        "observed_target_p": observed_target_p,
        "target_pointwise_empirical_p": point_empirical_p,
        "target_maxT_empirical_p": max_t_empirical_p,
        "fwer05_nominal_threshold": max(attainable) if attainable else None,
        "observed_permutation_min_p": float(min_p[observed_index]),
        "min_p_distribution": dict(sorted(min_p_counts.items(), key=lambda kv: float(kv[0]))),
        "max_t_adjusted_pvalues": max_t_adjusted,
    }


def write_tsv(path, fieldnames, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    out = args.output
    scripts_dir = out / "scripts"
    results_dir = out / "results"
    summary_dir = out / "summary"
    logs_dir = out / "logs"
    for directory in (scripts_dir, results_dir, summary_dir, logs_dir):
        directory.mkdir(parents=True, exist_ok=True)

    model_state = {
        name: {
            "name": name,
            "min_called_total": min_total,
            "min_called_per_group": min_group,
            "pvalues": [],
            "complete": 0,
            "group80": 0,
            "target_p": None,
        }
        for name, min_total, min_group in MODELS
    }
    primary_rows = []
    primary_patterns = set()
    source_value_counts = Counter()
    total_input = 0
    selected_input = 0
    target_raw_calls = None
    target_corrected_calls = None
    target_meta = None

    with args.input.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing_columns = [s for s in SAMPLES if s not in (reader.fieldnames or [])]
        if missing_columns:
            raise SystemExit("Missing sample columns: " + ",".join(missing_columns))

        for row in reader:
            total_input += 1
            if row.get("CHROM") not in CHROMS:
                continue
            if row.get("final_confidence") != "read_svgap_supported":
                continue
            selected_input += 1

            raw_calls = {}
            for sample in SAMPLES:
                raw_value = str(row[sample]).strip()
                source_value_counts[raw_value] += 1
                raw_calls[sample] = parse_call(raw_value)

            corrected_calls = dict(raw_calls)
            is_target = row.get("FinalSV_ID") == TARGET_ID
            if is_target:
                for sample in MANUAL_ZERO_SAMPLES:
                    corrected_calls[sample] = 0
                target_raw_calls = dict(raw_calls)
                target_corrected_calls = dict(corrected_calls)
                target_meta = {
                    "FinalSV_ID": row.get("FinalSV_ID"),
                    "CHROM": row.get("CHROM"),
                    "POS": row.get("POS"),
                    "END": row.get("END"),
                    "SVTYPE": row.get("SVTYPE"),
                    "SVLEN": row.get("SVLEN"),
                }

            stats = test_calls(corrected_calls)
            for name, min_total, min_group in MODELS:
                if passes_model(stats, min_total, min_group):
                    state = model_state[name]
                    state["pvalues"].append(stats["fisher_p"])
                    state["complete"] += int(stats["complete_separation"])
                    state["group80"] += int(stats["group_specific80"])
                    if is_target:
                        state["target_p"] = stats["fisher_p"]

            if passes_model(stats, 12, 5):
                output_row = {
                    "FinalSV_ID": row.get("FinalSV_ID"),
                    "CHROM": row.get("CHROM"),
                    "POS": row.get("POS"),
                    "END": row.get("END"),
                    "SVTYPE": row.get("SVTYPE"),
                    "SVLEN": row.get("SVLEN"),
                    "read_SV_ID": row.get("read_SV_ID", "."),
                    "svgap_ids": row.get("svgap_ids", "."),
                    "final_confidence": row.get("final_confidence", "."),
                    **stats,
                    "target_Chr23997": is_target,
                    "manual_IGV_correction_applied": is_target,
                }
                primary_rows.append(output_row)
                primary_patterns.add(calls_to_masks(corrected_calls))

    if target_raw_calls is None or target_corrected_calls is None:
        raise SystemExit("Target %s not found" % TARGET_ID)

    sensitivity_rows = []
    for name, min_total, min_group in MODELS:
        state = model_state[name]
        pvalues = state["pvalues"]
        qvalues = bh_adjust(pvalues)
        target_p = state["target_p"]
        target_rank_min = 1 + sum(p < target_p - 1e-15 for p in pvalues) if target_p is not None else None
        target_rank_max = sum(p <= target_p + 1e-15 for p in pvalues) if target_p is not None else None
        target_q = None
        if target_p is not None:
            candidates = [q for p, q in zip(pvalues, qvalues) if abs(p - target_p) < 1e-15]
            target_q = min(candidates) if candidates else None
        sensitivity_rows.append({
            "model": name,
            "min_called_total": min_total,
            "min_called_per_group": min_group,
            "tested_sites": len(pvalues),
            "p_lt_0.05": sum(p < 0.05 for p in pvalues),
            "p_lt_0.01": sum(p < 0.01 for p in pvalues),
            "complete_separation_sites": state["complete"],
            "group_specific80_sites": state["group80"],
            "minimum_p": min(pvalues) if pvalues else None,
            "minimum_BH_q": min(qvalues) if qvalues else None,
            "BH_q_le_0.05": sum(q <= 0.05 for q in qvalues),
            "BH_q_le_0.10": sum(q <= 0.10 for q in qvalues),
            "bonferroni_alpha_0.05": 0.05 / len(pvalues) if pvalues else None,
            "target_p": target_p,
            "target_BH_q": target_q,
            "target_rank_min": target_rank_min,
            "target_rank_max": target_rank_max,
        })

    primary_p = [row["fisher_p"] for row in primary_rows]
    primary_q = bh_adjust(primary_p)
    primary_min_p = min(primary_p)
    permutation = run_permutation_max_t(
        sorted(primary_patterns),
        calls_to_masks(target_corrected_calls),
        primary_p,
    )
    primary_m = len(primary_rows)
    bonf_alpha = 0.05 / primary_m
    for row, qvalue, max_t_p in zip(primary_rows, primary_q, permutation["max_t_adjusted_pvalues"]):
        row["BH_q"] = qvalue
        row["bonferroni_adjusted_p"] = min(1.0, row["fisher_p"] * primary_m)
        row["maxT_empirical_adjusted_p"] = max_t_p

    primary_rows.sort(key=lambda row: (row["fisher_p"], -row["abs_frequency_difference"], row["CHROM"], int(row["POS"])))
    for rank, row in enumerate(primary_rows, start=1):
        row["rank_by_nominal_p"] = rank

    target_raw_stats = test_calls(target_raw_calls)
    target_corrected_stats = test_calls(target_corrected_calls)
    target_primary = next(row for row in primary_rows if row["FinalSV_ID"] == TARGET_ID)
    primary_sensitivity = next(row for row in sensitivity_rows if row["model"] == "primary_callrate_ge80")

    primary_fields = [
        "rank_by_nominal_p", "FinalSV_ID", "CHROM", "POS", "END", "SVTYPE", "SVLEN",
        "read_SV_ID", "svgap_ids", "final_confidence",
        "del_present", "del_absent", "del_called", "intact_present", "intact_absent", "intact_called",
        "called_total", "missing_total", "del_freq", "intact_freq", "abs_frequency_difference", "direction",
        "fisher_p", "neglog10p", "BH_q", "bonferroni_adjusted_p", "maxT_empirical_adjusted_p",
        "complete_separation", "group_specific80", "target_Chr23997", "manual_IGV_correction_applied",
    ]
    write_tsv(results_dir / "missing_aware_primary_PAV_trait_cosegregation.tsv", primary_fields, primary_rows)

    call_rows = []
    for sample in SAMPLES:
        call_rows.append({
            "sample": sample,
            "phenotype_group": "DEL_spineless" if sample in DEL_GROUP else "INTACT_spiny_like",
            "raw_call": "NA" if target_raw_calls[sample] is None else target_raw_calls[sample],
            "IGV_corrected_call": "NA" if target_corrected_calls[sample] is None else target_corrected_calls[sample],
            "changed": target_raw_calls[sample] != target_corrected_calls[sample],
            "note": "documented false-positive raw DEL call; manual IGV/depth supports intact" if sample in MANUAL_ZERO_SAMPLES else ".",
        })
    write_tsv(
        results_dir / "Chr23997_genotypes_raw_vs_IGV_corrected.tsv",
        ["sample", "phenotype_group", "raw_call", "IGV_corrected_call", "changed", "note"],
        call_rows,
    )

    association_rows = []
    for label, stats in (("raw_matrix", target_raw_stats), ("manual_IGV_corrected", target_corrected_stats)):
        association_rows.append({
            "analysis": label,
            **target_meta,
            **stats,
            "BH_q_primary_family": target_primary["BH_q"] if label == "manual_IGV_corrected" else "NA",
            "bonferroni_adjusted_p_primary_family": target_primary["bonferroni_adjusted_p"] if label == "manual_IGV_corrected" else "NA",
            "maxT_empirical_adjusted_p_primary_family": target_primary["maxT_empirical_adjusted_p"] if label == "manual_IGV_corrected" else "NA",
        })
    association_fields = [
        "analysis", "FinalSV_ID", "CHROM", "POS", "END", "SVTYPE", "SVLEN",
        "del_present", "del_absent", "del_called", "intact_present", "intact_absent", "intact_called",
        "called_total", "missing_total", "del_freq", "intact_freq", "abs_frequency_difference", "direction",
        "fisher_p", "neglog10p", "complete_separation", "group_specific80",
        "BH_q_primary_family", "bonferroni_adjusted_p_primary_family", "maxT_empirical_adjusted_p_primary_family",
    ]
    write_tsv(results_dir / "Chr23997_association_raw_vs_IGV_corrected.tsv", association_fields, association_rows)

    write_tsv(
        summary_dir / "sensitivity_analysis.tsv",
        list(sensitivity_rows[0].keys()),
        sensitivity_rows,
    )

    permutation_rows = [{"min_nominal_p": p, "permutation_count": n}
                        for p, n in permutation["min_p_distribution"].items()]
    write_tsv(summary_dir / "permutation_minP_distribution.tsv",
              ["min_nominal_p", "permutation_count"], permutation_rows)

    summary = {
        "analysis": "missing_aware_PAV_trait_cosegregation",
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "source_matrix": str(args.input),
        "source_sha256": sha256_file(args.input),
        "total_input_rows": total_input,
        "read_svgap_supported_Chr1_Chr8_rows": selected_input,
        "genotype_cells_in_selected_rows": selected_input * len(SAMPLES),
        "source_call_0": source_value_counts["0"] + source_value_counts["0.0"],
        "source_call_1": source_value_counts["1"] + source_value_counts["1.0"],
        "source_call_missing": selected_input * len(SAMPLES) - source_value_counts["0"] - source_value_counts["0.0"] - source_value_counts["1"] - source_value_counts["1.0"],
        "source_missing_fraction": (selected_input * len(SAMPLES) - source_value_counts["0"] - source_value_counts["0.0"] - source_value_counts["1"] - source_value_counts["1.0"]) / float(selected_input * len(SAMPLES)),
        "primary_filter": "called_total>=12/14 and called_per_group>=5",
        "primary_tested_sites": primary_m,
        "primary_unique_genotype_patterns": len(primary_patterns),
        "primary_minimum_p": primary_min_p,
        "primary_sites_at_minimum_p": sum(abs(p - primary_min_p) < 1e-15 for p in primary_p),
        "primary_minimum_BH_q": min(primary_q),
        "primary_BH_q_le_0.05": sum(q <= 0.05 for q in primary_q),
        "primary_BH_q_le_0.10": sum(q <= 0.10 for q in primary_q),
        "primary_bonferroni_alpha_0.05": bonf_alpha,
        "primary_bonferroni_neglog10_threshold": -math.log10(bonf_alpha),
        "primary_bonferroni_significant_sites": sum(p <= bonf_alpha for p in primary_p),
        "primary_maxT_FWER_significant_sites": sum(row["maxT_empirical_adjusted_p"] <= 0.05 for row in primary_rows),
        "permutation_count": permutation["permutations"],
        "permutation_target_pointwise_empirical_p": permutation["target_pointwise_empirical_p"],
        "permutation_target_maxT_empirical_p": permutation["target_maxT_empirical_p"],
        "permutation_FWER05_nominal_threshold": permutation["fwer05_nominal_threshold"],
        "Chr23997_raw_p": target_raw_stats["fisher_p"],
        "Chr23997_corrected_p": target_corrected_stats["fisher_p"],
        "Chr23997_corrected_BH_q": target_primary["BH_q"],
        "Chr23997_corrected_bonferroni_adjusted_p": target_primary["bonferroni_adjusted_p"],
        "Chr23997_corrected_maxT_empirical_adjusted_p": target_primary["maxT_empirical_adjusted_p"],
        "Chr23997_rank_by_nominal_p": target_primary["rank_by_nominal_p"],
        "Chr23997_tie_aware_rank_min": primary_sensitivity["target_rank_min"],
        "Chr23997_tie_aware_rank_max": primary_sensitivity["target_rank_max"],
        "Chr23997_manual_zero_samples": sorted(MANUAL_ZERO_SAMPLES),
        "warning_phylogeny": "Unrestricted species-label permutations assume exchangeability and are diagnostic, not a phylogeny-aware formal test.",
        "warning_posthoc": "The missing-aware filter and manual target correction are post hoc reanalysis decisions and require transparent reporting.",
    }
    with (summary_dir / "reanalysis_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    with (summary_dir / "reanalysis_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        for key, value in summary.items():
            if isinstance(value, list):
                value = ",".join(str(x) for x in value)
            writer.writerow([key, value if value is not None else "NA"])

    readme = """# Missing-aware PAV trait-cosegregation reanalysis

Primary analysis:
- Source: high-confidence read+SVGAP-supported PAVs on Chr1-Chr8.
- Phenotype contrast: 7 DEL/spineless versus 7 INTACT/spiny-like genomes.
- Calls: 1=present, 0=absent, NA/.=missing (never converted to absence).
- Call filter: at least 12/14 samples called and at least 5 called per phenotype group.
- Test: two-sided Fisher exact test on called samples.
- Corrections: BH-FDR, Bonferroni, and exhaustive 3,432-label maxT permutation.
- Chr23997: raw result is retained separately; the corrected family replaces the three documented false-positive calls in genome_436, genome_457 and genome_474 with intact (0).

Sensitivity analyses use complete cases (14/14) and a relaxed >=10/14 filter.

Important limits:
- This is a species-level comparative screen, not an individual-level GWAS.
- Unrestricted label permutation is not phylogeny-aware and is reported as a diagnostic sensitivity analysis.
- Filters and manual correction are post hoc and must be reported transparently.
"""
    (out / "README.md").write_text(readme, encoding="utf-8")
    (summary_dir / "COMPLETED").write_text(datetime.now().isoformat(timespec="seconds") + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
