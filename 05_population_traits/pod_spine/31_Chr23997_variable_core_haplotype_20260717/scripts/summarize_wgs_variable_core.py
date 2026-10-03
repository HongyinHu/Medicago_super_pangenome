#!/usr/bin/env python3
"""Summarize direct read-based variable core-loss calls across depth cutoffs."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def count_core_states(rows, min_flank_depth):
    """Count high-confidence core loss/intact calls; phenotype 1 is spiny."""
    counts = {"spiny_loss": 0, "spiny_intact": 0, "spineless_loss": 0, "spineless_intact": 0}
    for row in rows:
        if row.get("evidence_class") != "HIGH_CONFIDENCE" or row.get("GT") not in {"0/0", "0/1", "1/1"}:
            continue
        if float(row["flank_mean_depth"]) < min_flank_depth:
            continue
        prefix = "spiny" if row["phenotype"] == "1" else "spineless"
        suffix = "loss" if row["GT"] != "0/0" else "intact"
        counts[f"{prefix}_{suffix}"] += 1
    return counts


def fisher_two_sided(a, b, c, d):
    row1, col1, total = a + b, a + c, a + b + c + d

    def prob(x):
        return math.comb(col1, x) * math.comb(total - col1, row1 - x) / math.comb(total, row1)

    observed = prob(a)
    lower = max(0, row1 - (total - col1))
    upper = min(row1, col1)
    return sum(prob(x) for x in range(lower, upper + 1) if prob(x) <= observed + 1e-12)


def odds_ratio(a, b, c, d):
    if b * c == 0:
        return "inf" if a * d else "NA"
    return f"{a * d / (b * c):.6g}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calls", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--depth-cutoffs", default="8,10,15,20")
    args = parser.parse_args()

    with open(args.calls, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    output = []
    for cutoff in [float(value) for value in args.depth_cutoffs.split(",")]:
        counts = count_core_states(rows, cutoff)
        a, b = counts["spiny_loss"], counts["spiny_intact"]
        c, d = counts["spineless_loss"], counts["spineless_intact"]
        output.append({
            "min_flank_depth": f"{cutoff:g}",
            "n_called": a + b + c + d,
            **counts,
            "odds_ratio_spiny_core_loss": odds_ratio(a, b, c, d),
            "fisher_two_sided_p": f"{fisher_two_sided(a, b, c, d):.8g}",
            "definition": "GT from split-read/core-coverage core-loss caller; exact DEL length is not required.",
        })
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output)


if __name__ == "__main__":
    main()
