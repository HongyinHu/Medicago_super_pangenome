#!/usr/bin/env python3
import csv
import re
import sys
from collections import Counter
from pathlib import Path


def classify_state(state):
    normalized = state.strip().lower()
    if normalized in {"spiral", "weak spiral"}:
        return 1
    if normalized == "non spiral":
        return 0
    raise ValueError(f"Unrecognized pod spiral state: {state!r}")


def extract_average_coil(raw):
    match = re.search(r"average_coil=~?([0-9]+(?:\.[0-9]+)?)", raw)
    if not match:
        raise ValueError(f"Cannot parse average_coil from: {raw!r}")
    return float(match.group(1))


def load_table(path):
    with open(path, encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise RuntimeError(f"No rows found in {path}")
    return rows


def prepare(master_path, audit_path, outdir):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    master_rows = load_table(master_path)
    master = {row["Sample_ID"]: row for row in master_rows}
    if len(master) != len(master_rows):
        raise RuntimeError("Duplicate Sample_ID values in master phenotype table")

    audit_rows = load_table(audit_path)
    sample_ids = [row["Sample_ID"] for row in audit_rows]
    if len(sample_ids) != len(set(sample_ids)):
        raise RuntimeError("Duplicate Sample_ID values in 144-sample audit table")

    missing = [sample_id for sample_id in sample_ids if sample_id not in master]
    if missing:
        raise RuntimeError(f"Samples missing from master phenotype table: {missing}")

    selected = []
    for sample_id in sample_ids:
        row = master[sample_id]
        average_coil = extract_average_coil(row["Pod_spiral_raw"])
        binary = int(average_coil >= 1.0)
        categorical_binary = classify_state(row["Pod_spiral"])
        if binary != categorical_binary:
            raise RuntimeError(
                f"Numeric/categorical phenotype mismatch for {sample_id}: "
                f"average_coil={average_coil}, Pod_spiral={row['Pod_spiral']!r}"
            )
        selected.append((sample_id, binary, average_coil, row))

    counts = Counter(binary for _, binary, _, _ in selected)
    if len(selected) != 144 or counts != Counter({1: 103, 0: 41}):
        raise RuntimeError(
            f"Unexpected coil>=1 composition: n={len(selected)}, counts={dict(counts)}"
        )

    with open(out / "samples_144.txt", "w", encoding="utf-8") as handle:
        for sample_id, _, _, _ in selected:
            handle.write(sample_id + "\n")

    with open(out / "samples_144.keep", "w", encoding="utf-8") as handle:
        for sample_id, _, _, _ in selected:
            handle.write(f"{sample_id}\t{sample_id}\n")

    with open(out / "coil_ge1.plink.pheno", "w", encoding="utf-8") as handle:
        for sample_id, binary, _, _ in selected:
            handle.write(f"{sample_id}\t{sample_id}\t{binary + 1}\n")

    fieldnames = [
        "FID",
        "IID",
        "coil_ge1_binary_spiral1_nonspiral0",
        "PLINK_case_control_spiral2_nonspiral1",
        "average_coil",
        "Pod_spiral",
        "Pod_spiral_raw",
        "Latin_name",
        "Section",
        "Trait_priority",
        "Pod_trait_source",
    ]
    with open(out / "coil_ge1_phenotype.tsv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for sample_id, binary, average_coil, row in selected:
            writer.writerow(
                {
                    "FID": sample_id,
                    "IID": sample_id,
                    "coil_ge1_binary_spiral1_nonspiral0": binary,
                    "PLINK_case_control_spiral2_nonspiral1": binary + 1,
                    "average_coil": average_coil,
                    "Pod_spiral": row["Pod_spiral"],
                    "Pod_spiral_raw": row.get("Pod_spiral_raw", ""),
                    "Latin_name": row.get("Latin_name", ""),
                    "Section": row.get("Section", ""),
                    "Trait_priority": row.get("Trait_priority", ""),
                    "Pod_trait_source": row.get("Pod_trait_source", ""),
                }
            )

    label_counts = Counter(
        row["Pod_spiral"].strip().lower() for _, _, _, row in selected
    )
    with open(out / "coil_ge1_sample_counts.tsv", "w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"input_samples\t{len(selected)}\n")
        handle.write(f"spiral_case_1\t{counts[1]}\n")
        handle.write(f"nonspiral_control_0\t{counts[0]}\n")
        handle.write(f"original_spiral\t{label_counts['spiral']}\n")
        handle.write(f"original_weak_spiral\t{label_counts['weak spiral']}\n")
        handle.write(f"original_non_spiral\t{label_counts['non spiral']}\n")
        handle.write("threshold\taverage_coil>=1 is spiral\n")

    return len(selected), counts[1], counts[0]


def main():
    if len(sys.argv) != 4:
        raise SystemExit(
            "Usage: prepare_coil_ge1.py master_phenotype.tsv phenotype_144.tsv outdir"
        )
    n, cases, controls = prepare(sys.argv[1], sys.argv[2], sys.argv[3])
    print(f"coil>=1 phenotype prepared: n={n}, spiral={cases}, non-spiral={controls}")


if __name__ == "__main__":
    main()
