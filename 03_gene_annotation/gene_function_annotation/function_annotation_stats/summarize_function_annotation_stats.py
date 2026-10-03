#!/usr/bin/env python3
"""Summarize functional annotation gene-list coverage for multiple genomes."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


DEFAULT_SAMPLES = [
    "genome_474",
    "genome_A17",
    "genome_Mpo",
    "genome_Msa",
    "genome_R108",
]

ANNOTATIONS = [
    ("Pfam", "pfam_anno", "{sample}.pfam.genelist"),
    ("Interproscan", "InterproScan_anno", "{sample}.interproscan.genelist"),
    ("KEGG", "KEGG_anno", "{sample}.KEGG.genelist"),
    ("NR", "NR_anno", "{sample}.NR.genelist"),
    ("Swiss-Prot", "SwissPort_anno", "{sample}.swissport.genelist"),
    ("KOG", "KOG_anno", "{sample}.KOG.genelist"),
    ("eggNOG", "eggNOG_anno", "{sample}.eggNOG.genelist"),
    ("TrEMBL", "TrEMBL_anno", "{sample}.TrEMBL.genelist"),
    ("GO", "GO_anno", "{sample}.GO.genelist"),
]

SPECIES_LABELS = {
    "genome_474": "Medicago carstiensis",
    "genome_A17": "genome_A17",
    "genome_Mpo": "genome_Mpo",
    "genome_Msa": "genome_Msa",
    "genome_R108": "genome_R108",
}


def read_gene_set(path: Path) -> set[str]:
    genes: set[str] = set()
    if not path.exists():
        return genes
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            genes.add(line.split()[0])
    return genes


def count_fasta_records(path: Path) -> int:
    count = 0
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def parse_samples(text: str | None) -> list[str]:
    if not text:
        return DEFAULT_SAMPLES
    return [item.strip() for item in text.split(",") if item.strip()]


def discover_samples(output_dir: Path, requested_samples: list[str]) -> list[str]:
    if requested_samples:
        return requested_samples

    samples: list[str] = []
    interpro_dir = output_dir / "InterproScan_anno"
    if not interpro_dir.exists():
        return samples
    for sample_dir in sorted(interpro_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        pep = sample_dir / f"{sample_dir.name}.pep.format"
        if pep.exists():
            samples.append(sample_dir.name)
    return samples


def collect(output_dir: Path, samples: list[str]) -> tuple[list[dict[str, object]], list[dict[str, str]]]:
    rows: list[dict[str, object]] = []
    missing: list[dict[str, str]] = []
    for sample in discover_samples(output_dir, samples):
        total_file = output_dir / "InterproScan_anno" / sample / f"{sample}.pep.format"
        if total_file.exists():
            total = count_fasta_records(total_file)
        else:
            total = 0
            missing.append({"sample_id": sample, "annotation": "Total_gene", "path": str(total_file)})

        row: dict[str, object] = {
            "sample_id": sample,
            "Species name": SPECIES_LABELS.get(sample, sample),
            "Total_gene Gene number": total,
            "Total_gene Percent (%)": 100.0,
            "total_gene_file": str(total_file),
        }

        annotated_genes: set[str] = set()
        for label, folder, pattern in ANNOTATIONS:
            gene_file = output_dir / folder / pattern.format(sample=sample)
            genes = read_gene_set(gene_file)
            if not gene_file.exists():
                missing.append({"sample_id": sample, "annotation": label, "path": str(gene_file)})
            annotated_genes.update(genes)
            count = len(genes)
            row[f"{label} Gene number"] = count
            row[f"{label} Percent (%)"] = count / total * 100 if total else 0.0
            row[f"{label} file"] = str(gene_file)

        annotated_file = output_dir / f"{sample}.all.anno.genelist"
        if annotated_file.exists():
            annotated_genes = read_gene_set(annotated_file)
        annotated = len(annotated_genes)
        not_annotated = max(total - annotated, 0)
        row["Not annotated Gene number"] = not_annotated
        row["Not annotated Percent (%)"] = not_annotated / total * 100 if total else 0.0
        row["all_annotated_file"] = str(annotated_file)
        rows.append(row)
    return rows, missing


def table_columns() -> list[str]:
    columns = ["Species name"]
    for label in ["Total_gene", *[item[0] for item in ANNOTATIONS], "Not annotated"]:
        columns += [f"{label} Gene number", f"{label} Percent (%)"]
    return columns


def write_table(rows: list[dict[str, object]], path: Path) -> None:
    columns = table_columns()
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            out: dict[str, object] = {}
            for column in columns:
                value = row[column]
                if column.endswith("Percent (%)"):
                    out[column] = f"{float(value):.2f}"
                else:
                    out[column] = value
            writer.writerow(out)


def write_detail(rows: list[dict[str, object]], path: Path) -> None:
    columns = ["sample_id", *table_columns(), "total_gene_file", "all_annotated_file"]
    for label, _folder, _pattern in ANNOTATIONS:
        columns.append(f"{label} file")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            out = row.copy()
            for key, value in list(out.items()):
                if key.endswith("Percent (%)"):
                    out[key] = f"{float(value):.4f}"
            writer.writerow({column: out.get(column, "") for column in columns})


def write_missing(missing: list[dict[str, str]], path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id", "annotation", "path"])
        writer.writeheader()
        writer.writerows(missing)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--table", type=Path, required=True)
    parser.add_argument("--detail", type=Path, required=True)
    parser.add_argument(
        "--missing-report",
        type=Path,
        help="Optional CSV report of missing protein/genelist inputs.",
    )
    parser.add_argument(
        "--samples",
        help="Comma-separated sample IDs. Defaults to the five target genomes.",
    )
    parser.add_argument(
        "--require-complete",
        action="store_true",
        help="Exit non-zero if any expected input file is missing.",
    )
    args = parser.parse_args()

    rows, missing = collect(args.output_dir, parse_samples(args.samples))
    write_table(rows, args.table)
    write_detail(rows, args.detail)
    if args.missing_report:
        write_missing(missing, args.missing_report)
    if args.require_complete and missing:
        for item in missing:
            print(
                f"missing {item['sample_id']} {item['annotation']}: {item['path']}",
                file=sys.stderr,
            )
        raise SystemExit(1)


if __name__ == "__main__":
    main()
