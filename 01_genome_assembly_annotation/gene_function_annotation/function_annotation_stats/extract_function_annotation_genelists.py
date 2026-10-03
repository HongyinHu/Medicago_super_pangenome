#!/usr/bin/env python3
"""Build gene-list files and final coverage tables from functional annotation outputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


SAMPLES = ["genome_474", "genome_A17", "genome_Mpo", "genome_Msa", "genome_R108"]

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


def read_genes(path: Path) -> set[str]:
    genes: set[str] = set()
    if not path.exists():
        return genes
    with path.open() as handle:
        for line in handle:
            if line.strip():
                genes.add(line.split()[0])
    return genes


def write_genes(path: Path, genes: set[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for gene in sorted(genes):
            print(gene, file=handle)


def fasta_count(path: Path) -> int:
    with path.open() as handle:
        return sum(1 for line in handle if line.startswith(">"))


def genes_from_first_column(path: Path, skip_hash: bool = False) -> set[str]:
    genes: set[str] = set()
    if not path.exists():
        return genes
    with path.open() as handle:
        for line in handle:
            if skip_hash and line.startswith("#"):
                continue
            fields = line.strip().split()
            if fields:
                genes.add(fields[0])
    return genes


def parse_emapper(sample: str, output_dir: Path) -> tuple[set[str], set[str], set[str]]:
    path = output_dir / "eggNOG_anno" / sample / f"{sample}.emapper.annotations"
    eggnog_genes: set[str] = set()
    kegg_genes: set[str] = set()
    go_genes: set[str] = set()
    kegg_rows: set[tuple[str, str]] = set()
    go_rows: set[tuple[str, str]] = set()
    if not path.exists():
        return eggnog_genes, kegg_genes, go_genes

    header: dict[str, int] = {}
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                names = line.lstrip("#").split("\t")
                header = {name: index for index, name in enumerate(names)}
                continue

            fields = line.split("\t")
            if not fields:
                continue
            query = fields[0]
            seed = get_field(fields, header, "seed_ortholog")
            if seed and seed != "-":
                eggnog_genes.add(query)

            for ko in split_values(get_field(fields, header, "KEGG_ko")):
                ko = ko.replace("ko:", "")
                if ko.startswith("K"):
                    kegg_genes.add(query)
                    kegg_rows.add((query, ko))

            for go in split_values(get_field(fields, header, "GOs")):
                if go.startswith("GO:"):
                    go_genes.add(query)
                    go_rows.add((query, go))

    write_pairs(output_dir / "KEGG_anno" / f"{sample}.emapper.annotations.tsv.KEGG_Knum.txt", kegg_rows)
    write_pairs(output_dir / "GO_anno" / f"{sample}.emapper.annotations.tsv.GO.txt", go_rows)
    return eggnog_genes, kegg_genes, go_genes


def get_field(fields: list[str], header: dict[str, int], name: str) -> str:
    index = header.get(name)
    if index is None or index >= len(fields):
        return ""
    return fields[index]


def split_values(value: str) -> list[str]:
    if not value or value == "-":
        return []
    return [item.strip() for item in value.replace(";", ",").split(",") if item.strip()]


def write_pairs(path: Path, rows: set[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for gene, value in sorted(rows):
            print(f"{gene}\t{value}", file=handle)


def build_genelists(output_dir: Path, samples: list[str]) -> None:
    for sample in samples:
        write_genes(
            output_dir / "pfam_anno" / sample / f"{sample}.pfam.genelist",
            genes_from_first_column(output_dir / "pfam_anno" / sample / f"{sample}.pep.format.pfam.out", skip_hash=True),
        )
        write_genes(
            output_dir / "InterproScan_anno" / sample / f"{sample}.interproscan.genelist",
            genes_from_first_column(
                output_dir / "InterproScan_anno" / sample / f"{sample}.pep.format.interproscan.out"
            ),
        )
        for folder, suffix in [
            ("NR_anno", "NR"),
            ("SwissPort_anno", "swissport"),
            ("KOG_anno", "KOG"),
            ("TrEMBL_anno", "TrEMBL"),
        ]:
            write_genes(
                output_dir / folder / sample / f"{sample}.{suffix}.genelist",
                genes_from_first_column(output_dir / folder / sample / f"{sample}.pep.format_blastp.tab"),
            )

        eggnog_genes, kegg_genes, go_genes = parse_emapper(sample, output_dir)
        write_genes(output_dir / "eggNOG_anno" / sample / f"{sample}.eggNOG.genelist", eggnog_genes)
        write_genes(output_dir / "KEGG_anno" / f"{sample}.KEGG.genelist", kegg_genes)
        write_genes(output_dir / "GO_anno" / f"{sample}.GO.genelist", go_genes)

        all_genes: set[str] = set()
        for _label, folder, pattern in ANNOTATIONS:
            path = output_dir / folder / pattern.format(sample=sample)
            if not path.exists():
                path = output_dir / folder / sample / pattern.format(sample=sample)
            all_genes.update(read_genes(path))
        write_genes(output_dir / f"{sample}.all.anno.genelist", all_genes)


def collect_rows(output_dir: Path, samples: list[str]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for sample in samples:
        total_file = output_dir / "InterproScan_anno" / sample / f"{sample}.pep.format"
        total = fasta_count(total_file)
        row: dict[str, object] = {
            "Species name": sample,
            "Total_gene Gene number": total,
            "Total_gene Percent (%)": 100.0,
        }
        for label, folder, pattern in ANNOTATIONS:
            path = output_dir / folder / pattern.format(sample=sample)
            if not path.exists():
                path = output_dir / folder / sample / pattern.format(sample=sample)
            count = len(read_genes(path))
            row[f"{label} Gene number"] = count
            row[f"{label} Percent (%)"] = count / total * 100 if total else 0.0
        annotated = len(read_genes(output_dir / f"{sample}.all.anno.genelist"))
        not_annotated = max(total - annotated, 0)
        row["Not annotated Gene number"] = not_annotated
        row["Not annotated Percent (%)"] = not_annotated / total * 100 if total else 0.0
        rows.append(row)
    return rows


def write_table(rows: list[dict[str, object]], path: Path) -> None:
    columns = ["Species name"]
    for label in ["Total_gene", *[item[0] for item in ANNOTATIONS], "Not annotated"]:
        columns.extend([f"{label} Gene number", f"{label} Percent (%)"])
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            out = {}
            for column in columns:
                value = row[column]
                out[column] = f"{float(value):.2f}" if column.endswith("Percent (%)") else value
            writer.writerow(out)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--samples", default=",".join(SAMPLES))
    parser.add_argument("--table", type=Path, required=True)
    args = parser.parse_args()

    samples = [item.strip() for item in args.samples.split(",") if item.strip()]
    build_genelists(args.output_dir, samples)
    write_table(collect_rows(args.output_dir, samples), args.table)


if __name__ == "__main__":
    main()
