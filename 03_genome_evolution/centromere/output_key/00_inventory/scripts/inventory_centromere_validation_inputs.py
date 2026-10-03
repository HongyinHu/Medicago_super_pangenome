#!/usr/bin/env python3
"""Inventory inputs for multi-species Mpo centromere-fate validation."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_fai(path: Path) -> list[tuple[str, int]]:
    if not path.exists():
        return []
    rows: list[tuple[str, int]] = []
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            chrom, size, *_ = line.rstrip("\n").split("\t")
            rows.append((chrom, int(size)))
    return rows


def count_primary_chroms(fai_rows: list[tuple[str, int]]) -> int:
    chr_rows = [x for x in fai_rows if x[0].startswith("Chr")]
    return len(chr_rows) if chr_rows else len(fai_rows)


def find_species_files(root: Path, species: str, keywords: list[str]) -> list[str]:
    out: list[str] = []
    token = f"genome_{species}"
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        path_str = str(path)
        name = path.name
        if token not in path_str:
            continue
        if all(keyword.lower() in name.lower() for keyword in keywords):
            out.append(str(path.relative_to(root)))
    return sorted(out)


def find_cen_candidates(root: Path, species: str) -> list[str]:
    token = f"genome_{species}"
    patterns = [
        "cen_core.bed",
        "functional_centromere.final.bed",
        "primary_functional_centromere.bed",
        "CENH3.domain.bed",
        "CENH3peak",
        "CENH3.cluster",
        "CENH3_signal.cluster.filtered.bed",
        "CENH3.bigwig.top4.clusters.bed",
    ]
    out: list[str] = []
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        path_str = str(path)
        if token not in path_str:
            continue
        name = path.name
        if any(pattern in name for pattern in patterns):
            out.append(str(path.relative_to(root)))
    return sorted(out)


def directory_exists(root: Path, names: list[str]) -> list[str]:
    out: list[str] = []
    for path in root.rglob("*"):
        if not path.is_dir():
            continue
        path_str = str(path)
        if all(name in path_str for name in names):
            out.append(str(path.relative_to(root)))
    return sorted(out)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default="path/to/project/10.centromere_analysis")
    args = parser.parse_args()

    root = Path(args.project_root)
    out_dir = root / "output_key/00_inventory"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "scripts").mkdir(parents=True, exist_ok=True)

    genome_dir = root / "output_all/00_genome"
    gff_dir = root / "output_all/01_gff"

    species: list[str] = []
    for fa in sorted(genome_dir.glob("genome_*.fa")):
        species.append(fa.stem.replace("genome_", ""))

    rows: list[dict[str, object]] = []
    for sp in species:
        fa = genome_dir / f"genome_{sp}.fa"
        fai = genome_dir / f"genome_{sp}.fa.fai"
        gff = gff_dir / f"genome_{sp}.gff3"
        fai_rows = read_fai(fai)
        cen_files = find_cen_candidates(root, sp)
        mpo_jcvi_dirs = []
        mpo_cent_dirs = []
        if sp != "Mpo":
            mpo_jcvi_dirs = directory_exists(root / "output_synthen/01_JCVI", [sp, "Mpo"])
            mpo_cent_dirs = directory_exists(root / "output_synthen/02_cent_synteny", [sp, "Mpo"])
        row = {
            "species": sp,
            "genome_fa": str(fa.relative_to(root)) if fa.exists() else "",
            "genome_fai": str(fai.relative_to(root)) if fai.exists() else "",
            "gff3": str(gff.relative_to(root)) if gff.exists() else "",
            "primary_chromosome_count": count_primary_chroms(fai_rows),
            "total_sequences_in_fai": len(fai_rows),
            "cen_candidate_file_count": len(cen_files),
            "cen_candidate_files": ";".join(cen_files),
            "has_pairwise_JCVI_with_Mpo": "yes" if mpo_jcvi_dirs or sp == "Mpo" else "no",
            "pairwise_JCVI_dirs_with_Mpo": ";".join(mpo_jcvi_dirs),
            "has_cent_synteny_with_Mpo": "yes" if mpo_cent_dirs or sp == "Mpo" else "no",
            "cent_synteny_dirs_with_Mpo": ";".join(mpo_cent_dirs),
        }
        if sp == "Mpo":
            row["role"] = "x7_target"
        elif row["primary_chromosome_count"] == 8:
            row["role"] = "x8_auxiliary_reference"
        else:
            row["role"] = "other_or_unresolved"
        rows.append(row)

    out_tsv = out_dir / "centromere_validation_input_inventory.tsv"
    with out_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    summary = out_dir / "centromere_validation_input_inventory.summary.md"
    with summary.open("w") as handle:
        handle.write("# Centromere Validation Input Inventory\n\n")
        handle.write(f"- Project root: `{root}`\n")
        handle.write(f"- Species detected: {', '.join(species)}\n")
        handle.write("- x7 target: " + ", ".join(row["species"] for row in rows if row["role"] == "x7_target") + "\n")
        handle.write("- x8 auxiliary references: " + ", ".join(row["species"] for row in rows if row["role"] == "x8_auxiliary_reference") + "\n\n")
        handle.write("## Run Readiness\n\n")
        for row in rows:
            missing = []
            if not row["genome_fai"]:
                missing.append("FAI")
            if not row["gff3"]:
                missing.append("GFF3")
            if int(row["cen_candidate_file_count"]) == 0:
                missing.append("CEN candidates")
            if row["species"] != "Mpo" and row["has_pairwise_JCVI_with_Mpo"] != "yes":
                missing.append("pairwise JCVI with Mpo")
            status = "ready-ish" if not missing else "incomplete: " + ", ".join(missing)
            handle.write(f"- {row['species']}: {status}\n")

    print(out_tsv)
    print(summary)


if __name__ == "__main__":
    main()
