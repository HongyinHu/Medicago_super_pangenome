#!/usr/bin/env python3
"""Project active CENH3 cores from x=8 species to R108 ancestral CEN1-CEN8."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_bed(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open() as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.strip().split()
            if len(fields) < 3:
                continue
            chrom, start, end = fields[:3]
            rows.append(
                {
                    "chr": chrom,
                    "start": int(float(start)),
                    "end": int(float(end)),
                    "name": fields[3] if len(fields) > 3 else "",
                }
            )
    return rows


def read_blocks(path: Path, min_len: int) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            q_len = int(float(row.get("query_len", 0)))
            if q_len < min_len:
                continue
            rows.append(
                {
                    "query_chr": row["query_chr"],
                    "query_start": int(float(row["query_start"])),
                    "query_end": int(float(row["query_end"])),
                    "target_chr": row["target_chr"],
                    "target_start": int(float(row["target_start"])),
                    "target_end": int(float(row["target_end"])),
                    "query_len": q_len,
                    "target_len": int(float(row.get("target_len", q_len))),
                }
            )
    return rows


def ov(a: int, b: int, c: int, d: int) -> int:
    return max(0, min(b, d) - max(a, c))


def cen_id_from_chr(chrom: str) -> str:
    return "CEN" + chrom.replace("Chr", "")


def nearest_cen(chrom: str, start: int, end: int, cens: list[dict[str, object]]) -> tuple[int, int, dict[str, object]] | None:
    best = None
    for cen in cens:
        if cen["chr"] != chrom:
            continue
        overlap = ov(start, end, int(cen["start"]), int(cen["end"]))
        if overlap:
            dist = 0
        elif end < int(cen["start"]):
            dist = int(cen["start"]) - end
        else:
            dist = start - int(cen["end"])
        item = (dist, -overlap, cen)
        if best is None or item < best:
            best = item
    return best


def project_source_to_r108(
    species: str,
    source_cens: list[dict[str, object]],
    r108_cens: list[dict[str, object]],
    blocks: list[dict[str, object]],
    direction: str,
    flank_bp: int,
    distance_cutoff: int,
    bridge: str,
) -> list[dict[str, object]]:
    if direction not in {"query_to_target", "target_to_query"}:
        raise ValueError(direction)
    rows: list[dict[str, object]] = []
    for cen in sorted(source_cens, key=lambda x: int(str(x["chr"]).replace("Chr", ""))):
        source_cen_id = cen_id_from_chr(str(cen["chr"]))
        source_start = max(0, int(cen["start"]) - flank_bp)
        source_end = int(cen["end"]) + flank_bp
        if direction == "query_to_target":
            selected = [
                b
                for b in blocks
                if b["query_chr"] == cen["chr"] and ov(source_start, source_end, int(b["query_start"]), int(b["query_end"])) > 0
            ]
            group_chr = "target_chr"
            proj_start = "target_start"
            proj_end = "target_end"
            weight = "query_len"
        else:
            selected = [
                b
                for b in blocks
                if b["target_chr"] == cen["chr"] and ov(source_start, source_end, int(b["target_start"]), int(b["target_end"])) > 0
            ]
            group_chr = "query_chr"
            proj_start = "query_start"
            proj_end = "query_end"
            weight = "query_len"

        by_chr: dict[str, list[dict[str, object]]] = defaultdict(list)
        for block in selected:
            by_chr[str(block[group_chr])].append(block)
        ranked = sorted(by_chr.items(), key=lambda kv: sum(int(x[weight]) for x in kv[1]), reverse=True)

        row: dict[str, object] = {
            "species": species,
            "bridge": bridge,
            "source_cen_id_by_chr": source_cen_id,
            "source_chr": cen["chr"],
            "source_start": cen["start"],
            "source_end": cen["end"],
            "source_query_window_start": source_start,
            "source_query_window_end": source_end,
            "flank_bp": flank_bp,
            "block_count": len(selected),
            "R108_candidate_chr_summary": ";".join(
                f"{chrom}:{sum(int(x[weight]) for x in xs)}bp/{len(xs)}blocks" for chrom, xs in ranked[:6]
            ),
            "best_R108_chr": "",
            "best_R108_projected_start": "",
            "best_R108_projected_end": "",
            "nearest_R108_CEN": "",
            "nearest_R108_CEN_id": "",
            "distance_to_nearest_R108_CEN_bp": "",
            "classification": "unresolved",
        }
        if ranked:
            best_chr, best_blocks = ranked[0]
            coords: list[int] = []
            for block in best_blocks:
                coords.extend([int(block[proj_start]), int(block[proj_end])])
            start, end = min(coords), max(coords)
            row.update(
                {
                    "best_R108_chr": best_chr,
                    "best_R108_projected_start": start,
                    "best_R108_projected_end": end,
                }
            )
            near = nearest_cen(best_chr, start, end, r108_cens)
            if near:
                dist, neg_overlap, rcen = near
                rcen_id = cen_id_from_chr(str(rcen["chr"]))
                row.update(
                    {
                        "nearest_R108_CEN": f"{rcen_id}|{rcen['chr']}:{rcen['start']}-{rcen['end']}",
                        "nearest_R108_CEN_id": rcen_id,
                        "distance_to_nearest_R108_CEN_bp": dist,
                    }
                )
                if dist <= distance_cutoff:
                    row["classification"] = "active_core_supports_R108_ancestral_CEN"
                else:
                    row["classification"] = "active_core_projects_to_R108_same_chr_distal"
        rows.append(row)
    return rows


def self_r108_rows(r108_cens: list[dict[str, object]], flank_bp: int) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for cen in sorted(r108_cens, key=lambda x: int(str(x["chr"]).replace("Chr", ""))):
        cid = cen_id_from_chr(str(cen["chr"]))
        rows.append(
            {
                "species": "R108",
                "bridge": "self",
                "source_cen_id_by_chr": cid,
                "source_chr": cen["chr"],
                "source_start": cen["start"],
                "source_end": cen["end"],
                "source_query_window_start": max(0, int(cen["start"]) - flank_bp),
                "source_query_window_end": int(cen["end"]) + flank_bp,
                "flank_bp": flank_bp,
                "block_count": "",
                "R108_candidate_chr_summary": f"{cen['chr']}:{int(cen['end']) - int(cen['start'])}bp/self",
                "best_R108_chr": cen["chr"],
                "best_R108_projected_start": cen["start"],
                "best_R108_projected_end": cen["end"],
                "nearest_R108_CEN": f"{cid}|{cen['chr']}:{cen['start']}-{cen['end']}",
                "nearest_R108_CEN_id": cid,
                "distance_to_nearest_R108_CEN_bp": 0,
                "classification": "self_active_R108_CEN",
            }
        )
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "species",
        "bridge",
        "source_cen_id_by_chr",
        "source_chr",
        "source_start",
        "source_end",
        "source_query_window_start",
        "source_query_window_end",
        "flank_bp",
        "block_count",
        "R108_candidate_chr_summary",
        "best_R108_chr",
        "best_R108_projected_start",
        "best_R108_projected_end",
        "nearest_R108_CEN",
        "nearest_R108_CEN_id",
        "distance_to_nearest_R108_CEN_bp",
        "classification",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def summarize_presence(rows: list[dict[str, object]], out_path: Path) -> None:
    by_species: dict[str, set[str]] = defaultdict(set)
    by_species_detail: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        cid = str(row.get("nearest_R108_CEN_id", ""))
        if cid and str(row["classification"]) in {
            "active_core_supports_R108_ancestral_CEN",
            "self_active_R108_CEN",
        }:
            by_species[str(row["species"])].add(cid)
            by_species_detail[str(row["species"])].append(f"{row['source_cen_id_by_chr']}->{cid}")
    summary_rows: list[dict[str, object]] = []
    all_cens = [f"CEN{i}" for i in range(1, 9)]
    for species in sorted({str(row["species"]) for row in rows}):
        present = by_species.get(species, set())
        summary_rows.append(
            {
                "species": species,
                "supported_R108_CEN_count": len(present),
                "supported_R108_CENs": ",".join(c for c in all_cens if c in present),
                "CEN5_supported": "yes" if "CEN5" in present else "no",
                "missing_R108_CENs": ",".join(c for c in all_cens if c not in present),
                "detail": ";".join(by_species_detail.get(species, [])),
            }
        )
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["species", "supported_R108_CEN_count", "supported_R108_CENs", "CEN5_supported", "missing_R108_CENs", "detail"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(summary_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default="path/to/project/10.centromere_analysis")
    parser.add_argument("--flank-bp", type=int, default=3_000_000)
    parser.add_argument("--min-block-len", type=int, default=1000)
    parser.add_argument("--distance-cutoff", type=int, default=3_000_000)
    args = parser.parse_args()

    root = Path(args.project_root)
    out_dir = root / "output_key/03_project_other_x8_CENs_to_Mpo"
    suffix = "core_only" if args.flank_bp == 0 else f"pm{args.flank_bp // 1_000_000}Mb"

    r108_cens = read_bed(root / "output_synthen/01_JCVI/genome_A17_genome_R108/genome_R108.cen_core.bed")
    rows: list[dict[str, object]] = []
    rows.extend(self_r108_rows(r108_cens, args.flank_bp))

    jobs = [
        {
            "species": "A17",
            "source_cen": root / "output_synthen/01_JCVI/genome_A17_genome_R108/genome_A17.cen_core.bed",
            "blocks": root / "output_synthen/02_cent_synteny/genome_A17_genome_R108/genome_A17_vs_genome_R108.rawPAF.blocks.tsv",
            "direction": "query_to_target",
            "bridge": "A17_to_R108_direct",
        },
        {
            "species": "Msa",
            "source_cen": root / "output_synthen/01_JCVI/R108_Msa_synteny_lz10_run/genome_Msa.cen_core.bed",
            "blocks": root / "output_synthen/02_cent_synteny/genome_R108_genome_Msa/genome_R108_vs_genome_Msa.rawPAF.blocks.tsv",
            "direction": "target_to_query",
            "bridge": "Msa_to_R108_reverse_rawPAF",
        },
        {
            "species": "474",
            "source_cen": root / "output_synthen/01_JCVI/genome_R108_genome_474_synteny_lz10_run/genome_474.cen_core.bed",
            "blocks": root / "output_synthen/02_cent_synteny/genome_R108_genome_474/genome_R108_vs_genome_474.rawPAF.blocks.tsv",
            "direction": "target_to_query",
            "bridge": "474_to_R108_reverse_rawPAF",
        },
    ]
    for job in jobs:
        missing = [str(job[key]) for key in ("source_cen", "blocks") if not Path(job[key]).exists()]
        if missing:
            for path in missing:
                rows.append(
                    {
                        "species": job["species"],
                        "bridge": job["bridge"],
                        "source_cen_id_by_chr": "",
                        "source_chr": "",
                        "source_start": "",
                        "source_end": "",
                        "source_query_window_start": "",
                        "source_query_window_end": "",
                        "flank_bp": args.flank_bp,
                        "block_count": "",
                        "R108_candidate_chr_summary": "",
                        "best_R108_chr": "",
                        "best_R108_projected_start": "",
                        "best_R108_projected_end": "",
                        "nearest_R108_CEN": "",
                        "nearest_R108_CEN_id": "",
                        "distance_to_nearest_R108_CEN_bp": "",
                        "classification": f"missing_input:{path}",
                    }
                )
            continue
        source_cens = read_bed(Path(job["source_cen"]))
        blocks = read_blocks(Path(job["blocks"]), args.min_block_len)
        rows.extend(
            project_source_to_r108(
                str(job["species"]),
                source_cens,
                r108_cens,
                blocks,
                str(job["direction"]),
                args.flank_bp,
                args.distance_cutoff,
                str(job["bridge"]),
            )
        )

    out_tsv = out_dir / f"x8_active_CENH3_cores_projected_to_R108_ancestral_CENs.{suffix}.tsv"
    summary_tsv = out_dir / f"x8_active_CENH3_cores_projected_to_R108_ancestral_CENs.{suffix}.summary.tsv"
    write_tsv(out_tsv, rows)
    summarize_presence(rows, summary_tsv)

    print(out_tsv)
    print(summary_tsv)


if __name__ == "__main__":
    main()
