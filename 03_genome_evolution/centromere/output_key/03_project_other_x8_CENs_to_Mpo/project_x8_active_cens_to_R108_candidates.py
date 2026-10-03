#!/usr/bin/env python3
"""Long-format x=8 active CENH3 core projections to R108 CEN candidates."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_bed(path: Path) -> list[dict[str, object]]:
    rows = []
    with path.open() as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.strip().split()
            rows.append({"chr": fields[0], "start": int(float(fields[1])), "end": int(float(fields[2]))})
    return rows


def read_blocks(path: Path, min_len: int) -> list[dict[str, object]]:
    rows = []
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


def cen_id(chrom: str) -> str:
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


def project_candidates(
    species: str,
    bridge: str,
    source_cens: list[dict[str, object]],
    r108_cens: list[dict[str, object]],
    blocks: list[dict[str, object]],
    direction: str,
    flank_bp: int,
    distance_cutoff: int,
) -> list[dict[str, object]]:
    rows = []
    for source in sorted(source_cens, key=lambda x: int(str(x["chr"]).replace("Chr", ""))):
        source_id = cen_id(str(source["chr"]))
        w_start = max(0, int(source["start"]) - flank_bp)
        w_end = int(source["end"]) + flank_bp
        if direction == "query_to_target":
            selected = [
                b
                for b in blocks
                if b["query_chr"] == source["chr"] and ov(w_start, w_end, int(b["query_start"]), int(b["query_end"])) > 0
            ]
            group_chr = "target_chr"
            proj_start = "target_start"
            proj_end = "target_end"
            length_field = "query_len"
        elif direction == "target_to_query":
            selected = [
                b
                for b in blocks
                if b["target_chr"] == source["chr"] and ov(w_start, w_end, int(b["target_start"]), int(b["target_end"])) > 0
            ]
            group_chr = "query_chr"
            proj_start = "query_start"
            proj_end = "query_end"
            length_field = "query_len"
        else:
            raise ValueError(direction)

        by_chr: dict[str, list[dict[str, object]]] = defaultdict(list)
        for block in selected:
            by_chr[str(block[group_chr])].append(block)
        ranked = sorted(by_chr.items(), key=lambda kv: sum(int(x[length_field]) for x in kv[1]), reverse=True)
        for rank, (r108_chr, xs) in enumerate(ranked, start=1):
            coords = []
            for block in xs:
                coords.extend([int(block[proj_start]), int(block[proj_end])])
            start, end = min(coords), max(coords)
            near = nearest_cen(r108_chr, start, end, r108_cens)
            nearest_id = ""
            nearest_label = ""
            distance = ""
            overlap_bp = ""
            classification = "projects_to_R108_no_near_CEN"
            if near:
                dist, neg_overlap, cen = near
                nearest_id = cen_id(str(cen["chr"]))
                nearest_label = f"{nearest_id}|{cen['chr']}:{cen['start']}-{cen['end']}"
                distance = dist
                overlap_bp = -neg_overlap
                if dist <= distance_cutoff:
                    classification = "supports_R108_CEN_candidate"
                else:
                    classification = "same_chr_distal_to_R108_CEN"
            rows.append(
                {
                    "species": species,
                    "bridge": bridge,
                    "source_cen_id_by_chr": source_id,
                    "source_chr": source["chr"],
                    "source_start": source["start"],
                    "source_end": source["end"],
                    "flank_bp": flank_bp,
                    "R108_chr": r108_chr,
                    "R108_chr_rank_by_block_length": rank,
                    "block_count": len(xs),
                    "sum_query_len": sum(int(x[length_field]) for x in xs),
                    "R108_projected_start": start,
                    "R108_projected_end": end,
                    "nearest_R108_CEN_id": nearest_id,
                    "nearest_R108_CEN": nearest_label,
                    "distance_to_nearest_R108_CEN_bp": distance,
                    "overlap_with_nearest_R108_CEN_bp": overlap_bp,
                    "classification": classification,
                }
            )
    return rows


def r108_self_candidates(r108_cens: list[dict[str, object]], flank_bp: int) -> list[dict[str, object]]:
    rows = []
    for cen in sorted(r108_cens, key=lambda x: int(str(x["chr"]).replace("Chr", ""))):
        cid = cen_id(str(cen["chr"]))
        rows.append(
            {
                "species": "R108",
                "bridge": "self",
                "source_cen_id_by_chr": cid,
                "source_chr": cen["chr"],
                "source_start": cen["start"],
                "source_end": cen["end"],
                "flank_bp": flank_bp,
                "R108_chr": cen["chr"],
                "R108_chr_rank_by_block_length": 1,
                "block_count": "",
                "sum_query_len": int(cen["end"]) - int(cen["start"]),
                "R108_projected_start": cen["start"],
                "R108_projected_end": cen["end"],
                "nearest_R108_CEN_id": cid,
                "nearest_R108_CEN": f"{cid}|{cen['chr']}:{cen['start']}-{cen['end']}",
                "distance_to_nearest_R108_CEN_bp": 0,
                "overlap_with_nearest_R108_CEN_bp": int(cen["end"]) - int(cen["start"]),
                "classification": "self_active_R108_CEN",
            }
        )
    return rows


def summarize(rows: list[dict[str, object]], min_support_len: int, min_support_blocks: int) -> list[dict[str, object]]:
    all_cens = [f"CEN{i}" for i in range(1, 9)]
    supported: dict[str, set[str]] = defaultdict(set)
    details: dict[str, list[str]] = defaultdict(list)
    best_cen5: dict[str, dict[str, object]] = {}

    for row in rows:
        species = str(row["species"])
        cid = str(row["nearest_R108_CEN_id"])
        if not cid:
            continue
        if row["classification"] not in {"supports_R108_CEN_candidate", "self_active_R108_CEN"}:
            continue
        block_count_raw = row["block_count"]
        block_count = 999 if block_count_raw == "" else int(block_count_raw)
        support_len = int(row["sum_query_len"])
        if support_len < min_support_len or block_count < min_support_blocks:
            continue
        supported[species].add(cid)
        details[species].append(f"{row['source_cen_id_by_chr']}->{cid}:{support_len}bp/{block_count}blocks/rank{row['R108_chr_rank_by_block_length']}")
        if cid == "CEN5":
            current = best_cen5.get(species)
            key = (support_len, block_count, -int(row["R108_chr_rank_by_block_length"]))
            if current is None or key > current["_key"]:
                best = dict(row)
                best["_key"] = key
                best_cen5[species] = best

    out = []
    for species in sorted({str(row["species"]) for row in rows}):
        present = supported.get(species, set())
        cen5 = best_cen5.get(species)
        out.append(
            {
                "species": species,
                "supported_R108_CEN_count": len(present),
                "supported_R108_CENs": ",".join(c for c in all_cens if c in present),
                "missing_R108_CENs": ",".join(c for c in all_cens if c not in present),
                "CEN5_supported": "yes" if "CEN5" in present else "no",
                "best_CEN5_source_core": cen5["source_cen_id_by_chr"] if cen5 else "",
                "best_CEN5_support": f"{cen5['sum_query_len']}bp/{cen5['block_count']}blocks/rank{cen5['R108_chr_rank_by_block_length']}" if cen5 else "",
                "detail": ";".join(details.get(species, [])),
            }
        )
    return out


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    rows = [{k: v for k, v in row.items() if k != "_key"} for row in rows]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default="path/to/project/10.centromere_analysis")
    parser.add_argument("--flank-bp", type=int, default=0)
    parser.add_argument("--min-block-len", type=int, default=1000)
    parser.add_argument("--distance-cutoff", type=int, default=3_000_000)
    parser.add_argument("--min-support-len", type=int, default=5000)
    parser.add_argument("--min-support-blocks", type=int, default=2)
    args = parser.parse_args()

    root = Path(args.project_root)
    out_dir = root / "output_key/03_project_other_x8_CENs_to_Mpo"
    suffix = "core_only" if args.flank_bp == 0 else f"pm{args.flank_bp // 1_000_000}Mb"

    r108_cens = read_bed(root / "output_synthen/01_JCVI/genome_A17_genome_R108/genome_R108.cen_core.bed")
    rows: list[dict[str, object]] = r108_self_candidates(r108_cens, args.flank_bp)

    jobs = [
        ("A17", "A17_to_R108_direct", root / "output_synthen/01_JCVI/genome_A17_genome_R108/genome_A17.cen_core.bed", root / "output_synthen/02_cent_synteny/genome_A17_genome_R108/genome_A17_vs_genome_R108.rawPAF.blocks.tsv", "query_to_target"),
        ("Msa", "Msa_to_R108_reverse_rawPAF", root / "output_synthen/01_JCVI/R108_Msa_synteny_lz10_run/genome_Msa.cen_core.bed", root / "output_synthen/02_cent_synteny/genome_R108_genome_Msa/genome_R108_vs_genome_Msa.rawPAF.blocks.tsv", "target_to_query"),
        ("474", "474_to_R108_reverse_rawPAF", root / "output_synthen/01_JCVI/genome_R108_genome_474_synteny_lz10_run/genome_474.cen_core.bed", root / "output_synthen/02_cent_synteny/genome_R108_genome_474/genome_R108_vs_genome_474.rawPAF.blocks.tsv", "target_to_query"),
    ]
    for species, bridge, cen_path, block_path, direction in jobs:
        if not cen_path.exists() or not block_path.exists():
            rows.append(
                {
                    "species": species,
                    "bridge": bridge,
                    "source_cen_id_by_chr": "",
                    "source_chr": "",
                    "source_start": "",
                    "source_end": "",
                    "flank_bp": args.flank_bp,
                    "R108_chr": "",
                    "R108_chr_rank_by_block_length": "",
                    "block_count": "",
                    "sum_query_len": "",
                    "R108_projected_start": "",
                    "R108_projected_end": "",
                    "nearest_R108_CEN_id": "",
                    "nearest_R108_CEN": "",
                    "distance_to_nearest_R108_CEN_bp": "",
                    "overlap_with_nearest_R108_CEN_bp": "",
                    "classification": "missing_input",
                }
            )
            continue
        rows.extend(
            project_candidates(
                species,
                bridge,
                read_bed(cen_path),
                r108_cens,
                read_blocks(block_path, args.min_block_len),
                direction,
                args.flank_bp,
                args.distance_cutoff,
            )
        )

    detail_path = out_dir / f"x8_active_CENH3_cores_to_R108_candidate_groups.{suffix}.tsv"
    summary_path = out_dir / f"x8_active_CENH3_cores_to_R108_candidate_groups.{suffix}.summary.tsv"
    write_tsv(detail_path, rows)
    write_tsv(summary_path, summarize(rows, args.min_support_len, args.min_support_blocks))
    print(detail_path)
    print(summary_path)


if __name__ == "__main__":
    main()
