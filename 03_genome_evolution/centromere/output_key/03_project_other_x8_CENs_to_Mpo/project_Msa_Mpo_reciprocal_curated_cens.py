#!/usr/bin/env python3
"""Project curated Msa CEN domains to Mpo and Mpo active CENs back to Msa."""

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
            p = raw.strip().split()
            rows.append({"chr": p[0], "start": int(float(p[1])), "end": int(float(p[2]))})
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
                }
            )
    return rows


def ov(a: int, b: int, c: int, d: int) -> int:
    return max(0, min(b, d) - max(a, c))


def cen_id(chrom: str) -> str:
    return "CEN" + chrom.replace("Chr", "")


def active_overlap(chrom: str, start: int, end: int, active: list[dict[str, object]]) -> list[tuple[dict[str, object], int]]:
    hits = []
    for cen in active:
        if cen["chr"] != chrom:
            continue
        x = ov(start, end, int(cen["start"]), int(cen["end"]))
        if x:
            hits.append((cen, x))
    hits.sort(key=lambda x: x[1], reverse=True)
    return hits


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


def msa_to_mpo(msa_cens, mpo_cens, blocks, flank_bp):
    rows = []
    for cen in sorted(msa_cens, key=lambda x: int(str(x["chr"]).replace("Chr", ""))):
        cid = cen_id(str(cen["chr"]))
        q_start = max(0, int(cen["start"]) - flank_bp)
        q_end = int(cen["end"]) + flank_bp
        selected = [b for b in blocks if b["query_chr"] == cen["chr"] and ov(q_start, q_end, int(b["query_start"]), int(b["query_end"])) > 0]
        by_target = defaultdict(list)
        for block in selected:
            by_target[str(block["target_chr"])].append(block)
        ranked = sorted(by_target.items(), key=lambda kv: sum(int(x["query_len"]) for x in kv[1]), reverse=True)
        row = {
            "centromere_id": cid,
            "Msa_chr": cen["chr"],
            "Msa_start": cen["start"],
            "Msa_end": cen["end"],
            "flank_bp": flank_bp,
            "block_count": len(selected),
            "Mpo_candidate_chr_summary": ";".join(f"{chrom}:{sum(int(x['query_len']) for x in xs)}bp/{len(xs)}blocks" for chrom, xs in ranked[:8]),
            "Mpo_chr": "",
            "Mpo_projected_start": "",
            "Mpo_projected_end": "",
            "overlaps_Mpo_CENH3_core": "no",
            "Mpo_CENH3_core": "",
            "classification": "unresolved",
        }
        if ranked:
            best_chr, best_blocks = ranked[0]
            coords = []
            for block in best_blocks:
                coords.extend([int(block["target_start"]), int(block["target_end"])])
            start, end = min(coords), max(coords)
            row.update({"Mpo_chr": best_chr, "Mpo_projected_start": start, "Mpo_projected_end": end})
            hits = active_overlap(best_chr, start, end, mpo_cens)
            if hits:
                hit, overlap = hits[0]
                row.update(
                    {
                        "overlaps_Mpo_CENH3_core": "yes",
                        "Mpo_CENH3_core": f"{hit['chr']}:{hit['start']}-{hit['end']}",
                        "classification": "projects_across_active_Mpo_CEN",
                    }
                )
            else:
                row["classification"] = "projects_without_active_Mpo_CEN"
        rows.append(row)
    return rows


def mpo_to_msa(mpo_cens, msa_cens, blocks, flank_bp):
    rows = []
    for cen in sorted(mpo_cens, key=lambda x: int(str(x["chr"]).replace("Chr", ""))):
        q_start = max(0, int(cen["start"]) - flank_bp)
        q_end = int(cen["end"]) + flank_bp
        selected = [b for b in blocks if b["target_chr"] == cen["chr"] and ov(q_start, q_end, int(b["target_start"]), int(b["target_end"])) > 0]
        by_query = defaultdict(list)
        for block in selected:
            by_query[str(block["query_chr"])].append(block)
        ranked = sorted(by_query.items(), key=lambda kv: sum(int(x["query_len"]) for x in kv[1]), reverse=True)
        row = {
            "Mpo_CENH3_core": f"{cen['chr']}:{cen['start']}-{cen['end']}",
            "Mpo_chr": cen["chr"],
            "Mpo_start": cen["start"],
            "Mpo_end": cen["end"],
            "flank_bp": flank_bp,
            "block_count": len(selected),
            "Msa_candidate_chr_summary": ";".join(f"{chrom}:{sum(int(x['query_len']) for x in xs)}bp/{len(xs)}blocks" for chrom, xs in ranked[:8]),
            "best_Msa_chr": "",
            "best_Msa_projected_start": "",
            "best_Msa_projected_end": "",
            "nearest_Msa_CEN": "",
            "nearest_Msa_CEN_id": "",
            "distance_to_nearest_Msa_CEN_bp": "",
            "classification": "unresolved",
        }
        if ranked:
            best_chr, best_blocks = ranked[0]
            coords = []
            for block in best_blocks:
                coords.extend([int(block["query_start"]), int(block["query_end"])])
            start, end = min(coords), max(coords)
            row.update({"best_Msa_chr": best_chr, "best_Msa_projected_start": start, "best_Msa_projected_end": end})
            near = nearest_cen(best_chr, start, end, msa_cens)
            if near:
                dist, _neg_overlap, mc = near
                mid = cen_id(str(mc["chr"]))
                row.update(
                    {
                        "nearest_Msa_CEN": f"{mid}|{mc['chr']}:{mc['start']}-{mc['end']}",
                        "nearest_Msa_CEN_id": mid,
                        "distance_to_nearest_Msa_CEN_bp": dist,
                        "classification": "reciprocal_support_near_Msa_CEN" if dist <= flank_bp else "reciprocal_support_same_chr_distal",
                    }
                )
        rows.append(row)
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default="path/to/project/10.centromere_analysis")
    parser.add_argument("--flank-bp", type=int, default=3_000_000)
    parser.add_argument("--min-block-len", type=int, default=1000)
    args = parser.parse_args()
    root = Path(args.project_root)
    suffix = "core_only" if args.flank_bp == 0 else f"pm{args.flank_bp // 1_000_000}Mb"

    msa_cens = read_bed(root / "output_synthen/01_JCVI/R108_Msa_synteny_lz10_run/genome_Msa.cen_core.bed")
    mpo_cens = read_bed(root / "output_synthen/01_JCVI/R108_Mpo_synteny_lz10_run/genome_Mpo.cen_core.bed")
    blocks = read_blocks(root / "output_synthen/02_cent_synteny/genome_Msa_genome_Mpo/genome_Msa_vs_genome_Mpo.rawPAF.blocks.tsv", args.min_block_len)

    out1 = root / f"output_key/03_project_other_x8_CENs_to_Mpo/Msa_curated_CEN1_CEN8_sequence_projection_to_Mpo.{suffix}.tsv"
    out2 = root / f"output_key/04_reciprocal_project_Mpo_CENs_to_x8/Mpo_active_CENs_back_projected_to_Msa_curated.{suffix}.tsv"
    write_tsv(out1, msa_to_mpo(msa_cens, mpo_cens, blocks, args.flank_bp))
    write_tsv(out2, mpo_to_msa(mpo_cens, msa_cens, blocks, args.flank_bp))
    print(out1)
    print(out2)


if __name__ == "__main__":
    main()
