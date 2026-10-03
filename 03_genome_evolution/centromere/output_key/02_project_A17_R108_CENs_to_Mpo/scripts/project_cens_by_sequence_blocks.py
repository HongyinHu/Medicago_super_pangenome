#!/usr/bin/env python3
import csv
import json
import os
from collections import defaultdict


def abspath(root, path):
    return path if os.path.isabs(path) else os.path.join(root, path)


def read_bed(path):
    rows = []
    with open(path) as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            rows.append({
                "chr": parts[0],
                "start": int(float(parts[1])),
                "end": int(float(parts[2])),
                "name": parts[3] if len(parts) > 3 else f"{parts[0]}:{parts[1]}-{parts[2]}",
            })
    return rows


def read_blocks(path, min_len):
    rows = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            q_len = int(float(row.get("query_len", 0)))
            if q_len < min_len:
                continue
            rows.append({
                "query_chr": row["query_chr"],
                "query_start": int(float(row["query_start"])),
                "query_end": int(float(row["query_end"])),
                "target_chr": row["target_chr"],
                "target_start": int(float(row["target_start"])),
                "target_end": int(float(row["target_end"])),
                "strand": row.get("strand", ""),
                "query_len": q_len,
                "target_len": int(float(row.get("target_len", q_len))),
                "block_id": row.get("block_id", ""),
            })
    return rows


def overlap_len(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def nearest_active(target_chr, start, end, active_cens):
    best = None
    for cen in active_cens:
        if cen["chr"] != target_chr:
            continue
        ov = overlap_len(start, end, cen["start"], cen["end"])
        if ov:
            dist = 0
        elif end < cen["start"]:
            dist = cen["start"] - end
        else:
            dist = start - cen["end"]
        item = (dist, ov, cen)
        if best is None or item < best:
            best = item
    return best


def write_tsv(path, rows, fields):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    config_path = os.environ.get("CONFIG", "config.yaml")
    with open(config_path) as handle:
        cfg = json.load(handle)
    root = cfg["project_root"]
    run_root = cfg["run_root"]
    flank_bp = int(os.environ.get("QUERY_FLANK_BP", str(cfg.get("flank_bp", 3000000))))
    min_block_len = int(os.environ.get("MIN_BLOCK_LEN", "1000"))
    active_distance_bp = int(os.environ.get("ACTIVE_DISTANCE_BP", str(flank_bp)))

    r108_cens = read_bed(abspath(root, cfg["species"]["genome_R108"]["cen_bed"]))
    mpo_cens = read_bed(abspath(root, cfg["species"]["genome_Mpo"]["cen_bed"]))
    blocks = read_blocks(abspath(root, cfg["sequence_blocks"]["R108_to_Mpo"]), min_block_len)

    rows = []
    detail = []
    for cen in sorted(r108_cens, key=lambda x: int(x["chr"].replace("Chr", ""))):
        cen_id = "CEN" + cen["chr"].replace("Chr", "")
        q_start = max(0, cen["start"] - flank_bp)
        q_end = cen["end"] + flank_bp
        selected = [
            b for b in blocks
            if b["query_chr"] == cen["chr"] and overlap_len(q_start, q_end, b["query_start"], b["query_end"]) > 0
        ]
        by_target = defaultdict(list)
        for b in selected:
            by_target[b["target_chr"]].append(b)
        ranked = sorted(
            by_target.items(),
            key=lambda kv: sum(x["query_len"] for x in kv[1]),
            reverse=True,
        )
        row = {
            "ancestral_centromere_id": cen_id,
            "R108_chr": cen["chr"],
            "R108_cen_start": cen["start"],
            "R108_cen_end": cen["end"],
            "R108_query_window_start": q_start,
            "R108_query_window_end": q_end,
            "block_count": len(selected),
            "candidate_target_summary": ";".join(
                f"{chrom}:{sum(x['query_len'] for x in xs)}bp/{len(xs)}blocks" for chrom, xs in ranked[:6]
            ),
            "Mpo_chr": "",
            "Mpo_projected_start": "",
            "Mpo_projected_end": "",
            "Mpo_projected_span_bp": "",
            "Mpo_nearest_active_CENH3_core": "",
            "distance_to_nearest_active_CENH3_bp": "",
            "overlaps_Mpo_CENH3_core": "no",
            "Mpo_CENH3_overlap_bp": 0,
            "classification": "deleted_or_unresolved",
            "interpretation": "No sequence blocks projected this R108 centromere window to Mpo.",
        }
        if ranked:
            best_chr, best_blocks = ranked[0]
            coords = []
            for b in best_blocks:
                coords.extend([b["target_start"], b["target_end"]])
                detail.append({
                    "ancestral_centromere_id": cen_id,
                    **b,
                })
            start, end = min(coords), max(coords)
            near = nearest_active(best_chr, start, end, mpo_cens)
            row.update({
                "Mpo_chr": best_chr,
                "Mpo_projected_start": start,
                "Mpo_projected_end": end,
                "Mpo_projected_span_bp": end - start,
                "classification": "projected_without_active_CENH3",
                "interpretation": "Sequence blocks project this R108 centromere window to Mpo but do not overlap the selected active CENH3 core.",
            })
            if near:
                dist, ov, hit = near
                row.update({
                    "Mpo_nearest_active_CENH3_core": f"{hit['chr']}:{hit['start']}-{hit['end']}",
                    "distance_to_nearest_active_CENH3_bp": dist,
                    "Mpo_CENH3_overlap_bp": ov,
                    "overlaps_Mpo_CENH3_core": "yes" if ov > 0 else "no",
                })
                if ov > 0:
                    row["classification"] = "conserved_active_centromere"
                    row["interpretation"] = "Sequence blocks project this R108 centromere window onto an active Mpo CENH3-defined centromere."
                elif dist <= active_distance_bp:
                    row["classification"] = "near_active_centromere"
                    row["interpretation"] = "Sequence blocks project this R108 centromere window near an active Mpo CENH3-defined centromere; finer anchor and signal checks are needed."
            if len(ranked) > 1:
                row["interpretation"] += " Multiple target chromosomes contain supporting blocks, consistent with rearrangement or repetitive centromeric sequence."
        rows.append(row)

    fields = [
        "ancestral_centromere_id", "R108_chr", "R108_cen_start", "R108_cen_end",
        "R108_query_window_start", "R108_query_window_end", "block_count",
        "candidate_target_summary", "Mpo_chr", "Mpo_projected_start", "Mpo_projected_end",
        "Mpo_projected_span_bp", "Mpo_nearest_active_CENH3_core",
        "distance_to_nearest_active_CENH3_bp", "overlaps_Mpo_CENH3_core",
        "Mpo_CENH3_overlap_bp", "classification", "interpretation",
    ]
    write_tsv(os.path.join(run_root, "02_project_A17_R108_CENs_to_Mpo/results/R108_CEN1_CEN8_sequence_block_projection_to_Mpo.tsv"), rows, fields)
    write_tsv(os.path.join(run_root, "02_project_A17_R108_CENs_to_Mpo/results/R108_CEN_sequence_block_projection_detail.tsv"), detail, list(detail[0].keys()) if detail else ["ancestral_centromere_id"])
    print(f"Wrote {len(rows)} sequence-block projection rows")


if __name__ == "__main__":
    main()
