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
            p = raw.rstrip("\n").split("\t")
            rows.append({"chr": p[0], "start": int(float(p[1])), "end": int(float(p[2])), "name": p[3] if len(p) > 3 else ""})
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
                "query_len": q_len,
                "target_len": int(float(row.get("target_len", q_len))),
            })
    return rows


def ov(a, b, c, d):
    return max(0, min(b, d) - max(a, c))


def nearest_r108_cen(query_chr, start, end, r108_cens):
    best = None
    for cen in r108_cens:
        if cen["chr"] != query_chr:
            continue
        overlap = ov(start, end, cen["start"], cen["end"])
        if overlap:
            dist = 0
        elif end < cen["start"]:
            dist = cen["start"] - end
        else:
            dist = start - cen["end"]
        item = (dist, -overlap, cen)
        if best is None or item < best:
            best = item
    return best


def write_tsv(path, rows, fields):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    config_path = os.environ.get("CONFIG", "config.yaml")
    with open(config_path) as handle:
        cfg = json.load(handle)
    root = cfg["project_root"]
    run_root = cfg["run_root"]
    flank_bp = int(os.environ.get("MPO_CEN_FLANK_BP", str(cfg.get("flank_bp", 3000000))))
    min_block_len = int(os.environ.get("MIN_BLOCK_LEN", "1000"))

    mpo_cens = read_bed(abspath(root, cfg["species"]["genome_Mpo"]["cen_bed"]))
    r108_cens = read_bed(abspath(root, cfg["species"]["genome_R108"]["cen_bed"]))
    blocks = read_blocks(abspath(root, cfg["sequence_blocks"]["R108_to_Mpo"]), min_block_len)

    rows = []
    for mpo in sorted(mpo_cens, key=lambda x: int(x["chr"].replace("Chr", ""))):
        w_start = max(0, mpo["start"] - flank_bp)
        w_end = mpo["end"] + flank_bp
        selected = [b for b in blocks if b["target_chr"] == mpo["chr"] and ov(w_start, w_end, b["target_start"], b["target_end"]) > 0]
        by_query = defaultdict(list)
        for b in selected:
            by_query[b["query_chr"]].append(b)
        ranked = sorted(by_query.items(), key=lambda kv: sum(x["query_len"] for x in kv[1]), reverse=True)
        row = {
            "Mpo_CENH3_core": f"{mpo['chr']}:{mpo['start']}-{mpo['end']}",
            "Mpo_chr": mpo["chr"],
            "Mpo_start": mpo["start"],
            "Mpo_end": mpo["end"],
            "Mpo_query_window_start": w_start,
            "Mpo_query_window_end": w_end,
            "block_count": len(selected),
            "R108_candidate_chr_summary": ";".join(f"{chrom}:{sum(x['query_len'] for x in xs)}bp/{len(xs)}blocks" for chrom, xs in ranked[:6]),
            "best_R108_chr": "",
            "best_R108_projected_start": "",
            "best_R108_projected_end": "",
            "nearest_R108_CEN": "",
            "distance_to_nearest_R108_CEN_bp": "",
            "classification": "unresolved",
        }
        if ranked:
            best_chr, best_blocks = ranked[0]
            coords = []
            for b in best_blocks:
                coords.extend([b["query_start"], b["query_end"]])
            start, end = min(coords), max(coords)
            near = nearest_r108_cen(best_chr, start, end, r108_cens)
            row.update({
                "best_R108_chr": best_chr,
                "best_R108_projected_start": start,
                "best_R108_projected_end": end,
            })
            if near:
                dist, neg_overlap, cen = near
                cen_id = "CEN" + cen["chr"].replace("Chr", "")
                row.update({
                    "nearest_R108_CEN": f"{cen_id}|{cen['chr']}:{cen['start']}-{cen['end']}",
                    "distance_to_nearest_R108_CEN_bp": dist,
                    "classification": "reciprocal_support_near_ancestral_CEN" if dist <= flank_bp else "reciprocal_support_same_chr_distal",
                })
        rows.append(row)

    fields = [
        "Mpo_CENH3_core", "Mpo_chr", "Mpo_start", "Mpo_end", "Mpo_query_window_start", "Mpo_query_window_end",
        "block_count", "R108_candidate_chr_summary", "best_R108_chr", "best_R108_projected_start", "best_R108_projected_end",
        "nearest_R108_CEN", "distance_to_nearest_R108_CEN_bp", "classification",
    ]
    write_tsv(os.path.join(run_root, "04_reciprocal_project_Mpo_CENs_to_x8/Mpo_active_CENs_back_projected_to_R108.tsv"), rows, fields)
    print(f"Wrote {len(rows)} reciprocal rows")


if __name__ == "__main__":
    main()
