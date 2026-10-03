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
            })
    return rows


def ov(a, b, c, d):
    return max(0, min(b, d) - max(a, c))


def nearest_cen(chrom, start, end, cens):
    best = None
    for cen in cens:
        if cen["chr"] != chrom:
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
    with open(os.environ.get("CONFIG", "config.yaml")) as handle:
        cfg = json.load(handle)
    root = cfg["project_root"]
    run_root = cfg["run_root"]
    flank_bp = int(os.environ.get("QUERY_FLANK_BP", str(cfg.get("flank_bp", 3000000))))
    min_block_len = int(os.environ.get("MIN_BLOCK_LEN", "1000"))

    a17_cens = read_bed(abspath(root, cfg["species"]["genome_A17"]["cen_bed"]))
    r108_cens = read_bed(abspath(root, cfg["species"]["genome_R108"]["cen_bed"]))
    blocks = read_blocks(abspath(root, "output_synthen/02_cent_synteny/genome_A17_genome_R108/genome_A17_vs_genome_R108.rawPAF.blocks.tsv"), min_block_len)

    rows = []
    detail = []
    for cen in sorted(a17_cens, key=lambda x: int(x["chr"].replace("Chr", ""))):
        cen_id = "CEN" + cen["chr"].replace("Chr", "")
        q_start = max(0, cen["start"] - flank_bp)
        q_end = cen["end"] + flank_bp
        selected = [b for b in blocks if b["query_chr"] == cen["chr"] and ov(q_start, q_end, b["query_start"], b["query_end"]) > 0]
        by_target = defaultdict(list)
        for b in selected:
            by_target[b["target_chr"]].append(b)
        ranked = sorted(by_target.items(), key=lambda kv: sum(x["query_len"] for x in kv[1]), reverse=True)
        row = {
            "ancestral_centromere_id": cen_id,
            "A17_chr": cen["chr"],
            "A17_start": cen["start"],
            "A17_end": cen["end"],
            "A17_query_window_start": q_start,
            "A17_query_window_end": q_end,
            "block_count": len(selected),
            "R108_candidate_chr_summary": ";".join(f"{chrom}:{sum(x['query_len'] for x in xs)}bp/{len(xs)}blocks" for chrom, xs in ranked[:6]),
            "R108_chr": "",
            "R108_projected_start": "",
            "R108_projected_end": "",
            "nearest_R108_CEN": "",
            "distance_to_nearest_R108_CEN_bp": "",
            "classification": "unresolved",
        }
        if ranked:
            best_chr, best_blocks = ranked[0]
            coords = []
            for b in best_blocks:
                coords.extend([b["target_start"], b["target_end"]])
                detail.append({"ancestral_centromere_id": cen_id, **b})
            start, end = min(coords), max(coords)
            near = nearest_cen(best_chr, start, end, r108_cens)
            row.update({"R108_chr": best_chr, "R108_projected_start": start, "R108_projected_end": end})
            if near:
                dist, neg_overlap, rcen = near
                rcen_id = "CEN" + rcen["chr"].replace("Chr", "")
                row.update({
                    "nearest_R108_CEN": f"{rcen_id}|{rcen['chr']}:{rcen['start']}-{rcen['end']}",
                    "distance_to_nearest_R108_CEN_bp": dist,
                    "classification": "supports_matching_R108_CEN" if rcen_id == cen_id and dist <= flank_bp else "projects_to_nonmatching_or_distal_R108_CEN",
                })
        rows.append(row)

    out_dir = os.path.join(run_root, "03_project_other_x8_CENs_to_Mpo")
    suffix = "core_only" if flank_bp == 0 else "pm3Mb"
    fields = [
        "ancestral_centromere_id", "A17_chr", "A17_start", "A17_end", "A17_query_window_start", "A17_query_window_end",
        "block_count", "R108_candidate_chr_summary", "R108_chr", "R108_projected_start", "R108_projected_end",
        "nearest_R108_CEN", "distance_to_nearest_R108_CEN_bp", "classification",
    ]
    write_tsv(os.path.join(out_dir, f"A17_CEN1_CEN8_sequence_projection_to_R108.{suffix}.tsv"), rows, fields)
    write_tsv(os.path.join(out_dir, f"A17_CEN_sequence_projection_to_R108_detail.{suffix}.tsv"), detail, list(detail[0].keys()) if detail else ["ancestral_centromere_id"])
    print(f"Wrote {len(rows)} A17->R108 rows ({suffix})")


if __name__ == "__main__":
    main()
