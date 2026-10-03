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


def active_overlap(chrom, start, end, active):
    hits = []
    for cen in active:
        if cen["chr"] != chrom:
            continue
        x = ov(start, end, cen["start"], cen["end"])
        if x:
            hits.append((cen, x))
    hits.sort(key=lambda x: x[1], reverse=True)
    return hits


def write_tsv(path, rows, fields):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def summarize_source_to_target(source_cens, target_cens, blocks, flank_bp, source_label, target_label):
    rows = []
    for cen in sorted(source_cens, key=lambda x: int(x["chr"].replace("Chr", ""))):
        cen_id = "CEN" + cen["chr"].replace("Chr", "")
        q_start = max(0, cen["start"] - flank_bp)
        q_end = cen["end"] + flank_bp
        selected = [b for b in blocks if b["query_chr"] == cen["chr"] and ov(q_start, q_end, b["query_start"], b["query_end"]) > 0]
        by_target = defaultdict(list)
        for b in selected:
            by_target[b["target_chr"]].append(b)
        ranked = sorted(by_target.items(), key=lambda kv: sum(x["query_len"] for x in kv[1]), reverse=True)
        row = {
            "centromere_id": cen_id,
            f"{source_label}_chr": cen["chr"],
            f"{source_label}_start": cen["start"],
            f"{source_label}_end": cen["end"],
            "block_count": len(selected),
            f"{target_label}_candidate_chr_summary": ";".join(f"{chrom}:{sum(x['query_len'] for x in xs)}bp/{len(xs)}blocks" for chrom, xs in ranked[:6]),
            f"{target_label}_chr": "",
            f"{target_label}_projected_start": "",
            f"{target_label}_projected_end": "",
            f"overlaps_{target_label}_CENH3_core": "no",
            f"{target_label}_CENH3_core": "",
            "classification": "unresolved",
        }
        if ranked:
            best_chr, best_blocks = ranked[0]
            coords = []
            for b in best_blocks:
                coords.extend([b["target_start"], b["target_end"]])
            start, end = min(coords), max(coords)
            row.update({f"{target_label}_chr": best_chr, f"{target_label}_projected_start": start, f"{target_label}_projected_end": end})
            hits = active_overlap(best_chr, start, end, target_cens)
            if hits:
                hit, _x = hits[0]
                row.update({
                    f"overlaps_{target_label}_CENH3_core": "yes",
                    f"{target_label}_CENH3_core": f"{hit['chr']}:{hit['start']}-{hit['end']}",
                    "classification": "projects_across_active_target_CEN",
                })
            else:
                row["classification"] = "projects_without_active_target_CEN"
        rows.append(row)
    return rows


def summarize_mpo_to_msa(mpo_cens, msa_cens, blocks, flank_bp):
    rows = []
    for mpo in sorted(mpo_cens, key=lambda x: int(x["chr"].replace("Chr", ""))):
        start = max(0, mpo["start"] - flank_bp)
        end = mpo["end"] + flank_bp
        selected = [b for b in blocks if b["target_chr"] == mpo["chr"] and ov(start, end, b["target_start"], b["target_end"]) > 0]
        by_query = defaultdict(list)
        for b in selected:
            by_query[b["query_chr"]].append(b)
        ranked = sorted(by_query.items(), key=lambda kv: sum(x["query_len"] for x in kv[1]), reverse=True)
        row = {
            "Mpo_CENH3_core": f"{mpo['chr']}:{mpo['start']}-{mpo['end']}",
            "block_count": len(selected),
            "Msa_candidate_chr_summary": ";".join(f"{chrom}:{sum(x['query_len'] for x in xs)}bp/{len(xs)}blocks" for chrom, xs in ranked[:6]),
            "best_Msa_chr": "",
            "best_Msa_projected_start": "",
            "best_Msa_projected_end": "",
            "nearest_Msa_CEN": "",
            "distance_to_nearest_Msa_CEN_bp": "",
            "classification": "unresolved",
        }
        if ranked:
            best_chr, best_blocks = ranked[0]
            coords = []
            for b in best_blocks:
                coords.extend([b["query_start"], b["query_end"]])
            q_start, q_end = min(coords), max(coords)
            near = nearest_cen(best_chr, q_start, q_end, msa_cens)
            row.update({"best_Msa_chr": best_chr, "best_Msa_projected_start": q_start, "best_Msa_projected_end": q_end})
            if near:
                dist, _neg_ov, cen = near
                cen_id = "CEN" + cen["chr"].replace("Chr", "")
                row.update({
                    "nearest_Msa_CEN": f"{cen_id}|{cen['chr']}:{cen['start']}-{cen['end']}",
                    "distance_to_nearest_Msa_CEN_bp": dist,
                    "classification": "reciprocal_support_near_Msa_CEN" if dist <= flank_bp else "reciprocal_support_same_chr_distal",
                })
        rows.append(row)
    return rows


def main():
    with open(os.environ.get("CONFIG", "config.yaml")) as handle:
        cfg = json.load(handle)
    root = cfg["project_root"]
    run_root = cfg["run_root"]
    flank_bp = int(os.environ.get("QUERY_FLANK_BP", str(cfg.get("flank_bp", 3000000))))
    min_len = int(os.environ.get("MIN_BLOCK_LEN", "1000"))
    msa_cens = read_bed(abspath(root, cfg["species"]["genome_Msa"]["cen_bed"]))
    mpo_cens = read_bed(abspath(root, cfg["species"]["genome_Mpo"]["cen_bed"]))
    blocks = read_blocks(abspath(root, cfg["sequence_blocks"]["Msa_to_Mpo"]), min_len)

    out_dir = os.path.join(run_root, "03_project_other_x8_CENs_to_Mpo")
    recip_dir = os.path.join(run_root, "04_reciprocal_project_Mpo_CENs_to_x8")
    suffix = "core_only" if flank_bp == 0 else "pm3Mb"

    msa_to_mpo = summarize_source_to_target(msa_cens, mpo_cens, blocks, flank_bp, "Msa", "Mpo")
    write_tsv(os.path.join(out_dir, f"Msa_CEN1_CEN8_sequence_projection_to_Mpo.{suffix}.tsv"), msa_to_mpo, [
        "centromere_id", "Msa_chr", "Msa_start", "Msa_end", "block_count", "Mpo_candidate_chr_summary",
        "Mpo_chr", "Mpo_projected_start", "Mpo_projected_end", "overlaps_Mpo_CENH3_core", "Mpo_CENH3_core", "classification",
    ])

    mpo_to_msa = summarize_mpo_to_msa(mpo_cens, msa_cens, blocks, flank_bp)
    write_tsv(os.path.join(recip_dir, f"Mpo_active_CENs_back_projected_to_Msa.{suffix}.tsv"), mpo_to_msa, [
        "Mpo_CENH3_core", "block_count", "Msa_candidate_chr_summary", "best_Msa_chr", "best_Msa_projected_start",
        "best_Msa_projected_end", "nearest_Msa_CEN", "distance_to_nearest_Msa_CEN_bp", "classification",
    ])
    print(f"Wrote Msa<->Mpo projection tables ({suffix})")


if __name__ == "__main__":
    main()
