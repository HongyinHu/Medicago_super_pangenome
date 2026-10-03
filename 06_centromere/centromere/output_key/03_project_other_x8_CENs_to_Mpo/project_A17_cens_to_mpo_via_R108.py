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


def read_gene_coord(path):
    by_chr = defaultdict(list)
    with open(path) as handle:
        for raw in handle:
            if not raw.strip():
                continue
            p = raw.rstrip("\n").split("\t")
            gene = {"gene": p[0], "chr": p[1], "start": int(p[2]), "end": int(p[3]), "strand": p[4]}
            by_chr[gene["chr"]].append(gene)
    for chrom in by_chr:
        by_chr[chrom].sort(key=lambda x: (x["start"], x["end"]))
    return by_chr


def read_links(path, src_col="sa_gene", dst_col="sm_gene"):
    links = defaultdict(list)
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            src = row[src_col]
            links[src].append({
                "src_gene": src,
                "dst_gene": row[dst_col],
                "src_chr": row["sa_chr"],
                "src_start": int(float(row["sa_start"])),
                "src_end": int(float(row["sa_end"])),
                "dst_chr": row["sm_chr"],
                "dst_start": int(float(row["sm_start"])),
                "dst_end": int(float(row["sm_end"])),
                "block_id": row.get("block_id", ""),
            })
    return links


def flank_genes(by_chr, chrom, start, end, flank_bp, n):
    genes = by_chr.get(chrom, [])
    left = [g for g in genes if g["end"] < start and g["end"] >= start - flank_bp]
    right = [g for g in genes if g["start"] > end and g["start"] <= end + flank_bp]
    return sorted(left, key=lambda g: start - g["end"])[:n], sorted(right, key=lambda g: g["start"] - end)[:n]


def overlap(a, b, c, d):
    return max(0, min(b, d) - max(a, c))


def active_hits(chrom, start, end, active):
    hits = []
    for c in active:
        if c["chr"] != chrom:
            continue
        ov = overlap(start, end, c["start"], c["end"])
        if ov:
            hits.append((c, ov))
    hits.sort(key=lambda x: x[1], reverse=True)
    return hits


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
    flank_bp = int(os.environ.get("FLANK_BP", str(cfg.get("flank_bp", 3000000))))
    n_anchor = int(os.environ.get("ANCHORS_PER_SIDE", str(cfg.get("anchors_per_side", 12))))

    a17_cens = read_bed(abspath(root, cfg["species"]["genome_A17"]["cen_bed"]))
    mpo_cens = read_bed(abspath(root, cfg["species"]["genome_Mpo"]["cen_bed"]))
    a17_by_chr = read_gene_coord(abspath(root, cfg["species"]["genome_A17"]["gene_coord"]))
    a17_to_r108 = read_links(abspath(root, cfg["links"]["A17_to_R108"]))
    r108_to_mpo = read_links(abspath(root, cfg["links"]["R108_to_Mpo"]))

    rows = []
    detail = []
    for cen in sorted(a17_cens, key=lambda x: int(x["chr"].replace("Chr", ""))):
        cen_id = "CEN" + cen["chr"].replace("Chr", "")
        left, right = flank_genes(a17_by_chr, cen["chr"], cen["start"], cen["end"], flank_bp, n_anchor)
        paths = []
        for side, genes in [("left", left), ("right", right)]:
            for g in genes:
                r108_links = a17_to_r108.get(g["gene"], [])
                if not r108_links:
                    detail.append({"ancestral_centromere_id": cen_id, "side": side, "A17_gene": g["gene"], "R108_gene": "", "Mpo_gene": "", "Mpo_chr": "", "Mpo_start": "", "Mpo_end": ""})
                for rlink in r108_links:
                    mpo_links = r108_to_mpo.get(rlink["dst_gene"], [])
                    if not mpo_links:
                        detail.append({"ancestral_centromere_id": cen_id, "side": side, "A17_gene": g["gene"], "R108_gene": rlink["dst_gene"], "Mpo_gene": "", "Mpo_chr": "", "Mpo_start": "", "Mpo_end": ""})
                    for mlink in mpo_links:
                        rec = {
                            "ancestral_centromere_id": cen_id,
                            "side": side,
                            "A17_gene": g["gene"],
                            "R108_gene": rlink["dst_gene"],
                            "Mpo_gene": mlink["dst_gene"],
                            "Mpo_chr": mlink["dst_chr"],
                            "Mpo_start": mlink["dst_start"],
                            "Mpo_end": mlink["dst_end"],
                        }
                        detail.append(rec)
                        paths.append(rec)

        by_chr = defaultdict(list)
        for p in paths:
            by_chr[p["Mpo_chr"]].append(p)
        ranked = sorted(by_chr.items(), key=lambda kv: len(kv[1]), reverse=True)
        row = {
            "ancestral_centromere_id": cen_id,
            "A17_chr": cen["chr"],
            "A17_start": cen["start"],
            "A17_end": cen["end"],
            "left_anchor_count": len(left),
            "right_anchor_count": len(right),
            "mapped_path_count": len(paths),
            "candidate_chr_summary": ";".join(f"{chrom}:{len(items)}" for chrom, items in ranked[:6]),
            "Mpo_chr": "",
            "Mpo_projected_start": "",
            "Mpo_projected_end": "",
            "overlaps_Mpo_CENH3_core": "no",
            "Mpo_CENH3_core": "",
            "Mpo_CENH3_overlap_bp": 0,
            "classification": "deleted_or_unresolved",
            "interpretation": "No chained A17->R108->Mpo anchors were found.",
        }
        if ranked:
            best_chr, best = ranked[0]
            coords = []
            for p in best:
                coords.extend([p["Mpo_start"], p["Mpo_end"]])
            start, end = min(coords), max(coords)
            row.update({"Mpo_chr": best_chr, "Mpo_projected_start": start, "Mpo_projected_end": end})
            hits = active_hits(best_chr, start, end, mpo_cens)
            if hits:
                hit, ov = hits[0]
                row.update({
                    "overlaps_Mpo_CENH3_core": "yes",
                    "Mpo_CENH3_core": f"{hit['chr']}:{hit['start']}-{hit['end']}",
                    "Mpo_CENH3_overlap_bp": ov,
                    "classification": "projected_near_or_across_active_centromere",
                    "interpretation": "Chained A17->R108->Mpo anchors project near or across a selected active Mpo CENH3 core.",
                })
            else:
                row.update({
                    "classification": "projected_without_active_CENH3",
                    "interpretation": "Chained A17->R108->Mpo anchors project to Mpo but do not overlap a selected active CENH3 core.",
                })
            if len(ranked) > 1:
                row["interpretation"] += " Multiple Mpo chromosomes carry projected anchors, suggesting rearrangement or ambiguity."
        rows.append(row)

    out_dir = os.path.join(run_root, "03_project_other_x8_CENs_to_Mpo")
    write_tsv(os.path.join(out_dir, "A17_CEN1_CEN8_projected_to_Mpo_via_R108.tsv"), rows, [
        "ancestral_centromere_id", "A17_chr", "A17_start", "A17_end", "left_anchor_count", "right_anchor_count",
        "mapped_path_count", "candidate_chr_summary", "Mpo_chr", "Mpo_projected_start", "Mpo_projected_end",
        "overlaps_Mpo_CENH3_core", "Mpo_CENH3_core", "Mpo_CENH3_overlap_bp", "classification", "interpretation",
    ])
    write_tsv(os.path.join(out_dir, "A17_CEN_anchor_chain_to_Mpo_detail.tsv"), detail, [
        "ancestral_centromere_id", "side", "A17_gene", "R108_gene", "Mpo_gene", "Mpo_chr", "Mpo_start", "Mpo_end",
    ])
    print(f"Wrote {len(rows)} A17 projection rows")


if __name__ == "__main__":
    main()
