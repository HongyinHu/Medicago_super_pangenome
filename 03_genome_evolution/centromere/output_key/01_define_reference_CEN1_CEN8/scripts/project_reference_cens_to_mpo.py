#!/usr/bin/env python3
import csv
import json
import os
from collections import defaultdict


def abspath(root, path):
    return path if os.path.isabs(path) else os.path.join(root, path)


def read_bed(path):
    rows = []
    with open(path, newline="") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            chrom, start, end = parts[0], int(float(parts[1])), int(float(parts[2]))
            name = parts[3] if len(parts) > 3 else f"{chrom}:{start}-{end}"
            rows.append({"chr": chrom, "start": start, "end": end, "name": name, "raw": parts})
    return rows


def read_gene_coord(path):
    genes = []
    with open(path, newline="") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            gene, chrom, start, end, strand = parts[:5]
            genes.append({
                "gene": gene,
                "chr": chrom,
                "start": int(float(start)),
                "end": int(float(end)),
                "strand": strand,
                "mid": (int(float(start)) + int(float(end))) // 2,
            })
    by_chr = defaultdict(list)
    by_gene = {}
    for gene in genes:
        by_chr[gene["chr"]].append(gene)
        by_gene[gene["gene"]] = gene
    for chrom in by_chr:
        by_chr[chrom].sort(key=lambda x: (x["start"], x["end"]))
    return by_chr, by_gene


def read_gene_links(path):
    links = defaultdict(list)
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            src_gene = row.get("sa_gene")
            dst_gene = row.get("sm_gene")
            if not src_gene or not dst_gene:
                continue
            rec = {
                "src_gene": src_gene,
                "dst_gene": dst_gene,
                "src_chr": row["sa_chr"],
                "src_start": int(float(row["sa_start"])),
                "src_end": int(float(row["sa_end"])),
                "dst_chr": row["sm_chr"],
                "dst_start": int(float(row["sm_start"])),
                "dst_end": int(float(row["sm_end"])),
                "block_id": row.get("block_id", ""),
            }
            links[src_gene].append(rec)
    return links


def nearest_flank_genes(by_chr, chrom, start, end, flank_bp, n):
    genes = by_chr.get(chrom, [])
    left = [g for g in genes if g["end"] < start and g["end"] >= start - flank_bp]
    right = [g for g in genes if g["start"] > end and g["start"] <= end + flank_bp]
    left = sorted(left, key=lambda g: start - g["end"])[:n]
    right = sorted(right, key=lambda g: g["start"] - end)[:n]
    return left, right


def active_overlap(projected_chr, projected_start, projected_end, active_cens):
    hits = []
    for cen in active_cens:
        if cen["chr"] != projected_chr:
            continue
        ov = max(0, min(projected_end, cen["end"]) - max(projected_start, cen["start"]))
        if ov > 0:
            hits.append((cen, ov))
    hits.sort(key=lambda x: x[1], reverse=True)
    return hits


def summarize_projection(anchor_records):
    mapped = [r for r in anchor_records if r["mapped"]]
    if not mapped:
        return None
    by_chr = defaultdict(list)
    for rec in mapped:
        for link in rec["links"]:
            by_chr[link["dst_chr"]].append(link)
    ranked = sorted(by_chr.items(), key=lambda kv: (len(kv[1]), sum(abs(x["dst_end"] - x["dst_start"]) for x in kv[1])), reverse=True)
    best_chr, best_links = ranked[0]
    coords = []
    for link in best_links:
        coords.extend([link["dst_start"], link["dst_end"]])
    return {
        "chr": best_chr,
        "start": min(coords),
        "end": max(coords),
        "mapped_anchor_count": len(best_links),
        "candidate_chr_count": len(ranked),
        "candidate_chr_summary": ";".join(f"{chrom}:{len(links)}" for chrom, links in ranked[:5]),
    }


def write_tsv(path, rows, fieldnames):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    config_path = os.environ.get("CONFIG", "config.yaml")
    with open(config_path) as handle:
        cfg = json.load(handle)
    root = cfg["project_root"]
    run_root = cfg["run_root"]
    flank_bp = int(cfg.get("flank_bp", 3000000))
    n_anchor = int(cfg.get("anchors_per_side", 12))

    r108_cens = read_bed(abspath(root, cfg["species"]["genome_R108"]["cen_bed"]))
    a17_cens = read_bed(abspath(root, cfg["species"]["genome_A17"]["cen_bed"]))
    mpo_cens = read_bed(abspath(root, cfg["species"]["genome_Mpo"]["cen_bed"]))
    r108_by_chr, _ = read_gene_coord(abspath(root, cfg["species"]["genome_R108"]["gene_coord"]))
    r108_to_mpo = read_gene_links(abspath(root, cfg["links"]["R108_to_Mpo"]))

    ref_rows = []
    projection_rows = []
    anchor_rows = []

    a17_by_chr = {c["chr"]: c for c in a17_cens}

    for cen in sorted(r108_cens, key=lambda x: int(x["chr"].replace("Chr", "")) if x["chr"].replace("Chr", "").isdigit() else x["chr"]):
        cen_id = "CEN" + cen["chr"].replace("Chr", "")
        left_genes, right_genes = nearest_flank_genes(r108_by_chr, cen["chr"], cen["start"], cen["end"], flank_bp, n_anchor)
        left_ids = [g["gene"] for g in left_genes]
        right_ids = [g["gene"] for g in right_genes]

        ref_rows.append({
            "ancestral_centromere_id": cen_id,
            "reference_species": "genome_R108",
            "reference_chr": cen["chr"],
            "reference_start": cen["start"],
            "reference_end": cen["end"],
            "reference_size_bp": cen["end"] - cen["start"],
            "A17_chr": a17_by_chr.get(cen["chr"], {}).get("chr", ""),
            "A17_start": a17_by_chr.get(cen["chr"], {}).get("start", ""),
            "A17_end": a17_by_chr.get(cen["chr"], {}).get("end", ""),
            "left_anchor_genes": ",".join(left_ids),
            "right_anchor_genes": ",".join(right_ids),
            "left_anchor_count": len(left_ids),
            "right_anchor_count": len(right_ids),
        })

        anchor_records = []
        for side, genes in [("left", left_genes), ("right", right_genes)]:
            for gene in genes:
                gene_links = r108_to_mpo.get(gene["gene"], [])
                anchor_records.append({"side": side, "gene": gene, "links": gene_links, "mapped": bool(gene_links)})
                if gene_links:
                    for link in gene_links:
                        anchor_rows.append({
                            "ancestral_centromere_id": cen_id,
                            "side": side,
                            "R108_gene": gene["gene"],
                            "R108_chr": gene["chr"],
                            "R108_start": gene["start"],
                            "R108_end": gene["end"],
                            "Mpo_gene": link["dst_gene"],
                            "Mpo_chr": link["dst_chr"],
                            "Mpo_start": link["dst_start"],
                            "Mpo_end": link["dst_end"],
                            "block_id": link["block_id"],
                        })
                else:
                    anchor_rows.append({
                        "ancestral_centromere_id": cen_id,
                        "side": side,
                        "R108_gene": gene["gene"],
                        "R108_chr": gene["chr"],
                        "R108_start": gene["start"],
                        "R108_end": gene["end"],
                        "Mpo_gene": "",
                        "Mpo_chr": "",
                        "Mpo_start": "",
                        "Mpo_end": "",
                        "block_id": "",
                    })

        proj = summarize_projection(anchor_records)
        row = {
            "ancestral_centromere_id": cen_id,
            "R108_chr": cen["chr"],
            "R108_start": cen["start"],
            "R108_end": cen["end"],
            "flank_bp": flank_bp,
            "left_anchor_count": len(left_ids),
            "right_anchor_count": len(right_ids),
            "mapped_anchor_count": 0,
            "Mpo_chr": "",
            "Mpo_projected_start": "",
            "Mpo_projected_end": "",
            "candidate_chr_count": 0,
            "candidate_chr_summary": "",
            "overlaps_Mpo_CENH3_core": "no",
            "Mpo_CENH3_core": "",
            "Mpo_CENH3_overlap_bp": 0,
            "classification": "deleted_or_unresolved",
            "interpretation": "No flanking anchors could be confidently projected to Mpo.",
        }
        if proj:
            row.update({
                "mapped_anchor_count": proj["mapped_anchor_count"],
                "Mpo_chr": proj["chr"],
                "Mpo_projected_start": proj["start"],
                "Mpo_projected_end": proj["end"],
                "candidate_chr_count": proj["candidate_chr_count"],
                "candidate_chr_summary": proj["candidate_chr_summary"],
            })
            hits = active_overlap(proj["chr"], proj["start"], proj["end"], mpo_cens)
            if hits:
                hit, ov = hits[0]
                row.update({
                    "overlaps_Mpo_CENH3_core": "yes",
                    "Mpo_CENH3_core": f"{hit['chr']}:{hit['start']}-{hit['end']}",
                    "Mpo_CENH3_overlap_bp": ov,
                    "classification": "conserved_active_centromere",
                    "interpretation": "Projected flanking anchors span an active Mpo CENH3-defined centromere.",
                })
            else:
                row.update({
                    "classification": "projected_without_active_CENH3",
                    "interpretation": "Projected flanking anchors define an Mpo interval that does not overlap the selected active CENH3 core; repeat and coverage controls are required before calling an inactive/paleo-centromere.",
                })
            if proj["candidate_chr_count"] > 1:
                row["interpretation"] += " Multiple Mpo chromosomes carry projected anchors, suggesting rearrangement or anchor ambiguity."
        projection_rows.append(row)

    write_tsv(
        os.path.join(run_root, "01_define_reference_CEN1_CEN8/results/R108_reference_CEN1_CEN8.tsv"),
        ref_rows,
        ["ancestral_centromere_id", "reference_species", "reference_chr", "reference_start", "reference_end", "reference_size_bp", "A17_chr", "A17_start", "A17_end", "left_anchor_genes", "right_anchor_genes", "left_anchor_count", "right_anchor_count"],
    )
    fields = ["ancestral_centromere_id", "R108_chr", "R108_start", "R108_end", "flank_bp", "left_anchor_count", "right_anchor_count", "mapped_anchor_count", "Mpo_chr", "Mpo_projected_start", "Mpo_projected_end", "candidate_chr_count", "candidate_chr_summary", "overlaps_Mpo_CENH3_core", "Mpo_CENH3_core", "Mpo_CENH3_overlap_bp", "classification", "interpretation"]
    write_tsv(os.path.join(run_root, "02_project_A17_R108_CENs_to_Mpo/results/R108_CEN1_CEN8_projected_to_Mpo.tsv"), projection_rows, fields)
    write_tsv(os.path.join(run_root, "02_project_A17_R108_CENs_to_Mpo/results/Mpo_projected_ancestral_centromere_fate.tsv"), projection_rows, fields)
    write_tsv(
        os.path.join(run_root, "02_project_A17_R108_CENs_to_Mpo/results/R108_CEN_flanking_anchor_projection_detail.tsv"),
        anchor_rows,
        ["ancestral_centromere_id", "side", "R108_gene", "R108_chr", "R108_start", "R108_end", "Mpo_gene", "Mpo_chr", "Mpo_start", "Mpo_end", "block_id"],
    )

    print(f"Wrote {len(ref_rows)} reference centromeres")
    print(f"Wrote {len(projection_rows)} projected centromere rows")


if __name__ == "__main__":
    main()
