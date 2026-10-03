#!/usr/bin/env python3
import csv
import os
from collections import defaultdict

RUN_ROOT = "path/to/project/10.centromere_analysis/output_key"
DETAIL = os.path.join(RUN_ROOT, "02_project_A17_R108_CENs_to_Mpo/results/R108_CEN_sequence_block_projection_detail.tsv")
OUT = os.path.join(RUN_ROOT, "05_fusion_chr_CEN_fate/R108_CEN_projection_block_clusters.tsv")
MERGE_GAP = int(os.environ.get("MERGE_GAP_BP", "100000"))


def main():
    groups = defaultdict(list)
    with open(DETAIL, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            groups[(row["ancestral_centromere_id"], row["target_chr"])].append({
                "target_start": int(row["target_start"]),
                "target_end": int(row["target_end"]),
                "query_chr": row["query_chr"],
                "query_start": int(row["query_start"]),
                "query_end": int(row["query_end"]),
                "query_len": int(row["query_len"]),
            })

    rows = []
    for (cen_id, chrom), blocks in groups.items():
        blocks.sort(key=lambda x: (x["target_start"], x["target_end"]))
        current = None
        for b in blocks:
            if current is None or b["target_start"] > current["cluster_end"] + MERGE_GAP:
                if current is not None:
                    rows.append(current)
                current = {
                    "ancestral_centromere_id": cen_id,
                    "Mpo_chr": chrom,
                    "cluster_start": b["target_start"],
                    "cluster_end": b["target_end"],
                    "block_count": 1,
                    "sum_query_len": b["query_len"],
                    "R108_query_chr": b["query_chr"],
                    "R108_query_min": b["query_start"],
                    "R108_query_max": b["query_end"],
                }
            else:
                current["cluster_end"] = max(current["cluster_end"], b["target_end"])
                current["block_count"] += 1
                current["sum_query_len"] += b["query_len"]
                current["R108_query_min"] = min(current["R108_query_min"], b["query_start"])
                current["R108_query_max"] = max(current["R108_query_max"], b["query_end"])
        if current is not None:
            rows.append(current)

    rows.sort(key=lambda r: (r["ancestral_centromere_id"], -r["sum_query_len"]))
    fields = ["ancestral_centromere_id", "Mpo_chr", "cluster_start", "cluster_end", "block_count", "sum_query_len", "R108_query_chr", "R108_query_min", "R108_query_max"]
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(OUT)


if __name__ == "__main__":
    main()
