#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict
from pathlib import Path


THRESHOLDS = [
    ("permissive", 500, 0.70),
    ("standard", 1000, 0.75),
    ("stringent", 2000, 0.85),
]
REFERENCE_ORDER = ["genome_R108", "genome_A17", "genome_Msa", "genome_474"]
DISPLAY = {
    "genome_R108": "Mtru_R108",
    "genome_A17": "Mtru_A17",
    "genome_Msa": "Msat_cae",
    "genome_474": "Mcar",
}
RELIC_ORDER = ["Mpo_Chr3_relic", "Mpo_Chr5_relic"]


def read_tsv(path):
    with Path(path).open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def interval_distance(a0, a1, b0, b1):
    if overlap(a0, a1, b0, b1) > 0:
        return 0
    if a1 <= b0:
        return b0 - a1
    return a0 - b1


def merge_intervals(intervals, gap=0):
    if not intervals:
        return []
    ordered = sorted((min(a, b), max(a, b)) for a, b in intervals)
    merged = [list(ordered[0])]
    for start, end in ordered[1:]:
        if start <= merged[-1][1] + gap:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(a, b) for a, b in merged]


def union_bp(intervals):
    return sum(end - start for start, end in merge_intervals(intervals))


def overlap_union_bp(intervals, feature_start, feature_end):
    clipped = []
    for start, end in intervals:
        left = max(start, feature_start)
        right = min(end, feature_end)
        if right > left:
            clipped.append((left, right))
    return union_bp(clipped)


def parse_tags(parts):
    tags = {}
    for item in parts[12:]:
        fields = item.split(":", 2)
        if len(fields) == 3:
            tags[fields[0]] = fields[2]
    return tags


def parse_paf(path, direction, target_species):
    rows = []
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            tags = parse_tags(parts)
            qlen = int(parts[1])
            nmatch = int(parts[9])
            alen = int(parts[10])
            rows.append(
                {
                    "direction": direction,
                    "target_species": target_species,
                    "query_id": parts[0],
                    "query_length_bp": qlen,
                    "query_start_bp": int(parts[2]),
                    "query_end_bp": int(parts[3]),
                    "strand": parts[4],
                    "target_chromosome": parts[5],
                    "target_length_bp": int(parts[6]),
                    "target_start_bp": int(parts[7]),
                    "target_end_bp": int(parts[8]),
                    "matching_bp": nmatch,
                    "alignment_block_bp": alen,
                    "identity": nmatch / max(1, alen),
                    "mapq": int(parts[11]),
                    "alignment_type": tags.get("tp", ""),
                    "source_paf": Path(path).name,
                    "source_line": line_number,
                }
            )
    return rows


def cluster_alignments(rows, gap_bp, threshold_name, min_alignment_bp, min_identity):
    kept = [
        row
        for row in rows
        if row["alignment_block_bp"] >= min_alignment_bp
        and row["identity"] >= min_identity
    ]
    grouped = defaultdict(list)
    for row in kept:
        grouped[
            (
                row["direction"],
                row["target_species"],
                row["query_id"],
                row["target_chromosome"],
            )
        ].append(row)

    clusters = []
    cluster_members = []
    for key, group in grouped.items():
        group = sorted(
            group,
            key=lambda row: (
                row["target_start_bp"],
                row["target_end_bp"],
                -row["alignment_block_bp"],
            ),
        )
        current = []
        current_end = None
        grouped_clusters = []
        for row in group:
            if current and row["target_start_bp"] > current_end + gap_bp:
                grouped_clusters.append(current)
                current = []
                current_end = None
            current.append(row)
            current_end = max(current_end or row["target_end_bp"], row["target_end_bp"])
        if current:
            grouped_clusters.append(current)

        for index, members in enumerate(grouped_clusters, 1):
            target_intervals = [
                (row["target_start_bp"], row["target_end_bp"]) for row in members
            ]
            query_intervals = [
                (row["query_start_bp"], row["query_end_bp"]) for row in members
            ]
            total_alignment = sum(row["alignment_block_bp"] for row in members)
            total_matches = sum(row["matching_bp"] for row in members)
            query_length = members[0]["query_length_bp"]
            cluster_id = "{}|{}|{}|{}|{}|{}".format(
                threshold_name, key[0], key[1], key[2], key[3], index
            )
            cluster = {
                "threshold": threshold_name,
                "min_alignment_bp": min_alignment_bp,
                "min_identity": min_identity,
                "cluster_id": cluster_id,
                "direction": key[0],
                "target_species": key[1],
                "query_id": key[2],
                "target_chromosome": key[3],
                "target_start_bp": min(start for start, _ in target_intervals),
                "target_end_bp": max(end for _, end in target_intervals),
                "target_span_bp": max(end for _, end in target_intervals)
                - min(start for start, _ in target_intervals),
                "target_union_bp": union_bp(target_intervals),
                "query_length_bp": query_length,
                "query_union_bp": union_bp(query_intervals),
                "query_coverage_pct": union_bp(query_intervals)
                / max(1, query_length)
                * 100.0,
                "alignment_count": len(members),
                "primary_alignment_count": sum(
                    row["alignment_type"] == "P" for row in members
                ),
                "alignment_sum_bp": total_alignment,
                "weighted_identity": total_matches / max(1, total_alignment),
                "max_identity": max(row["identity"] for row in members),
                "max_mapq": max(row["mapq"] for row in members),
                "_target_intervals": target_intervals,
            }
            clusters.append(cluster)
            for row in members:
                member = dict(row)
                member["threshold"] = threshold_name
                member["cluster_id"] = cluster_id
                cluster_members.append(member)
    return clusters, cluster_members


def feature_overlap(cluster, feature):
    if cluster["target_chromosome"] != feature["chromosome"]:
        return 0
    return overlap_union_bp(
        cluster["_target_intervals"], feature["start_bp"], feature["end_bp"]
    )


def feature_distance(cluster, feature):
    if cluster["target_chromosome"] != feature["chromosome"]:
        return ""
    return min(
        interval_distance(start, end, feature["start_bp"], feature["end_bp"])
        for start, end in cluster["_target_intervals"]
    )


def choose_cluster(clusters, feature):
    supporting = []
    for cluster in clusters:
        ov = feature_overlap(cluster, feature)
        if ov > 0:
            supporting.append((ov, cluster))
    if supporting:
        supporting.sort(
            key=lambda item: (
                item[0],
                item[1]["query_union_bp"],
                item[1]["weighted_identity"],
            ),
            reverse=True,
        )
        return supporting[0][1], supporting[0][0], "overlaps_target_feature"
    if not clusters:
        return None, 0, "no_alignment_cluster"
    best = max(
        clusters,
        key=lambda cluster: (
            cluster["query_union_bp"] * cluster["weighted_identity"],
            cluster["target_union_bp"],
            cluster["max_identity"],
        ),
    )
    return best, 0, "best_genome_wide_non_overlapping_cluster"


def genome_wide_rank(clusters, selected):
    if selected is None:
        return ""
    ranked = sorted(
        clusters,
        key=lambda cluster: (
            cluster["query_union_bp"] * cluster["weighted_identity"],
            cluster["target_union_bp"],
            cluster["max_identity"],
        ),
        reverse=True,
    )
    for rank, cluster in enumerate(ranked, 1):
        if cluster["cluster_id"] == selected["cluster_id"]:
            return rank
    return ""


def public_cluster(cluster):
    if cluster is None:
        return {}
    return {key: value for key, value in cluster.items() if not key.startswith("_")}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True)
    args = parser.parse_args()
    run = Path(args.run)
    results = run / "results"

    manifest_rows = read_tsv(results / "input_intervals.tsv")
    features = {}
    for row in manifest_rows:
        row["start_bp"] = int(row["start_bp"])
        row["end_bp"] = int(row["end_bp"])
        features[row["feature_id"]] = row

    all_alignments = []
    all_alignments.extend(
        parse_paf(
            run / "paf" / "x8_CEN5_cores_to_Mpo.asm20.paf",
            "reference_CEN5_to_Mpol",
            "genome_Mpo",
        )
    )
    for species in REFERENCE_ORDER:
        all_alignments.extend(
            parse_paf(
                run / "paf" / ("Mpo_relics_to_{}.asm20.paf".format(species)),
                "Mpol_relic_to_reference",
                species,
            )
        )

    alignment_fields = [
        "direction",
        "target_species",
        "query_id",
        "query_length_bp",
        "query_start_bp",
        "query_end_bp",
        "strand",
        "target_chromosome",
        "target_length_bp",
        "target_start_bp",
        "target_end_bp",
        "matching_bp",
        "alignment_block_bp",
        "identity",
        "mapq",
        "alignment_type",
        "source_paf",
        "source_line",
    ]
    write_tsv(results / "all_raw_alignments.tsv", all_alignments, alignment_fields)

    all_clusters = []
    all_cluster_members = []
    threshold_clusters = {}
    threshold_members = {}
    for threshold_name, min_alignment, min_identity in THRESHOLDS:
        clusters, members = cluster_alignments(
            all_alignments,
            100000,
            threshold_name,
            min_alignment,
            min_identity,
        )
        threshold_clusters[threshold_name] = clusters
        threshold_members[threshold_name] = members
        all_clusters.extend(clusters)
        all_cluster_members.extend(members)

    cluster_fields = [
        "threshold",
        "min_alignment_bp",
        "min_identity",
        "cluster_id",
        "direction",
        "target_species",
        "query_id",
        "target_chromosome",
        "target_start_bp",
        "target_end_bp",
        "target_span_bp",
        "target_union_bp",
        "query_length_bp",
        "query_union_bp",
        "query_coverage_pct",
        "alignment_count",
        "primary_alignment_count",
        "alignment_sum_bp",
        "weighted_identity",
        "max_identity",
        "max_mapq",
    ]
    write_tsv(
        results / "all_alignment_clusters.tsv",
        [public_cluster(row) for row in all_clusters],
        cluster_fields,
    )

    sensitivity = []
    selected_by_threshold = {}
    for threshold_name, _, _ in THRESHOLDS:
        clusters = threshold_clusters[threshold_name]
        for species in REFERENCE_ORDER:
            forward_query = "{}_CEN5_core".format(species)
            forward_pool = [
                cluster
                for cluster in clusters
                if cluster["direction"] == "reference_CEN5_to_Mpol"
                and cluster["query_id"] == forward_query
                and cluster["target_species"] == "genome_Mpo"
            ]
            active_cen5 = features[forward_query]
            for relic_id in RELIC_ORDER:
                relic = features[relic_id]
                reverse_pool = [
                    cluster
                    for cluster in clusters
                    if cluster["direction"] == "Mpol_relic_to_reference"
                    and cluster["query_id"] == relic_id
                    and cluster["target_species"] == species
                ]
                f_cluster, f_overlap, f_status = choose_cluster(forward_pool, relic)
                r_cluster, r_overlap, r_status = choose_cluster(
                    reverse_pool, active_cen5
                )
                if f_overlap > 0 and r_overlap > 0:
                    reciprocal = "bidirectional_support"
                elif f_overlap > 0:
                    reciprocal = "forward_only"
                elif r_overlap > 0:
                    reciprocal = "reverse_only"
                else:
                    reciprocal = "no_CEN5_overlap"
                key = (threshold_name, species, relic_id)
                selected_by_threshold[key] = (f_cluster, r_cluster)
                row = {
                    "threshold": threshold_name,
                    "reference_species": species,
                    "display_species": DISPLAY[species],
                    "relic_id": relic_id,
                    "forward_status": f_status,
                    "forward_overlap_Mpol_relic_bp": f_overlap,
                    "reverse_status": r_status,
                    "reverse_overlap_reference_CEN5_bp": r_overlap,
                    "reciprocal_support": reciprocal,
                    "forward_genome_wide_rank": genome_wide_rank(
                        forward_pool, f_cluster
                    ),
                    "reverse_genome_wide_rank": genome_wide_rank(
                        reverse_pool, r_cluster
                    ),
                }
                for prefix, cluster in [
                    ("forward", f_cluster),
                    ("reverse", r_cluster),
                ]:
                    if cluster is None:
                        continue
                    row[prefix + "_cluster_id"] = cluster["cluster_id"]
                    row[prefix + "_target_chromosome"] = cluster[
                        "target_chromosome"
                    ]
                    row[prefix + "_target_start_bp"] = cluster["target_start_bp"]
                    row[prefix + "_target_end_bp"] = cluster["target_end_bp"]
                    row[prefix + "_target_union_bp"] = cluster["target_union_bp"]
                    row[prefix + "_query_union_bp"] = cluster["query_union_bp"]
                    row[prefix + "_query_coverage_pct"] = round(
                        cluster["query_coverage_pct"], 5
                    )
                    row[prefix + "_weighted_identity"] = round(
                        cluster["weighted_identity"], 6
                    )
                    row[prefix + "_max_identity"] = round(
                        cluster["max_identity"], 6
                    )
                    row[prefix + "_alignment_count"] = cluster["alignment_count"]
                sensitivity.append(row)

    sensitivity_fields = [
        "threshold",
        "reference_species",
        "display_species",
        "relic_id",
        "forward_status",
        "forward_overlap_Mpol_relic_bp",
        "forward_genome_wide_rank",
        "forward_cluster_id",
        "forward_target_chromosome",
        "forward_target_start_bp",
        "forward_target_end_bp",
        "forward_target_union_bp",
        "forward_query_union_bp",
        "forward_query_coverage_pct",
        "forward_weighted_identity",
        "forward_max_identity",
        "forward_alignment_count",
        "reverse_status",
        "reverse_overlap_reference_CEN5_bp",
        "reverse_genome_wide_rank",
        "reverse_cluster_id",
        "reverse_target_chromosome",
        "reverse_target_start_bp",
        "reverse_target_end_bp",
        "reverse_target_union_bp",
        "reverse_query_union_bp",
        "reverse_query_coverage_pct",
        "reverse_weighted_identity",
        "reverse_max_identity",
        "reverse_alignment_count",
        "reciprocal_support",
    ]
    write_tsv(
        results / "reciprocal_support_threshold_sensitivity.tsv",
        sensitivity,
        sensitivity_fields,
    )

    standard_rows = [
        row for row in sensitivity if row["threshold"] == "standard"
    ]
    for row in standard_rows:
        species = row["reference_species"]
        relic_id = row["relic_id"]
        robust = sum(
            other["reciprocal_support"] == "bidirectional_support"
            for other in sensitivity
            if other["reference_species"] == species
            and other["relic_id"] == relic_id
        )
        row["bidirectional_threshold_count"] = robust
        reverse_cluster = selected_by_threshold[
            ("standard", species, relic_id)
        ][1]
        active_cen5 = features["{}_CEN5_core".format(species)]
        row["reference_CEN5_chromosome"] = active_cen5["chromosome"]
        row["reference_CEN5_start_bp"] = active_cen5["start_bp"]
        row["reference_CEN5_end_bp"] = active_cen5["end_bp"]
        relic = features[relic_id]
        row["Mpol_relic_chromosome"] = relic["chromosome"]
        row["Mpol_relic_start_bp"] = relic["start_bp"]
        row["Mpol_relic_end_bp"] = relic["end_bp"]
        row["reverse_distance_to_reference_CEN5_bp"] = (
            feature_distance(reverse_cluster, active_cen5)
            if reverse_cluster is not None
            else ""
        )

    standard_fields = [
        "reference_species",
        "display_species",
        "relic_id",
        "Mpol_relic_chromosome",
        "Mpol_relic_start_bp",
        "Mpol_relic_end_bp",
        "reference_CEN5_chromosome",
        "reference_CEN5_start_bp",
        "reference_CEN5_end_bp",
        "forward_status",
        "forward_overlap_Mpol_relic_bp",
        "forward_genome_wide_rank",
        "forward_target_chromosome",
        "forward_target_start_bp",
        "forward_target_end_bp",
        "forward_target_union_bp",
        "forward_query_union_bp",
        "forward_query_coverage_pct",
        "forward_weighted_identity",
        "forward_max_identity",
        "forward_alignment_count",
        "reverse_status",
        "reverse_overlap_reference_CEN5_bp",
        "reverse_genome_wide_rank",
        "reverse_distance_to_reference_CEN5_bp",
        "reverse_target_chromosome",
        "reverse_target_start_bp",
        "reverse_target_end_bp",
        "reverse_target_union_bp",
        "reverse_query_union_bp",
        "reverse_query_coverage_pct",
        "reverse_weighted_identity",
        "reverse_max_identity",
        "reverse_alignment_count",
        "reciprocal_support",
        "bidirectional_threshold_count",
    ]
    write_tsv(
        results / "reciprocal_CEN5_relic_locations.standard.tsv",
        standard_rows,
        standard_fields,
    )

    global_best_rows = []
    standard_clusters = threshold_clusters["standard"]
    mpo_active = [
        feature
        for feature in features.values()
        if feature["species"] == "genome_Mpo"
        and feature["feature_type"] == "active_CENH3"
    ]
    for species in REFERENCE_ORDER:
        forward_query = "{}_CEN5_core".format(species)
        forward_pool = [
            cluster
            for cluster in standard_clusters
            if cluster["direction"] == "reference_CEN5_to_Mpol"
            and cluster["query_id"] == forward_query
            and cluster["target_species"] == "genome_Mpo"
        ]
        if forward_pool:
            best = max(
                forward_pool,
                key=lambda cluster: (
                    cluster["query_union_bp"]
                    * cluster["weighted_identity"],
                    cluster["target_union_bp"],
                    cluster["max_identity"],
                ),
            )
            row = {
                "direction": "reference_CEN5_to_Mpol",
                "reference_species": species,
                "display_species": DISPLAY[species],
                "query_id": forward_query,
                "target_species": "genome_Mpo",
                "cluster_id": best["cluster_id"],
                "target_chromosome": best["target_chromosome"],
                "target_start_bp": best["target_start_bp"],
                "target_end_bp": best["target_end_bp"],
                "target_union_bp": best["target_union_bp"],
                "query_union_bp": best["query_union_bp"],
                "query_coverage_pct": round(best["query_coverage_pct"], 5),
                "weighted_identity": round(best["weighted_identity"], 6),
                "max_identity": round(best["max_identity"], 6),
                "alignment_count": best["alignment_count"],
                "overlap_Mpo_Chr3_relic_bp": feature_overlap(
                    best, features["Mpo_Chr3_relic"]
                ),
                "overlap_Mpo_Chr5_relic_bp": feature_overlap(
                    best, features["Mpo_Chr5_relic"]
                ),
                "overlap_target_active_CENH3_bp": sum(
                    feature_overlap(best, feature) for feature in mpo_active
                ),
            }
            global_best_rows.append(row)
        for relic_id in RELIC_ORDER:
            reverse_pool = [
                cluster
                for cluster in standard_clusters
                if cluster["direction"] == "Mpol_relic_to_reference"
                and cluster["query_id"] == relic_id
                and cluster["target_species"] == species
            ]
            if not reverse_pool:
                continue
            best = max(
                reverse_pool,
                key=lambda cluster: (
                    cluster["query_union_bp"]
                    * cluster["weighted_identity"],
                    cluster["target_union_bp"],
                    cluster["max_identity"],
                ),
            )
            active_cen5 = features[forward_query]
            global_best_rows.append(
                {
                    "direction": "Mpol_relic_to_reference",
                    "reference_species": species,
                    "display_species": DISPLAY[species],
                    "query_id": relic_id,
                    "target_species": species,
                    "cluster_id": best["cluster_id"],
                    "target_chromosome": best["target_chromosome"],
                    "target_start_bp": best["target_start_bp"],
                    "target_end_bp": best["target_end_bp"],
                    "target_union_bp": best["target_union_bp"],
                    "query_union_bp": best["query_union_bp"],
                    "query_coverage_pct": round(best["query_coverage_pct"], 5),
                    "weighted_identity": round(best["weighted_identity"], 6),
                    "max_identity": round(best["max_identity"], 6),
                    "alignment_count": best["alignment_count"],
                    "overlap_reference_CEN5_bp": feature_overlap(
                        best, active_cen5
                    ),
                    "distance_to_reference_CEN5_bp": feature_distance(
                        best, active_cen5
                    ),
                }
            )
    write_tsv(
        results / "genome_wide_best_clusters.standard.tsv",
        global_best_rows,
        [
            "direction",
            "reference_species",
            "display_species",
            "query_id",
            "target_species",
            "cluster_id",
            "target_chromosome",
            "target_start_bp",
            "target_end_bp",
            "target_union_bp",
            "query_union_bp",
            "query_coverage_pct",
            "weighted_identity",
            "max_identity",
            "alignment_count",
            "overlap_Mpo_Chr3_relic_bp",
            "overlap_Mpo_Chr5_relic_bp",
            "overlap_target_active_CENH3_bp",
            "overlap_reference_CEN5_bp",
            "distance_to_reference_CEN5_bp",
        ],
    )

    selected_cluster_ids = set()
    for row in standard_rows:
        for direction in ["forward", "reverse"]:
            cluster_id = row.get(direction + "_cluster_id", "")
            if cluster_id:
                selected_cluster_ids.add(cluster_id)
    selected_members = [
        row
        for row in threshold_members["standard"]
        if row["cluster_id"] in selected_cluster_ids
    ]
    write_tsv(
        results / "selected_cluster_alignment_blocks.standard.tsv",
        selected_members,
        ["threshold", "cluster_id"] + alignment_fields,
    )

    context_alignments = []
    for species in REFERENCE_ORDER:
        context_alignments.extend(
            parse_paf(
                run
                / "paf"
                / ("Mpo_relic_context_pm2Mb_to_{}.asm20.paf".format(species)),
                "Mpol_relic_context_to_reference",
                species,
            )
        )
    context_clusters, _ = cluster_alignments(
        context_alignments, 500000, "context", 5000, 0.80
    )
    context_rows = []
    for species in REFERENCE_ORDER:
        active_cen5 = features["{}_CEN5_core".format(species)]
        for relic_id in RELIC_ORDER:
            query_id = relic_id + "_pm2Mb"
            pool = [
                cluster
                for cluster in context_clusters
                if cluster["target_species"] == species
                and cluster["query_id"] == query_id
            ]
            if pool:
                chosen = max(
                    pool,
                    key=lambda cluster: (
                        cluster["query_union_bp"]
                        * cluster["weighted_identity"],
                        cluster["target_union_bp"],
                        cluster["max_identity"],
                    ),
                )
                ov = feature_overlap(chosen, active_cen5)
                status = (
                    "best_context_cluster_overlaps_reference_CEN5"
                    if ov > 0
                    else "best_context_cluster_non_overlapping"
                )
            else:
                chosen = None
                ov = 0
                status = "no_alignment_cluster"
            row = {
                "reference_species": species,
                "display_species": DISPLAY[species],
                "relic_id": relic_id,
                "context_status": status,
                "context_overlap_reference_CEN5_bp": ov,
            }
            if chosen is not None:
                row.update(
                    {
                        "context_target_chromosome": chosen["target_chromosome"],
                        "context_target_start_bp": chosen["target_start_bp"],
                        "context_target_end_bp": chosen["target_end_bp"],
                        "context_target_union_bp": chosen["target_union_bp"],
                        "context_query_union_bp": chosen["query_union_bp"],
                        "context_query_coverage_pct": round(
                            chosen["query_coverage_pct"], 5
                        ),
                        "context_weighted_identity": round(
                            chosen["weighted_identity"], 6
                        ),
                        "context_max_identity": round(
                            chosen["max_identity"], 6
                        ),
                        "context_alignment_count": chosen["alignment_count"],
                        "context_distance_to_reference_CEN5_bp": feature_distance(
                            chosen, active_cen5
                        ),
                    }
                )
            context_rows.append(row)
    write_tsv(
        results / "relic_pm2Mb_context_locations.tsv",
        context_rows,
        [
            "reference_species",
            "display_species",
            "relic_id",
            "context_status",
            "context_overlap_reference_CEN5_bp",
            "context_distance_to_reference_CEN5_bp",
            "context_target_chromosome",
            "context_target_start_bp",
            "context_target_end_bp",
            "context_target_union_bp",
            "context_query_union_bp",
            "context_query_coverage_pct",
            "context_weighted_identity",
            "context_max_identity",
            "context_alignment_count",
        ],
    )

    qc_rows = [
        {
            "metric": "raw_alignment_rows",
            "value": len(all_alignments),
        },
        {
            "metric": "standard_clusters",
            "value": len(threshold_clusters["standard"]),
        },
        {
            "metric": "standard_reference_relic_pairs",
            "value": len(standard_rows),
        },
        {
            "metric": "standard_bidirectional_pairs",
            "value": sum(
                row["reciprocal_support"] == "bidirectional_support"
                for row in standard_rows
            ),
        },
        {
            "metric": "standard_forward_only_pairs",
            "value": sum(
                row["reciprocal_support"] == "forward_only"
                for row in standard_rows
            ),
        },
        {
            "metric": "standard_reverse_only_pairs",
            "value": sum(
                row["reciprocal_support"] == "reverse_only"
                for row in standard_rows
            ),
        },
        {
            "metric": "standard_no_CEN5_overlap_pairs",
            "value": sum(
                row["reciprocal_support"] == "no_CEN5_overlap"
                for row in standard_rows
            ),
        },
    ]
    write_tsv(results / "QC_summary.tsv", qc_rows, ["metric", "value"])


if __name__ == "__main__":
    main()
