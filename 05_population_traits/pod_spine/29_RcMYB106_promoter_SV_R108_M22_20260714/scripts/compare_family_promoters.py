#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze_rcmyb106_promoter as core


PAIRS = [
    ("closest_broad_coverage", "R10827621.1", "Chr24685.1", 74, 74, False),
    ("secondary_broad_coverage", "R10819971.1", "Chr17985.1", 60, 98, False),
    ("domain_only_related_3", "R10830912.1", "Chr27578.1", 45, 45, False),
    ("domain_only_Rc10_RBH", "R10824448.1", "Chr21982.1", 40, 37, True),
]


def merge_coverage(intervals, region_start, region_end):
    clipped = []
    for start, end in intervals:
        start = max(start, region_start)
        end = min(end, region_end)
        if end > start:
            clipped.append((start, end))
    if not clipped:
        return 0
    clipped.sort()
    total = 0
    current_start, current_end = clipped[0]
    for start, end in clipped[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            total += current_end - current_start
            current_start, current_end = start, end
    return total + current_end - current_start


def gene_id(records, protein_id):
    header = records[protein_id][0]
    return core.gene_from_header(header, protein_id)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--r108-genome", required=True)
    parser.add_argument("--r108-gff", required=True)
    parser.add_argument("--r108-pep", required=True)
    parser.add_argument("--m22-genome", required=True)
    parser.add_argument("--m22-gff", required=True)
    parser.add_argument("--m22-pep", required=True)
    parser.add_argument("--minimap2", required=True)
    args = parser.parse_args()

    out = Path(args.outdir)
    promoter_dir = out / "02_promoter" / "family_screen"
    alignment_dir = out / "03_alignment" / "family_screen"
    summary_dir = out / "04_summary"
    log_path = out / "logs" / "family_screen.log"
    for directory in (promoter_dir, alignment_dir, summary_dir):
        directory.mkdir(parents=True, exist_ok=True)

    r108_records = core.read_fasta(args.r108_pep)
    m22_records = core.read_fasta(args.m22_pep)
    rows = []
    event_rows = []
    coordinate_rows = []

    for label, r108_protein, m22_protein, r108_rc_qcov, m22_rc_qcov, rc10_rbh in PAIRS:
        r108_gene_id = gene_id(r108_records, r108_protein)
        m22_gene_id = gene_id(m22_records, m22_protein)
        r108_gene = core.parse_simple_gff_gene(args.r108_gff, r108_gene_id)
        m22_gene = core.parse_simple_gff_gene(args.m22_gff, m22_gene_id)
        pair_promoter_dir = promoter_dir / label
        pair_promoter_dir.mkdir(parents=True, exist_ok=True)
        r108_windows = core.extract_windows(Path(args.r108_genome), r108_gene, pair_promoter_dir, "R108")
        m22_windows = core.extract_windows(Path(args.m22_genome), m22_gene, pair_promoter_dir, "M22")

        for species, gene, windows in (("R108", r108_gene, r108_windows), ("M22", m22_gene, m22_windows)):
            info = windows["promoter_10kb"]
            coordinate_rows.append({
                "pair": label, "species": species, "gene_id": gene["gene_id"],
                "chrom": gene["seqid"], "gene_start": gene["start"], "gene_end": gene["end"],
                "strand": gene["strand"], "promoter_start": info["coords"][0],
                "promoter_end": info["coords"][1], "promoter_length": info["length"],
            })

        paf = alignment_dir / (label + ".R108_vs_M22.promoter_10kb.asm20.paf")
        core.run([
            args.minimap2, "-x", "asm20", "-c", "--cs=long", "-N", "20",
            str(r108_windows["promoter_10kb"]["path"]),
            str(m22_windows["promoter_10kb"]["path"]),
        ], log_path, paf)
        alignments, events = core.parse_paf(paf, "R108", "M22", 10000)
        plus = [aln for aln in alignments if aln["strand"] == "+"]
        target_intervals = [(aln["tstart"], aln["tend"]) for aln in plus]
        query_intervals = [(aln["qstart"], aln["qend"]) for aln in plus]
        target_coverage = merge_coverage(target_intervals, 0, 10000)
        query_coverage = merge_coverage(query_intervals, 0, 10000)
        analogous_start = 10000 - 4607
        analogous_end = 10000 - 254
        analogous_coverage = merge_coverage(target_intervals, analogous_start, analogous_end)
        deletions = [e for e in events if e["event"] == "M22_deletion_relative_to_R108"]
        large_deletions = [e for e in deletions if e["length"] >= 1000]
        overlapping_large = [
            e for e in large_deletions
            if e["target_end0"] > analogous_start and e["target_start0"] < analogous_end
        ]
        max_deletion = max([e["length"] for e in deletions] or [0])
        rows.append({
            "pair": label, "r108_gene": r108_gene_id, "r108_protein": r108_protein,
            "m22_gene": m22_gene_id, "m22_protein": m22_protein,
            "r108_Rc10_query_coverage_pct": r108_rc_qcov,
            "m22_Rc10_query_coverage_pct": m22_rc_qcov,
            "strict_Rc10_reciprocal_best_hit": rc10_rbh,
            "alignment_count": len(plus), "r108_promoter_coverage_bp": target_coverage,
            "r108_promoter_coverage_pct": round(target_coverage / 100.0, 2),
            "m22_promoter_coverage_bp": query_coverage,
            "m22_promoter_coverage_pct": round(query_coverage / 100.0, 2),
            "analogous_interval_coverage_bp": analogous_coverage,
            "analogous_interval_coverage_pct": round(100.0 * analogous_coverage / (analogous_end - analogous_start), 2),
            "max_M22_deletion_relative_to_R108_bp": max_deletion,
            "deletions_ge1kb": len(large_deletions),
            "analogous_overlapping_deletions_ge1kb": len(overlapping_large),
        })
        for event in events:
            event = dict(event)
            event["pair"] = label
            event_rows.append(event)

    core.write_tsv(alignment_dir / "MYB106_family_promoter_comparison.tsv", rows, [
        "pair", "r108_gene", "r108_protein", "m22_gene", "m22_protein",
        "r108_Rc10_query_coverage_pct", "m22_Rc10_query_coverage_pct", "strict_Rc10_reciprocal_best_hit",
        "alignment_count",
        "r108_promoter_coverage_bp", "r108_promoter_coverage_pct", "m22_promoter_coverage_bp",
        "m22_promoter_coverage_pct", "analogous_interval_coverage_bp", "analogous_interval_coverage_pct",
        "max_M22_deletion_relative_to_R108_bp", "deletions_ge1kb", "analogous_overlapping_deletions_ge1kb",
    ])
    core.write_tsv(alignment_dir / "MYB106_family_promoter_indels_ge50bp.tsv", event_rows, [
        "pair", "source", "event", "length", "target_start0", "target_end0", "query_start0",
        "query_end0", "target_upstream_start", "target_upstream_end",
    ])
    core.write_tsv(promoter_dir / "MYB106_family_promoter_coordinates.tsv", coordinate_rows, [
        "pair", "species", "gene_id", "chrom", "gene_start", "gene_end", "strand",
        "promoter_start", "promoter_end", "promoter_length",
    ])

    lines = [
        "# MYB106 family promoter screen",
        "",
        "Four R108/M22 candidate pairs were compared with minimap2 `asm20`.",
        "The castor analogous interval is defined as -4607 to -254 bp relative to TSS.",
        "",
        "| Pair | R108 | M22 | Rc10 query coverage (R108/M22) | Analogous coverage | Max M22 deletion | >=1 kb deletion |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            "| {pair} | {r108_gene} | {m22_gene} | {r108_Rc10_query_coverage_pct}%/{m22_Rc10_query_coverage_pct}% | {analogous_interval_coverage_pct}% | "
            "{max_M22_deletion_relative_to_R108_bp} bp | {deletions_ge1kb} |".format(**row)
        )
    (summary_dir / "MYB106_family_promoter_screen.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
