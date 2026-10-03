#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path("path/to/project/10.centromere_analysis/output_key/10_jcvi_mpo_synteny")


def load_gene_chr(bed: Path) -> dict[str, tuple[str, int, int]]:
    out: dict[str, tuple[str, int, int]] = {}
    with bed.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, gene, *_ = line.rstrip("\n").split("\t")
            out[gene] = (chrom, int(start), int(end))
    return out


def summarize_pair(pair_dir: Path, qname: str, sname: str, simple_name: str) -> list[dict[str, object]]:
    qmap = load_gene_chr(pair_dir / f"{qname}.bed")
    smap = load_gene_chr(pair_dir / f"{sname}.bed")
    stats: dict[tuple[str, str], dict[str, object]] = {}
    with (pair_dir / simple_name).open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            q1, q2, s1, s2, nanchors, orientation = line.rstrip("\n").split("\t")[:6]
            if q1 not in qmap or q2 not in qmap or s1 not in smap or s2 not in smap:
                continue
            qchr = qmap[q1][0]
            schr = smap[s1][0]
            key = (qchr, schr)
            rec = stats.setdefault(
                key,
                {
                    "comparison": f"{qname}_vs_{sname}",
                    "query_chr": qchr,
                    "target_chr": schr,
                    "block_count": 0,
                    "sum_anchor_genes": 0,
                    "plus_blocks": 0,
                    "minus_blocks": 0,
                    "query_span_start": min(qmap[q1][1], qmap[q2][1]),
                    "query_span_end": max(qmap[q1][2], qmap[q2][2]),
                    "target_span_start": min(smap[s1][1], smap[s2][1]),
                    "target_span_end": max(smap[s1][2], smap[s2][2]),
                },
            )
            rec["block_count"] = int(rec["block_count"]) + 1
            rec["sum_anchor_genes"] = int(rec["sum_anchor_genes"]) + int(nanchors)
            if orientation == "+":
                rec["plus_blocks"] = int(rec["plus_blocks"]) + 1
            else:
                rec["minus_blocks"] = int(rec["minus_blocks"]) + 1
            rec["query_span_start"] = min(int(rec["query_span_start"]), qmap[q1][1], qmap[q2][1])
            rec["query_span_end"] = max(int(rec["query_span_end"]), qmap[q1][2], qmap[q2][2])
            rec["target_span_start"] = min(int(rec["target_span_start"]), smap[s1][1], smap[s2][1])
            rec["target_span_end"] = max(int(rec["target_span_end"]), smap[s1][2], smap[s2][2])
    rows = list(stats.values())
    rows.sort(key=lambda r: (r["comparison"], r["target_chr"], -int(r["sum_anchor_genes"])))
    return rows


def main() -> None:
    global_dir = ROOT / "global_R108_Mpo_Msa"
    rows = []
    rows += summarize_pair(global_dir, "genome_R108", "genome_Mpo", "genome_R108.genome_Mpo.anchors.simple")
    rows += summarize_pair(global_dir, "genome_Msa", "genome_Mpo", "genome_Msa.genome_Mpo.anchors.simple")
    a17_dir = ROOT / "pairwise_A17_Mpo"
    p474_dir = ROOT / "pairwise_474_Mpo"
    if (a17_dir / "genome_A17.genome_Mpo.anchors.simple").exists():
        rows += summarize_pair(a17_dir, "genome_A17", "genome_Mpo", "genome_A17.genome_Mpo.anchors.simple")
    if (p474_dir / "genome_474.genome_Mpo.anchors.simple").exists():
        rows += summarize_pair(p474_dir, "genome_474", "genome_Mpo", "genome_474.genome_Mpo.anchors.simple")

    out = ROOT / "Mpo_centric_JCVI_chr_pair_summary.tsv"
    fields = [
        "comparison",
        "query_chr",
        "target_chr",
        "block_count",
        "sum_anchor_genes",
        "plus_blocks",
        "minus_blocks",
        "query_span_start",
        "query_span_end",
        "target_span_start",
        "target_span_end",
    ]
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    focus = [r for r in rows if r["target_chr"] == "Chr4"]
    focus.sort(key=lambda r: (r["comparison"], -int(r["sum_anchor_genes"])))
    out_focus = ROOT / "Mpo_Chr4_JCVI_chr_pair_summary.tsv"
    with out_focus.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(focus)

    for p in [out, out_focus]:
        print(p)


if __name__ == "__main__":
    main()
