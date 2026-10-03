#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_interval(text: str) -> tuple[int, int] | None:
    if not text or "-" not in text:
        return None
    a, b = text.split("-", 1)
    return int(a), int(b)


def overlap(a: tuple[int, int], b: tuple[int, int]) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def span(iv: tuple[int, int]) -> int:
    return max(0, iv[1] - iv[0])


def union(intervals: list[tuple[int, int]]) -> tuple[int, int]:
    return min(x[0] for x in intervals), max(x[1] for x in intervals)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", required=True)
    parser.add_argument("--results", required=True)
    args = parser.parse_args()

    species = args.species
    results = Path(args.results)
    rows = read_tsv(results / f"{species}_CENH3_recall_consensus.tsv")
    curated = []

    for r in rows:
        raw = parse_interval(r["raw_domain"])
        relaxed = parse_interval(r["relaxed_domain"])
        q20 = parse_interval(r["q20_domain"])
        used = []
        note = []

        if raw and relaxed and overlap(raw, relaxed) / max(1, min(span(raw), span(relaxed))) >= 0.5:
            core = union([raw, relaxed])
            used = [raw, relaxed]
            note.append("raw_relaxed_concordant")
        elif raw and relaxed:
            core = union([raw, relaxed])
            used = [raw, relaxed]
            note.append("raw_relaxed_broad_or_split")
        elif raw:
            core = raw
            used = [raw]
            note.append("raw_only")
        elif relaxed:
            core = relaxed
            used = [relaxed]
            note.append("relaxed_only")
        elif q20:
            core = q20
            used = [q20]
            note.append("q20_only")
        else:
            continue

        q20_status = "missing"
        q20_overlap_bp = 0
        if q20:
            q20_overlap_bp = overlap(core, q20)
            q20_overlap_frac_vs_q20 = q20_overlap_bp / max(1, span(q20))
            q20_overlap_frac_vs_core = q20_overlap_bp / max(1, span(core))
            if q20_overlap_frac_vs_q20 >= 0.25 or q20_overlap_frac_vs_core >= 0.25:
                used.append(q20)
                q20_status = "concordant_or_partly_overlapping"
                note.append("q20_included")
            else:
                q20_status = "discordant_excluded"
                note.append("q20_discordant_excluded")

        start, end = union(used)
        curated.append(
            {
                "chrom": r["chrom"],
                "curated_start": start,
                "curated_end": end,
                "curated_span_bp": end - start,
                "used_modes": ",".join(
                    mode
                    for mode, iv in [("raw", raw), ("relaxed", relaxed), ("q20", q20)]
                    if iv and any(iv == x for x in used)
                ),
                "q20_status": q20_status,
                "q20_overlap_with_raw_relaxed_bp": q20_overlap_bp,
                "raw_domain": r["raw_domain"],
                "relaxed_domain": r["relaxed_domain"],
                "q20_domain": r["q20_domain"],
                "old_core_overlap_bp": r.get("overlap_old_core_bp", ""),
                "old_core_names": r.get("old_core_names", ""),
                "curation_note": ";".join(note),
            }
        )

    fields = [
        "chrom",
        "curated_start",
        "curated_end",
        "curated_span_bp",
        "used_modes",
        "q20_status",
        "q20_overlap_with_raw_relaxed_bp",
        "raw_domain",
        "relaxed_domain",
        "q20_domain",
        "old_core_overlap_bp",
        "old_core_names",
        "curation_note",
    ]
    write_tsv(results / f"{species}_CENH3_recall_curated_consensus.tsv", curated, fields)

    bed = results / f"{species}.CENH3.functional_centromere.recall.high_confidence.bed"
    with bed.open("w") as fh:
        for r in curated:
            fh.write(
                "\t".join(
                    [
                        str(r["chrom"]),
                        str(r["curated_start"]),
                        str(r["curated_end"]),
                        f"{r['chrom']}_CENH3_high_confidence",
                        str(r["used_modes"]),
                        str(r["q20_status"]),
                    ]
                )
                + "\n"
            )


if __name__ == "__main__":
    main()
