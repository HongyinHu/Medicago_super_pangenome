#!/usr/bin/env python3
from __future__ import annotations

from collections import Counter
from pathlib import Path
import csv

OUT = Path("path/to/project/N_4.pod_spiny/05_direction_free_sv_cosegregation_20260701")


def rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def count_by(rows_: list[dict[str, str]], *cols: str) -> Counter:
    c = Counter()
    for row in rows_:
        key = tuple(row.get(col, "") for col in cols)
        c[key] += 1
    return c


def print_counter(title: str, counter: Counter) -> None:
    print(title)
    for key, value in sorted(counter.items()):
        label = "\t".join(key) if isinstance(key, tuple) else str(key)
        print(f"{label}\t{value}")


def scan_logs() -> None:
    patterns = (
        "Traceback",
        "No space left",
        "Disk quota exceeded",
        "ENOSPC",
        "ERROR",
        "error:",
        "CalledProcessError",
    )
    hits = []
    for path in sorted((OUT / "logs").glob("*")):
        if not path.is_file():
            continue
        text = path.read_text(errors="ignore")
        for pattern in patterns:
            if pattern in text:
                hits.append((path.name, pattern))
                break
    print("log_issue_files")
    if hits:
        for name, pattern in hits:
            print(f"{name}\t{pattern}")
    else:
        print("NONE")


def main() -> None:
    summary = rows(OUT / "results" / "candidate_strict_flankQC_summary.tsv")
    selected = rows(OUT / "results" / "direction_free_igv_selected.tsv")
    manifest = rows(OUT / "igv_reports_html" / "igv_reports_manifest.tsv")
    html = sorted((OUT / "igv_reports_html").glob("*.html"))
    bams = sorted((OUT / "bams").glob("*/*.bam"))

    print(f"run_dir\t{OUT}")
    print(f"direction_free_igv_selected_rows\t{len(selected)}")
    print(f"flankQC_summary_rows\t{len(summary)}")
    print(f"igv_manifest_rows\t{len(manifest)}")
    print(f"html_files\t{len(html)}")
    print(f"bam_window_files\t{len(bams)}")
    print_counter("selected_by_dataset_ref", count_by(selected, "dataset", "ref"))
    print_counter("selected_by_best_direction", count_by(selected, "best_direction"))
    print_counter("flankQC_by_ref_verdict", count_by(summary, "ref", "strict_flankQC_verdict"))
    print_counter("igv_by_ref_verdict", count_by(manifest, "ref", "verdict"))
    scan_logs()
    print("important_paths")
    for rel in (
        "results/direction_free_all_candidates.tsv",
        "results/direction_free_exact_core.tsv",
        "results/direction_free_near_or_exact.tsv",
        "results/direction_free_igv_selected.tsv",
        "results/candidate_strict_flankQC_summary.tsv",
        "results/candidate_strict_pass.tsv",
        "results/candidate_manual_review.tsv",
        "results/candidate_fail_conflict.tsv",
        "results/per_sample_event_flankQC.tsv",
        "igv_reports_html/index.html",
        "igv_reports_html/igv_reports_manifest.tsv",
    ):
        path = OUT / rel
        print(f"{rel}\t{path.exists()}\t{path}")


if __name__ == "__main__":
    main()
