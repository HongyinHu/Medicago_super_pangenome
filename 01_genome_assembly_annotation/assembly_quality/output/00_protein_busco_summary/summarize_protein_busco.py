#!/usr/bin/env python3
"""Collect one selected protein-mode BUSCO result per top-level genome directory."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
from collections import Counter
from datetime import datetime
from pathlib import Path


FIELDS = [
    "species_id",
    "display_name",
    "group_code",
    "status",
    "complete_pct",
    "single_copy_pct",
    "duplicated_pct",
    "fragmented_pct",
    "missing_pct",
    "complete_single_n",
    "complete_duplicated_n",
    "fragmented_n",
    "missing_n",
    "busco_n",
    "protein_sequences",
    "lineage",
    "lineage_creation_date",
    "busco_version",
    "input_file",
    "current_input_file",
    "input_exists",
    "selected_summary_json",
    "selected_summary_txt",
    "candidate_count",
    "selection_note",
]


def display_name(species_id: str) -> tuple[str, str]:
    group, _, suffix = species_id.partition("_genome_")
    if not suffix:
        return species_id, ""
    if group == "P" and suffix.startswith("Medicago_"):
        return suffix.replace("Medicago_", "Medicago ", 1).replace("_", " "), group
    return suffix, group


def fasta_sequence_count(path: Path | None) -> int | str:
    if path is None or not path.is_file():
        return ""
    count = 0
    with path.open("rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                count += 1
    return count


def full_table_counts(summary_json: Path) -> Counter:
    table = summary_json.parent / "full_table.tsv"
    states: dict[str, str] = {}
    if not table.is_file():
        return Counter()
    with table.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                states.setdefault(fields[0], fields[1])
    return Counter(states.values())


def candidate_priority(path: Path) -> tuple[int, float]:
    text = path.as_posix().lower()
    if "busco_assess_finally" in text:
        rank = 30
    elif "/annotation_evalate/" in text:
        rank = 20
    else:
        rank = 10
    return rank, path.stat().st_mtime


def load_candidates(species_dir: Path) -> list[dict]:
    candidates = []
    for path in species_dir.rglob("short_summary.json"):
        text = path.as_posix().lower()
        if "busco_downloads" in text or "/assembly_evaluate/" in text or "/assessbly_evaluate/" in text:
            continue
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        params = payload.get("parameters", {})
        if params.get("mode") != "proteins":
            continue
        candidates.append({"path": path, "payload": payload})
    return candidates


def relative(path: str | Path, root: Path) -> str:
    path = Path(path)
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return str(path)


def locate_current_input(species_dir: Path, recorded_input: str) -> Path | None:
    basename = Path(recorded_input).name if recorded_input else ""
    candidates = []
    if recorded_input:
        candidates.append(Path(recorded_input))
    if basename:
        candidates.extend(
            [
                species_dir / "annotation_evalate" / basename,
                species_dir / "data" / basename,
                species_dir / basename,
            ]
        )
    for path in candidates:
        if path.is_file():
            return path
    return None


def row_for_candidate(species_dir: Path, root: Path, item: dict, candidate_count: int) -> dict:
    path = item["path"]
    payload = item["payload"]
    params = payload.get("parameters", {})
    lineage = payload.get("lineage_dataset", {})
    versions = payload.get("versions", {})
    results = payload.get("results", {})
    counts = full_table_counts(path)
    input_file = str(params.get("in", ""))
    current_input = locate_current_input(species_dir, input_file)
    name, group = display_name(species_dir.name)
    summary_txt = path.with_suffix(".txt")
    note = "selected by priority: finally > annotation_evalate > other protein-mode result; newest breaks ties"
    return {
        "species_id": species_dir.name,
        "display_name": name,
        "group_code": group,
        "status": "OK",
        "complete_pct": results.get("Complete", ""),
        "single_copy_pct": results.get("Single copy", ""),
        "duplicated_pct": results.get("Multi copy", ""),
        "fragmented_pct": results.get("Fragmented", ""),
        "missing_pct": results.get("Missing", ""),
        "complete_single_n": counts.get("Complete", ""),
        "complete_duplicated_n": counts.get("Duplicated", ""),
        "fragmented_n": counts.get("Fragmented", ""),
        "missing_n": counts.get("Missing", ""),
        "busco_n": results.get("n_markers", lineage.get("number_of_buscos", "")),
        "protein_sequences": fasta_sequence_count(current_input),
        "lineage": lineage.get("name", ""),
        "lineage_creation_date": lineage.get("creation_date", ""),
        "busco_version": versions.get("busco", ""),
        "input_file": input_file,
        "current_input_file": relative(current_input, root) if current_input else "",
        "input_exists": current_input is not None,
        "selected_summary_json": relative(path, root),
        "selected_summary_txt": relative(summary_txt, root) if summary_txt.is_file() else "",
        "candidate_count": candidate_count,
        "selection_note": note if candidate_count > 1 else "only protein-mode candidate",
    }


def missing_row(species_dir: Path) -> dict:
    name, group = display_name(species_dir.name)
    row = {field: "" for field in FIELDS}
    row.update(
        species_id=species_dir.name,
        display_name=name,
        group_code=group,
        status="MISSING_PROTEIN_BUSCO",
        candidate_count=0,
        selection_note="no protein-mode short_summary.json found",
    )
    return row


def write_delimited(path: Path, rows: list[dict], delimiter: str) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter=delimiter, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def fmt(value) -> str:
    return "NA" if value == "" else str(value)


def write_markdown(path: Path, rows: list[dict]) -> None:
    lines = [
        "# Protein-coding gene BUSCO summary",
        "",
        "| Species/material | Status | C (%) | S (%) | D (%) | F (%) | M (%) | Proteins | BUSCO set |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            "| {display_name} | {status} | {c} | {s} | {d} | {f} | {m} | {proteins} | {n} |".format(
                display_name=row["display_name"],
                status=row["status"],
                c=fmt(row["complete_pct"]),
                s=fmt(row["single_copy_pct"]),
                d=fmt(row["duplicated_pct"]),
                f=fmt(row["fragmented_pct"]),
                m=fmt(row["missing_pct"]),
                proteins=fmt(row["protein_sequences"]),
                n=fmt(row["busco_n"]),
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_candidates(path: Path, root: Path, records: list[tuple[Path, dict, bool]]) -> None:
    fields = ["species_id", "selected", "priority", "mtime", "complete_pct", "input_file", "summary_json"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for species_dir, item, selected in records:
            candidate = item["path"]
            results = item["payload"].get("results", {})
            params = item["payload"].get("parameters", {})
            writer.writerow(
                {
                    "species_id": species_dir.name,
                    "selected": selected,
                    "priority": candidate_priority(candidate)[0],
                    "mtime": datetime.fromtimestamp(candidate.stat().st_mtime).isoformat(timespec="seconds"),
                    "complete_pct": results.get("Complete", ""),
                    "input_file": params.get("in", ""),
                    "summary_json": relative(candidate, root),
                }
            )


def write_manifest(output: Path) -> None:
    entries = []
    for path in sorted(output.rglob("*")):
        if not path.is_file() or path.name == "manifest.sha256":
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        entries.append(f"{digest}  {path.relative_to(output).as_posix()}")
    (output / "manifest.sha256").write_text("\n".join(entries) + "\n", encoding="ascii")


def validate_rows(rows: list[dict]) -> None:
    seen = set()
    for row in rows:
        species_id = row["species_id"]
        if species_id in seen:
            raise ValueError(f"duplicate species row: {species_id}")
        seen.add(species_id)
        if row["status"] != "OK":
            continue
        exact_total = sum(
            int(row[field])
            for field in ("complete_single_n", "complete_duplicated_n", "fragmented_n", "missing_n")
        )
        if exact_total != int(row["busco_n"]):
            raise ValueError(f"BUSCO exact counts do not sum to n for {species_id}: {exact_total}")
        percent_total = float(row["complete_pct"]) + float(row["fragmented_pct"]) + float(row["missing_pct"])
        if abs(percent_total - 100.0) > 0.11:
            raise ValueError(f"BUSCO percentages do not sum to 100 for {species_id}: {percent_total}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--output-name", default="00_protein_busco_summary")
    args = parser.parse_args()
    root = args.root.resolve()
    output = root / args.output_name
    if not root.is_dir():
        raise SystemExit(f"root does not exist: {root}")
    if output.exists():
        raise SystemExit(f"output already exists; refusing to replace: {output}")

    species_dirs = sorted(
        path for path in root.iterdir() if path.is_dir() and path.name.startswith(("N_genome_", "P_genome_", "X_genome_"))
    )
    rows = []
    candidate_records = []
    selected_items = []
    for species_dir in species_dirs:
        candidates = load_candidates(species_dir)
        if not candidates:
            rows.append(missing_row(species_dir))
            continue
        selected = max(candidates, key=lambda item: candidate_priority(item["path"]))
        rows.append(row_for_candidate(species_dir, root, selected, len(candidates)))
        selected_items.append((species_dir, selected))
        for item in sorted(candidates, key=lambda x: candidate_priority(x["path"]), reverse=True):
            candidate_records.append((species_dir, item, item is selected))

    validate_rows(rows)

    output.mkdir()
    raw_dir = output / "selected_short_summaries"
    raw_dir.mkdir()
    write_delimited(output / "protein_busco_summary.tsv", rows, "\t")
    write_delimited(output / "protein_busco_summary.csv", rows, ",")
    write_markdown(output / "protein_busco_summary.md", rows)
    write_candidates(output / "all_protein_busco_candidates.tsv", root, candidate_records)
    for species_dir, item in selected_items:
        json_path = item["path"]
        shutil.copy2(json_path, raw_dir / f"{species_dir.name}.short_summary.json")
        txt_path = json_path.with_suffix(".txt")
        if txt_path.is_file():
            shutil.copy2(txt_path, raw_dir / f"{species_dir.name}.short_summary.txt")
    shutil.copy2(Path(__file__), output / "summarize_protein_busco.py")

    ok = sum(row["status"] == "OK" for row in rows)
    missing = len(rows) - ok
    accessible_inputs = sum(row["status"] == "OK" and bool(row["input_exists"]) for row in rows)
    duplicated_candidates = sum(int(row["candidate_count"] or 0) > 1 for row in rows)
    readme = f"""# 蛋白编码基因 BUSCO 结果整理

- 生成时间：{datetime.now().astimezone().isoformat(timespec='seconds')}
- 扫描目录：`{root}`
- 物种/材料目录：{len(rows)}
- 已找到蛋白模式 BUSCO：{ok}
- 缺少蛋白模式 BUSCO：{missing}
- 存在多个蛋白候选结果：{duplicated_candidates}
- 当前仍可访问的原始蛋白输入：{accessible_inputs}/{ok}

## 文件说明

- `protein_busco_summary.tsv`：主汇总表，适合 Linux/R/Python。
- `protein_busco_summary.csv`：同内容 CSV。
- `protein_busco_summary.md`：便于直接查看的简表。
- `all_protein_busco_candidates.tsv`：全部蛋白模式候选及最终选择，便于审计重复结果。
- `selected_short_summaries/`：每个物种最终采用的 BUSCO 原始短摘要副本。
- `manifest.sha256`：整理目录内文件校验值。

## 选择规则

仅纳入 BUSCO JSON 中 `mode=proteins` 的结果，并排除 `assembly_evaluate`/`assessbly_evaluate` 下的组装评估。
同一物种出现多个蛋白结果时，优先级为 `busco_assess_finally` > `annotation_evalate` > 其他蛋白模式结果；同级取修改时间较新者。所有候选均记录在审计表中。
`input_file` 保留 BUSCO JSON 中记录的历史绝对路径；若项目迁移导致该路径失效，`current_input_file` 给出当前结果树中可访问的对应蛋白文件。
"""
    (output / "README.md").write_text(readme, encoding="utf-8")
    write_manifest(output)
    print(f"OUTPUT={output}")
    print(f"SPECIES={len(rows)} OK={ok} MISSING={missing} MULTI_CANDIDATE={duplicated_candidates}")
    print(f"VALIDATION_OK unique_species={len(rows)} accessible_inputs={accessible_inputs}/{ok} exact_count_checks={ok} percentage_checks={ok}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
