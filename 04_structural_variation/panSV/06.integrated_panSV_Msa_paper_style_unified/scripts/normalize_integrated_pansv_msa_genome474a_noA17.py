#!/usr/bin/env python3
"""Normalize Msa single-reference integrated panSV tables.

Rules:
  - Use genome_Msa as the only reference coordinate system.
  - Use genome_G474a as the representative assembly for final genome_474.
  - Exclude genome_G474b and genome_A17 from final SVGAP support.
  - Keep only read-mapping SV events as final panSV events.
  - Treat overlapping SVGAP events as evidence only; never overwrite PAV.
  - Drop all SVGAP-only events because they have no read support.
  - Normalize SVGAP chromosome names from Msa.ChrN to ChrN.
"""

import argparse
import csv
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Set


FINAL_SAMPLES = [
    "genome_395",
    "genome_410",
    "genome_436",
    "genome_454",
    "genome_457",
    "genome_461",
    "genome_468",
    "genome_472",
    "genome_474",
    "genome_482",
    "genome_M22",
    "genome_M46",
    "genome_Mar",
    "genome_Mpo",
    "genome_Mru",
    "genome_Msa",
    "genome_R108",
    "genome_ZM4",
]

DROP_LABELS = {"genome_A17", "A17", "genome_G474b", "G474b"}
MAP_LABELS = {"genome_G474a": "genome_474", "G474a": "genome_474"}
DROP_SUBSTRINGS = ("genome_A17", "A17", "genome_G474b", "G474b")


def split_labels(value):
    if not value or value == ".":
        return []
    labels = []
    for token in value.replace(";", ",").split(","):
        token = token.strip()
        if token and token != ".":
            labels.append(token)
    return labels


def normalize_labels(value):
    out = []
    seen = set()
    for label in split_labels(value):
        if label in DROP_LABELS:
            continue
        label = MAP_LABELS.get(label, label)
        if label not in FINAL_SAMPLES:
            continue
        if label not in seen:
            out.append(label)
            seen.add(label)
    return out


def clean_provenance_ids(value):
    if not value or value == ".":
        return "."
    out = []
    seen = set()
    for token in split_labels(value):
        if any(x in token for x in DROP_SUBSTRINGS):
            continue
        token = token.replace("genome_G474a", "genome_474")
        token = token.replace("G474a", "genome_474")
        if token and token not in seen:
            out.append(token)
            seen.add(token)
    return ",".join(out) if out else "."


def clean_note(notes):
    notes = notes.replace("svgap_support_samples_normalized_G474a_noA17", "svgap_support_samples_normalized")
    notes = notes.replace("svgap_support_excluded_after_A17_G474b_removal", "svgap_support_excluded_after_sample_policy")
    return notes


def norm_chrom(chrom):
    return chrom[4:] if chrom.startswith("Msa.") else chrom


def append_note(notes, note):
    if not notes or notes == ".":
        return note
    if note in notes.split(";"):
        return notes
    return notes + ";" + note


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input-dir", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--run-dir", required=True, type=Path)
    return p.parse_args()


def main():
    args = parse_args()
    in_results = args.input_dir / "results"
    out_results = args.output_dir / "results"
    out_summary = args.output_dir / "summary"
    out_results.mkdir(parents=True, exist_ok=True)
    out_summary.mkdir(parents=True, exist_ok=True)

    events_in = in_results / "Msa_single_ref.final_panSV.events.tsv"
    pav_in = in_results / "Msa_single_ref.final_panSV.PAV.matrix.tsv"
    stats_in = in_results / "Msa_single_ref.final_panSV.stats.tsv"
    evidence_in = in_results / "Msa_single_ref.svgap_readbased.overlap_evidence.tsv"

    prefix = "paper_style_unified_panSV.Msa_ref"
    events_out = out_results / f"{prefix}.events.tsv"
    pav_out = out_results / f"{prefix}.PAV.matrix.tsv"
    stats_out = out_results / f"{prefix}.stats.tsv"
    counts_out = out_summary / f"{prefix}.counts.tsv"
    method_out = out_summary / f"{prefix}.method.tsv"
    done_out = out_summary / f"{prefix}.done"
    evidence_out = out_results / f"{prefix}.svgap_readbased.overlap_evidence.tsv"
    dropped_out = out_summary / f"{prefix}.dropped_excluded_samples.tsv"

    row_info = {}
    counts = Counter()
    dropped_counts = Counter()

    with events_in.open(newline="") as fin, events_out.open("w", newline="") as fout, dropped_out.open("w", newline="") as fd:
        reader = csv.DictReader(fin, delimiter="\t")
        fieldnames = reader.fieldnames or []
        writer = csv.DictWriter(fout, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        drop_writer = csv.writer(fd, delimiter="\t", lineterminator="\n")
        drop_writer.writerow(["FinalSV_ID", "old_source_methods", "reason", "old_svgap_support_samples", "old_svgap_ids"])

        for row in reader:
            old_source = row["source_methods"]
            norm_support = normalize_labels(row.get("svgap_support_samples", "."))
            keep = True
            new_source = old_source
            new_conf = row["final_confidence"]
            note = row.get("notes", ".")

            if old_source == "SVGAP":
                keep = False
                dropped_counts["SVGAP_only_removed_no_read_support"] += 1
                drop_writer.writerow([
                    row["FinalSV_ID"],
                    old_source,
                    "svgap_only_removed_no_read_mapping_support",
                    row.get("svgap_support_samples", "."),
                    row.get("svgap_ids", "."),
                ])
            elif old_source == "read_based+SVGAP":
                if norm_support:
                    new_source = "read_based+SVGAP"
                    new_conf = "read_svgap_supported"
                    note = append_note(note, "svgap_evidence_only_pav_from_read_mapping")
                else:
                    new_source = "read_based"
                    new_conf = "read_based"
                    note = append_note(note, "svgap_support_excluded_after_sample_policy")
                    row["svgap_ids"] = "."
                    row["svgap_source_id"] = "."
                    row["svgap_source_file"] = "."
                    dropped_counts["read_svgap_demoted_to_read_based"] += 1

            if not keep:
                row_info[row["FinalSV_ID"]] = {"keep": False}
                continue

            row["source_methods"] = new_source
            row["reference"] = "genome_Msa"
            row["CHROM"] = norm_chrom(row["CHROM"])
            row["read_SV_ID"] = clean_provenance_ids(row.get("read_SV_ID", "."))
            row["svgap_support_count"] = str(len(norm_support)) if new_source != "read_based" else "0"
            row["svgap_support_samples"] = ",".join(norm_support) if norm_support and new_source != "read_based" else "."
            if new_source == "read_based":
                row["svgap_ids"] = "."
                row["svgap_source_id"] = "."
                row["svgap_source_file"] = "."
            else:
                row["svgap_ids"] = clean_provenance_ids(row.get("svgap_ids", "."))
                row["svgap_source_id"] = clean_provenance_ids(row.get("svgap_source_id", "."))
                row["svgap_source_file"] = clean_provenance_ids(row.get("svgap_source_file", "."))
            row["final_confidence"] = new_conf
            row["notes"] = clean_note(note)

            row_info[row["FinalSV_ID"]] = {
                "keep": True,
                "source_methods": new_source,
                "final_confidence": new_conf,
                "norm_support": set(norm_support),
                "chrom": row["CHROM"],
            }
            counts["events_retained"] += 1
            counts[f"source_methods:{new_source}"] += 1
            counts[f"final_confidence:{new_conf}"] += 1
            counts[f"SVTYPE:{row['SVTYPE']}"] += 1
            writer.writerow(row)

    with pav_in.open(newline="") as fin, pav_out.open("w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        fieldnames = reader.fieldnames or []
        missing = [s for s in FINAL_SAMPLES if s not in fieldnames]
        if missing:
            raise SystemExit(f"Missing expected sample columns in PAV: {missing}")
        writer = csv.DictWriter(fout, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            info = row_info.get(row["FinalSV_ID"])
            if not info or not info.get("keep"):
                continue
            row["source_methods"] = str(info["source_methods"])
            row["reference"] = "genome_Msa"
            row["CHROM"] = norm_chrom(row["CHROM"])
            if "read_SV_ID" in row:
                row["read_SV_ID"] = clean_provenance_ids(row["read_SV_ID"])
            row["final_confidence"] = str(info["final_confidence"])
            support = info["norm_support"]
            if row["source_methods"] == "SVGAP":
                for sample in FINAL_SAMPLES:
                    row[sample] = "1" if sample in support else "0"
            elif row["source_methods"] == "read_based+SVGAP":
                pass
            writer.writerow(row)

    with evidence_in.open(newline="") as fin, evidence_out.open("w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        fieldnames = reader.fieldnames or []
        writer = csv.DictWriter(fout, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        kept_read_ids = {
            row_id for row_id, info in row_info.items()
            if info.get("keep") and info.get("source_methods") == "read_based+SVGAP"
        }
        # evidence uses read_SV_ID, not FinalSV_ID. Keep rows unless they are
        # explicitly from excluded A17/G474b source IDs; this is supporting
        # evidence only, while final PAV is controlled by normalized event rows.
        excluded_evidence = 0
        kept_evidence = 0
        for row in reader:
            source_id = row.get("svgap_source_id", "")
            if "A17" in source_id or "G474b" in source_id:
                excluded_evidence += 1
                continue
            if row.get("chrom", "").startswith("Msa."):
                row["chrom"] = norm_chrom(row["chrom"])
            if "read_SV_ID" in row:
                row["read_SV_ID"] = clean_provenance_ids(row["read_SV_ID"])
            if "svgap_source_id" in row:
                row["svgap_source_id"] = clean_provenance_ids(row["svgap_source_id"])
            for key in row:
                row[key] = row[key].replace("genome_G474a", "genome_474").replace("G474a", "genome_474")
            writer.writerow(row)
            kept_evidence += 1
        counts["evidence_rows_retained"] = kept_evidence
        dropped_counts["evidence_rows_A17_or_G474b"] = excluded_evidence

    # Recompute sample count summaries from final PAV.
    sample_cell_counts = {s: Counter() for s in FINAL_SAMPLES}
    with pav_out.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            for sample in FINAL_SAMPLES:
                sample_cell_counts[sample][row[sample]] += 1

    with stats_out.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["section", "key", "count"])
        w.writerow(["input", "source_dir", str(args.input_dir)])
        w.writerow(["parameters", "reference", "genome_Msa"])
        w.writerow(["parameters", "event_policy", "read_mapping_primary_svgap_evidence_only_no_svgap_only"])
        w.writerow(["parameters", "genome_474_policy", "use_genome_G474a_only_as_genome_474"])
        w.writerow(["parameters", "excluded_samples", "genome_A17,genome_G474b"])
        w.writerow(["parameters", "chromosome_name_policy", "strip_Msa_prefix"])
        for key, value in sorted(counts.items()):
            if ":" in key:
                section, k = key.split(":", 1)
                w.writerow([section, k, value])
            else:
                w.writerow(["integrated", key, value])
        for key, value in sorted(dropped_counts.items()):
            w.writerow(["dropped_or_demoted", key, value])
        for sample in FINAL_SAMPLES:
            c = sample_cell_counts[sample]
            w.writerow(["sample_called_0", sample, c["0"]])
            w.writerow(["sample_called_1", sample, c["1"]])
            w.writerow(["sample_called_NA", sample, c["NA"] + c["."]])

    with counts_out.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["section", "key", "count"])
        for key, value in sorted(counts.items()):
            if ":" in key:
                section, k = key.split(":", 1)
                w.writerow([section, k, value])
            else:
                w.writerow(["rows", key, value])
        for key, value in sorted(dropped_counts.items()):
            w.writerow(["dropped_or_demoted", key, value])

    with method_out.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["metric", "value"])
        w.writerow(["method", "Msa_single_ref_read_mapping_primary_no_svgap_only_genome474a_noA17"])
        w.writerow(["reference", "genome_Msa"])
        w.writerow(["input_single_ref", str(args.input_dir)])
        w.writerow(["genome_474_policy", "map genome_G474a to genome_474; ignore genome_G474b"])
        w.writerow(["excluded_samples", "genome_A17,genome_G474b"])
        w.writerow(["event_policy", "read-mapping SV is the representative event; one or more overlapping SVGAP events are retained only as assembly-support evidence"])
        w.writerow(["PAV_rule", "read_mapping PAV retained as the final presence/absence state; overlapping SVGAP is evidence only and is not ORed into PAV; SVGAP-only events are removed"])
        w.writerow(["chromosome_name_policy", "Msa.ChrN normalized to ChrN"])
        w.writerow(["created_at", datetime.now().strftime("%Y-%m-%d %H:%M:%S")])
        w.writerow(["run_dir", str(args.run_dir)])

    done_out.write_text("done\n" + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
