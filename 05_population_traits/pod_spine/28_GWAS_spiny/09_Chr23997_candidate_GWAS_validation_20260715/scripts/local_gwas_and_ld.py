#!/usr/bin/env python3
"""Collect Chr23997 +/- window association results and un-imputed target-marker LD."""

import argparse
import csv
import gzip
import math
import subprocess
from pathlib import Path

import numpy as np


TARGET_CHROM = "Chr4"
TARGET_POS = 88980831
TARGET_END = 88981052


def normalize_chrom(value):
    text = str(value).strip().lower()
    if text.startswith("chr"):
        text = text[3:]
    return text.lstrip("0") or "0"


def r_squared(left, right):
    pairs = [(float(x), float(y)) for x, y in zip(left, right) if x is not None and y is not None]
    if len(pairs) < 3:
        return math.nan
    x = np.asarray([item[0] for item in pairs], dtype=float)
    y = np.asarray([item[1] for item in pairs], dtype=float)
    if np.std(x) == 0.0 or np.std(y) == 0.0:
        return math.nan
    value = float(np.corrcoef(x, y)[0, 1] ** 2)
    if abs(value - 1.0) < 1e-12:
        return 1.0
    return min(1.0, max(0.0, value))


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def write_tsv(path, fieldnames, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_association(path, label, start, end):
    rows = []
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if normalize_chrom(row["chr"]) != "4":
                continue
            try:
                position = int(float(row["pos"]))
                pvalue = float(row["p_wald"])
            except (TypeError, ValueError):
                continue
            if start <= position <= end and 0.0 < pvalue <= 1.0:
                rows.append({
                    "class": row.get("class", label), "chr": "Chr4", "pos": position,
                    "end": row.get("end", position), "id": row.get("id", ""),
                    "pvalue": pvalue, "minus_log10_p": -math.log10(pvalue),
                    "source": label,
                })
    return rows


def candidate_row(path):
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["callset"] == "strict" and row["model"] == "Fisher_exact_two_sided":
                pvalue = float(row["pvalue"])
                return {
                    "class": "Candidate_SV", "chr": TARGET_CHROM, "pos": TARGET_POS,
                    "end": TARGET_END, "id": "Chr23997_intron2_DEL221",
                    "pvalue": pvalue, "minus_log10_p": -math.log10(pvalue),
                    "source": "targeted_direct_evidence",
                }
    raise ValueError("strict Fisher candidate result is absent")


def vcf_samples(path):
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    raise ValueError("VCF header is missing sample columns: {}".format(path))


def vcf_chromosome(tabix, path):
    completed = subprocess.run([tabix, "-l", path], check=True, text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for contig in completed.stdout.splitlines():
        if normalize_chrom(contig) == "4":
            return contig
    raise ValueError("{} has no chromosome 4 contig".format(path))


def dosage(gt):
    allele = gt.split(":", 1)[0]
    if allele in {".", "./.", ".|."}:
        return None
    tokens = allele.replace("|", "/").split("/")
    if len(tokens) != 2 or any(token not in {"0", "1"} for token in tokens):
        return None
    return int(tokens[0]) + int(tokens[1])


def local_ld_from_vcf(tabix, path, label, target_by_sample, start, end):
    samples = vcf_samples(path)
    shared = [(index, sample) for index, sample in enumerate(samples) if sample in target_by_sample]
    if not shared:
        raise ValueError("{} has no samples shared with target callset".format(path))
    region = "{}:{}-{}".format(vcf_chromosome(tabix, path), start, end)
    completed = subprocess.run([tabix, path, region], check=True,
                               text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    records = []
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 10:
            continue
        marker_values = [dosage(value) for value in fields[9:]]
        target = [target_by_sample[sample] for _, sample in shared]
        marker = [marker_values[index] if index < len(marker_values) else None for index, _ in shared]
        paired = sum(x is not None and y is not None for x, y in zip(target, marker))
        rsq = r_squared(target, marker)
        records.append({
            "class": label, "chr": fields[0], "pos": int(fields[1]), "id": fields[2], "alt": fields[4],
            "n_pairwise_called": paired, "r2_with_Chr23997_221bp_PAV": "NA" if math.isnan(rsq) else "{:.10g}".format(rsq),
        })
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strict-target", required=True)
    parser.add_argument("--candidate-association", required=True)
    parser.add_argument("--snp-indel-association", required=True)
    parser.add_argument("--sv-association", required=True)
    parser.add_argument("--snp-vcf", required=True)
    parser.add_argument("--indel-vcf", required=True)
    parser.add_argument("--sv-vcf", required=True)
    parser.add_argument("--tabix", required=True)
    parser.add_argument("--window", type=int, default=100000)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    start = TARGET_POS - args.window
    end = TARGET_END + args.window
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    association = []
    association.extend(load_association(args.snp_indel_association, "SNP_INDEL_GEMMA", start, end))
    association.extend(load_association(args.sv_association, "SV_GEMMA", start, end))
    association.append(candidate_row(args.candidate_association))
    association.sort(key=lambda row: (row["pos"], row["class"], row["id"]))
    fields = ["class", "chr", "pos", "end", "id", "pvalue", "minus_log10_p", "source"]
    write_tsv(outdir / "Chr23997.local_association.tsv", fields, association)

    strict_rows = read_strict = []
    with open(args.strict_target, "r", newline="") as handle:
        strict_rows = list(csv.DictReader(handle, delimiter="\t"))
    target_by_sample = {
        row["sample"]: (None if row["genotype"] == "NA" else int(row["genotype"])) for row in strict_rows
    }
    ld_rows = []
    ld_rows.extend(local_ld_from_vcf(args.tabix, args.snp_vcf, "SNP", target_by_sample, start, end))
    ld_rows.extend(local_ld_from_vcf(args.tabix, args.indel_vcf, "INDEL", target_by_sample, start, end))
    ld_rows.extend(local_ld_from_vcf(args.tabix, args.sv_vcf, "SV", target_by_sample, start, end))
    ld_rows.sort(key=lambda row: (normalize_chrom(row["chr"]), row["pos"], row["id"]))
    write_tsv(outdir / "Chr23997.local_LD.tsv",
              ["class", "chr", "pos", "id", "alt", "n_pairwise_called", "r2_with_Chr23997_221bp_PAV"], ld_rows)

    tag_candidates = [row for row in association if row["class"] in {"SNP", "INDEL"}]
    if tag_candidates:
        best = min(tag_candidates, key=lambda row: row["pvalue"])
        best["selection"] = "lowest_nominal_p_in_Chr23997_100kb_window"
        write_tsv(outdir / "Chr23997.local_best_SNP_INDEL_tag.tsv", fields + ["selection"], [best])
    else:
        write_tsv(outdir / "Chr23997.local_best_SNP_INDEL_tag.tsv", fields + ["selection"], [])


if __name__ == "__main__":
    main()
