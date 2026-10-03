#!/usr/bin/env python3
"""Small helpers for audited SHP/AG-like functional alleles."""

import csv
import gzip
import math
import re
import argparse
from collections import defaultdict
from pathlib import Path


FUNCTIONAL_LABELS = {
    ("Chr3", 53196520, "G", "A"): "MsaSHP_P201S",
    ("Chr3", 53196478, "C", "T"): "MsaSHP_A215T",
    ("Chr4", 91177192, "A", "G"): "AGlike_S138G",
}

FUNCTIONAL_METADATA = {
    "MsaSHP_P201S": {"gene": "Chr14120", "gene_class": "MsaSHP", "aa_change": "P201S"},
    "MsaSHP_A215T": {"gene": "Chr14120", "gene_class": "MsaSHP", "aa_change": "A215T"},
    "AGlike_S138G": {"gene": "Chr24229", "gene_class": "AG-like", "aa_change": "S138G"},
}


def alt_dosage(gt, allele_index):
    if not gt or gt in {".", "./.", ".|."}:
        return None
    alleles = gt.replace("|", "/").split("/")
    if any(allele == "." for allele in alleles):
        return None
    return sum(int(allele) == allele_index for allele in alleles)


def functional_label(chrom, pos, ref, alt):
    return FUNCTIONAL_LABELS.get((chrom, int(pos), ref.upper(), alt.upper()))


def fisher_exact_two_sided(a, b, c, d):
    row1, row2, col1 = a + b, c + d, a + c
    total = row1 + row2
    if total == 0:
        return math.nan

    def choose(n, k):
        if k < 0 or k > n:
            return 0
        return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

    def probability(x):
        return choose(row1, x) * choose(row2, col1 - x) / choose(total, col1)

    observed = probability(a)
    lower = max(0, col1 - row2)
    upper = min(row1, col1)
    return sum(probability(x) for x in range(lower, upper + 1) if probability(x) <= observed + 1e-12)


def _open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, encoding="utf-8")


def collect_functional_rows(vcf_path, phenotypes, analysis_ids):
    rows = []
    sample_ids = None
    with _open_text(vcf_path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                sample_ids = line.rstrip("\n").split("\t")[9:]
                continue
            if not line.startswith("Chr"):
                continue
            if sample_ids is None:
                raise ValueError("VCF sample header is missing")
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
            site = functional_label(chrom, pos, ref, alt)
            if site is None:
                continue
            format_keys = fields[8].split(":")
            metadata = FUNCTIONAL_METADATA[site]
            for sample_id, value in zip(sample_ids, fields[9:]):
                values = dict(zip(format_keys, value.split(":")))
                genotype = values.get("GT", ".")
                phenotype = phenotypes.get(sample_id, {})
                rows.append(
                    {
                        "sample_id": sample_id,
                        "site": site,
                        "gene": metadata["gene"],
                        "gene_class": metadata["gene_class"],
                        "aa_change": metadata["aa_change"],
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "gt": genotype,
                        "dosage": alt_dosage(genotype, 1),
                        "dp": values.get("DP", ""),
                        "ad": values.get("AD", ""),
                        "phenotype_available": int(sample_id in phenotypes),
                        "coil_ge1": phenotype.get("coil_ge1", ""),
                        "average_coil": phenotype.get("average_coil", ""),
                        "latin_name": phenotype.get("latin_name", ""),
                        "analysis_sample": int(sample_id in analysis_ids),
                    }
                )
    return rows


def load_phenotypes(path):
    phenotypes = {}
    with open(path, encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sample_id = row.get("Sample_ID") or row.get("IID")
            if not sample_id:
                raise ValueError("Phenotype table lacks Sample_ID or IID")
            match = re.search(r"average_coil=~?([0-9]+(?:\.[0-9]+)?)", row["Pod_spiral_raw"])
            if not match:
                raise ValueError(f"Missing average_coil for {sample_id}")
            average_coil = float(match.group(1))
            phenotypes[sample_id] = {
                "coil_ge1": int(average_coil >= 1.0),
                "average_coil": average_coil,
                "latin_name": row.get("Latin_name", ""),
            }
    return phenotypes


def summarize_functional_rows(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["site"]].append(row)
    summaries = []
    for site, site_rows in sorted(grouped.items()):
        called = [row for row in site_rows if row["dosage"] is not None]
        phenotyped = [row for row in called if row["phenotype_available"]]
        analysis = [row for row in phenotyped if row["analysis_sample"]]
        alt_case = sum(row["dosage"] > 0 and row["coil_ge1"] == 1 for row in phenotyped)
        alt_control = sum(row["dosage"] > 0 and row["coil_ge1"] == 0 for row in phenotyped)
        ref_case = sum(row["dosage"] == 0 and row["coil_ge1"] == 1 for row in phenotyped)
        ref_control = sum(row["dosage"] == 0 and row["coil_ge1"] == 0 for row in phenotyped)
        alt_allele_count = sum(row["dosage"] for row in called)
        total_alleles = 2 * len(called)
        alt_frequency = alt_allele_count / total_alleles if total_alleles else math.nan
        summaries.append(
            {
                "site": site,
                "gene": site_rows[0].get("gene", ""),
                "gene_class": site_rows[0].get("gene_class", ""),
                "aa_change": site_rows[0].get("aa_change", ""),
                "chrom": site_rows[0].get("chrom", ""),
                "pos": site_rows[0].get("pos", ""),
                "ref": site_rows[0].get("ref", ""),
                "alt": site_rows[0].get("alt", ""),
                "vcf_samples": len(site_rows),
                "called_samples": len(called),
                "missing_samples": len(site_rows) - len(called),
                "phenotyped_called_samples": len(phenotyped),
                "analysis_called_samples": len(analysis),
                "alt_allele_count": alt_allele_count,
                "alt_frequency": alt_frequency,
                "maf": min(alt_frequency, 1 - alt_frequency) if math.isfinite(alt_frequency) else math.nan,
                "alt_case_carriers": alt_case,
                "alt_control_carriers": alt_control,
                "ref_case_samples": ref_case,
                "ref_control_samples": ref_control,
                "raw_fisher_p": fisher_exact_two_sided(alt_case, alt_control, ref_case, ref_control),
            }
        )
    return summaries


def _read_keep_ids(path):
    with open(path, encoding="utf-8") as handle:
        return {line.split()[1] for line in handle if line.split()}


def _write_tsv(path, rows):
    if not rows:
        raise ValueError(f"No rows to write: {path}")
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("audit", nargs="?")
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--phenotype", required=True, type=Path)
    parser.add_argument("--analysis-keep", required=True, type=Path)
    parser.add_argument("--out-genotypes", required=True, type=Path)
    parser.add_argument("--out-summary", required=True, type=Path)
    parser.add_argument("--expected-vcf-samples", required=True, type=int)
    args = parser.parse_args()

    phenotypes = load_phenotypes(args.phenotype)
    analysis_ids = _read_keep_ids(args.analysis_keep)
    rows = collect_functional_rows(args.vcf, phenotypes, analysis_ids)
    if not rows:
        raise ValueError("No predefined functional allele was found in the split VCF")
    sites = {row["site"] for row in rows}
    expected_sites = set(FUNCTIONAL_METADATA)
    if sites != expected_sites:
        raise ValueError(f"Functional sites in VCF differ from expectation: {sites} vs {expected_sites}")
    observed_sizes = {site: sum(row["site"] == site for row in rows) for site in sites}
    if set(observed_sizes.values()) != {args.expected_vcf_samples}:
        raise ValueError(f"Unexpected VCF sample counts by site: {observed_sizes}")

    args.out_genotypes.parent.mkdir(parents=True, exist_ok=True)
    args.out_summary.parent.mkdir(parents=True, exist_ok=True)
    _write_tsv(args.out_genotypes, rows)
    _write_tsv(args.out_summary, summarize_functional_rows(rows))
    print(f"Wrote {len(rows)} functional genotype rows and {len(sites)} site summaries")


if __name__ == "__main__":
    main()
