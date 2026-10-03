#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

N4 = Path("path/to/project/N_4.pod_spiny")
MAIN = Path("path/to/project/N_3.call_SV/01.read_based_dualref_hifi")
RUN = N4 / "05_direction_free_sv_cosegregation_20260701"
OUTDIR = N4 / "05_direction_free_sv_cosegregation_20260701" / "summary"
TABIX = Path("path/to/home/anaconda3/bin/tabix")

GENE = "Chr23997"
REF = "Msa"
CALLERS = ["pbsv", "sniffles2", "cutesv"]


def parse_attrs(value: str) -> dict[str, str]:
    out = {}
    for item in value.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
    return out


def gff_records() -> list[dict[str, str]]:
    gff = N4 / "00_data/4.reference_anno/genome_Msa.gff"
    rows = []
    with gff.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attrs = parse_attrs(parts[8])
            if attrs.get("ID", "").startswith(GENE) or attrs.get("Parent", "").startswith(GENE):
                rows.append(
                    {
                        "chrom": parts[0],
                        "source": parts[1],
                        "type": parts[2],
                        "start": int(parts[3]),
                        "end": int(parts[4]),
                        "strand": parts[6],
                        "attrs": parts[8],
                    }
                )
    return rows


def sample_rows() -> list[dict[str, str]]:
    path = RUN / "config/samples_extended.tsv"
    if path.exists():
        with path.open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    return []


def parse_info(info: str) -> dict[str, str]:
    out = {}
    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
        else:
            out[item] = "1"
    return out


def svlen_value(info: dict[str, str], pos: int, end: int) -> int:
    for key in ("SVLEN", "SV_LENGTH", "END"):
        if key in info:
            try:
                if key == "END":
                    return abs(int(info[key]) - pos)
                vals = re.split(r",", info[key])
                return max(abs(int(float(v))) for v in vals if v not in {"", "."})
            except Exception:
                pass
    return abs(end - pos)


def vcf_records(path: Path, region: str) -> list[dict[str, str]]:
    if not path.exists():
        return []
    try:
        proc = subprocess.run(
            [str(TABIX), str(path), region],
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except Exception:
        return []
    rows = []
    if proc.returncode not in (0, 1):
        return []
    for raw in proc.stdout.splitlines():
        if not raw or raw.startswith("#"):
            continue
        parts = raw.split("\t")
        if len(parts) < 8:
            continue
        chrom = parts[0]
        pos = int(parts[1])
        info = parse_info(parts[7])
        end = int(info.get("END", pos + 1))
        svtype = info.get("SVTYPE", parts[4].strip("<>"))
        rows.append(
            {
                "chrom": chrom,
                "pos": pos,
                "end": end,
                "id": parts[2],
                "ref": parts[3],
                "alt": parts[4],
                "filter": parts[6],
                "svtype": svtype,
                "svlen": svlen_value(info, pos, end),
                "info": parts[7],
            }
        )
    return rows


def overlaps(a_start: int, a_end: int, b_start: int, b_end: int, pad: int = 0) -> bool:
    return max(a_start, b_start - pad) <= min(a_end, b_end + pad)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    grecs = gff_records()
    gene = next(r for r in grecs if r["type"] == "gene")
    chrom = gene["chrom"]
    gstart = gene["start"]
    gend = gene["end"]
    window_start = max(1, gstart - 10000)
    window_end = gend + 10000
    region = f"{chrom}:{window_start}-{window_end}"

    # User-prior targeted interval from manual IGV screenshots.
    manual_start = 88978818
    manual_end = 88982780
    manual_region = f"{chrom}:{manual_start}-{manual_end}"

    samples = sample_rows()
    sample_order = [r["sample"] for r in samples]
    role_by_sample = {r["sample"]: r.get("sample_role", "") for r in samples}
    pheno_by_sample = {}
    for r in samples:
        if r.get("spiny_primary") == "1":
            pheno_by_sample[r["sample"]] = "spiny"
        elif r.get("spiny_primary") == "0":
            pheno_by_sample[r["sample"]] = "spineless"
        else:
            pheno_by_sample[r["sample"]] = r.get("sample_role", "")

    all_records = []
    for sample in sample_order:
        for caller in CALLERS:
            vcf = MAIN / "03_per_sample" / REF / caller / f"{sample}.{caller}.vcf.gz"
            for rec in vcf_records(vcf, manual_region):
                rec.update(
                    {
                        "sample": sample,
                        "caller": caller,
                        "sample_role": role_by_sample.get(sample, ""),
                        "phenotype": pheno_by_sample.get(sample, ""),
                        "overlap_gene": overlaps(rec["pos"], rec["end"], gstart, gend),
                        "overlap_manual_window": overlaps(rec["pos"], rec["end"], manual_start, manual_end),
                        "near_manual_del": rec["svtype"] == "DEL"
                        and 50 <= rec["svlen"] <= 2000
                        and overlaps(rec["pos"], rec["end"], manual_start, manual_end, pad=300),
                    }
                )
                all_records.append(rec)

    records_path = OUTDIR / "Chr23997_caller_records_20260701.tsv"
    with records_path.open("w", newline="") as out:
        fieldnames = [
            "sample",
            "phenotype",
            "sample_role",
            "caller",
            "chrom",
            "pos",
            "end",
            "svtype",
            "svlen",
            "filter",
            "ref",
            "near_manual_del",
            "overlap_gene",
            "overlap_manual_window",
            "id",
            "alt",
            "info",
        ]
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_records)

    support = defaultdict(set)
    per_sample = []
    for sample in sample_order:
        row = {
            "sample": sample,
            "phenotype": pheno_by_sample.get(sample, ""),
            "sample_role": role_by_sample.get(sample, ""),
        }
        for caller in CALLERS:
            hits = [
                r
                for r in all_records
                if r["sample"] == sample and r["caller"] == caller and r["near_manual_del"]
            ]
            row[caller] = len(hits)
            if hits:
                support[sample].add(caller)
        row["n_callers_support"] = len(support[sample])
        row["support_callers"] = ",".join(sorted(support[sample])) if support[sample] else "NA"
        per_sample.append(row)

    support_path = OUTDIR / "Chr23997_per_sample_caller_support_20260701.tsv"
    with support_path.open("w", newline="") as out:
        fieldnames = [
            "sample",
            "phenotype",
            "sample_role",
            "pbsv",
            "sniffles2",
            "cutesv",
            "n_callers_support",
            "support_callers",
        ]
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(per_sample)

    print(f"gene\t{GENE}\t{chrom}:{gstart}-{gend}\tstrand={gene['strand']}")
    print(f"manual_window\t{manual_region}")
    print(f"records_file\t{records_path}")
    print(f"support_file\t{support_path}")
    print("support_by_phenotype_and_callers")
    c = Counter((r["phenotype"], r["n_callers_support"], r["support_callers"]) for r in per_sample)
    for key, value in sorted(c.items()):
        print("\t".join(map(str, key)), value, sep="\t")
    print("per_sample_support")
    for r in per_sample:
        print(
            r["sample"],
            r["phenotype"],
            r["sample_role"],
            r["n_callers_support"],
            r["support_callers"],
            sep="\t",
        )


if __name__ == "__main__":
    main()
