#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path


ROOT = Path("path/to/project/N_5.pod_coil/02_SNP_INDEL_SV_GWAS_20260713")
STAGE = Path("path/to/project/N_4.pod_spiny/28_GWAS_spiny/05_sv_genotyping")
MANIFEST = STAGE / "summary" / "sv_analysis_manifest.tsv"
BCFTOOLS = "path/to/home/anaconda3/envs/panpop/bin/bcftools"
REPORT = ROOT / "03_sv_gwas" / "summary" / "existing_sv_genotype_validation.tsv"
SUMMARY = ROOT / "03_sv_gwas" / "summary" / "existing_sv_genotype_validation.done"


def command_ok(args: list[str]) -> tuple[bool, str]:
    proc = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return proc.returncode == 0, proc.stderr.strip()


def validate(row: dict[str, str]) -> dict[str, str | int]:
    sid = row["sample_id"]
    vcf = STAGE / "joint_genotype_per_sample" / sid / f"{sid}-smoove.genotyped.vcf.gz"
    index = Path(f"{vcf}.csi")
    result: dict[str, str | int] = {
        "sample_id": sid,
        "vcf": str(vcf),
        "vcf_bytes": vcf.stat().st_size if vcf.exists() else 0,
        "index_present": int(index.is_file() and index.stat().st_size > 0),
        "gzip_ok": 0,
        "header_ok": 0,
        "sample_ok": 0,
        "record_count": 0,
        "status": "FAIL",
        "message": "",
    }
    if not vcf.is_file() or vcf.stat().st_size == 0:
        result["message"] = "missing_or_empty_vcf"
        return result
    if not result["index_present"]:
        result["message"] = "missing_or_empty_csi"
        return result

    try:
        with gzip.open(vcf, "rb") as handle:
            while handle.read(8 * 1024 * 1024):
                pass
        result["gzip_ok"] = 1
    except Exception as exc:  # pragma: no cover - reports external file corruption
        result["message"] = f"gzip:{exc}"
        return result

    ok, err = command_ok([BCFTOOLS, "view", "-h", str(vcf)])
    result["header_ok"] = int(ok)
    if not ok:
        result["message"] = f"header:{err}"
        return result

    proc = subprocess.run(
        [BCFTOOLS, "query", "-l", str(vcf)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    samples = [x for x in proc.stdout.splitlines() if x]
    result["sample_ok"] = int(proc.returncode == 0 and samples == [sid])
    if not result["sample_ok"]:
        result["message"] = f"sample_columns:{','.join(samples)};stderr:{proc.stderr.strip()}"
        return result

    proc = subprocess.run(
        [BCFTOOLS, "index", "-n", str(vcf)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    try:
        count = int(proc.stdout.strip())
    except ValueError:
        count = 0
    result["record_count"] = count
    if proc.returncode != 0 or count <= 0:
        result["message"] = f"index_count:{proc.stdout.strip()};stderr:{proc.stderr.strip()}"
        return result

    result["status"] = "PASS"
    result["message"] = "core_smoove_genotype_valid;duphold_optional_step_failed"
    return result


def main() -> int:
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    with MANIFEST.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 143:
        raise RuntimeError(f"Expected 143 manifest rows, found {len(rows)}")

    results = []
    with ThreadPoolExecutor(max_workers=12) as pool:
        futures = {pool.submit(validate, row): row["sample_id"] for row in rows}
        for future in as_completed(futures):
            results.append(future.result())
    results.sort(key=lambda x: str(x["sample_id"]))

    columns = [
        "sample_id",
        "vcf",
        "vcf_bytes",
        "index_present",
        "gzip_ok",
        "header_ok",
        "sample_ok",
        "record_count",
        "status",
        "message",
    ]
    temp_report = ROOT / "tmp" / f"existing_sv_genotype_validation.{os.getpid()}.tsv"
    temp_report.parent.mkdir(parents=True, exist_ok=True)
    with temp_report.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)
    os.replace(temp_report, REPORT)

    passed = [row for row in results if row["status"] == "PASS"]
    failed = [row for row in results if row["status"] != "PASS"]
    if failed:
        print(f"Validation failed for {len(failed)} of {len(results)} samples", file=sys.stderr)
        for row in failed[:20]:
            print(f"{row['sample_id']}\t{row['message']}", file=sys.stderr)
        return 1

    timestamp = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z")
    for row in passed:
        sid = str(row["sample_id"])
        done = STAGE / "joint_genotype_per_sample" / sid / f"{sid}.done"
        done.write_text(
            f"completed\t{timestamp}\n"
            f"host\tvalidated_existing_output\n"
            f"sample\t{sid}\n"
            "note\tcore_smoove_genotype_valid;duphold_optional_step_skipped_due_to_abi_error\n"
        )

    counts = sorted({int(row["record_count"]) for row in passed})
    SUMMARY.write_text(
        f"completed\t{timestamp}\n"
        f"samples_validated\t{len(passed)}\n"
        f"record_count_min\t{min(counts)}\n"
        f"record_count_max\t{max(counts)}\n"
        f"report\t{REPORT}\n"
    )
    print(f"Validated {len(passed)} existing SV genotype VCFs")
    print(f"Record count range: {min(counts)}-{max(counts)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
