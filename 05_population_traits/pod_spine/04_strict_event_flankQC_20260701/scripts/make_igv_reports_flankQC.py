#!/usr/bin/env python3
"""Build IGV HTML reports for strict flank-QC pod-spiny SV candidates."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


DEFAULT_RUN = Path(
    "path/to/project/"
    "N_4.pod_spiny/04_strict_event_flankQC_20260701"
)
DEFAULT_N4 = Path(
    "path/to/project/N_4.pod_spiny"
)
DEFAULT_MAIN = Path(
    "path/to/project/"
    "N_3.call_SV/01.read_based_dualref_hifi"
)
DEFAULT_SAMTOOLS = Path("path/to/home/anaconda3/bin/samtools")
DEFAULT_CREATE_REPORT = Path("path/to/home/anaconda3/bin/create_report")

REF_TO_GFF = {
    "Msa": "genome_Msa.gff",
    "R108": "genome_R108.gff",
}


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fai(path: Path) -> dict[str, int]:
    out: dict[str, int] = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            out[parts[0]] = int(parts[1])
    return out


def parse_locus(locus: str) -> tuple[str, int, int]:
    chrom, span = locus.split(":", 1)
    start_s, end_s = span.split("-", 1)
    return chrom, int(start_s), int(end_s)


def ensure_reference_links(run_dir: Path, n4_dir: Path) -> None:
    ref_dir = run_dir / "references"
    ref_dir.mkdir(parents=True, exist_ok=True)
    src_dir = n4_dir / "02_IGV_validation_29SV" / "references"
    for ref in ("Msa", "R108"):
        for suffix in (".candidate_chroms.fa", ".candidate_chroms.fa.fai"):
            src = src_dir / f"{ref}{suffix}"
            dst = ref_dir / f"{ref}{suffix}"
            if dst.exists():
                continue
            if not src.exists():
                raise FileNotFoundError(src)
            os.symlink(src, dst)


def sample_track_info(row: dict[str, str], order: int) -> dict[str, str]:
    sample = row["sample"]
    role = row.get("sample_role", "")
    primary = row.get("spiny_primary", "NA")
    quality = row.get("assembly_quality", "")
    if role.startswith("core_spiny"):
        group = "spiny"
        color = "rgb(216,60,60)"
    elif role.startswith("core_spineless"):
        group = "spineless"
        color = "rgb(75,120,210)"
    elif role == "sensitivity_weak_spiny":
        group = "weak_spiny"
        color = "rgb(230,150,40)"
    elif role == "sensitivity_spiny_lowQ":
        group = "spiny_lowQ"
        color = "rgb(230,120,40)"
    elif role == "extra_accession":
        group = "extra_spiny"
        color = "rgb(190,70,70)"
    elif role == "ref_self":
        group = "ref_self"
        color = "rgb(120,120,120)"
    elif role.startswith("exclude"):
        group = "excluded_lowQ"
        color = "rgb(150,150,150)"
    elif primary == "1":
        group = "spiny"
        color = "rgb(216,60,60)"
    elif primary == "0":
        group = "spineless"
        color = "rgb(75,120,210)"
    else:
        group = "ambiguous"
        color = "rgb(230,150,40)"
    label = safe_name(f"{order:02d}_{group}_{sample}")
    return {
        "sample": sample,
        "role": role,
        "quality": quality,
        "group": group,
        "label": label,
        "color": color,
    }


def write_marker_bed(path: Path, chrom: str, start: int, end: int, svtype: str, svlen: str) -> None:
    left = min(start, end)
    right = max(start, end)
    with path.open("w") as out:
        out.write(
            f"{chrom}\t{max(left - 1, 0)}\t{right}\t{svtype}_len{svlen}_Area\t0\t.\t"
            f"{max(left - 1, 0)}\t{right}\t170,170,170\n"
        )
        s0 = max(start - 1, 0)
        out.write(f"{chrom}\t{s0}\t{s0 + 1}\tSTART\t1000\t.\t{s0}\t{s0 + 1}\t255,0,0\n")
        if abs(end - start) > 5:
            e0 = max(end - 1, 0)
            out.write(f"{chrom}\t{e0}\t{e0 + 1}\tEND\t1000\t.\t{e0}\t{e0 + 1}\t0,0,255\n")


def write_single_vcf(path: Path, contig_lengths: dict[str, int], row: dict[str, str]) -> None:
    chrom = row["chrom"]
    start = int(row["start"])
    end = int(row["end"])
    svtype = row["SVTYPE"] or "SV"
    svlen = row["SVLEN"] or str(abs(end - start))
    vid = safe_name(
        f"{row['candidate_id']}|{row['MetaSV_ID']}|{svtype}|len{svlen}|{chrom}:{start}-{end}"
    )
    info = (
        f"END={end};SVTYPE={svtype};SVLEN={svlen};"
        f"METASV={row['MetaSV_ID']};ASSOC={row.get('association_direction', '')};"
        f"GENE={row.get('best_gene_id', '')};FLANKQC={row.get('strict_flankQC_verdict', '')}"
    )
    with path.open("w") as out:
        out.write("##fileformat=VCFv4.2\n")
        for contig, length in contig_lengths.items():
            out.write(f"##contig=<ID={contig},length={length}>\n")
        out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n')
        out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n')
        out.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">\n')
        out.write('##INFO=<ID=METASV,Number=1,Type=String,Description="Meta pan-SV ID">\n')
        out.write('##INFO=<ID=ASSOC,Number=1,Type=String,Description="Phenotype association direction">\n')
        out.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Candidate gene">\n')
        out.write('##INFO=<ID=FLANKQC,Number=1,Type=String,Description="Strict flank QC verdict">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        out.write(f"{chrom}\t{start}\t{vid}\tN\t<{svtype}>\t.\tPASS\t{info}\n")


def interval_overlap(start_a: int, end_a: int, start_b: int, end_b: int) -> bool:
    return max(start_a, start_b) <= min(end_a, end_b)


def write_gene_subset(gff: Path, out_path: Path, chrom: str, start: int, end: int) -> int:
    count = 0
    with gff.open() as inp, out_path.open("w") as out:
        out.write("##gff-version 3\n")
        for raw in inp:
            if not raw.strip() or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[0] != chrom:
                continue
            f_start = int(parts[3])
            f_end = int(parts[4])
            if not interval_overlap(start, end, f_start, f_end):
                continue
            out.write(raw)
            count += 1
    return count


def slice_one_bam(
    samtools: Path,
    in_bam: Path,
    out_bam: Path,
    region: str,
    log_path: Path,
) -> tuple[str, bool]:
    out_bam.parent.mkdir(parents=True, exist_ok=True)
    if out_bam.exists() and Path(str(out_bam) + ".bai").exists() and out_bam.stat().st_size > 0:
        return str(out_bam), False
    tmp_bam = out_bam.with_suffix(out_bam.suffix + ".tmp")
    if tmp_bam.exists():
        tmp_bam.unlink()
    with log_path.open("a") as log:
        log.write("+ " + " ".join(map(str, [samtools, "view", "-b", in_bam, region, "-o", tmp_bam])) + "\n")
        log.flush()
        subprocess.run([str(samtools), "view", "-b", str(in_bam), region, "-o", str(tmp_bam)], check=True)
        log.write("+ " + " ".join(map(str, [samtools, "index", tmp_bam])) + "\n")
        log.flush()
        subprocess.run([str(samtools), "index", str(tmp_bam)], check=True)
    os.replace(tmp_bam, out_bam)
    tmp_bai = Path(str(tmp_bam) + ".bai")
    if tmp_bai.exists():
        os.replace(tmp_bai, Path(str(out_bam) + ".bai"))
    elif not Path(str(out_bam) + ".bai").exists():
        with log_path.open("a") as log:
            subprocess.run([str(samtools), "index", str(out_bam)], check=True, stdout=log, stderr=log)
    return str(out_bam), True


def build_track_config(
    track_json: Path,
    gene_gff: Path | None,
    ref: str,
    tracks: list[dict[str, str]],
) -> None:
    config = []
    if gene_gff is not None and gene_gff.exists() and gene_gff.stat().st_size > 20:
        config.append(
            {
                "name": f"{ref} gene model",
                "url": str(gene_gff),
                "format": "gff3",
                "type": "annotation",
                "displayMode": "EXPANDED",
                "height": 150,
                "color": "rgb(40,40,40)",
            }
        )
    for track in tracks:
        config.append(
            {
                "name": track["name"],
                "url": track["bam"],
                "indexURL": track["bai"],
                "format": "bam",
                "type": "alignment",
                "height": 90,
                "color": track["color"],
                "showCoverage": True,
            }
        )
    with track_json.open("w") as out:
        json.dump(config, out, indent=2)


def run_create_report(cmd: list[str], log_path: Path) -> None:
    with log_path.open("a") as log:
        log.write("+ " + subprocess.list2cmdline(cmd) + "\n")
        log.flush()
        subprocess.run(cmd, check=True, stdout=log, stderr=log)


def write_index_html(path: Path, rows: list[dict[str, str]]) -> None:
    by_verdict = {"strict_pass": [], "needs_manual_review_missing_or_ambiguous": [], "fail_conflict": []}
    for row in rows:
        by_verdict.setdefault(row["verdict"], []).append(row)
    with path.open("w") as out:
        out.write("<!doctype html><html><head><meta charset='utf-8'>")
        out.write("<title>04 strict flankQC IGV reports</title>")
        out.write("<style>body{font-family:Arial,sans-serif;margin:24px;}")
        out.write("table{border-collapse:collapse;margin-bottom:28px;}td,th{border:1px solid #ccc;padding:4px 8px;}")
        out.write(".pass{background:#e8f5e9}.manual{background:#fff8e1}.fail{background:#ffebee}</style>")
        out.write("</head><body><h1>04 strict flankQC IGV reports</h1>")
        out.write("<p>Track order: core spiny, core spineless, then sensitivity/excluded/reference samples.</p>")
        for verdict, group_rows in by_verdict.items():
            if not group_rows:
                continue
            css = "pass" if verdict == "strict_pass" else "manual" if "manual" in verdict else "fail"
            out.write(f"<h2>{verdict} ({len(group_rows)})</h2><table>")
            out.write("<tr><th>#</th><th>candidate</th><th>locus</th><th>SV</th><th>gene</th><th>HTML</th></tr>")
            for i, row in enumerate(group_rows, 1):
                href = Path(row["html"]).name
                out.write(
                    f"<tr class='{css}'><td>{i}</td><td>{row['candidate_id']}</td>"
                    f"<td>{row['locus']}</td><td>{row['SVTYPE']} len={row['SVLEN']}</td>"
                    f"<td>{row.get('best_gene_id', '')}</td><td><a href='{href}'>open</a></td></tr>"
                )
            out.write("</table>")
        out.write("</body></html>\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN)
    parser.add_argument("--n4-dir", type=Path, default=DEFAULT_N4)
    parser.add_argument("--main-sv-dir", type=Path, default=DEFAULT_MAIN)
    parser.add_argument("--samtools", type=Path, default=DEFAULT_SAMTOOLS)
    parser.add_argument("--create-report", type=Path, default=DEFAULT_CREATE_REPORT)
    parser.add_argument("--slice-flank", type=int, default=10000)
    parser.add_argument("--report-flanking", type=int, default=2000)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--report-jobs", type=int, default=2)
    parser.add_argument("--only-ref", default="Msa", help="Ref to draw; use ALL for all refs in summary.")
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()

    run_dir = args.run_dir
    out_root = run_dir / "igv_reports_html"
    work_root = out_root / "report_inputs"
    bams_root = run_dir / "bams"
    logs_dir = run_dir / "logs"
    out_root.mkdir(parents=True, exist_ok=True)
    work_root.mkdir(parents=True, exist_ok=True)
    bams_root.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    ensure_reference_links(run_dir, args.n4_dir)

    detail_rows = {row["candidate_id"]: row for row in load_tsv(run_dir / "config" / "candidates_29.tsv")}
    summary_rows = load_tsv(run_dir / "results" / "candidate_strict_flankQC_summary.tsv")
    samples = load_tsv(run_dir / "config" / "samples_extended.tsv")
    selected = []
    only_ref = args.only_ref.strip()
    for row in summary_rows:
        cid = row["candidate_id"]
        if cid not in detail_rows:
            continue
        full = {**detail_rows[cid], **row}
        if only_ref != "ALL" and full["ref"].strip() != only_ref:
            continue
        selected.append(full)
    selected.sort(
        key=lambda row: (
            {"strict_pass": 0, "needs_manual_review_missing_or_ambiguous": 1, "fail_conflict": 2}.get(
                row.get("strict_flankQC_verdict", ""), 9
            ),
            row["candidate_id"],
        )
    )
    if args.limit:
        selected = selected[: args.limit]
    if not selected:
        raise SystemExit("No candidates selected")

    gff_dir = args.n4_dir / "00_data" / "4.reference_anno"
    slice_jobs = []
    report_rows = []
    commands = []
    command_list = out_root / "igv_report_commands.list"
    slice_log = logs_dir / "make_igv_reports_flankQC.slice.log"
    report_log = logs_dir / "make_igv_reports_flankQC.create_report.log"
    slice_log.write_text("")
    report_log.write_text("")

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as executor:
        for row in selected:
            ref = row["ref"]
            chrom = row["chrom"]
            start = int(row["start"])
            end = int(row["end"])
            view_start = max(1, start - args.slice_flank)
            view_end = end + args.slice_flank
            region = f"{chrom}:{view_start}-{view_end}"
            cid = row["candidate_id"]
            for order, sample_row in enumerate(samples, 1):
                track = sample_track_info(sample_row, order)
                in_bam = args.main_sv_dir / "03_per_sample" / ref / "bam" / f"{track['sample']}.sorted.bam"
                out_bam = bams_root / cid / f"{track['label']}.bam"
                if not in_bam.exists():
                    raise FileNotFoundError(in_bam)
                slice_jobs.append(
                    executor.submit(slice_one_bam, args.samtools, in_bam, out_bam, region, slice_log)
                )
        done = 0
        for future in as_completed(slice_jobs):
            future.result()
            done += 1
            if done % 50 == 0:
                print(f"sliced {done}/{len(slice_jobs)} BAM windows", flush=True)

    for row in selected:
        ref = row["ref"]
        chrom = row["chrom"]
        start = int(row["start"])
        end = int(row["end"])
        svtype = row["SVTYPE"]
        svlen = row["SVLEN"]
        verdict = row.get("strict_flankQC_verdict", "")
        cid = row["candidate_id"]
        tag = safe_name(f"{verdict}__{cid}_{svtype}_len{svlen}_{chrom}_{start}_{end}")
        cand_work = work_root / tag
        cand_work.mkdir(parents=True, exist_ok=True)
        ref_fasta = run_dir / "references" / f"{ref}.candidate_chroms.fa"
        contig_lengths = read_fai(Path(str(ref_fasta) + ".fai"))
        vcf = cand_work / f"{tag}.vcf"
        marker = cand_work / f"{tag}.mark.bed"
        gene_subset = cand_work / f"{tag}.gene_model.gff3"
        track_json = cand_work / f"{tag}.tracks.json"
        html = out_root / f"{tag}.html"
        write_single_vcf(vcf, contig_lengths, row)
        write_marker_bed(marker, chrom, start, end, svtype, svlen)
        gff = gff_dir / REF_TO_GFF.get(ref, f"genome_{ref}.gff")
        gene_track = None
        if gff.exists():
            if write_gene_subset(gff, gene_subset, chrom, max(1, start - args.slice_flank), end + args.slice_flank):
                gene_track = gene_subset
        tracks = []
        for order, sample_row in enumerate(samples, 1):
            track = sample_track_info(sample_row, order)
            bam = bams_root / cid / f"{track['label']}.bam"
            tracks.append(
                {
                    "name": track["label"],
                    "bam": str(bam),
                    "bai": str(Path(str(bam) + ".bai")),
                    "color": track["color"],
                }
            )
        build_track_config(track_json, gene_track, ref, tracks)
        title = f"{cid} {row['MetaSV_ID']} {svtype} len={svlen} {chrom}:{start}-{end} {verdict}"
        cmd = [
            str(args.create_report),
            str(vcf),
            "--fasta",
            str(ref_fasta),
            "--flanking",
            str(args.report_flanking),
            "--tracks",
            str(marker),
            "--roi",
            str(marker),
            "--track-config",
            str(track_json),
            "--title",
            title,
            "--output",
            str(html),
        ]
        commands.append(cmd)
        report_rows.append(
            {
                "candidate_id": cid,
                "MetaSV_ID": row["MetaSV_ID"],
                "ref": ref,
                "locus": f"{chrom}:{start}-{end}",
                "SVTYPE": svtype,
                "SVLEN": svlen,
                "best_gene_id": row.get("best_gene_id", ""),
                "verdict": verdict,
                "html": str(html),
                "vcf": str(vcf),
                "marker_bed": str(marker),
                "track_config": str(track_json),
            }
        )

    with command_list.open("w") as out:
        for cmd in commands:
            out.write(subprocess.list2cmdline(cmd) + "\n")

    with ThreadPoolExecutor(max_workers=max(1, args.report_jobs)) as executor:
        futures = [executor.submit(run_create_report, cmd, report_log) for cmd in commands]
        for future in as_completed(futures):
            future.result()

    manifest = out_root / "igv_reports_manifest.tsv"
    fieldnames = [
        "candidate_id",
        "MetaSV_ID",
        "ref",
        "locus",
        "SVTYPE",
        "SVLEN",
        "best_gene_id",
        "verdict",
        "html",
        "vcf",
        "marker_bed",
        "track_config",
    ]
    with manifest.open("w", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(report_rows)
    write_index_html(out_root / "index.html", report_rows)
    print(f"selected_candidates={len(selected)}")
    print(f"sliced_bam_windows={len(slice_jobs)}")
    print(f"html_reports={len(report_rows)}")
    print(f"manifest={manifest}")
    print(f"index={out_root / 'index.html'}")


if __name__ == "__main__":
    main()
