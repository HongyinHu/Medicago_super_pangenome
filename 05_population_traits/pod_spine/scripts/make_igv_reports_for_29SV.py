#!/usr/bin/env python3
"""Build igv-reports HTML pages for pod-spiny SV validation candidates."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
from pathlib import Path


DEFAULT_PACKAGE = Path(
    "path/to/project/"
    "N_4.pod_spiny/02_IGV_validation_29SV"
)
DEFAULT_CREATE_REPORT = Path("path/to/home/anaconda3/bin/create_report")
REF_TO_GFF_NAME = {
    "Msa": "genome_Msa.gff",
    "R108": "genome_R108.gff",
}


def read_fai(fai_path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with fai_path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            lengths[parts[0]] = int(parts[1])
    return lengths


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def write_marker_bed(path: Path, chrom: str, start_1based: int, end_1based: int, svtype: str, svlen: str) -> None:
    left = min(start_1based, end_1based)
    right = max(start_1based, end_1based)
    if right <= left:
        right = left + 1

    def one_base_marker(pos_1based: int) -> tuple[int, int]:
        start0 = max(pos_1based - 1, 0)
        return start0, start0 + 1

    with path.open("w") as out:
        out.write(f"{chrom}\t{max(left - 1, 0)}\t{right}\t{svtype}_len{svlen}_Area\t0\t.\t"
                  f"{max(left - 1, 0)}\t{right}\t170,170,170\n")
        s0, e0 = one_base_marker(start_1based)
        out.write(f"{chrom}\t{s0}\t{e0}\tSTART\t1000\t.\t{s0}\t{e0}\t255,0,0\n")
        if abs(end_1based - start_1based) > 5:
            s1, e1 = one_base_marker(end_1based)
            out.write(f"{chrom}\t{s1}\t{e1}\tEND\t1000\t.\t{s1}\t{e1}\t0,0,255\n")


def write_single_vcf(path: Path, contig_lengths: dict[str, int], row: dict[str, str]) -> None:
    chrom = row["chrom"]
    pos = int(row["start"])
    end = int(row["end"])
    svtype = row["SVTYPE"] or "SV"
    svlen = row["SVLEN"] or str(abs(end - pos))
    alt = f"<{svtype}>"
    info = (
        f"END={end};SVTYPE={svtype};SVLEN={svlen};"
        f"METASV={row['MetaSV_ID']};ASSOC={row['association_direction']};"
        f"GENE={row.get('best_gene_id', '')}"
    )
    vid = safe_name(f"{row['candidate_id']}|{row['MetaSV_ID']}|{svtype}|len{svlen}|{chrom}:{pos}-{end}")
    with path.open("w") as out:
        out.write("##fileformat=VCFv4.2\n")
        for contig, length in contig_lengths.items():
            out.write(f"##contig=<ID={contig},length={length}>\n")
        out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n')
        out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n')
        out.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">\n')
        out.write('##INFO=<ID=METASV,Number=1,Type=String,Description="Meta pan-SV ID">\n')
        out.write('##INFO=<ID=ASSOC,Number=1,Type=String,Description="Phenotype association direction">\n')
        out.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Nearest or overlapping gene">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        out.write(f"{chrom}\t{pos}\t{vid}\tN\t{alt}\t.\tPASS\t{info}\n")


def build_track_config(
    candidate_dir: Path,
    out_json: Path,
    samples: list[dict[str, str]],
    gene_model_gff: Path | None = None,
    ref_name: str = "",
) -> None:
    tracks = []
    if gene_model_gff is not None and gene_model_gff.exists():
        tracks.append(
            {
                "name": f"{ref_name} gene model",
                "url": str(gene_model_gff),
                "format": "gff3",
                "type": "annotation",
                "displayMode": "EXPANDED",
                "height": 150,
                "color": "rgb(40,40,40)",
            }
        )
    for sample in samples:
        order = int(sample["order"])
        phenotype = sample["phenotype"]
        sample_id = sample["sample"]
        prefix = f"{order:02d}_{phenotype}_{sample_id}"
        bam = candidate_dir / f"{prefix}.bam"
        bai = candidate_dir / f"{prefix}.bam.bai"
        if not bam.exists():
            raise FileNotFoundError(f"Missing BAM: {bam}")
        if not bai.exists():
            raise FileNotFoundError(f"Missing BAI: {bai}")
        color = "rgb(216,60,60)" if phenotype == "spiny" else "rgb(75,120,210)"
        tracks.append(
            {
                "name": prefix,
                "url": str(bam),
                "indexURL": str(bai),
                "format": "bam",
                "type": "alignment",
                "height": 90,
                "color": color,
                "showCoverage": True,
            }
        )
    with out_json.open("w") as out:
        json.dump(tracks, out, indent=2)


def parse_gff_attrs(attr_text: str) -> dict[str, str]:
    attrs = {}
    for part in attr_text.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
        elif " " in part:
            key, value = part.split(" ", 1)
            value = value.strip('"')
        else:
            continue
        attrs[key] = value
    return attrs


def interval_overlap_bp(start1: int, end1: int, start2: int, end2: int) -> int:
    left = max(start1, start2)
    right = min(end1, end2)
    return max(0, right - left + 1)


def load_gene_models(gff_path: Path) -> dict[str, list[dict[str, object]]]:
    genes_by_chrom: dict[str, list[dict[str, object]]] = {}
    transcript_to_gene: dict[str, str] = {}
    feature_rows: list[tuple[str, str, int, int, str, dict[str, str]]] = []

    with gff_path.open() as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _source, ftype, start, end, _score, strand, _phase, attrs_text = parts
            attrs = parse_gff_attrs(attrs_text)
            start_i, end_i = int(start), int(end)
            feature_rows.append((chrom, ftype, start_i, end_i, strand, attrs))
            if ftype == "gene" and "ID" in attrs:
                gene = {
                    "id": attrs["ID"],
                    "chrom": chrom,
                    "start": start_i,
                    "end": end_i,
                    "strand": strand,
                    "exons": [],
                    "cds": [],
                }
                genes_by_chrom.setdefault(chrom, []).append(gene)
            elif ftype in {"mRNA", "transcript"} and "ID" in attrs and "Parent" in attrs:
                transcript_to_gene[attrs["ID"]] = attrs["Parent"].split(",")[0]

    gene_lookup = {
        gene["id"]: gene
        for chrom_genes in genes_by_chrom.values()
        for gene in chrom_genes
    }

    for chrom, ftype, start_i, end_i, _strand, attrs in feature_rows:
        if ftype not in {"exon", "CDS"}:
            continue
        parents = attrs.get("Parent", "").split(",")
        for parent in parents:
            gene_id = transcript_to_gene.get(parent, parent)
            gene = gene_lookup.get(gene_id)
            if gene is None:
                continue
            if ftype == "exon":
                gene["exons"].append((start_i, end_i, attrs.get("ID", "")))
            elif ftype == "CDS":
                gene["cds"].append((start_i, end_i, attrs.get("ID", "")))

    for chrom in genes_by_chrom:
        genes_by_chrom[chrom].sort(key=lambda item: (int(item["start"]), int(item["end"])))
    return genes_by_chrom


def write_gene_model_subset(gff_path: Path, out_path: Path, chrom: str, start: int, end: int) -> int:
    count = 0
    with gff_path.open() as inp, out_path.open("w") as out:
        out.write("##gff-version 3\n")
        for raw in inp:
            if not raw.strip() or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[0] != chrom:
                continue
            f_start, f_end = int(parts[3]), int(parts[4])
            if interval_overlap_bp(start, end, f_start, f_end) <= 0:
                continue
            out.write(raw)
            count += 1
    return count


def classify_gene_context(
    genes_by_chrom: dict[str, list[dict[str, object]]],
    chrom: str,
    start: int,
    end: int,
    upstream_window: int,
) -> dict[str, str]:
    genes = genes_by_chrom.get(chrom, [])
    event_len = max(1, end - start + 1)
    overlapping = []
    for gene in genes:
        gene_start, gene_end = int(gene["start"]), int(gene["end"])
        gene_bp = interval_overlap_bp(start, end, gene_start, gene_end)
        if gene_bp <= 0:
            continue
        exon_bp = sum(interval_overlap_bp(start, end, a, b) for a, b, _fid in gene["exons"])
        cds_bp = sum(interval_overlap_bp(start, end, a, b) for a, b, _fid in gene["cds"])
        intron_bp = max(0, gene_bp - exon_bp)
        if cds_bp and intron_bp:
            context = "CDS+intron_boundary"
        elif cds_bp:
            context = "CDS"
        elif exon_bp and intron_bp:
            context = "exon+intron_boundary"
        elif exon_bp:
            context = "exon_nonCDS"
        else:
            context = "intron"
        overlapping.append((gene, context, gene_bp, exon_bp, cds_bp, intron_bp))

    if overlapping:
        priority = {
            "CDS": 0,
            "CDS+intron_boundary": 1,
            "exon_nonCDS": 2,
            "exon+intron_boundary": 3,
            "intron": 4,
        }
        overlapping.sort(key=lambda item: (priority.get(item[1], 99), -item[2]))
        gene, context, gene_bp, exon_bp, cds_bp, intron_bp = overlapping[0]
        all_gene_ids = ",".join(str(item[0]["id"]) for item in overlapping)
        return {
            "gene_model_context": context,
            "context_gene_id": str(gene["id"]),
            "context_gene_strand": str(gene["strand"]),
            "context_gene_locus": f'{gene["chrom"]}:{gene["start"]}-{gene["end"]}',
            "context_gene_overlap_bp": str(gene_bp),
            "context_exon_overlap_bp": str(exon_bp),
            "context_CDS_overlap_bp": str(cds_bp),
            "context_intron_overlap_bp": str(intron_bp),
            "context_distance_bp": "0",
            "nearest_gene_relation": "overlap_gene",
            "all_overlapping_genes": all_gene_ids,
            "event_overlap_fraction_in_gene": f"{gene_bp / event_len:.4f}",
        }

    nearest = None
    for gene in genes:
        gene_start, gene_end = int(gene["start"]), int(gene["end"])
        if end < gene_start:
            distance = gene_start - end
            side = "before_gene"
        elif start > gene_end:
            distance = start - gene_end
            side = "after_gene"
        else:
            continue
        if nearest is None or distance < nearest[0]:
            nearest = (distance, side, gene)

    if nearest is None:
        return {
            "gene_model_context": "no_gene_on_contig",
            "context_gene_id": "",
            "context_gene_strand": "",
            "context_gene_locus": "",
            "context_gene_overlap_bp": "0",
            "context_exon_overlap_bp": "0",
            "context_CDS_overlap_bp": "0",
            "context_intron_overlap_bp": "0",
            "context_distance_bp": "",
            "nearest_gene_relation": "",
            "all_overlapping_genes": "",
            "event_overlap_fraction_in_gene": "0.0000",
        }

    distance, side, gene = nearest
    strand = str(gene["strand"])
    if strand == "+":
        relation = "upstream" if side == "before_gene" else "downstream"
    elif strand == "-":
        relation = "upstream" if side == "after_gene" else "downstream"
    else:
        relation = side
    context = relation if distance <= upstream_window else f"intergenic_{relation}"
    return {
        "gene_model_context": context,
        "context_gene_id": str(gene["id"]),
        "context_gene_strand": strand,
        "context_gene_locus": f'{gene["chrom"]}:{gene["start"]}-{gene["end"]}',
        "context_gene_overlap_bp": "0",
        "context_exon_overlap_bp": "0",
        "context_CDS_overlap_bp": "0",
        "context_intron_overlap_bp": "0",
        "context_distance_bp": str(distance),
        "nearest_gene_relation": relation,
        "all_overlapping_genes": "",
        "event_overlap_fraction_in_gene": "0.0000",
    }


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def html_escape(value: object) -> str:
    text = "" if value is None else str(value)
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-dir", type=Path, default=DEFAULT_PACKAGE)
    parser.add_argument("--create-report", type=Path, default=DEFAULT_CREATE_REPORT)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--limit", type=int, default=0, help="Only process first N candidates; 0 means all.")
    parser.add_argument("--no-run", action="store_true", help="Only create inputs and command list.")
    parser.add_argument("--standalone", action="store_true", help="Embed igv.js as well as data in each HTML.")
    parser.add_argument("--flanking", type=int, default=2000)
    parser.add_argument("--gene-model-dir", type=Path, default=None)
    parser.add_argument("--gene-context-window", type=int, default=10000)
    args = parser.parse_args()

    pkg = args.package_dir
    spiny_root = pkg.parent
    gene_model_dir = args.gene_model_dir or (spiny_root / "00_data" / "4.reference_anno")
    candidates = load_tsv(pkg / "config" / "candidates_29.tsv")
    samples = load_tsv(pkg / "config" / "samples_main.tsv")
    support_path = pkg / "candidate_read_support_summary.tsv"
    support_by_id = {}
    if support_path.exists():
        support_by_id = {row["candidate_id"]: row for row in load_tsv(support_path)}
    if args.limit:
        candidates = candidates[: args.limit]

    refs_needed = sorted({row["ref"] for row in candidates})
    gff_by_ref: dict[str, Path] = {}
    genes_by_ref: dict[str, dict[str, list[dict[str, object]]]] = {}
    for ref in refs_needed:
        gff_name = REF_TO_GFF_NAME.get(ref, f"genome_{ref}.gff")
        gff_path = gene_model_dir / gff_name
        if gff_path.exists():
            gff_by_ref[ref] = gff_path
            genes_by_ref[ref] = load_gene_models(gff_path)

    out_root = pkg / "igv_reports_html"
    work_root = out_root / "report_inputs"
    out_root.mkdir(parents=True, exist_ok=True)
    work_root.mkdir(parents=True, exist_ok=True)

    commands = []
    index_rows = []
    context_rows = []
    manifest_path = out_root / "igv_reports_manifest.tsv"
    with manifest_path.open("w") as manifest:
        manifest.write(
            "candidate_id\tMetaSV_ID\tref\tlocus\tSVTYPE\tSVLEN\tassociation_direction\t"
            "html\tvcf\tmarker_bed\ttrack_config\n"
        )
        for row in candidates:
            candidate_id = row["candidate_id"]
            ref = row["ref"]
            ref_fasta = pkg / "references" / f"{ref}.candidate_chroms.fa"
            fai = Path(str(ref_fasta) + ".fai")
            if not ref_fasta.exists() or not fai.exists():
                raise FileNotFoundError(f"Missing reference or FAI for {ref}: {ref_fasta}")
            contig_lengths = read_fai(fai)

            candidate_dir = pkg / "bams" / candidate_id
            if not candidate_dir.exists():
                raise FileNotFoundError(f"Missing candidate BAM directory: {candidate_dir}")

            tag = safe_name(f"{candidate_id}_{row['SVTYPE']}_len{row['SVLEN']}_{row['chrom']}_{row['start']}_{row['end']}")
            cand_work = work_root / tag
            cand_work.mkdir(parents=True, exist_ok=True)
            vcf = cand_work / f"{tag}.vcf"
            marker = cand_work / f"{tag}.mark.bed"
            gene_subset = cand_work / f"{tag}.gene_model.gff3"
            track_json = cand_work / f"{tag}.tracks.json"
            html = out_root / f"{tag}.html"

            write_single_vcf(vcf, contig_lengths, row)
            write_marker_bed(marker, row["chrom"], int(row["start"]), int(row["end"]), row["SVTYPE"], row["SVLEN"])
            gene_context = {
                "gene_model_context": "gene_model_not_found",
                "context_gene_id": "",
                "context_gene_strand": "",
                "context_gene_locus": "",
                "context_gene_overlap_bp": "0",
                "context_exon_overlap_bp": "0",
                "context_CDS_overlap_bp": "0",
                "context_intron_overlap_bp": "0",
                "context_distance_bp": "",
                "nearest_gene_relation": "",
                "all_overlapping_genes": "",
                "event_overlap_fraction_in_gene": "0.0000",
            }
            gene_subset_for_track: Path | None = None
            if ref in gff_by_ref:
                view_start = max(1, int(row["start"]) - args.gene_context_window)
                view_end = int(row["end"]) + args.gene_context_window
                subset_count = write_gene_model_subset(gff_by_ref[ref], gene_subset, row["chrom"], view_start, view_end)
                if subset_count:
                    gene_subset_for_track = gene_subset
                gene_context = classify_gene_context(
                    genes_by_ref[ref],
                    row["chrom"],
                    int(row["start"]),
                    int(row["end"]),
                    args.gene_context_window,
                )
            build_track_config(candidate_dir, track_json, samples, gene_subset_for_track, ref)
            context_rows.append({**row, **gene_context})

            title = f"{candidate_id} {row['MetaSV_ID']} {row['SVTYPE']} len={row['SVLEN']} {row['chrom']}:{row['start']}-{row['end']}"
            cmd = [
                str(args.create_report),
                str(vcf),
                "--fasta",
                str(ref_fasta),
                "--flanking",
                str(args.flanking),
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
            if args.standalone:
                cmd.insert(-2, "--standalone")
            commands.append(cmd)
            manifest.write(
                f"{candidate_id}\t{row['MetaSV_ID']}\t{ref}\t{row['chrom']}:{row['start']}-{row['end']}\t"
                f"{row['SVTYPE']}\t{row['SVLEN']}\t{row['association_direction']}\t"
                f"{html}\t{vcf}\t{marker}\t{track_json}\n"
            )
            support = support_by_id.get(candidate_id, {})
            index_rows.append(
                {
                    "candidate_id": candidate_id,
                    "MetaSV_ID": row["MetaSV_ID"],
                    "ref": ref,
                    "locus": f"{row['chrom']}:{row['start']}-{row['end']}",
                    "SVTYPE": row["SVTYPE"],
                    "SVLEN": row["SVLEN"],
                    "association_direction": row["association_direction"],
                    "best_gene_id": row.get("best_gene_id", ""),
                    "gene_model_context": gene_context.get("gene_model_context", ""),
                    "context_gene_id": gene_context.get("context_gene_id", ""),
                    "context_gene_locus": gene_context.get("context_gene_locus", ""),
                    "context_distance_bp": gene_context.get("context_distance_bp", ""),
                    "support_verdict": support.get("support_verdict", ""),
                    "expected_presence_supported": support.get("expected_presence_supported", ""),
                    "expected_absence_supported": support.get("expected_absence_supported", ""),
                    "html": html,
                }
            )

    cmd_list = out_root / "igv_report_commands.list"
    with cmd_list.open("w") as out:
        for cmd in commands:
            out.write(" ".join(subprocess.list2cmdline([part]) for part in cmd) + "\n")

    context_path = out_root / "candidate_gene_model_context.tsv"
    context_cols = [
        "candidate_id",
        "MetaSV_ID",
        "ref",
        "chrom",
        "start",
        "end",
        "SVTYPE",
        "SVLEN",
        "association_direction",
        "best_gene_id",
        "gene_model_context",
        "context_gene_id",
        "context_gene_strand",
        "context_gene_locus",
        "context_gene_overlap_bp",
        "context_exon_overlap_bp",
        "context_CDS_overlap_bp",
        "context_intron_overlap_bp",
        "context_distance_bp",
        "nearest_gene_relation",
        "all_overlapping_genes",
        "event_overlap_fraction_in_gene",
    ]
    with context_path.open("w", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=context_cols, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(context_rows)

    index_path = out_root / "index.html"
    with index_path.open("w") as out:
        out.write(
            "<!doctype html><html><head><meta charset='utf-8'>"
            "<title>Pod spiny SV IGV reports</title>"
            "<style>body{font-family:Arial,sans-serif;margin:24px}"
            "table{border-collapse:collapse;font-size:13px}"
            "th,td{border:1px solid #ddd;padding:6px 8px;vertical-align:top}"
            "th{background:#f3f5f7} .spiny_present{color:#a30000;font-weight:600}"
            ".IGV_priority_confirm_likely{background:#e7f5e7}"
            ".IGV_manual_review_possible{background:#fff8d8}"
            ".IGV_suspect_low_expected_support,.IGV_suspect_unexpected_absent_group_support{background:#fde8e8}"
            "</style></head><body>\n"
            "<h1>Pod spiny SV IGV reports</h1>\n"
            f"<p>Total reports: {len(index_rows)}. Open a report link to inspect 15 sliced BAM tracks "
            "with the SV interval shaded and breakpoints marked.</p>\n"
            "<table><thead><tr>"
        )
        cols = [
            "report",
            "candidate_id",
            "MetaSV_ID",
            "ref",
            "locus",
            "SVTYPE",
            "SVLEN",
            "association_direction",
            "best_gene_id",
            "gene_model_context",
            "context_gene_id",
            "context_gene_locus",
            "context_distance_bp",
            "support_verdict",
            "expected_presence_supported",
            "expected_absence_supported",
        ]
        for col in cols:
            out.write(f"<th>{html_escape(col)}</th>")
        out.write("</tr></thead><tbody>\n")
        for row in index_rows:
            rel_html = os.path.relpath(row["html"], out_root)
            out.write(f"<tr class='{html_escape(row['support_verdict'])}'>")
            out.write(f"<td><a href='{html_escape(rel_html)}'>open</a></td>")
            for col in cols[1:]:
                cls = " class='spiny_present'" if col == "association_direction" and row[col] == "spiny_present" else ""
                out.write(f"<td{cls}>{html_escape(row.get(col, ''))}</td>")
            out.write("</tr>\n")
        out.write("</tbody></table></body></html>\n")

    if args.no_run:
        print(f"Prepared {len(commands)} IGV report commands")
        print(f"Command list: {cmd_list}")
        print(f"Manifest: {manifest_path}")
        print(f"Gene model context: {context_path}")
        print(f"Index: {index_path}")
        return

    failures = []
    if args.jobs <= 1:
        for cmd in commands:
            print("RUN", " ".join(subprocess.list2cmdline([part]) for part in cmd), flush=True)
            result = subprocess.run(cmd)
            if result.returncode != 0:
                failures.append((cmd, result.returncode))
    else:
        from concurrent.futures import ThreadPoolExecutor, as_completed

        def run_one(cmd: list[str]) -> tuple[list[str], int]:
            result = subprocess.run(cmd)
            return cmd, result.returncode

        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(run_one, cmd) for cmd in commands]
            for future in as_completed(futures):
                cmd, code = future.result()
                print(("OK" if code == 0 else f"FAIL={code}"), cmd[-1], flush=True)
                if code != 0:
                    failures.append((cmd, code))

    if failures:
        print(f"Failed {len(failures)} of {len(commands)} commands", flush=True)
        for cmd, code in failures:
            print(f"FAIL {code}: {' '.join(subprocess.list2cmdline([part]) for part in cmd)}", flush=True)
        raise SystemExit(1)
    print(f"Generated {len(commands)} IGV HTML reports in {out_root}")
    print(f"Index: {index_path}")


if __name__ == "__main__":
    main()
