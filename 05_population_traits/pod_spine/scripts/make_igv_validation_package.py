#!/usr/bin/env python3
import csv
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


SPINY = ["genome_410", "genome_474", "genome_Mpo", "genome_R108"]
SPINELESS = [
    "genome_395",
    "genome_454",
    "genome_461",
    "genome_468",
    "genome_472",
    "genome_482",
    "genome_M22",
    "genome_Mar",
    "genome_Mru",
    "genome_Msa",
    "genome_ZM4",
]
SAMPLES = [(s, "spiny") for s in SPINY] + [(s, "spineless") for s in SPINELESS]


def run(cmd, log=None):
    if log:
        log.write("+ " + " ".join(map(str, cmd)) + "\n")
        log.flush()
    subprocess.run(cmd, check=True)


def cigar_deletion_support(cigar, ref_start_1, sv_start, sv_end):
    ref_pos = ref_start_1
    best = 0
    for length_s, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        length = int(length_s)
        if op in "M=X":
            ref_pos += length
        elif op == "D":
            d_start = ref_pos
            d_end = ref_pos + length - 1
            overlap = max(0, min(d_end, sv_end) - max(d_start, sv_start) + 1)
            if overlap > best:
                best = overlap
            ref_pos += length
        elif op in "N":
            ref_pos += length
        elif op in "ISH":
            pass
    return best


def softclip_near_breakpoint(cigar, ref_start_1, breakpoint, window=100):
    ref_len = 0
    ops = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    for length_s, op in ops:
        length = int(length_s)
        if op in "MDN=X":
            ref_len += length
    ref_end = ref_start_1 + max(0, ref_len - 1)
    left_clip = int(ops[0][0]) if ops and ops[0][1] == "S" else 0
    right_clip = int(ops[-1][0]) if ops and ops[-1][1] == "S" else 0
    return (left_clip >= 30 and abs(ref_start_1 - breakpoint) <= window) or (
        right_clip >= 30 and abs(ref_end - breakpoint) <= window
    )


def parse_depth_mean(samtools, bam, chrom, start, end):
    if start > end:
        return 0.0
    region = f"{chrom}:{start}-{end}"
    proc = subprocess.run(
        [samtools, "depth", "-aa", "-r", region, bam],
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    )
    vals = []
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3:
            vals.append(int(parts[2]))
    return round(sum(vals) / len(vals), 3) if vals else 0.0


def support_metrics(samtools, bam, chrom, start, end, svtype):
    event_len = max(1, end - start + 1)
    flank = min(5000, max(1000, event_len * 10))
    left_start = max(1, start - flank)
    left_end = max(1, start - 1)
    right_start = end + 1
    right_end = end + flank
    region = f"{chrom}:{max(1, start - 200)}-{end + 200}"
    proc = subprocess.run(
        [samtools, "view", bam, region],
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    )
    total = 0
    del_support = 0
    clip_support = 0
    min_del_overlap = min(event_len * 0.5, max(30, event_len * 0.8))
    for line in proc.stdout.splitlines():
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 6:
            continue
        flag = int(parts[1])
        if flag & 0x904:
            continue
        total += 1
        pos = int(parts[3])
        cigar = parts[5]
        overlap = cigar_deletion_support(cigar, pos, start, end)
        if svtype == "DEL" and overlap >= min_del_overlap:
            del_support += 1
        if softclip_near_breakpoint(cigar, pos, start) or softclip_near_breakpoint(cigar, pos, end):
            clip_support += 1
    left_cov = parse_depth_mean(samtools, bam, chrom, left_start, left_end)
    event_cov = parse_depth_mean(samtools, bam, chrom, start, end)
    right_cov = parse_depth_mean(samtools, bam, chrom, right_start, right_end)
    flank_cov = round((left_cov + right_cov) / 2, 3)
    cov_ratio = round(event_cov / flank_cov, 3) if flank_cov else "NA"
    return {
        "reads_in_window": total,
        "del_cigar_support_reads": del_support,
        "softclip_near_breakpoint_reads": clip_support,
        "left_cov_mean": left_cov,
        "event_cov_mean": event_cov,
        "right_cov_mean": right_cov,
        "flank_cov_mean": flank_cov,
        "event_to_flank_cov_ratio": cov_ratio,
    }


def support_call(metrics, svtype, phenotype, direction):
    if svtype != "DEL":
        return "manual_review"
    del_reads = int(metrics["del_cigar_support_reads"])
    cov_ratio = metrics["event_to_flank_cov_ratio"]
    ratio_ok = isinstance(cov_ratio, float) and cov_ratio <= 0.65
    present = del_reads >= 3 or ratio_ok
    if direction == "spiny_present":
        expected_present = phenotype == "spiny"
    elif direction == "spiny_absent":
        expected_present = phenotype == "spineless"
    else:
        expected_present = None
    if expected_present is None:
        return "manual_review"
    if expected_present and present:
        return "supports_expected_presence"
    if expected_present and not present:
        return "lacks_expected_presence_support"
    if (not expected_present) and present:
        return "unexpected_presence_support"
    return "supports_expected_absence"


def read_candidates(path):
    out = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("RefA_CHROM") and row.get("RefA_CHROM") not in {".", "NA", ""}:
                ref = "Msa"
                chrom = row["RefA_CHROM"]
                start = int(row["RefA_POS"])
                end = int(row["RefA_END"])
            else:
                ref = "R108"
                chrom = row["RefB_CHROM"]
                start = int(row["RefB_POS"])
                end = int(row["RefB_END"])
            out.append({**row, "igv_ref": ref, "igv_chrom": chrom, "igv_start": start, "igv_end": end})
    return out


def make_candidate_id(rank, metasv):
    return f"{int(rank):02d}_{metasv}"


def main():
    if len(sys.argv) != 8:
        raise SystemExit(
            "usage: make_igv_validation_package.py <candidate_tsv> <outdir> <main_sv_dir> "
            "<refa_fa> <refb_fa> <samtools> <flank>"
        )
    candidate_tsv, outdir, main_sv_dir, refa_fa, refb_fa, samtools, flank_s = sys.argv[1:]
    outdir = Path(outdir)
    flank = int(flank_s)
    outdir.mkdir(parents=True, exist_ok=True)
    for sub in ["config", "bams", "references", "igv_batch_templates", "logs", "snapshots"]:
        (outdir / sub).mkdir(parents=True, exist_ok=True)
    candidates = read_candidates(candidate_tsv)
    refs = {"Msa": refa_fa, "R108": refb_fa}
    ref_chroms = defaultdict(set)
    for c in candidates:
        ref_chroms[c["igv_ref"]].add(c["igv_chrom"])

    with (outdir / "logs" / "make_package.log").open("w") as log:
        for ref, chroms in ref_chroms.items():
            ref_out = outdir / "references" / f"{ref}.candidate_chroms.fa"
            if not ref_out.exists():
                with ref_out.open("wb") as out_f:
                    for chrom in sorted(chroms):
                        log.write(f"+ {samtools} faidx {refs[ref]} {chrom} >> {ref_out}\n")
                        log.flush()
                        seq = subprocess.run([samtools, "faidx", refs[ref], chrom], check=True, stdout=subprocess.PIPE)
                        out_f.write(seq.stdout)
                run([samtools, "faidx", str(ref_out)], log=log)

        candidate_rows = []
        support_rows = []
        for c in candidates:
            cid = make_candidate_id(c["final_validation_rank"], c["MetaSV_ID"])
            chrom = c["igv_chrom"]
            start = c["igv_start"]
            end = c["igv_end"]
            view_start = max(1, start - flank)
            view_end = end + flank
            ref = c["igv_ref"]
            candidate_dir = outdir / "bams" / cid
            candidate_dir.mkdir(parents=True, exist_ok=True)
            candidate_rows.append(
                {
                    "candidate_id": cid,
                    "MetaSV_ID": c["MetaSV_ID"],
                    "ref": ref,
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "view_region": f"{chrom}:{view_start}-{view_end}",
                    "SVTYPE": c["SVTYPE"],
                    "SVLEN": c["SVLEN"],
                    "association_direction": c["association_direction"],
                    "best_gene_id": c.get("best_gene_id", ""),
                    "candidate_locus": c.get("candidate_locus", ""),
                }
            )
            batch = outdir / "igv_batch_templates" / f"{cid}.igv.batch.template.txt"
            with batch.open("w") as b:
                b.write("new\n")
                b.write(f"genome {{PACKAGE_DIR}}/references/{ref}.candidate_chroms.fa\n")
                b.write("maxPanelHeight 900\n")
                b.write(f"snapshotDirectory {{PACKAGE_DIR}}/snapshots/{cid}\n")
                b.write("squish\n")
                b.write(f"goto {chrom}:{view_start}-{view_end}\n")
                for i, (sample, phenotype) in enumerate(SAMPLES, 1):
                    in_bam = Path(main_sv_dir) / "03_per_sample" / ref / "bam" / f"{sample}.sorted.bam"
                    out_bam = candidate_dir / f"{i:02d}_{phenotype}_{sample}.bam"
                    if not out_bam.exists() or not Path(str(out_bam) + ".bai").exists():
                        run([samtools, "view", "-b", str(in_bam), f"{chrom}:{view_start}-{view_end}", "-o", str(out_bam)], log=log)
                        run([samtools, "index", str(out_bam)], log=log)
                    b.write(f"load {{PACKAGE_DIR}}/bams/{cid}/{out_bam.name}\n")
                    metrics = support_metrics(samtools, str(in_bam), chrom, start, end, c["SVTYPE"])
                    support_rows.append(
                        {
                            "candidate_id": cid,
                            "MetaSV_ID": c["MetaSV_ID"],
                            "ref": ref,
                            "sample": sample,
                            "phenotype": phenotype,
                            "expected_pattern": c["association_direction"],
                            "chrom": chrom,
                            "start": start,
                            "end": end,
                            "SVTYPE": c["SVTYPE"],
                            "SVLEN": c["SVLEN"],
                            **metrics,
                            "support_call": support_call(metrics, c["SVTYPE"], phenotype, c["association_direction"]),
                        }
                    )
                b.write(f"goto {chrom}:{view_start}-{view_end}\n")
                b.write(f"snapshot {cid}.png\n")
                b.write("exit\n")

        with (outdir / "config" / "candidates_29.tsv").open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(candidate_rows[0].keys()), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(candidate_rows)
        with (outdir / "config" / "samples_main.tsv").open("w", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(["order", "sample", "phenotype"])
            for i, (sample, phenotype) in enumerate(SAMPLES, 1):
                w.writerow([i, sample, phenotype])
        with (outdir / "read_support_metrics.tsv").open("w", newline="") as f:
            fieldnames = list(support_rows[0].keys())
            w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(support_rows)
        with (outdir / "run_all_igv_batches.template.txt").open("w") as f:
            for c in candidate_rows:
                cid = c["candidate_id"]
                f.write(f"# {cid} {c['candidate_locus']} {c['association_direction']} {c['best_gene_id']}\n")
                f.write(f"batch {{PACKAGE_DIR}}/igv_batch_templates/{cid}.igv.batch.txt\n")
        print(f"PACKAGE\t{outdir}")
        print(f"CANDIDATES\t{len(candidates)}")
        print(f"READ_SUPPORT_ROWS\t{len(support_rows)}")


if __name__ == "__main__":
    main()
