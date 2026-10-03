#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import subprocess
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
CALLERS = ("pbsv", "sniffles2", "cutesv")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def as_int(value: object, default: int = 0) -> int:
    try:
        if value in (None, "", "NA", "."):
            return default
        return int(float(str(value)))
    except Exception:
        return default


def parse_info(info: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            key, val = item.split("=", 1)
            parsed[key] = val
        else:
            parsed[item] = "1"
    return parsed


def run_capture(cmd: list[str]) -> str:
    return subprocess.check_output(cmd, text=True, stderr=subprocess.DEVNULL)


def region(chrom: str, start: int, end: int) -> str:
    start = max(1, start)
    end = max(start, end)
    return f"{chrom}:{start}-{end}"


def depth_stats(samtools: str, bam: Path, chrom: str, start: int, end: int) -> dict[str, object]:
    if start > end:
        return {"mean": "NA", "covered_prop": "NA", "bases": 0, "covered_bases": 0}
    reg = region(chrom, start, end)
    n_bases = end - start + 1
    try:
        out = run_capture([samtools, "depth", "-aa", "-r", reg, str(bam)])
    except subprocess.CalledProcessError:
        return {"mean": "NA", "covered_prop": "NA", "bases": n_bases, "covered_bases": 0}
    total = 0
    covered = 0
    observed = 0
    for line in out.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        depth = as_int(parts[2])
        total += depth
        observed += 1
        if depth > 0:
            covered += 1
    # samtools depth -aa may omit invalid regions; keep denominator as requested span.
    denom = n_bases if n_bases > 0 else max(1, observed)
    return {
        "mean": round(total / denom, 4),
        "covered_prop": round(covered / denom, 4),
        "bases": n_bases,
        "covered_bases": covered,
    }


def cigar_deletions(cigar: str, ref_start: int) -> list[tuple[int, int, int]]:
    ref_pos = ref_start
    dels: list[tuple[int, int, int]] = []
    for length_s, op in CIGAR_RE.findall(cigar):
        length = int(length_s)
        if op in ("M", "=", "X"):
            ref_pos += length
        elif op == "D":
            del_start = ref_pos
            del_end = ref_pos + length - 1
            dels.append((del_start, del_end, length))
            ref_pos += length
        elif op in ("N",):
            ref_pos += length
        elif op in ("I", "S", "H", "P"):
            continue
    return dels


def cigar_query_span(cigar: str) -> int:
    span = 0
    for length_s, op in CIGAR_RE.findall(cigar):
        if op in ("M", "=", "X", "I", "S"):
            span += int(length_s)
    return span


def overlap_len(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def read_support_stats(samtools: str, bam: Path, chrom: str, start: int, end: int, pad: int) -> dict[str, object]:
    reg = region(chrom, start - pad, end + pad)
    target_len = max(1, end - start + 1)
    try:
        out = run_capture([samtools, "view", str(bam), reg])
    except subprocess.CalledProcessError:
        return {
            "reads_in_window": 0,
            "del_cigar_support_reads": 0,
            "max_del_overlap": 0,
            "max_del_len": 0,
            "softclip_near_breakpoint_reads": 0,
            "spanning_reference_reads": 0,
        }
    reads = 0
    del_support = 0
    max_del_overlap = 0
    max_del_len = 0
    softclip = 0
    ref_span = 0
    seen_names = set()
    for line in out.splitlines():
        if not line or line.startswith("@"):
            continue
        parts = line.split("\t")
        if len(parts) < 6:
            continue
        qname = parts[0]
        flag = as_int(parts[1])
        if flag & 0x900:
            continue
        rstart = as_int(parts[3])
        cigar = parts[5]
        if cigar == "*":
            continue
        reads += 1
        if qname not in seen_names:
            seen_names.add(qname)
        has_support = False
        for ds, de, dlen in cigar_deletions(cigar, rstart):
            ov = overlap_len(ds, de, start, end)
            max_del_overlap = max(max_del_overlap, ov)
            max_del_len = max(max_del_len, dlen)
            if ov >= min(30, target_len * 0.5) or (dlen >= 50 and ov >= target_len * 0.25):
                has_support = True
        if has_support:
            del_support += 1
        # Simple breakpoint soft-clip proxy.
        ops = CIGAR_RE.findall(cigar)
        if ops:
            left_clip = int(ops[0][0]) if ops[0][1] == "S" else 0
            right_clip = int(ops[-1][0]) if ops[-1][1] == "S" else 0
            # Approximate reference end.
            ref_consumed = 0
            for length_s, op in ops:
                if op in ("M", "=", "X", "D", "N"):
                    ref_consumed += int(length_s)
            rend = rstart + ref_consumed - 1
            if (left_clip >= 20 and abs(rstart - start) <= pad) or (right_clip >= 20 and abs(rend - end) <= pad):
                softclip += 1
        ref_consumed_no_del = 0
        for length_s, op in CIGAR_RE.findall(cigar):
            if op in ("M", "=", "X", "N"):
                ref_consumed_no_del += int(length_s)
        rend = rstart + ref_consumed_no_del - 1
        if rstart <= start and rend >= end:
            ref_span += 1
    return {
        "reads_in_window": reads,
        "unique_read_names": len(seen_names),
        "del_cigar_support_reads": del_support,
        "max_del_overlap": max_del_overlap,
        "max_del_len": max_del_len,
        "softclip_near_breakpoint_reads": softclip,
        "spanning_reference_reads": ref_span,
    }


def open_vcf(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", errors="replace")
    return path.open(errors="replace")


def parse_vcf_lines(lines: list[str], chrom: str, qstart: int, qend: int) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for line in lines:
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 8:
            continue
        rchrom = parts[0]
        if rchrom != chrom:
            continue
        pos = as_int(parts[1])
        info = parse_info(parts[7])
        svtype = info.get("SVTYPE", "")
        alt = parts[4]
        if not svtype and alt.startswith("<") and alt.endswith(">"):
            svtype = alt.strip("<>")
        rend = as_int(info.get("END"), pos)
        svlen_raw = info.get("SVLEN", "")
        if "," in svlen_raw:
            svlen_raw = svlen_raw.split(",", 1)[0]
        svlen = abs(as_int(svlen_raw, rend - pos + 1))
        if rend < qstart or pos > qend:
            continue
        records.append({"id": parts[2], "chrom": rchrom, "pos": pos, "end": rend, "svtype": svtype, "svlen": svlen})
    return records


def vcf_records_near(tabix: str, vcf: Path, chrom: str, start: int, end: int, max_scan_bp: int = 2000) -> list[dict[str, object]]:
    if not vcf.exists():
        return []
    qstart = max(1, start - max_scan_bp)
    qend = end + max_scan_bp
    if Path(str(vcf) + ".tbi").exists():
        try:
            out = run_capture([tabix, str(vcf), region(chrom, qstart, qend)])
        except subprocess.CalledProcessError:
            out = ""
        return parse_vcf_lines(out.splitlines(), chrom, qstart, qend)

    # Last resort for unindexed tiny files only. Avoid full scans of large VCFs.
    if vcf.stat().st_size > 20_000_000:
        return []
    try:
        with open_vcf(vcf) as handle:
            return parse_vcf_lines(list(handle), chrom, qstart, qend)
    except Exception:
        return []


def match_vcf_event(records: list[dict[str, object]], svtype: str, start: int, end: int, svlen: int, max_dist: int, min_ro: float, len_ratio_min: float, len_ratio_max: float) -> list[str]:
    target_len = max(1, abs(svlen) if svlen else end - start + 1)
    matches: list[str] = []
    for rec in records:
        if str(rec.get("svtype", "")).upper() != svtype.upper():
            continue
        rstart = as_int(rec["pos"])
        rend = as_int(rec["end"], rstart)
        rlen = max(1, as_int(rec.get("svlen"), rend - rstart + 1))
        if abs(rstart - start) > max_dist and abs(rend - end) > max_dist:
            continue
        ov = overlap_len(start, end, rstart, rend)
        ro = min(ov / max(1, end - start + 1), ov / max(1, rend - rstart + 1))
        ratio = rlen / target_len
        if ro >= min_ro and len_ratio_min <= ratio <= len_ratio_max:
            matches.append(f"{rec.get('id','.')}:pos={rstart}-{rend}:len={rlen}:ro={ro:.3f}:ratio={ratio:.3f}")
    return matches


def classify_sample(row: dict[str, object], min_flank_depth: float, min_flank_cov: float) -> tuple[str, str]:
    left_mean = float(row["left_mean_depth"])
    right_mean = float(row["right_mean_depth"])
    body_mean = float(row["body_mean_depth"])
    left_cov = float(row["left_covered_prop"])
    right_cov = float(row["right_covered_prop"])
    flank_mean = (left_mean + right_mean) / 2 if left_mean >= 0 and right_mean >= 0 else -1
    min_side = min(left_mean, right_mean)
    max_side = max(left_mean, right_mean)
    flank_ratio = min_side / max_side if max_side > 0 else 0
    body_ratio = body_mean / flank_mean if flank_mean > 0 else math.inf
    del_reads = as_int(row["del_cigar_support_reads"])
    vcf_match_callers = as_int(row["vcf_match_callers"])
    ref_span = as_int(row["spanning_reference_reads"])
    if left_mean < min_flank_depth or right_mean < min_flank_depth or left_cov < min_flank_cov or right_cov < min_flank_cov:
        return "NA", "LOW_FLANK_COV"
    if flank_ratio < 0.25:
        return "NA", "IMBALANCED_FLANK"
    if (del_reads >= 3 and body_ratio <= 0.65) or (vcf_match_callers >= 1 and del_reads >= 2 and body_ratio <= 0.75):
        return "present", f"DEL_supported;del_reads={del_reads};vcf_callers={vcf_match_callers};body_flank_ratio={body_ratio:.3f}"
    if vcf_match_callers >= 2 and body_ratio <= 0.80:
        return "present", f"VCF_multi_caller_supported;vcf_callers={vcf_match_callers};body_flank_ratio={body_ratio:.3f}"
    if del_reads <= 1 and vcf_match_callers == 0 and body_ratio >= 0.75 and ref_span >= 2:
        return "absent", f"reference_like;ref_span={ref_span};body_flank_ratio={body_ratio:.3f}"
    return "ambiguous", f"borderline;del_reads={del_reads};vcf_callers={vcf_match_callers};ref_span={ref_span};body_flank_ratio={body_ratio:.3f}"


def sample_expected(sample: dict[str, str], direction: str) -> str:
    spiny = sample.get("spiny_primary", "NA")
    if spiny == "NA":
        return "NA"
    if direction == "spiny_present":
        return "present" if spiny == "1" else "absent"
    if direction == "spineless_present":
        return "present" if spiny == "0" else "absent"
    return "NA"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--samples", type=Path, required=True)
    ap.add_argument("--bam-root", type=Path, required=True)
    ap.add_argument("--vcf-root", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--samtools", default="path/to/home/anaconda3/bin/samtools")
    ap.add_argument("--tabix", default="path/to/home/anaconda3/bin/tabix")
    ap.add_argument("--flank", type=int, default=1000)
    ap.add_argument("--pad", type=int, default=300)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--min-flank-depth", type=float, default=5.0)
    ap.add_argument("--min-flank-covered", type=float, default=0.70)
    ap.add_argument("--max-breakpoint-dist", type=int, default=200)
    ap.add_argument("--min-reciprocal-overlap", type=float, default=0.50)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    candidates = read_tsv(args.candidates)
    samples = read_tsv(args.samples)
    sample_by_id = {r["sample"]: r for r in samples}

    jobs = []
    for cand in candidates:
        ref = cand["ref"]
        if ref not in {"Msa", "RefA"}:
            continue
        chrom = cand["chrom"]
        start = as_int(cand["start"])
        end = as_int(cand["end"])
        svtype = cand["SVTYPE"]
        svlen = abs(as_int(cand["SVLEN"], end - start + 1))
        for sample in samples:
            sid = sample["sample"]
            bam = args.bam_root / "Msa" / "bam" / f"{sid}.sorted.bam"
            if not bam.exists():
                continue
            jobs.append((cand, sample, bam, chrom, start, end, svtype, svlen))

    def process(job):
        cand, sample, bam, chrom, start, end, svtype, svlen = job
        sid = sample["sample"]
        left = depth_stats(args.samtools, bam, chrom, start - args.flank, start - 1)
        body = depth_stats(args.samtools, bam, chrom, start, end)
        right = depth_stats(args.samtools, bam, chrom, end + 1, end + args.flank)
        support = read_support_stats(args.samtools, bam, chrom, start, end, args.pad)
        vcf_hits = {}
        for caller in CALLERS:
            vcf = args.vcf_root / "Msa" / caller / f"{sid}.{caller}.vcf.gz"
            recs = vcf_records_near(args.tabix, vcf, chrom, start, end)
            hits = match_vcf_event(
                recs, svtype, start, end, svlen,
                max_dist=args.max_breakpoint_dist,
                min_ro=args.min_reciprocal_overlap,
                len_ratio_min=0.5,
                len_ratio_max=2.0,
            )
            if hits:
                vcf_hits[caller] = hits
        left_mean = left["mean"] if left["mean"] != "NA" else -1
        right_mean = right["mean"] if right["mean"] != "NA" else -1
        body_mean = body["mean"] if body["mean"] != "NA" else -1
        flank_mean = (float(left_mean) + float(right_mean)) / 2 if float(left_mean) >= 0 and float(right_mean) >= 0 else -1
        body_ratio = float(body_mean) / flank_mean if flank_mean > 0 else "NA"
        out = {
            "candidate_id": cand["candidate_id"],
            "MetaSV_ID": cand["MetaSV_ID"],
            "ref": cand["ref"],
            "chrom": chrom,
            "start": start,
            "end": end,
            "SVTYPE": svtype,
            "SVLEN": cand["SVLEN"],
            "association_direction": cand.get("association_direction", ""),
            "best_gene_id": cand.get("best_gene_id", ""),
            "sample": sid,
            "spiny_primary": sample.get("spiny_primary", "NA"),
            "sample_role": sample.get("sample_role", ""),
            "use_core_screen": sample.get("use_core_screen", "0"),
            "assembly_quality": sample.get("assembly_quality", ""),
            "left_mean_depth": left_mean,
            "left_covered_prop": left["covered_prop"],
            "body_mean_depth": body_mean,
            "body_covered_prop": body["covered_prop"],
            "right_mean_depth": right_mean,
            "right_covered_prop": right["covered_prop"],
            "flank_mean_depth": round(flank_mean, 4) if flank_mean >= 0 else "NA",
            "body_to_flank_depth_ratio": round(body_ratio, 4) if body_ratio != "NA" else "NA",
            **support,
            "vcf_match_callers": len(vcf_hits),
            "vcf_match_caller_names": ",".join(sorted(vcf_hits)) if vcf_hits else "NA",
            "vcf_match_details": ";".join(f"{caller}:" + "|".join(vals[:3]) for caller, vals in sorted(vcf_hits.items())) if vcf_hits else "NA",
        }
        call, reason = classify_sample(out, args.min_flank_depth, args.min_flank_covered)
        out["strict_event_call"] = call
        out["strict_event_reason"] = reason
        out["expected_call"] = sample_expected(sample, cand.get("association_direction", ""))
        if out["expected_call"] in {"present", "absent"} and call in {"present", "absent"}:
            out["matches_expected"] = "1" if out["expected_call"] == call else "0"
        else:
            out["matches_expected"] = "NA"
        return out

    per_sample_rows = []
    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures = [pool.submit(process, job) for job in jobs]
        for fut in as_completed(futures):
            per_sample_rows.append(fut.result())

    per_sample_rows.sort(key=lambda r: (as_int(str(r["candidate_id"]).split("_", 1)[0], 9999), r["sample"]))
    per_fields = [
        "candidate_id", "MetaSV_ID", "ref", "chrom", "start", "end", "SVTYPE", "SVLEN",
        "association_direction", "best_gene_id", "sample", "spiny_primary", "sample_role",
        "use_core_screen", "assembly_quality", "expected_call", "strict_event_call",
        "matches_expected", "strict_event_reason", "left_mean_depth", "left_covered_prop",
        "body_mean_depth", "body_covered_prop", "right_mean_depth", "right_covered_prop",
        "flank_mean_depth", "body_to_flank_depth_ratio", "reads_in_window", "unique_read_names",
        "del_cigar_support_reads", "max_del_overlap", "max_del_len",
        "softclip_near_breakpoint_reads", "spanning_reference_reads", "vcf_match_callers",
        "vcf_match_caller_names", "vcf_match_details",
    ]
    write_tsv(args.out_dir / "per_sample_event_flankQC.tsv", per_sample_rows, per_fields)

    by_candidate: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in per_sample_rows:
        by_candidate[str(row["candidate_id"])].append(row)
    summary_rows: list[dict[str, object]] = []
    for cid, rows in sorted(by_candidate.items(), key=lambda kv: as_int(kv[0].split("_", 1)[0], 9999)):
        core = [r for r in rows if str(r["use_core_screen"]) == "1" and str(r["sample_role"]) not in {"ref_self", "exclude"}]
        sens = [r for r in rows if str(r["use_core_screen"]) != "1" or str(r["sample_role"]) in {"ref_self", "exclude", "sensitivity"}]
        counts = Counter(str(r["strict_event_call"]) for r in core)
        exp_present = [r for r in core if r["expected_call"] == "present"]
        exp_absent = [r for r in core if r["expected_call"] == "absent"]
        present_ok = sum(r["strict_event_call"] == "present" for r in exp_present)
        absent_ok = sum(r["strict_event_call"] == "absent" for r in exp_absent)
        conflicts = [r for r in core if r["matches_expected"] == "0"]
        unknown = [r for r in core if r["strict_event_call"] in {"NA", "ambiguous"} or r["expected_call"] == "NA"]
        low_flank = [r for r in core if "LOW_FLANK" in str(r["strict_event_reason"])]
        if conflicts:
            verdict = "fail_conflict"
        elif unknown:
            verdict = "needs_manual_review_missing_or_ambiguous"
        elif present_ok == len(exp_present) and absent_ok == len(exp_absent):
            verdict = "strict_pass"
        else:
            verdict = "needs_manual_review"
        first = rows[0]
        summary_rows.append({
            "candidate_id": cid,
            "MetaSV_ID": first["MetaSV_ID"],
            "ref": first["ref"],
            "locus": f"{first['chrom']}:{first['start']}-{first['end']}",
            "SVTYPE": first["SVTYPE"],
            "SVLEN": first["SVLEN"],
            "association_direction": first["association_direction"],
            "best_gene_id": first["best_gene_id"],
            "core_expected_present_ok": f"{present_ok}/{len(exp_present)}",
            "core_expected_absent_ok": f"{absent_ok}/{len(exp_absent)}",
            "core_present_calls": counts["present"],
            "core_absent_calls": counts["absent"],
            "core_ambiguous_calls": counts["ambiguous"],
            "core_NA_calls": counts["NA"],
            "core_conflict_samples": ",".join(str(r["sample"]) + ":" + str(r["strict_event_call"]) for r in conflicts) if conflicts else "NA",
            "core_unknown_samples": ",".join(str(r["sample"]) + ":" + str(r["strict_event_call"]) for r in unknown) if unknown else "NA",
            "core_low_flank_samples": ",".join(str(r["sample"]) for r in low_flank) if low_flank else "NA",
            "sensitivity_calls": ";".join(str(r["sample"]) + ":" + str(r["strict_event_call"]) for r in sens) if sens else "NA",
            "strict_flankQC_verdict": verdict,
        })
    sum_fields = [
        "candidate_id", "MetaSV_ID", "ref", "locus", "SVTYPE", "SVLEN", "association_direction",
        "best_gene_id", "core_expected_present_ok", "core_expected_absent_ok",
        "core_present_calls", "core_absent_calls", "core_ambiguous_calls", "core_NA_calls",
        "core_conflict_samples", "core_unknown_samples", "core_low_flank_samples",
        "sensitivity_calls", "strict_flankQC_verdict",
    ]
    write_tsv(args.out_dir / "candidate_strict_flankQC_summary.tsv", summary_rows, sum_fields)
    write_tsv(args.out_dir / "candidate_strict_pass.tsv", [r for r in summary_rows if r["strict_flankQC_verdict"] == "strict_pass"], sum_fields)
    write_tsv(args.out_dir / "candidate_manual_review.tsv", [r for r in summary_rows if r["strict_flankQC_verdict"].startswith("needs")], sum_fields)
    write_tsv(args.out_dir / "candidate_fail_conflict.tsv", [r for r in summary_rows if r["strict_flankQC_verdict"] == "fail_conflict"], sum_fields)


if __name__ == "__main__":
    main()
