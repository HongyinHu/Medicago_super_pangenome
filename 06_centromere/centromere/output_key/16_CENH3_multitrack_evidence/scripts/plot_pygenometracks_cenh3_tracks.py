#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
import subprocess
from pathlib import Path


SPECIES = ["genome_Msa", "genome_474"]
BIN = 100_000
LOCAL_FLANK = 10_000_000


def read_bed(path: Path) -> list[dict]:
    rows = []
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            rows.append(
                {
                    "chrom": p[0],
                    "start": int(float(p[1])),
                    "end": int(float(p[2])),
                    "name": p[3] if len(p) > 3 else ".",
                    "monomer": p[4] if len(p) > 4 else motif_from_name(p[3] if len(p) > 3 else ""),
                }
            )
    return rows


def read_sizes(path: Path) -> dict[str, int]:
    sizes = {}
    with path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            c, s = line.rstrip("\n").split("\t")[:2]
            if c.startswith("Chr"):
                sizes[c] = int(float(s))
    return dict(sorted(sizes.items(), key=lambda kv: int(kv[0].replace("Chr", ""))))


def motif_from_name(name: str) -> str:
    m = re.search(r"TRASH_(\d+)bp", name)
    return m.group(1) if m else "sat"


def read_intervals_from_gff(path: Path, chroms: set[str], feature: str | None = None) -> dict[str, list[tuple[int, int]]]:
    out = {c: [] for c in chroms}
    with path.open() as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith("#"):
                continue
            p = raw.rstrip("\n").split("\t")
            if len(p) < 5 or p[0] not in chroms:
                continue
            if feature is not None and p[2] != feature:
                continue
            out[p[0]].append((int(p[3]) - 1, int(p[4])))
    return out


def read_sat_bed(path: Path, chroms: set[str]) -> dict[str, list[dict]]:
    out = {c: [] for c in chroms}
    for row in read_bed(path):
        if row["chrom"] in chroms:
            out[row["chrom"]].append(row)
    return out


def overlap(a: tuple[int, int], b: tuple[int, int]) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    intervals = sorted(intervals)
    if not intervals:
        return []
    out = [intervals[0]]
    for s, e in intervals[1:]:
        ls, le = out[-1]
        if s <= le:
            out[-1] = (ls, max(le, e))
        else:
            out.append((s, e))
    return out


def binned_coverage(intervals: list[tuple[int, int]], chrom: str, chrom_size: int, bin_size: int = BIN):
    n = math.ceil(chrom_size / bin_size)
    cov = [0] * n
    merged = merge_intervals([(max(0, s), min(chrom_size, e)) for s, e in intervals if e > 0 and s < chrom_size])
    for s, e in merged:
        b0 = max(0, s // bin_size)
        b1 = min(n - 1, (max(s, e - 1)) // bin_size)
        for b in range(b0, b1 + 1):
            bs = b * bin_size
            be = min((b + 1) * bin_size, chrom_size)
            cov[b] += overlap((bs, be), (s, e))
    for b, val in enumerate(cov):
        bs = b * bin_size
        be = min((b + 1) * bin_size, chrom_size)
        yield chrom, bs, be, val / max(1, be - bs)


def binned_counts(intervals: list[tuple[int, int]], chrom: str, chrom_size: int, bin_size: int = BIN):
    n = math.ceil(chrom_size / bin_size)
    counts = [0] * n
    for s, e in intervals:
        if e <= 0 or s >= chrom_size:
            continue
        ss, ee = max(0, s), min(chrom_size, e)
        b0 = max(0, ss // bin_size)
        b1 = min(n - 1, max(ss, ee - 1) // bin_size)
        for b in range(b0, b1 + 1):
            counts[b] += 1
    for b, val in enumerate(counts):
        bs = b * bin_size
        be = min((b + 1) * bin_size, chrom_size)
        yield chrom, bs, be, val


def write_bedgraph(path: Path, rows):
    with path.open("w") as out:
        for chrom, start, end, value in rows:
            out.write(f"{chrom}\t{start}\t{end}\t{value:.6g}\n")


def write_ab_bedgraph(path: Path, ab_tsv: Path, chroms: set[str]):
    with ab_tsv.open(newline="") as fh, path.open("w") as out:
        for row in csv.DictReader(fh, delimiter="\t"):
            chrom = row.get("chrom", "")
            if chrom not in chroms or not row.get("E1"):
                continue
            try:
                out.write(f"{chrom}\t{int(row['start'])}\t{int(row['end'])}\t{float(row['E1']):.6g}\n")
            except ValueError:
                pass


def region_satellite_families(rows: list[dict], start: int, end: int, n: int = 4) -> list[str]:
    cov = {}
    for r in rows:
        ov = overlap((start, end), (r["start"], r["end"]))
        if ov:
            cov[str(r["monomer"])] = cov.get(str(r["monomer"]), 0) + ov
    families = [k for k, _ in sorted(cov.items(), key=lambda kv: (-kv[1], kv[0]))[:n]]
    while len(families) < n:
        families.append("")
    return families


def write_region_beds(region_dir: Path, chrom: str, start: int, end: int, cen: dict, sat_rows: list[dict]):
    region_dir.mkdir(parents=True, exist_ok=True)
    chrom_bar = region_dir / "chrom_bar.bed"
    with chrom_bar.open("w") as out:
        out.write(f"{chrom}\t{start}\t{end}\t{chrom}\n")
    cen_bed = region_dir / "cenh3_domain.bed"
    with cen_bed.open("w") as out:
        out.write(f"{cen['chrom']}\t{cen['start']}\t{cen['end']}\tCENH3\n")
    sat_paths = []
    for fam in region_satellite_families(sat_rows, start, end):
        label = f"sat{fam}bp" if fam else "satNA"
        bed_path = region_dir / f"{label}.bed"
        with bed_path.open("w") as out:
            if fam:
                for r in sat_rows:
                    if str(r["monomer"]) == fam and r["end"] > start and r["start"] < end:
                        out.write(f"{r['chrom']}\t{max(start, r['start'])}\t{min(end, r['end'])}\t{label}\n")
        sat_paths.append((label, bed_path))
    return chrom_bar, cen_bed, sat_paths


def load_data(project: Path, species: str):
    recall = project / "output_key" / "15_recall_Msa_474_CENH3" / species
    hic = project / "output_all" / "08_hic_AB_TAD" / species
    repeat_root = project / "output_all" / "03_repeat"
    cen = read_bed(recall / "results" / f"{species}.CENH3.functional_centromere.recall.high_confidence.bed")
    cen = sorted(cen, key=lambda r: int(r["chrom"].replace("Chr", "")))
    chroms = {x["chrom"] for x in cen}
    sizes = read_sizes(hic / f"{species}.chrom.sizes")
    return {
        "recall": recall,
        "hic": hic,
        "cen": cen,
        "chroms": chroms,
        "sizes": sizes,
        "sat": read_sat_bed(repeat_root / f"{species}.TRASH" / f"{species}.TRASH_arrays.bed", chroms),
        "repeat": read_intervals_from_gff(repeat_root / f"{species}.EDTA" / f"{species}.fa.mod.EDTA.TEanno.gff3", chroms),
        "gene": read_intervals_from_gff(project / "output_all" / "01_gff" / f"{species}.gff3", chroms, feature="gene"),
    }


def prepare_species_inputs(project: Path, outdir: Path, species: str, data: dict):
    input_dir = outdir / species / "inputs"
    input_dir.mkdir(parents=True, exist_ok=True)
    repeat_bg = input_dir / f"{species}.repeat_density.100kb.bedGraph"
    gene_bg = input_dir / f"{species}.gene_count.100kb.bedGraph"
    ab_bg = input_dir / f"{species}.AB_E1.100kb.bedGraph"
    if not repeat_bg.exists():
        rows = []
        for chrom in sorted(data["chroms"], key=lambda x: int(x.replace("Chr", ""))):
            rows.extend(binned_coverage(data["repeat"].get(chrom, []), chrom, data["sizes"][chrom]))
        write_bedgraph(repeat_bg, rows)
    if not gene_bg.exists():
        rows = []
        for chrom in sorted(data["chroms"], key=lambda x: int(x.replace("Chr", ""))):
            rows.extend(binned_counts(data["gene"].get(chrom, []), chrom, data["sizes"][chrom]))
        write_bedgraph(gene_bg, rows)
    if not ab_bg.exists():
        write_ab_bedgraph(ab_bg, data["hic"] / f"{species}.AB.100kb.cis.vecs.tsv", data["chroms"])
    return {
        "repeat_bg": repeat_bg,
        "gene_bg": gene_bg,
        "ab_bg": ab_bg,
        "tad_score": data["hic"] / f"{species}.TAD.100kb_score.bedgraph",
        "cool": data["hic"] / f"{species}.100000.cool",
    }


def track_bedgraph(name: str, file: Path, title: str, height: float, color: str, extra: str = "") -> str:
    return f"""
[{name}]
file = {file}
title = {title}
height = {height}
color = {color}
min_value = auto
max_value = auto
nans_to_zeros = true
show_data_range = false
summary_method = mean
number_of_bins = 900
file_type = bedgraph
{extra}
"""


def track_bed(name: str, file: Path, title: str, height: float, color: str, extra: str = "") -> str:
    return f"""
[{name}]
file = {file}
title = {title}
height = {height}
color = {color}
border_color = none
display = collapsed
labels = false
file_type = bed
{extra}
"""


def write_tracks_ini(
    ini: Path,
    species: str,
    chrom: str,
    start: int,
    end: int,
    cen: dict,
    prepared: dict,
    region_files: tuple[Path, Path, list[tuple[str, Path]]],
    mode: str,
):
    recall = prepared["recall"] if "recall" in prepared else None
    depth = max(1_000_000, min((end - start) // 2, 30_000_000 if mode == "global" else 12_000_000))
    raw = recall / "bedgraph" / f"{species}.CENH3_vs_Input.raw.10k.smooth50k.CPM.log2.bedGraph"
    relaxed = recall / "bedgraph" / f"{species}.CENH3_vs_Input.relaxed.10k.smooth50k.CPM.log2.bedGraph"
    q20 = recall / "bedgraph" / f"{species}.CENH3_vs_Input.q20.10k.smooth50k.CPM.log2.bedGraph"
    chrom_bar, cen_bed, sat_paths = region_files
    contents = f"""
[CENH3 raw]
file = {raw}
title = CENH3
height = 1.4
color = #9b9b9b
alpha = 0.45
min_value = -1
max_value = 3.2
nans_to_zeros = true
show_data_range = false
type = line:0.45
summary_method = mean
number_of_bins = 1000
file_type = bedgraph

[CENH3 relaxed]
file = {relaxed}
title = 
height = 1.4
color = #e24a3b
alpha = 0.80
min_value = -1
max_value = 3.2
nans_to_zeros = true
show_data_range = false
type = line:0.65
summary_method = mean
number_of_bins = 1000
overlay_previous = share-y
file_type = bedgraph

[CENH3 q20]
file = {q20}
title = 
height = 1.4
color = #2f6f9f
alpha = 0.60
min_value = -1
max_value = 3.2
nans_to_zeros = true
show_data_range = false
type = line:0.55
summary_method = mean
number_of_bins = 1000
overlay_previous = share-y
file_type = bedgraph
"""
    for label, path in sat_paths:
        contents += track_bed(label, path, label if label != "satNA" else "sat", 0.25, "#8ea9cf")
    contents += track_bed("chrom_bar", chrom_bar, chrom, 0.35, "#bfbfbf")
    contents += track_bed("CENH3_domain", cen_bed, f"CENH3 {cen['start']/1e6:.2f}-{cen['end']/1e6:.2f} Mb", 0.28, "#8199bd")
    contents += """
[x-axis]
where = top
fontsize = 8
"""
    contents += track_bedgraph("repeat_density", prepared["repeat_bg"], "Repeat Density", 0.55, "#cfcfcf", "min_value = 0\nmax_value = 1\ntype = fill")
    contents += track_bedgraph("gene_density_fill", prepared["gene_bg"], "Gene Density", 0.55, "#d4d4d4", "min_value = 0\ntype = fill")
    contents += track_bedgraph(
        "gene_density_line",
        prepared["gene_bg"],
        "Gene Density",
        0.55,
        "#e07a73",
        "min_value = 0\ntype = line:0.45\noverlay_previous = share-y",
    )
    contents += track_bedgraph(
        "AB_E1",
        prepared["ab_bg"],
        "Eigv",
        0.65,
        "#c94737",
        "negative_color = #2f7f9f\ntype = fill",
    )
    contents += track_bedgraph("TAD_score", prepared["tad_score"], "TAD score", 0.45, "#777777", "type = line:0.35")
    contents += f"""
[Hi-C]
file = {prepared["cool"]}
title = TAD
height = 5.8
colormap = ['white', '#fff7bc', '#fec44f', '#d95f0e', '#7f2704']
transform = log1p
show_masked_bins = false
depth = {depth}
file_type = hic_matrix
"""
    ini.write_text(contents.strip() + "\n")


def run_pygenometracks(ini: Path, region: str, out_prefix: Path, formats: list[str], width: int, dpi: int):
    for ext in formats:
        out = Path(f"{out_prefix}.{ext}")
        cmd = [
            "pyGenomeTracks",
            "--tracks",
            str(ini),
            "--region",
            region,
            "--outFileName",
            str(out),
            "--width",
            str(width),
            "--trackLabelFraction",
            "0.16",
            "--trackLabelHAlign",
            "right",
            "--fontSize",
            "8",
            "--dpi",
            str(dpi),
        ]
        subprocess.run(cmd, check=True)


def draw_species(project: Path, outdir: Path, species: str, chrom_filter: str | None, views: list[str], formats: list[str], run: bool):
    data = load_data(project, species)
    prepared = prepare_species_inputs(project, outdir, species, data)
    prepared["recall"] = data["recall"]
    for cen in data["cen"]:
        chrom = cen["chrom"]
        if chrom_filter and chrom != chrom_filter:
            continue
        for view in views:
            if view == "global":
                start, end = 0, data["sizes"][chrom]
            else:
                mid = (cen["start"] + cen["end"]) // 2
                start, end = max(0, mid - LOCAL_FLANK), min(data["sizes"][chrom], mid + LOCAL_FLANK)
            region_dir = outdir / species / "tracks" / f"{chrom}.{view}"
            region_files = write_region_beds(region_dir, chrom, start, end, cen, data["sat"].get(chrom, []))
            ini = region_dir / f"{species}.{chrom}.{view}.tracks.ini"
            write_tracks_ini(ini, species, chrom, start, end, cen, prepared, region_files, view)
            out_prefix = outdir / species / f"{species}.{chrom}.pygt_{view}"
            region = f"{chrom}:{start}-{end}"
            if run:
                run_pygenometracks(ini, region, out_prefix, formats, width=36, dpi=450)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", default="path/to/project/10.centromere_analysis")
    ap.add_argument("--outdir", default="path/to/project/10.centromere_analysis/output_key/16_CENH3_multitrack_evidence/pygenometracks_style")
    ap.add_argument("--species", nargs="+", default=SPECIES)
    ap.add_argument("--chrom", default=None)
    ap.add_argument("--views", nargs="+", default=["local"])
    ap.add_argument("--formats", nargs="+", default=["png", "pdf", "svg"])
    ap.add_argument("--write-only", action="store_true")
    args = ap.parse_args()
    project = Path(args.project)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for species in args.species:
        (outdir / species).mkdir(parents=True, exist_ok=True)
        draw_species(project, outdir, species, args.chrom, args.views, args.formats, run=not args.write_only)


if __name__ == "__main__":
    main()
