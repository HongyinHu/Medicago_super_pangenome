#!/usr/bin/env python3
import argparse
import csv
import gzip
import os
import subprocess
from collections import Counter, defaultdict
from pathlib import Path


BIN_SIZE = 100000
DIST = 500
INS_DIST = 1000
MIN_RO = 0.5

FINAL_SAMPLES = [
    "genome_395", "genome_410", "genome_436", "genome_454", "genome_457", "genome_461",
    "genome_468", "genome_472", "genome_474", "genome_482", "genome_M22", "genome_M46",
    "genome_Mar", "genome_Mpo", "genome_Mru", "genome_Msa", "genome_R108", "genome_ZM4",
]

DROP_LABELS = {"genome_A17", "A17", "genome_G474b", "G474b", "genome_474b"}
MAP_LABELS = {"genome_G474a": "genome_474", "G474a": "genome_474", "genome_474a": "genome_474"}

INFO_DEFS = [
    '##INFO=<ID=FINAL_SV_ID,Number=1,Type=String,Description="Final per-species merged SV identifier">',
    '##INFO=<ID=SPECIES,Number=1,Type=String,Description="Species/sample ID for this per-species merged SV file">',
    '##INFO=<ID=READ_SV_ID,Number=1,Type=String,Description="Original per-species Jasmine read-mapping consensus SV identifier">',
    '##INFO=<ID=SOURCE_METHODS,Number=1,Type=String,Description="Evidence source category: read_based or read_based+SVGAP">',
    '##INFO=<ID=SVGAP_SUPPORT_COUNT,Number=1,Type=Integer,Description="Number of overlapping SVGAP/genome-synteny SV records supporting this read-based event">',
    '##INFO=<ID=SVGAP_SUPPORT_SAMPLES,Number=.,Type=String,Description="Species with overlapping SVGAP/genome-synteny evidence after sample-name normalization">',
    '##INFO=<ID=FINAL_CONFIDENCE,Number=1,Type=String,Description="Final confidence label for per-species merged SV event">',
]


def split_labels(value):
    if not value or value == ".":
        return []
    out = []
    for token in str(value).replace(";", ",").split(","):
        token = token.strip()
        if token and token != ".":
            out.append(token)
    return out


def normalize_sample(label):
    label = MAP_LABELS.get(label, label)
    if label in DROP_LABELS:
        return None
    if label in FINAL_SAMPLES:
        return label
    return None


def sample_from_query_token(token):
    if not token:
        return None
    base = token.split(".", 1)[0]
    base = base.split("_Chr", 1)[0]
    aliases = {
        "395": "genome_395", "410": "genome_410", "436": "genome_436",
        "454": "genome_454", "457": "genome_457", "461": "genome_461",
        "468": "genome_468", "472": "genome_472", "482": "genome_482",
        "G395": "genome_395", "G410": "genome_410", "G436": "genome_436",
        "G454": "genome_454", "G457": "genome_457", "G461": "genome_461",
        "G468": "genome_468", "G472": "genome_472", "G474a": "genome_474",
        "G474b": None, "G482": "genome_482",
        "A17": None, "M22": "genome_M22", "M46": "genome_M46",
        "Mar": "genome_Mar", "Mpo": "genome_Mpo", "Mru": "genome_Mru",
        "R108": "genome_R108", "ZM4": "genome_ZM4", "Msa": "genome_Msa",
    }
    if base in aliases:
        return aliases[base]
    if base.startswith("genome_"):
        return normalize_sample(base)
    return normalize_sample("genome_" + base)


def samples_from_svgap_tokens(value):
    out = set()
    for token in split_labels(value):
        sample = sample_from_query_token(token)
        if sample:
            out.add(sample)
    return out


def samples_from_merge_samples(value):
    out = set()
    for token in split_labels(value):
        sample = normalize_sample(token)
        if sample:
            out.add(sample)
    out.update(samples_from_svgap_tokens(value))
    return out


def safe_int(value, default=None):
    try:
        return int(float(str(value).split(",", 1)[0]))
    except Exception:
        return default


def norm_chrom(chrom):
    return chrom[4:] if chrom.startswith("Msa.") else chrom


def parse_info(info):
    out = {}
    if not info or info == ".":
        return out
    for item in info.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value
        elif item:
            out[item] = True
    return out


def clean_info_value(value):
    if value is None or value == "" or value == ".":
        return "."
    return str(value).replace(" ", "_").replace(";", ",")


def interval_len(start, end):
    return max(1, int(end) - int(start) + 1)


def reciprocal_overlap(a1, a2, b1, b2):
    ov = max(0, min(a2, b2) - max(a1, b1) + 1)
    if ov <= 0:
        return 0.0
    return min(float(ov) / interval_len(a1, a2), float(ov) / interval_len(b1, b2))


def length_ratio(a, b):
    a = abs(safe_int(a, 0) or 0)
    b = abs(safe_int(b, 0) or 0)
    if a == 0 or b == 0:
        return 1.0
    return min(float(a) / b, float(b) / a)


def bins_for(chrom, start, end):
    if start is None or end is None:
        return []
    start = max(0, int(start))
    end = max(start, int(end))
    return [(chrom, i) for i in range(start // BIN_SIZE, end // BIN_SIZE + 1)]


def match_reason(read_ev, svg_ev):
    if read_ev["SVTYPE"] != svg_ev["SVTYPE"]:
        return None
    if read_ev["CHROM"] != svg_ev["CHROM"]:
        return None
    if read_ev["SVTYPE"] == "INS":
        return "ins_position_distance" if abs(read_ev["POS"] - svg_ev["POS"]) <= INS_DIST else None
    if read_ev["SVTYPE"] == "TRA":
        return "tra_breakpoint_distance" if abs(read_ev["POS"] - svg_ev["POS"]) <= INS_DIST else None
    ro = reciprocal_overlap(read_ev["POS"], read_ev["END"], svg_ev["POS"], svg_ev["END"])
    if ro >= MIN_RO and length_ratio(read_ev["SVLEN"], svg_ev["SVLEN"]) >= 0.5:
        return "reciprocal_overlap"
    if abs(read_ev["POS"] - svg_ev["POS"]) <= DIST or abs(read_ev["END"] - svg_ev["END"]) <= DIST:
        if length_ratio(read_ev["SVLEN"], svg_ev["SVLEN"]) >= 0.5:
            return "breakpoint_distance"
    return None


def iter_vcf_records(vcf_gz):
    with gzip.open(vcf_gz, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                yield line, None
            else:
                yield None, line.rstrip("\n").split("\t")


def read_species_events(consensus_dir, species_list):
    events_by_species = {s: [] for s in species_list}
    index = {s: defaultdict(list) for s in species_list}
    by_read_id = {}
    paths = {}
    for species in species_list:
        path = Path(consensus_dir) / ("%s.Msa.consensus.flankQC.sorted.vcf.gz" % species)
        if not path.exists() or path.stat().st_size == 0:
            raise SystemExit("missing flankQC consensus VCF: %s" % path)
        paths[species] = path
        seq = 0
        for header, parts in iter_vcf_records(path):
            if header is not None:
                continue
            if len(parts) < 8:
                continue
            seq += 1
            info = parse_info(parts[7])
            svtype = info.get("SVTYPE", "NA")
            if svtype == "BND":
                svtype = "TRA"
            pos = safe_int(parts[1], 0) or 0
            end = safe_int(info.get("END", pos), pos) or pos
            if end < pos:
                pos, end = end, pos
            svlen = safe_int(info.get("SVLEN", interval_len(pos, end)), interval_len(pos, end))
            ev = {
                "species": species,
                "FinalSV_ID": "%s_Msa_mergedSV_%09d" % (species, seq),
                "read_SV_ID": parts[2],
                "CHROM": norm_chrom(parts[0]),
                "POS": pos,
                "END": end,
                "SVTYPE": svtype,
                "SVLEN": svlen if svlen is not None else interval_len(pos, end),
                "svgap_count": 0,
                "svgap_ids": [],
                "svgap_samples": set(),
                "svgap_source_ids": [],
            }
            events_by_species[species].append(ev)
            by_read_id[(species, parts[2])] = ev
            for key in bins_for(ev["CHROM"], ev["POS"], ev["END"]):
                index[species][key].append(ev)
    return paths, events_by_species, index, by_read_id


def iter_svgap_events(svg_dir):
    svg_dir = Path(svg_dir)
    files = [
        (svg_dir / "All.DELs.50bplarge.bed.combined.sorted.txt", "combined"),
        (svg_dir / "All.INTs.50bplarge.bed.combined.sorted.txt", "combined"),
        (svg_dir / "All.CNVs.bed", "raw_cnv"),
        (svg_dir / "All.INVs.bed", "raw_inv"),
        (svg_dir / "All.TRLs.bed", "raw_trl"),
    ]
    seq = 0
    for source_file, kind in files:
        if not source_file.exists():
            continue
        with open(source_file, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t", 14)
                try:
                    if kind == "combined":
                        if len(parts) < 8:
                            continue
                        chrom = norm_chrom(parts[0])
                        start = safe_int(parts[1])
                        end = safe_int(parts[2])
                        svtype = parts[4]
                        if svtype == "INT":
                            svtype = "INS"
                        svlen = abs(safe_int(parts[5], 0) or interval_len(start, end))
                        source_id = "%s:%s-%s:%s" % (chrom, start, end, parts[3])
                        support_samples = samples_from_merge_samples(parts[7])
                    elif kind == "raw_cnv":
                        if len(parts) < 8:
                            continue
                        chrom = norm_chrom(parts[0])
                        start = safe_int(parts[1])
                        end = safe_int(parts[2])
                        query = parts[3]
                        svtype = "DUP"
                        svlen = abs(safe_int(parts[7], 0) or interval_len(start, end))
                        source_id = "%s_%s_%s_%s" % (query, parts[5], parts[6], svlen)
                        sample = sample_from_query_token(query)
                        support_samples = {sample} if sample else set()
                    elif kind == "raw_inv":
                        if len(parts) < 8:
                            continue
                        chrom = norm_chrom(parts[0])
                        start = safe_int(parts[1])
                        end = safe_int(parts[2])
                        query = parts[3]
                        svtype = "INV"
                        svlen = abs(safe_int(parts[7], 0) or interval_len(start, end))
                        source_id = "%s_%s_%s_%s" % (query, parts[5], parts[6], svlen)
                        sample = sample_from_query_token(query)
                        support_samples = {sample} if sample else set()
                    elif kind == "raw_trl":
                        if len(parts) < 9:
                            continue
                        chrom = norm_chrom(parts[0])
                        start = safe_int(parts[1])
                        end = safe_int(parts[2])
                        query = parts[3]
                        svtype = "TRA"
                        svlen = max(abs(safe_int(parts[7], 0) or 0), abs(safe_int(parts[8], 0) or 0), interval_len(start, end))
                        source_id = "%s_%s_%s_%s_%s_%s" % (query, parts[5], parts[6], chrom, start, end)
                        sample = sample_from_query_token(query)
                        support_samples = {sample} if sample else set()
                    else:
                        continue
                except Exception:
                    continue
                if start is None or end is None or not support_samples:
                    continue
                if start > end:
                    start, end = end, start
                support_samples = {s for s in support_samples if s in FINAL_SAMPLES}
                if not support_samples:
                    continue
                seq += 1
                yield {
                    "svgap_id": "SVGAP_Msa_%s_%09d" % (svtype, seq),
                    "CHROM": chrom,
                    "POS": start,
                    "END": end,
                    "SVTYPE": svtype,
                    "SVLEN": svlen,
                    "source_file": str(source_file),
                    "source_id": source_id.replace("G474a", "genome_474"),
                    "support_samples": support_samples,
                }


def candidate_keys(svg_ev):
    if svg_ev["SVTYPE"] in ("INS", "TRA"):
        return bins_for(svg_ev["CHROM"], svg_ev["POS"] - INS_DIST, svg_ev["POS"] + INS_DIST)
    return bins_for(svg_ev["CHROM"], svg_ev["POS"], svg_ev["END"])


def add_svgap_support(events_by_species, index, out_dir, species_list, svg_dir):
    evidence_dir = Path(out_dir) / "evidence"
    evidence_dir.mkdir(parents=True, exist_ok=True)
    evidence_fields = [
        "species", "svgap_id", "read_SV_ID", "FinalSV_ID", "chrom",
        "svgap_start", "svgap_end", "read_start", "read_end", "svtype",
        "match_reason", "svgap_source_file", "svgap_source_id", "svgap_support_samples",
    ]
    writers = {}
    handles = {}
    for species in species_list:
        fh = open(evidence_dir / ("%s.Msa.mergedSV.svgap_overlap_evidence.tsv" % species), "w", newline="")
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=evidence_fields)
        writer.writeheader()
        handles[species] = fh
        writers[species] = writer

    seen = Counter()
    svgap_only = Counter()
    evidence_rows = Counter()
    try:
        for svg_ev in iter_svgap_events(svg_dir):
            keys = candidate_keys(svg_ev)
            for species in svg_ev["support_samples"]:
                if species not in index:
                    continue
                seen[species] += 1
                matched = False
                visited = set()
                for key in keys:
                    for read_ev in index[species].get(key, []):
                        rid = id(read_ev)
                        if rid in visited:
                            continue
                        visited.add(rid)
                        reason = match_reason(read_ev, svg_ev)
                        if not reason:
                            continue
                        matched = True
                        read_ev["svgap_count"] += 1
                        if len(read_ev["svgap_ids"]) < 100:
                            read_ev["svgap_ids"].append(svg_ev["svgap_id"])
                        read_ev["svgap_samples"].update(svg_ev["support_samples"])
                        if len(read_ev["svgap_source_ids"]) < 100:
                            read_ev["svgap_source_ids"].append(svg_ev["source_id"])
                        evidence_rows[species] += 1
                        writers[species].writerow({
                            "species": species,
                            "svgap_id": svg_ev["svgap_id"],
                            "read_SV_ID": read_ev["read_SV_ID"],
                            "FinalSV_ID": read_ev["FinalSV_ID"],
                            "chrom": read_ev["CHROM"],
                            "svgap_start": svg_ev["POS"],
                            "svgap_end": svg_ev["END"],
                            "read_start": read_ev["POS"],
                            "read_end": read_ev["END"],
                            "svtype": read_ev["SVTYPE"],
                            "match_reason": reason,
                            "svgap_source_file": svg_ev["source_file"],
                            "svgap_source_id": svg_ev["source_id"],
                            "svgap_support_samples": ",".join(sorted(svg_ev["support_samples"])),
                        })
                if not matched:
                    svgap_only[species] += 1
    finally:
        for fh in handles.values():
            fh.close()
    return seen, svgap_only, evidence_rows


def source_method(ev):
    return "read_based+SVGAP" if ev["svgap_count"] > 0 else "read_based"


def confidence(ev):
    return "read_based+svgap_evidence" if ev["svgap_count"] > 0 else "read_based"


def write_event_tables(out_dir, species_list, events_by_species):
    events_dir = Path(out_dir) / "events"
    events_dir.mkdir(parents=True, exist_ok=True)
    fields = [
        "SpeciesSV_ID", "source_methods", "reference", "species", "CHROM", "POS", "END",
        "SVTYPE", "SVLEN", "read_SV_ID", "svgap_support_count", "svgap_ids",
        "svgap_support_samples", "svgap_source_id", "final_confidence", "notes",
    ]
    for species in species_list:
        with open(events_dir / ("%s.Msa.mergedSV.events.tsv" % species), "w", newline="") as fh:
            writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
            writer.writeheader()
            for ev in events_by_species[species]:
                writer.writerow({
                    "SpeciesSV_ID": ev["FinalSV_ID"],
                    "source_methods": source_method(ev),
                    "reference": "genome_Msa",
                    "species": species,
                    "CHROM": ev["CHROM"],
                    "POS": ev["POS"],
                    "END": ev["END"],
                    "SVTYPE": ev["SVTYPE"],
                    "SVLEN": ev["SVLEN"],
                    "read_SV_ID": ev["read_SV_ID"],
                    "svgap_support_count": ev["svgap_count"],
                    "svgap_ids": ",".join(ev["svgap_ids"]) if ev["svgap_ids"] else ".",
                    "svgap_support_samples": ",".join(sorted(ev["svgap_samples"])) if ev["svgap_samples"] else ".",
                    "svgap_source_id": ",".join(ev["svgap_source_ids"]) if ev["svgap_source_ids"] else ".",
                    "final_confidence": confidence(ev),
                    "notes": "per_species_read_mapping_consensus;flankQC_pass;svgap_evidence_only;no_svgap_only",
                })


def annotate_species_vcfs(consensus_paths, events_by_species, out_dir, bgzip, tabix):
    vcf_dir = Path(out_dir) / "vcf"
    vcf_dir.mkdir(parents=True, exist_ok=True)
    by_read = {}
    for species, events in events_by_species.items():
        for ev in events:
            by_read[(species, ev["read_SV_ID"])] = ev

    for species, in_vcf in consensus_paths.items():
        out_vcf = vcf_dir / ("%s.Msa.mergedSV.read_mapping_flankQC.svgap_evidence.vcf.gz" % species)
        tmp = Path(str(out_vcf) + ".inprogress")
        for path in [out_vcf, Path(str(out_vcf) + ".tbi"), tmp]:
            if path.exists() or path.is_symlink():
                path.unlink()
        with open(tmp, "wb") as out_bin:
            proc = subprocess.Popen([bgzip, "-c"], stdin=subprocess.PIPE, stdout=out_bin)
            assert proc.stdin is not None
            inserted = False
            for header, parts in iter_vcf_records(in_vcf):
                if header is not None:
                    if header.startswith("##"):
                        proc.stdin.write(header.encode())
                        continue
                    if header.startswith("#CHROM"):
                        if not inserted:
                            for info in INFO_DEFS:
                                proc.stdin.write((info + "\n").encode())
                            inserted = True
                        proc.stdin.write(header.encode())
                    continue
                if len(parts) < 8:
                    continue
                read_id = parts[2]
                ev = by_read[(species, read_id)]
                parts[2] = ev["FinalSV_ID"]
                extra = [
                    ("FINAL_SV_ID", ev["FinalSV_ID"]),
                    ("SPECIES", species),
                    ("READ_SV_ID", read_id),
                    ("SOURCE_METHODS", source_method(ev)),
                    ("SVGAP_SUPPORT_COUNT", ev["svgap_count"]),
                    ("SVGAP_SUPPORT_SAMPLES", ",".join(sorted(ev["svgap_samples"])) if ev["svgap_samples"] else "."),
                    ("FINAL_CONFIDENCE", confidence(ev)),
                ]
                add = ";".join("%s=%s" % (k, clean_info_value(v)) for k, v in extra)
                parts[7] = add if parts[7] in ("", ".") else parts[7] + ";" + add
                proc.stdin.write(("\t".join(map(str, parts)) + "\n").encode())
            proc.stdin.close()
            ret = proc.wait()
            if ret != 0:
                raise SystemExit("bgzip failed for %s" % species)
        tmp.rename(out_vcf)
        subprocess.check_call([tabix, "-f", "-p", "vcf", str(out_vcf)])


def write_stats(out_dir, species_list, events_by_species, svgap_seen, svgap_only, evidence_rows):
    summary_dir = Path(out_dir) / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)
    counts_path = summary_dir / "species_level_mergedSV.counts_by_species_type.tsv"
    with open(counts_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["species", "SVTYPE", "source_methods", "count"])
        for species in species_list:
            counter = Counter((ev["SVTYPE"], source_method(ev)) for ev in events_by_species[species])
            for (svtype, method), count in sorted(counter.items()):
                writer.writerow([species, svtype, method, count])

    overview_path = summary_dir / "species_level_mergedSV.overview.tsv"
    with open(overview_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "species", "total", "DEL", "INS", "DUP", "INV", "TRA",
            "read_based_only", "read_based_plus_SVGAP", "svgap_events_seen",
            "svgap_only_removed", "svgap_overlap_evidence_rows",
        ])
        for species in species_list:
            svtype_counts = Counter(ev["SVTYPE"] for ev in events_by_species[species])
            method_counts = Counter(source_method(ev) for ev in events_by_species[species])
            writer.writerow([
                species,
                len(events_by_species[species]),
                svtype_counts.get("DEL", 0),
                svtype_counts.get("INS", 0),
                svtype_counts.get("DUP", 0),
                svtype_counts.get("INV", 0),
                svtype_counts.get("TRA", 0),
                method_counts.get("read_based", 0),
                method_counts.get("read_based+SVGAP", 0),
                svgap_seen.get(species, 0),
                svgap_only.get(species, 0),
                evidence_rows.get(species, 0),
            ])

    method_path = Path(out_dir) / "README.method.tsv"
    with open(method_path, "w") as w:
        w.write("item\tvalue\n")
        w.write("reference\tgenome_Msa\n")
        w.write("species_count\t18\n")
        w.write("input_per_species_read_mapping_sv\tpbsv+sniffles2+cuteSV Jasmine consensus with min_support=2\n")
        w.write("flank_coverage_filter\t500bp left/right flanks; min_each_flank_depth>=5; mean_flank_depth>=8; mosdepth\n")
        w.write("genome_synteny_svgap_policy\tSVGAP is overlap evidence only; SVGAP-only SVs are counted as removed and not retained\n")
        w.write("svtypes\tDEL;INS;DUP;INV;TRA\n")
        w.write("output_unit\tone merged SV VCF and one event table per species\n")

    manifest_path = Path(out_dir) / "MANIFEST.tsv"
    with open(manifest_path, "w") as w:
        w.write("path\tdescription\n")
        w.write("vcf/*.vcf.gz\tPer-species merged SV VCFs: read-mapping consensus + flankQC, annotated with SVGAP evidence; no SVGAP-only records\n")
        w.write("events/*.events.tsv\tPer-species event tables with read_SV_ID and SVGAP support columns\n")
        w.write("evidence/*.svgap_overlap_evidence.tsv\tPer-species SVGAP overlap evidence rows\n")
        w.write("summary/species_level_mergedSV.overview.tsv\tMain table for per-species SV number and type statistics\n")
        w.write("summary/species_level_mergedSV.counts_by_species_type.tsv\tLong-format counts by species, SVTYPE and source_methods\n")
        w.write("README.method.tsv\tMethod and policy summary\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species-list", required=True)
    ap.add_argument("--consensus-dir", required=True)
    ap.add_argument("--svgap-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--bgzip", required=True)
    ap.add_argument("--tabix", required=True)
    args = ap.parse_args()

    with open(args.species_list) as fh:
        species_list = [line.strip() for line in fh if line.strip()]
    if species_list != FINAL_SAMPLES:
        missing = set(FINAL_SAMPLES).difference(species_list)
        extra = set(species_list).difference(FINAL_SAMPLES)
        if missing or extra:
            raise SystemExit("unexpected species list; missing=%s extra=%s" % (sorted(missing), sorted(extra)))

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    consensus_paths, events_by_species, index, _ = read_species_events(args.consensus_dir, species_list)
    svgap_seen, svgap_only, evidence_rows = add_svgap_support(
        events_by_species, index, out_dir, species_list, args.svgap_dir
    )
    write_event_tables(out_dir, species_list, events_by_species)
    annotate_species_vcfs(consensus_paths, events_by_species, out_dir, args.bgzip, args.tabix)
    write_stats(out_dir, species_list, events_by_species, svgap_seen, svgap_only, evidence_rows)

    total = sum(len(events_by_species[s]) for s in species_list)
    supported = sum(1 for s in species_list for ev in events_by_species[s] if ev["svgap_count"] > 0)
    print("PER_SPECIES_MERGED_SV_DIR=%s" % out_dir)
    print("SPECIES=%d" % len(species_list))
    print("TOTAL_READ_MAPPING_FLANKQC_EVENTS=%d" % total)
    print("EVENTS_WITH_SVGAP_EVIDENCE=%d" % supported)
    print("OVERVIEW=%s" % (out_dir / "summary" / "species_level_mergedSV.overview.tsv"))


if __name__ == "__main__":
    main()
