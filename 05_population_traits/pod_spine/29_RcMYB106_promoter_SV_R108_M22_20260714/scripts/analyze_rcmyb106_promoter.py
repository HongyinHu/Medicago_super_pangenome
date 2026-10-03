#!/usr/bin/env python3
import argparse
import gzip
import html
import re
import shutil
import subprocess
from pathlib import Path


AT_MYB106 = (
    "MGRSPCCDKAGLKKGPWTPEEDQKLLAYIEEHGHGSWRSLPEKAGLQRCGKSCRLRWTNYLRPDIKRGKF"
    "TVQEEQTIIQLHALLGNRWSAIATHLPKRTDNEIKNYWNTHLKKRLIKMGIDPVTHKHKNETLSSSTGQS"
    "KNAATLSHMAQWESARLEAEARLARESKLLHLQHYQNNNNLNKSAAPQQHCFTQKTSTNWTKPNQGNGDQ"
    "QLESPTSTVTFSENLLMPLGIPTDSSRNRNNNNNESSAMIELAVSSSTSSDVSLVKEHEHDWIRQINCGS"
    "GGIGEGFTSLLIGDSVGRGLPTGKNEATAGVGNESEYNYYEDNKNYWNSILNLVDSSPSDSATMF"
)

AT_MYB16_ACCESSION = "NP_197035.1"
CASTOR_BASE = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/578/655/"
    "GCF_019578655.1_ASM1957865v1/"
)


def run(cmd, log_path=None, stdout_path=None):
    if log_path:
        with open(log_path, "a", encoding="utf-8") as log:
            log.write("CMD: " + " ".join(map(str, cmd)) + "\n")
    if stdout_path:
        with open(stdout_path, "w", encoding="utf-8") as out:
            result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, universal_newlines=True)
    else:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if log_path:
        with open(log_path, "a", encoding="utf-8") as log:
            if result.stdout and not stdout_path:
                log.write(result.stdout)
            if result.stderr:
                log.write(result.stderr)
    if result.returncode:
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(map(str, cmd))}")
    return result.stdout


def download(url, output, log_path):
    output = Path(output)
    if output.exists() and output.stat().st_size > 0:
        return
    run(["curl", "-fL", "--retry", "3", "--connect-timeout", "30", "-o", str(output), url], log_path)


def gunzip_file(source, output):
    output = Path(output)
    if output.exists() and output.stat().st_size > 0:
        return
    with gzip.open(source, "rb") as src, open(output, "wb") as dst:
        shutil.copyfileobj(src, dst)


def read_fasta(path):
    records = {}
    header = None
    chunks = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records[header.split()[0]] = (header, "".join(chunks))
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        records[header.split()[0]] = (header, "".join(chunks))
    return records


def write_fasta(path, records):
    with open(path, "w", encoding="utf-8") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                handle.write(seq[i:i + 60] + "\n")


def parse_attrs(text):
    attrs = {}
    for item in text.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs


def parse_gff(path):
    genes = {}
    cds_by_protein = {}
    chromosome_names = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9:
                continue
            seqid, source, feature, start, end, score, strand, phase, attr_text = fields
            attrs = parse_attrs(attr_text)
            start = int(start)
            end = int(end)
            if feature == "region" and "chromosome" in attrs:
                chromosome_names[attrs["chromosome"]] = seqid
            elif feature == "gene" and "ID" in attrs:
                genes[attrs["ID"]] = {
                    "gene_id": attrs["ID"], "seqid": seqid, "start": start,
                    "end": end, "strand": strand, "attrs": attrs,
                }
            elif feature == "CDS" and "protein_id" in attrs:
                protein_id = attrs["protein_id"].split(".")[0] + "." + attrs["protein_id"].split(".")[-1]
                cds_by_protein.setdefault(protein_id, []).append({
                    "seqid": seqid, "start": start, "end": end, "strand": strand,
                    "attrs": attrs,
                })
    return genes, cds_by_protein, chromosome_names


def parse_simple_gff_gene(path, gene_id):
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = parse_attrs(fields[8])
            if attrs.get("ID") == gene_id:
                return {
                    "gene_id": gene_id, "seqid": fields[0], "start": int(fields[3]),
                    "end": int(fields[4]), "strand": fields[6],
                }
    raise KeyError(f"Gene {gene_id} not found in {path}")


def blast(blastp, query, subject, output, log_path, max_hits=50):
    fields = "qseqid sseqid pident length qlen slen evalue bitscore qcovs"
    run([
        blastp, "-query", str(query), "-subject", str(subject), "-evalue", "1e-10",
        "-max_target_seqs", str(max_hits), "-outfmt", f"6 {fields}", "-out", str(output),
    ], log_path)


def read_blast(path):
    rows = []
    names = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore", "qcovs"]
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            values = line.rstrip().split("\t")
            row = dict(zip(names, values))
            for key in ("pident", "evalue", "bitscore", "qcovs"):
                row[key] = float(row[key])
            for key in ("length", "qlen", "slen"):
                row[key] = int(row[key])
            rows.append(row)
    return rows


def gene_from_header(header, seq_id):
    match = re.search(r"(?:^|\s)gene=([^\s]+)", header)
    if match:
        return match.group(1)
    return re.sub(r"\.\d+$", "", seq_id)


def fasta_index(path):
    index = {}
    with open(str(path) + ".fai", encoding="utf-8") as handle:
        for line in handle:
            name, length, offset, line_bases, line_width = line.rstrip().split("\t")[:5]
            index[name] = tuple(map(int, (length, offset, line_bases, line_width)))
    return index


def fetch_fasta_region(path, index, seqid, start, end):
    length, offset, line_bases, line_width = index[seqid]
    start = max(1, start)
    end = min(length, end)
    if end < start:
        return ""
    start0 = start - 1
    end0 = end
    byte_start = offset + (start0 // line_bases) * line_width + (start0 % line_bases)
    byte_end = offset + (end0 // line_bases) * line_width + (end0 % line_bases)
    with open(path, "rb") as handle:
        handle.seek(byte_start)
        data = handle.read(byte_end - byte_start + line_width)
    seq = re.sub(b"[\r\n]", b"", data).decode("ascii")
    return seq[:end0 - start0].upper()


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def extract_windows(genome, gene, output_dir, prefix):
    index = fasta_index(genome)
    chrom_len = index[gene["seqid"]][0]
    start, end, strand = gene["start"], gene["end"], gene["strand"]
    if strand == "+":
        promoter10_coords = (max(1, start - 10000), start - 1)
        promoter2_coords = (max(1, start - 2000), start - 1)
        locus_coords = (max(1, start - 10000), min(chrom_len, end + 2000))
    else:
        promoter10_coords = (end + 1, min(chrom_len, end + 10000))
        promoter2_coords = (end + 1, min(chrom_len, end + 2000))
        locus_coords = (max(1, start - 2000), min(chrom_len, end + 10000))

    outputs = {}
    for label, coords in (
        ("promoter_10kb", promoter10_coords),
        ("promoter_2kb", promoter2_coords),
        ("locus_up10kb_down2kb", locus_coords),
        ("gene", (start, end)),
    ):
        seq = fetch_fasta_region(genome, index, gene["seqid"], coords[0], coords[1])
        if strand == "-":
            seq = revcomp(seq)
        path = output_dir / f"{prefix}.{label}.fa"
        write_fasta(path, [(f"{prefix}|{gene['gene_id']}|{gene['seqid']}:{coords[0]}-{coords[1]}|{strand}", seq)])
        outputs[label] = {"path": path, "coords": coords, "length": len(seq)}
    return outputs


def parse_paf(path, target_prefix, query_prefix, promoter_len=None):
    alignments = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 12:
                continue
            tags = {}
            for tag in fields[12:]:
                parts = tag.split(":", 2)
                if len(parts) == 3:
                    tags[parts[0]] = parts[2]
            alignments.append({
                "qname": fields[0], "qlen": int(fields[1]), "qstart": int(fields[2]),
                "qend": int(fields[3]), "strand": fields[4], "tname": fields[5],
                "tlen": int(fields[6]), "tstart": int(fields[7]), "tend": int(fields[8]),
                "nmatch": int(fields[9]), "alen": int(fields[10]), "mapq": int(fields[11]),
                "cg": tags.get("cg", ""), "tp": tags.get("tp", ""),
            })

    events = []
    for aln in alignments:
        if aln["strand"] != "+" or not aln["cg"]:
            continue
        qpos, tpos = aln["qstart"], aln["tstart"]
        for length_text, op in re.findall(r"(\d+)([MIDNSHP=X])", aln["cg"]):
            length = int(length_text)
            if op in ("M", "=", "X"):
                qpos += length
                tpos += length
            elif op == "I":
                if length >= 50:
                    events.append({
                        "source": "CIGAR", "event": f"{query_prefix}_insertion_relative_to_{target_prefix}",
                        "length": length, "target_start0": tpos, "target_end0": tpos,
                        "query_start0": qpos, "query_end0": qpos + length,
                    })
                qpos += length
            elif op in ("D", "N"):
                if length >= 50:
                    events.append({
                        "source": "CIGAR", "event": f"{query_prefix}_deletion_relative_to_{target_prefix}",
                        "length": length, "target_start0": tpos, "target_end0": tpos + length,
                        "query_start0": qpos, "query_end0": qpos,
                    })
                tpos += length

    plus = sorted([a for a in alignments if a["strand"] == "+"], key=lambda x: (x["tstart"], x["qstart"]))
    for left, right in zip(plus, plus[1:]):
        target_gap = right["tstart"] - left["tend"]
        query_gap = right["qstart"] - left["qend"]
        if target_gap < 0 or query_gap < 0:
            continue
        delta = target_gap - query_gap
        if abs(delta) < 50:
            continue
        event = (
            f"{query_prefix}_deletion_relative_to_{target_prefix}"
            if delta > 0 else f"{query_prefix}_insertion_relative_to_{target_prefix}"
        )
        events.append({
            "source": "split_alignment_gap", "event": event, "length": abs(delta),
            "target_start0": left["tend"], "target_end0": right["tstart"],
            "query_start0": left["qend"], "query_end0": right["qstart"],
        })

    for event in events:
        if promoter_len is not None:
            event["target_upstream_start"] = event["target_start0"] - promoter_len
            event["target_upstream_end"] = event["target_end0"] - promoter_len
        else:
            event["target_upstream_start"] = "NA"
            event["target_upstream_end"] = "NA"
    return alignments, events


def write_tsv(path, rows, columns):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(column, "NA")) for column in columns) + "\n")


def make_svg(path, r108_gene, m22_gene, events):
    width, height = 1200, 310
    x0, x1 = 170, 1130
    scale = (x1 - x0) / 10000.0
    large_deletions = [e for e in events if "deletion" in e["event"] and e["length"] >= 1000]
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#222}.label{font-size:18px;font-weight:600}.tick{font-size:13px}.note{font-size:14px}</style>',
        '<text x="30" y="38" class="label">RcMYB106 ortholog promoter comparison</text>',
    ]
    for y, label, color in ((105, f"R108  {r108_gene}", "#14866d"), (205, f"M22  {m22_gene}", "#565b66")):
        parts.append(f'<text x="30" y="{y + 6}" class="label">{html.escape(label)}</text>')
        parts.append(f'<line x1="{x0}" y1="{y}" x2="{x1}" y2="{y}" stroke="{color}" stroke-width="10"/>')
        parts.append(f'<polygon points="{x1},{y - 10} {x1 + 20},{y} {x1},{y + 10}" fill="{color}"/>')
    for value in (-10000, -8000, -6000, -4000, -2000, 0):
        x = x0 + (value + 10000) * scale
        parts.append(f'<line x1="{x:.1f}" y1="75" x2="{x:.1f}" y2="230" stroke="#d4d7dc" stroke-width="1"/>')
        parts.append(f'<text x="{x:.1f}" y="258" text-anchor="middle" class="tick">{value}</text>')
    for event in large_deletions:
        start = int(event["target_start0"])
        end = int(event["target_end0"])
        x = x0 + start * scale
        w = max(2, (end - start) * scale)
        parts.append(f'<rect x="{x:.1f}" y="184" width="{w:.1f}" height="42" fill="#d64b3c" opacity="0.88"/>')
        parts.append(f'<line x1="{x:.1f}" y1="105" x2="{x:.1f}" y2="184" stroke="#d64b3c" stroke-dasharray="5,4"/>')
        parts.append(f'<line x1="{x + w:.1f}" y1="105" x2="{x + w:.1f}" y2="184" stroke="#d64b3c" stroke-dasharray="5,4"/>')
        parts.append(f'<text x="{x + w / 2:.1f}" y="174" text-anchor="middle" class="note">{event["length"]} bp</text>')
    if not large_deletions:
        parts.append('<text x="650" y="285" text-anchor="middle" class="note">No M22 deletion >=1 kb detected in the aligned 10 kb promoter.</text>')
    else:
        parts.append('<text x="650" y="290" text-anchor="middle" class="note">Red block: sequence absent from M22 relative to R108.</text>')
    parts.append('</svg>')
    Path(path).write_text("\n".join(parts) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--r108-genome", required=True)
    parser.add_argument("--r108-gff", required=True)
    parser.add_argument("--r108-pep", required=True)
    parser.add_argument("--m22-genome", required=True)
    parser.add_argument("--m22-gff", required=True)
    parser.add_argument("--m22-pep", required=True)
    parser.add_argument("--blastp", required=True)
    parser.add_argument("--minimap2", required=True)
    args = parser.parse_args()

    out = Path(args.outdir)
    inputs = out / "00_inputs"
    ortholog = out / "01_ortholog"
    promoter = out / "02_promoter"
    alignment = out / "03_alignment"
    summary = out / "04_summary"
    logs = out / "logs"
    for directory in (inputs, ortholog, promoter, alignment, summary, logs):
        directory.mkdir(parents=True, exist_ok=True)
    log_path = logs / "pipeline.log"

    query_fa = inputs / "Arabidopsis_MYB106_MYB16.pep.fa"
    if not query_fa.exists():
        at_myb16 = run([
            "curl", "-sL", "--get", "--data-urlencode", "db=protein",
            "--data-urlencode", f"id={AT_MYB16_ACCESSION}", "--data-urlencode", "rettype=fasta",
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        ], log_path)
        myb16_seq = "".join(line.strip() for line in at_myb16.splitlines() if not line.startswith(">"))
        write_fasta(query_fa, [("AtMYB106|NP_001326423.1", AT_MYB106), ("AtMYB16|NP_197035.1", myb16_seq)])

    castor_protein_gz = inputs / "GCF_019578655.1_ASM1957865v1_protein.faa.gz"
    castor_gff_gz = inputs / "GCF_019578655.1_ASM1957865v1_genomic.gff.gz"
    download(CASTOR_BASE + castor_protein_gz.name, castor_protein_gz, log_path)
    download(CASTOR_BASE + castor_gff_gz.name, castor_gff_gz, log_path)
    castor_protein = inputs / "GCF_019578655.1_ASM1957865v1_protein.faa"
    castor_gff = inputs / "GCF_019578655.1_ASM1957865v1_genomic.gff"
    gunzip_file(castor_protein_gz, castor_protein)
    gunzip_file(castor_gff_gz, castor_gff)

    at_vs_castor = ortholog / "Arabidopsis_MYB106_MYB16_vs_castor.blastp.tsv"
    blast(args.blastp, query_fa, castor_protein, at_vs_castor, log_path, 100)
    blast_rows = read_blast(at_vs_castor)
    _, castor_cds, chromosome_names = parse_gff(castor_gff)
    chrom10 = chromosome_names.get("10")
    if not chrom10:
        raise RuntimeError("Chromosome 10 was not found in the castor NCBI GFF")
    region_proteins = {
        protein_id for protein_id, entries in castor_cds.items()
        if any(
            entry["seqid"] == chrom10 and "MYB106" in entry["attrs"].get("product", "")
            for entry in entries
        )
    }
    candidates = [
        row for row in blast_rows
        if row["qseqid"].startswith("AtMYB106") and row["sseqid"] in region_proteins
    ]
    if not candidates:
        raise RuntimeError("No AtMYB106-like protein annotated as MYB106 was found on castor chromosome 10")
    rc_hit = max(candidates, key=lambda row: row["bitscore"])
    rc_protein_id = rc_hit["sseqid"]
    castor_records = read_fasta(castor_protein)
    rc_header, rc_seq = castor_records[rc_protein_id]
    rc_candidate_fa = ortholog / "RcMYB106_WT05_candidate.pep.fa"
    write_fasta(rc_candidate_fa, [(f"RcMYB106_WT05_candidate|{rc_protein_id}", rc_seq)])

    rc_region_rows = []
    for protein_id in sorted(region_proteins):
        relevant = [row for row in blast_rows if row["sseqid"] == protein_id]
        if not relevant:
            continue
        best = max(relevant, key=lambda row: row["bitscore"])
        entry = castor_cds[protein_id][0]
        rc_region_rows.append({
            "protein_id": protein_id, "seqid": entry["seqid"],
            "start": min(x["start"] for x in castor_cds[protein_id]),
            "end": max(x["end"] for x in castor_cds[protein_id]),
            "strand": entry["strand"], "query": best["qseqid"],
            "pident": best["pident"], "qcovs": best["qcovs"],
            "evalue": best["evalue"], "bitscore": best["bitscore"],
            "selected": protein_id == rc_protein_id,
        })
    write_tsv(ortholog / "castor_chr10_MYB106_candidates.tsv", rc_region_rows,
              ["protein_id", "seqid", "start", "end", "strand", "query", "pident", "qcovs", "evalue", "bitscore", "selected"])

    species = {
        "R108": {"pep": Path(args.r108_pep), "gff": Path(args.r108_gff), "genome": Path(args.r108_genome)},
        "M22": {"pep": Path(args.m22_pep), "gff": Path(args.m22_gff), "genome": Path(args.m22_genome)},
    }
    selected = {}
    ortholog_rows = []
    for name, data in species.items():
        forward_path = ortholog / f"RcMYB106_vs_{name}.blastp.tsv"
        blast(args.blastp, rc_candidate_fa, data["pep"], forward_path, log_path, 30)
        forward = read_blast(forward_path)
        if not forward:
            raise RuntimeError(f"No RcMYB106 hit found in {name}")
        top = forward[0]
        records = read_fasta(data["pep"])
        header, seq = records[top["sseqid"]]
        gene_id = gene_from_header(header, top["sseqid"])
        candidate_fa = ortholog / f"{name}_RcMYB106_candidate.pep.fa"
        write_fasta(candidate_fa, [(f"{name}|{top['sseqid']}|gene={gene_id}", seq)])
        reciprocal_path = ortholog / f"{name}_candidate_vs_castor.blastp.tsv"
        blast(args.blastp, candidate_fa, castor_protein, reciprocal_path, log_path, 20)
        reciprocal = read_blast(reciprocal_path)
        reciprocal_top = reciprocal[0]["sseqid"] if reciprocal else "NA"
        gene = parse_simple_gff_gene(data["gff"], gene_id)
        selected[name] = {
            "protein_id": top["sseqid"], "gene_id": gene_id, "gene": gene,
            "candidate_fa": candidate_fa, "seq": seq, "forward": top,
            "reciprocal_top": reciprocal_top,
        }
        ortholog_rows.append({
            "species": name, "gene_id": gene_id, "protein_id": top["sseqid"],
            "chrom": gene["seqid"], "start": gene["start"], "end": gene["end"], "strand": gene["strand"],
            "rc_pident": top["pident"], "rc_qcovs": top["qcovs"], "rc_evalue": top["evalue"],
            "rc_bitscore": top["bitscore"], "reciprocal_castor_top": reciprocal_top,
            "rc_reciprocal_best_hit": reciprocal_top == rc_protein_id,
        })

    for query_name, subject_name in (("R108", "M22"), ("M22", "R108")):
        cross_path = ortholog / f"{query_name}_candidate_vs_{subject_name}.blastp.tsv"
        blast(args.blastp, selected[query_name]["candidate_fa"], species[subject_name]["pep"], cross_path, log_path, 20)
        cross = read_blast(cross_path)
        selected[query_name]["cross_top"] = cross[0]["sseqid"] if cross else "NA"
    for row in ortholog_rows:
        other = "M22" if row["species"] == "R108" else "R108"
        row["cross_species_top"] = selected[row["species"]]["cross_top"]
        row["cross_reciprocal_best_hit"] = selected[row["species"]]["cross_top"] == selected[other]["protein_id"]
    write_tsv(ortholog / "RcMYB106_R108_M22_ortholog_summary.tsv", ortholog_rows,
              ["species", "gene_id", "protein_id", "chrom", "start", "end", "strand", "rc_pident", "rc_qcovs",
               "rc_evalue", "rc_bitscore", "reciprocal_castor_top", "rc_reciprocal_best_hit", "cross_species_top",
               "cross_reciprocal_best_hit"])
    write_fasta(ortholog / "RcMYB106_R108_M22_candidates.pep.fa", [
        (f"RcMYB106_WT05|{rc_protein_id}", rc_seq),
        (f"R108|{selected['R108']['protein_id']}|gene={selected['R108']['gene_id']}", selected["R108"]["seq"]),
        (f"M22|{selected['M22']['protein_id']}|gene={selected['M22']['gene_id']}", selected["M22"]["seq"]),
    ])

    windows = {}
    window_rows = []
    for name, data in species.items():
        windows[name] = extract_windows(data["genome"], selected[name]["gene"], promoter, name)
        for label, info in windows[name].items():
            window_rows.append({
                "species": name, "gene_id": selected[name]["gene_id"], "feature": label,
                "chrom": selected[name]["gene"]["seqid"], "genomic_start": info["coords"][0],
                "genomic_end": info["coords"][1], "strand": selected[name]["gene"]["strand"],
                "sequence_length": info["length"], "fasta": info["path"],
            })
    write_tsv(promoter / "RcMYB106_promoter_window_coordinates.tsv", window_rows,
              ["species", "gene_id", "feature", "chrom", "genomic_start", "genomic_end", "strand", "sequence_length", "fasta"])

    all_events = []
    alignment_rows = []
    for label in ("promoter_2kb", "promoter_10kb", "locus_up10kb_down2kb", "gene"):
        paf = alignment / f"R108_vs_M22.{label}.paf"
        run([
            args.minimap2, "-x", "asm20", "-c", "--cs=long", "-N", "20",
            str(windows["R108"][label]["path"]), str(windows["M22"][label]["path"]),
        ], log_path, paf)
        promoter_len = windows["R108"][label]["length"] if label.startswith("promoter") else None
        alignments, events = parse_paf(paf, "R108", "M22", promoter_len)
        for event in events:
            event["comparison"] = label
        all_events.extend(events)
        for aln in alignments:
            alignment_rows.append({
                "comparison": label, "query_length": aln["qlen"], "query_start0": aln["qstart"],
                "query_end0": aln["qend"], "strand": aln["strand"], "target_length": aln["tlen"],
                "target_start0": aln["tstart"], "target_end0": aln["tend"], "matches": aln["nmatch"],
                "alignment_block_length": aln["alen"], "identity_pct": round(100 * aln["nmatch"] / max(1, aln["alen"]), 3),
                "mapq": aln["mapq"], "cigar": aln["cg"],
            })
    write_tsv(alignment / "R108_vs_M22.alignment_summary.tsv", alignment_rows,
              ["comparison", "query_length", "query_start0", "query_end0", "strand", "target_length", "target_start0",
               "target_end0", "matches", "alignment_block_length", "identity_pct", "mapq", "cigar"])
    write_tsv(alignment / "R108_vs_M22.indels_ge50bp.tsv", all_events,
              ["comparison", "source", "event", "length", "target_start0", "target_end0", "query_start0", "query_end0",
               "target_upstream_start", "target_upstream_end"])

    promoter10_events = [e for e in all_events if e["comparison"] == "promoter_10kb"]
    large_deletions = [e for e in promoter10_events if e["event"] == "M22_deletion_relative_to_R108" and e["length"] >= 1000]
    analogous = [
        e for e in large_deletions
        if e["target_upstream_start"] < -254 and e["target_upstream_end"] > -4607
    ]
    make_svg(summary / "RcMYB106_R108_M22_promoter_comparison.svg",
             selected["R108"]["gene_id"], selected["M22"]["gene_id"], promoter10_events)

    reciprocal_ok = all(row["rc_reciprocal_best_hit"] and row["cross_reciprocal_best_hit"] for row in ortholog_rows)
    if large_deletions:
        conclusion = (
            f"Detected {len(large_deletions)} M22 deletion event(s) >=1 kb relative to R108 "
            "within the aligned 10 kb upstream promoter."
        )
    else:
        conclusion = "No M22 deletion >=1 kb relative to R108 was detected within the aligned 10 kb upstream promoter."
    if analogous:
        analogy = "At least one deletion overlaps the castor LYY9 analogous interval (-4607 to -254 bp)."
    else:
        analogy = "No detected >=1 kb deletion overlaps the castor LYY9 analogous interval (-4607 to -254 bp)."

    md = [
        "# RcMYB106 ortholog promoter comparison: R108 vs M22",
        "",
        "## Ortholog identification",
        "",
        f"- Castor WT05 candidate: `{rc_protein_id}` selected from chromosome 10 MYB106 annotations using AtMYB106 similarity.",
        f"- R108 candidate: `{selected['R108']['gene_id']}` / `{selected['R108']['protein_id']}`.",
        f"- M22 candidate: `{selected['M22']['gene_id']}` / `{selected['M22']['protein_id']}`.",
        f"- Reciprocal best-hit support across castor, R108 and M22: `{'PASS' if reciprocal_ok else 'CHECK'}`.",
        "",
        "## Promoter structural comparison",
        "",
        f"- {conclusion}",
        f"- {analogy}",
        "- Promoter coordinates are oriented by transcription direction; upstream positions are reported relative to TSS (0).",
        "- A negative result means no assembly-level deletion >=1 kb was found in this pair; it does not exclude smaller indels or cis-regulatory SNPs.",
        "",
        "## Main outputs",
        "",
        "- `01_ortholog/RcMYB106_R108_M22_ortholog_summary.tsv`",
        "- `02_promoter/RcMYB106_promoter_window_coordinates.tsv`",
        "- `03_alignment/R108_vs_M22.indels_ge50bp.tsv`",
        "- `03_alignment/R108_vs_M22.alignment_summary.tsv`",
        "- `04_summary/RcMYB106_R108_M22_promoter_comparison.svg`",
    ]
    (summary / "README.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md))


if __name__ == "__main__":
    main()
