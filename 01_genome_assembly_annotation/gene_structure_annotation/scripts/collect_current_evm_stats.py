#!/usr/bin/env python3
import csv
import datetime as dt
import glob
import os
import statistics
import sys


SPECIES = {
    "genome_395": ("Mlu_395", "Medicago_lupulina"),
    "genome_410": ("Mpr_410", "Medicago_praecox"),
    "genome_436": ("Mra_436", "Medicago_radiata"),
    "genome_454": ("Mla_454", "Medicago_lanigera"),
    "genome_457": ("Mma_457", "Medicago_marina"),
    "genome_461": ("Mse_461", "Medicago_secundiflora"),
    "genome_468": ("Mcr_468", "Medicago_cretacea"),
    "genome_472": ("Msu_472", "Medicago_suffruticosa"),
    "genome_474a": ("Mca_474a", "Medicago_carstiensis_haplotype_a"),
    "genome_474b": ("Mca_474b", "Medicago_carstiensis_haplotype_b"),
    "genome_482": ("Med_482", "Medicago_edgeworthii"),
    "genome_A17": ("Mtr_A17", "Medicago_truncatula_A17"),
    "genome_M22": ("Mor_M22", "Medicago_orbicularis"),
    "genome_M46": ("Mfi_M46", "Medicago_fischeriana"),
    "genome_Mar": ("Mar_100", "Medicago_archiducis-nicolai"),
    "genome_Mpo": ("Mpo_200", "Medicago_polymorpha"),
    "genome_Mru": ("Mru_300", "Medicago_ruthenica"),
    "genome_Msa1": ("Msa_T2T_haplotype_1", "Medicago_sativa_haplotype_1"),
    "genome_Msa2": ("Msa_T2T_haplotype_2", "Medicago_sativa_haplotype_2"),
    "genome_R108": ("Mtr_R108", "Medicago_truncatula_R108"),
    "genome_ZM4": ("Msa_ZM4", "Medicago_sativa_cv_Zhongmu-4"),
}

VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*")


def parse_fasta(path):
    records = []
    current_id = None
    chunks = []
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, "".join(chunks).upper()))
                current_id = line[1:].split()[0]
                chunks = []
            else:
                if current_id is None:
                    raise ValueError(f"sequence before first FASTA header: {path}")
                chunks.append(line)
    if current_id is not None:
        records.append((current_id, "".join(chunks).upper()))
    return records


def parse_attrs(text):
    attrs = {}
    for item in text.rstrip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs


def parse_gff(path):
    counts = {"gene": 0, "mRNA": 0, "CDS": 0, "exon": 0}
    gene_ids = set()
    mrna_ids = set()
    cds_parent_ids = set()
    gene_span_bp = 0
    cds_bp = 0
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            if feature not in counts:
                continue
            start, end = int(fields[3]), int(fields[4])
            attrs = parse_attrs(fields[8])
            counts[feature] += 1
            if feature == "gene":
                gene_span_bp += end - start + 1
                if attrs.get("ID"):
                    gene_ids.add(attrs["ID"])
            elif feature == "mRNA":
                if attrs.get("ID"):
                    mrna_ids.add(attrs["ID"])
            elif feature == "CDS":
                cds_bp += end - start + 1
                cds_parent_ids.update(x for x in attrs.get("Parent", "").split(",") if x)
    return counts, gene_ids, mrna_ids, cds_parent_ids, gene_span_bp, cds_bp


def n50(lengths):
    target = sum(lengths) / 2
    cumulative = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= target:
            return length
    return 0


def percentile_nearest_rank(lengths, fraction):
    ordered = sorted(lengths)
    if not ordered:
        return 0
    index = max(0, min(len(ordered) - 1, int(len(ordered) * fraction + 0.999999) - 1))
    return ordered[index]


def collect_sample(sample, pep_path, gff_path):
    records = parse_fasta(pep_path)
    protein_ids = [ident for ident, _ in records]
    protein_id_set = set(protein_ids)
    duplicate_ids = len(protein_ids) - len(protein_id_set)
    lengths = []
    terminal_stop = internal_stop = invalid_seq = empty_seq = 0
    for _, sequence in records:
        if sequence.endswith("*"):
            terminal_stop += 1
            core = sequence[:-1]
        else:
            core = sequence
        if "*" in core:
            internal_stop += 1
        if any(residue not in VALID_AA for residue in sequence):
            invalid_seq += 1
        aa_len = sum(1 for residue in core if residue != "*")
        lengths.append(aa_len)
        if aa_len == 0:
            empty_seq += 1

    counts, gene_ids, mrna_ids, cds_parent_ids, gene_span_bp, cds_bp = parse_gff(gff_path)
    sample_id, latin_name = SPECIES.get(sample, (sample, "unknown"))
    protein_count = len(records)
    gene_count = counts["gene"]
    mrna_count = counts["mRNA"]
    return {
        "assembly": sample,
        "sample_id": sample_id,
        "latin_name": latin_name,
        "genes": gene_count,
        "mrnas": mrna_count,
        "cds_features": counts["CDS"],
        "exons": counts["exon"],
        "proteins": protein_count,
        "total_aa": sum(lengths),
        "mean_aa": f"{statistics.mean(lengths):.2f}" if lengths else "0.00",
        "median_aa": f"{statistics.median(lengths):.1f}" if lengths else "0.0",
        "p05_aa": percentile_nearest_rank(lengths, 0.05),
        "p95_aa": percentile_nearest_rank(lengths, 0.95),
        "min_aa": min(lengths, default=0),
        "max_aa": max(lengths, default=0),
        "n50_aa": n50(lengths),
        "lt50_aa": sum(length < 50 for length in lengths),
        "lt100_aa": sum(length < 100 for length in lengths),
        "ge1000_aa": sum(length >= 1000 for length in lengths),
        "internal_stop": internal_stop,
        "terminal_stop": terminal_stop,
        "invalid_residue_seq": invalid_seq,
        "empty_seq": empty_seq,
        "duplicate_protein_ids": duplicate_ids,
        "protein_ids_missing_from_gff_mrna": len(protein_id_set - mrna_ids),
        "gff_mrna_ids_missing_from_protein": len(mrna_ids - protein_id_set),
        "cds_parent_ids_missing_from_gff_mrna": len(cds_parent_ids - mrna_ids),
        "mean_gene_span_bp": f"{gene_span_bp / gene_count:.2f}" if gene_count else "0.00",
        "mean_cds_bp_per_mrna": f"{cds_bp / mrna_count:.2f}" if mrna_count else "0.00",
        "genes_equal_proteins": "yes" if gene_count == protein_count else "no",
        "mrnas_equal_proteins": "yes" if mrna_count == protein_count else "no",
        "done_marker": "yes",
    }


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: collect_current_evm_stats.py BASE")
    base = os.path.abspath(sys.argv[1])
    out_dir = os.path.join(base, "09_summary_tables")
    os.makedirs(out_dir, exist_ok=True)
    final_glob = os.path.join(base, "05_EVM_integration", "work", "*", "02_final", "*.evm.pep.fa")
    rows = []
    for pep_path in sorted(glob.glob(final_glob)):
        filename = os.path.basename(pep_path)
        sample = filename[:-len(".evm.pep.fa")]
        gff_path = os.path.join(os.path.dirname(pep_path), f"{sample}.evm.gff3")
        done_path = os.path.join(base, "05_EVM_integration", "status", f"{sample}.done")
        if not os.path.isfile(gff_path):
            raise FileNotFoundError(gff_path)
        row = collect_sample(sample, pep_path, gff_path)
        row["done_marker"] = "yes" if os.path.isfile(done_path) else "no"
        rows.append(row)
    if not rows:
        raise RuntimeError(f"no protein files found under {final_glob}")

    summary_path = os.path.join(out_dir, "current_evm_protein_annotation_summary.tsv")
    with open(summary_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    status_path = os.path.join(out_dir, "current_annotation_stage_status.tsv")
    with open(status_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["stage", "state", "interpretation"])
        writer.writerow(["05_EVM_integration", "21_done" if len(rows) == 21 else f"{len(rows)}_protein_sets", "current structural protein-coding gene models"])
        writer.writerow(["07_filter_TE_and_QC", "empty" if not os.listdir(os.path.join(base, "07_filter_TE_and_QC")) else "nonempty", "post-EVM filtering not represented in this summary"])
        writer.writerow(["08_function_annotation", "empty" if not os.listdir(os.path.join(base, "08_function_annotation")) else "nonempty", "no database functional annotation results available here"])

    readme_path = os.path.join(out_dir, "README_current_evm_protein_annotation_stats.txt")
    with open(readme_path, "w", encoding="utf-8") as handle:
        handle.write("Current EVM protein-coding gene annotation statistics\n")
        handle.write(f"Generated: {dt.datetime.now().astimezone().isoformat(timespec='seconds')}\n")
        handle.write(f"Inputs: {base}/05_EVM_integration/work/*/02_final/*.evm.pep.fa and matching GFF3 files\n")
        handle.write(f"Assemblies or haplotypes: {len(rows)}\n")
        handle.write("Protein lengths exclude a terminal stop symbol. internal_stop counts proteins that still contain a stop after removal of a terminal stop.\n")
        handle.write("Important: 07_filter_TE_and_QC and 08_function_annotation are currently empty. This table summarizes current EVM structural protein sets, not post-TE/QC final sets or GO/KEGG/InterPro/eggNOG functional-annotation coverage.\n")
        handle.write("Mca 474a/474b and Msa1/Msa2 are paired haplotypes; A17/R108 are two M. truncatula accessions.\n")
    print(summary_path)
    print(status_path)
    print(readme_path)


if __name__ == "__main__":
    main()
