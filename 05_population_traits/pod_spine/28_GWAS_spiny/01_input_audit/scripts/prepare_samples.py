#!/usr/bin/env python3
import csv
import glob
import gzip
import os
import zipfile
import xml.etree.ElementTree as ET


GWAS = "path/to/project/N_4.pod_spiny/28_GWAS_spiny"
PHENO_XLSX = os.path.join(GWAS, "00_data/144_phenotype/重测序保留个体-144.xlsx")
SV_RUN = (
    "path/to/project/N_4.pod_spiny/"
    "09_short_read_SV_Msa_20260706/medicago144_3caller_20260709"
)
SV_MANIFEST = os.path.join(SV_RUN, "input/medicago144_sample_bams.tsv")
SV_GLOB = os.path.join(SV_RUN, "results/consensus_per_sample/*/*.vcf.gz")
SNP_VCF = (
    "path/to/project/38.medicago_resequence/"
    "2.call_SNP_new/07_filter_snp/DP_6-85_miss_0.2.all_samples.biallelic.SNP.vcf.gz"
)
INDEL_VCF = (
    "path/to/project/38.medicago_resequence/"
    "2.call_SNP_new/08_filter_indel/DP_6-85_miss_0.2.all_samples.biallelic.INDEL.vcf.gz"
)
OUTDIR = os.path.join(GWAS, "01_input_audit/summary")


def normalize_id(value):
    if value is None:
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip()


def vcf_samples(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    raise RuntimeError(f"No #CHROM header in {path}")


def load_xlsx_sheet(path, sheet_name):
    main_ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    rel_ns = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    pkg_rel_ns = "http://schemas.openxmlformats.org/package/2006/relationships"
    with zipfile.ZipFile(path) as archive:
        shared_strings = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in root.findall(f"{{{main_ns}}}si"):
                shared_strings.append(
                    "".join(node.text or "" for node in item.iter(f"{{{main_ns}}}t"))
                )

        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relationship_id = None
        for sheet in workbook.findall(f".//{{{main_ns}}}sheet"):
            if sheet.attrib.get("name") == sheet_name:
                relationship_id = sheet.attrib.get(f"{{{rel_ns}}}id")
                break
        if relationship_id is None:
            raise RuntimeError(f"Worksheet not found: {sheet_name}")

        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        target = None
        for relationship in relationships.findall(f"{{{pkg_rel_ns}}}Relationship"):
            if relationship.attrib.get("Id") == relationship_id:
                target = relationship.attrib["Target"]
                break
        if target is None:
            raise RuntimeError(f"Worksheet relationship not found: {relationship_id}")
        worksheet_path = "xl/" + target.lstrip("/")
        worksheet_path = os.path.normpath(worksheet_path).replace("\\", "/")
        root = ET.fromstring(archive.read(worksheet_path))

        parsed_rows = []
        for row in root.findall(f".//{{{main_ns}}}row"):
            cells = {}
            max_column = -1
            for cell in row.findall(f"{{{main_ns}}}c"):
                reference = cell.attrib.get("r", "")
                letters = "".join(char for char in reference if char.isalpha())
                column = 0
                for char in letters.upper():
                    column = column * 26 + ord(char) - ord("A") + 1
                column -= 1
                max_column = max(max_column, column)
                cell_type = cell.attrib.get("t")
                value_node = cell.find(f"{{{main_ns}}}v")
                if cell_type == "inlineStr":
                    inline = cell.find(f"{{{main_ns}}}is")
                    value = "" if inline is None else "".join(
                        node.text or "" for node in inline.iter(f"{{{main_ns}}}t")
                    )
                elif value_node is None:
                    value = None
                elif cell_type == "s":
                    value = shared_strings[int(value_node.text)]
                else:
                    raw = value_node.text
                    try:
                        number = float(raw)
                        value = int(number) if number.is_integer() else number
                    except (TypeError, ValueError):
                        value = raw
                cells[column] = value
            parsed_rows.append([cells.get(i) for i in range(max_column + 1)])
    return parsed_rows


os.makedirs(OUTDIR, exist_ok=True)

rows = iter(load_xlsx_sheet(PHENO_XLSX, "resequence_sample_pod_traits_me"))
headers = [str(value).strip() if value is not None else "" for value in next(rows)]
required = [
    "Sample_ID",
    "Latin_name",
    "Section",
    "Pod_spine",
    "Pod_spine_binary_spiny1_spineless0",
    "PLINK2_case_control_spiny2_spineless1_missing-9",
]
missing_headers = [name for name in required if name not in headers]
if missing_headers:
    raise RuntimeError("Missing phenotype columns: " + ",".join(missing_headers))

phenotypes = []
for values in rows:
    record = dict(zip(headers, values))
    sample_id = normalize_id(record["Sample_ID"])
    if not sample_id:
        continue
    record["Sample_ID"] = sample_id
    record["Pod_spine_binary_spiny1_spineless0"] = int(
        record["Pod_spine_binary_spiny1_spineless0"]
    )
    record["PLINK2_case_control_spiny2_spineless1_missing-9"] = int(
        record["PLINK2_case_control_spiny2_spineless1_missing-9"]
    )
    phenotypes.append(record)

sample_order = [row["Sample_ID"] for row in phenotypes]
if len(sample_order) != 144 or len(set(sample_order)) != 144:
    raise RuntimeError(
        f"Expected 144 unique phenotype samples, observed {len(sample_order)} rows and "
        f"{len(set(sample_order))} unique IDs"
    )

manifest = {}
with open(SV_MANIFEST, newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        manifest[normalize_id(row["sample_id"])] = row

snp_samples = vcf_samples(SNP_VCF)
indel_samples = vcf_samples(INDEL_VCF)
sv_samples = {
    os.path.basename(os.path.dirname(path)) for path in glob.glob(SV_GLOB)
}

sets = {
    "phenotype": set(sample_order),
    "SNP": set(snp_samples),
    "INDEL": set(indel_samples),
    "SV": sv_samples,
    "BAM_manifest": set(manifest),
}
for source, values in sets.items():
    missing = sorted(set(sample_order) - values)
    if source != "phenotype" and missing:
        raise RuntimeError(f"{source} is missing phenotype samples: {','.join(missing)}")

with open(os.path.join(OUTDIR, "sample_order.txt"), "w", encoding="utf-8") as handle:
    handle.write("\n".join(sample_order) + "\n")

with open(os.path.join(OUTDIR, "bam.list"), "w", encoding="utf-8") as handle:
    for sample_id in sample_order:
        bam = manifest[sample_id]["bam"]
        bai = manifest[sample_id]["bai"]
        if not os.path.isfile(bam) or not os.path.isfile(bai):
            raise RuntimeError(f"Missing BAM/BAI for {sample_id}: {bam} {bai}")
        handle.write(bam + "\n")

metadata_headers = [
    "Sample_ID",
    "Latin_name",
    "Section",
    "Pod_spine",
    "Pod_spine_binary_spiny1_spineless0",
    "PLINK2_case_control_spiny2_spineless1_missing-9",
    "bam",
    "bai",
]
with open(
    os.path.join(OUTDIR, "phenotype_144.tsv"), "w", newline="", encoding="utf-8"
) as handle:
    writer = csv.DictWriter(handle, fieldnames=metadata_headers, delimiter="\t")
    writer.writeheader()
    for row in phenotypes:
        output = {key: row.get(key, "") for key in metadata_headers}
        output["bam"] = manifest[row["Sample_ID"]]["bam"]
        output["bai"] = manifest[row["Sample_ID"]]["bai"]
        writer.writerow(output)

with open(
    os.path.join(OUTDIR, "pod_spine.plink.pheno"), "w", newline="", encoding="utf-8"
) as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["#FID", "IID", "POD_SPINE"])
    for row in phenotypes:
        sample_id = row["Sample_ID"]
        writer.writerow(
            [
                sample_id,
                sample_id,
                row["PLINK2_case_control_spiny2_spineless1_missing-9"],
            ]
        )

species_counts = {}
for row in phenotypes:
    key = (str(row.get("Latin_name", "")), str(row.get("Section", "")))
    counts = species_counts.setdefault(key, [0, 0, 0])
    phenotype = row["Pod_spine_binary_spiny1_spineless0"]
    counts[0] += 1
    counts[1] += phenotype == 0
    counts[2] += phenotype == 1

with open(
    os.path.join(OUTDIR, "phenotype_counts_by_species.tsv"),
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["Latin_name", "Section", "N", "spineless_0", "spiny_1"])
    for (latin_name, section), counts in sorted(species_counts.items()):
        writer.writerow([latin_name, section] + counts)

with open(
    os.path.join(OUTDIR, "sample_overlap.tsv"), "w", newline="", encoding="utf-8"
) as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["source", "total_samples", "phenotype_overlap", "extra_samples"])
    for source, values in sets.items():
        writer.writerow(
            [source, len(values), len(values & set(sample_order)), len(values - set(sample_order))]
        )

print(f"Prepared {len(sample_order)} samples in {OUTDIR}")
print("Phenotype counts: spineless=95, spiny=49")
