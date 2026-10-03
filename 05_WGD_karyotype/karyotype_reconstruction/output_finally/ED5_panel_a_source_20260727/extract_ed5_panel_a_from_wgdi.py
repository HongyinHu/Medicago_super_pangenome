#!/usr/bin/env python3
"""Extract ED5 panel-a source data from the underlying WGDI analysis outputs.

This script deliberately reads WGDI ``km_result.txt`` tables and their matching
``*.lens``, ``*.gff`` and ``total.conf`` files.  It does not inspect, rasterize,
or measure any PDF/image geometry.
"""

from __future__ import print_function

import argparse
import csv
import hashlib
import os
import re
import sys


ANALYSIS_ROOT_DEFAULT = (
    "path/to/project/"
    "37.karyotype_reconstruction"
)

# Species order follows the simplified phylogeny used in ED5a.
SPECIES = [
    ("Malbus", "Mal_out", 8, "genome_Mal"),
    ("Marc", "Mar_100", 8, "genome_Mar"),
    ("Mrut", "Mru_300", 8, "genome_Mru"),
    ("Mlan", "Mla_454", 8, "genome_454"),
    ("Mcar", "Mca_474", 8, "genome_474"),
    ("Mcre", "Mcr_468", 8, "genome_468"),
    ("Msat_cae", "Msa_T2T", 8, "genome_Msa"),
    ("Msat_zm4", "Msa_zm4", 8, "genome_ZM4"),
    ("Mmar", "Mma_457", 8, "genome_457"),
    ("Mpra", "Mpr_410", 7, "genome_410"),
    ("Mtru_R108", "Mtr_500", 8, "genome_R108"),
    # The 2026-06-18 rerun is the completed, corrected Mpol analysis.
    ("Mpol", "Mpo_200", 7, "genome_Mpo_rerun_20260618"),
    ("Morb", "Mor_22", 8, "genome_M22"),
    ("Msec", "Mse_461", 8, "genome_461"),
    ("Mlup", "Mlu_395", 8, "genome_395"),
    ("Msuf", "Msu_472", 8, "genome_472"),
    ("Medg", "Med_482", 8, "genome_482"),
    ("Mfis", "Mfi_46", 8, "genome_M46"),
    ("Mrad", "Mra_436", 8, "genome_436"),
]

WGDI_COLOR_TO_AMK = {
    "royalblue": ("AMK1", "#4169E1", "#3A6EA5"),
    "red": ("AMK2", "#FF0000", "#C44E52"),
    "#99cc00": ("AMK3", "#99CC00", "#6F74B8"),
    "deepskyblue": ("AMK4", "#00BFFF", "#4FB6C2"),
    "#339966": ("AMK5", "#339966", "#2F8F83"),
    "#ffcc00": ("AMK6", "#FFCC00", "#D99058"),
    "fuchsia": ("AMK7", "#FF00FF", "#B15AA3"),
    "#aa6e20": ("AMK8", "#AA6E20", "#8C6D31"),
}


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path, fieldnames, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def chromosome_key(value):
    text = str(value)
    numbers = re.findall(r"\d+", text)
    if numbers:
        return (0, int(numbers[-1]), text)
    return (1, 0, text)


def parse_section(path, wanted_section):
    section = None
    result = {}
    with open(path, "r") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith(";"):
                continue
            if line.startswith("[") and line.endswith("]"):
                section = line[1:-1].strip()
                continue
            if section == wanted_section and "=" in line:
                key, value = line.split("=", 1)
                result[key.strip()] = value.strip()
    return result


def resolve_input(output_dir, configured_name):
    if os.path.isabs(configured_name):
        path = configured_name
    else:
        path = os.path.join(output_dir, configured_name)
    path = os.path.abspath(path)
    if os.path.exists(path):
        return path, "configured_path"

    # Two datasets (Mcar and Mtru_R108) had their original WGDI inputs moved
    # into data/<genome>/old when newer annotations were staged in June 2026.
    # Their output-directory symlinks now dangle, but the exact same-named
    # analysis inputs remain archived.  Use those archived originals rather
    # than mixing newer annotations with the historical km_result.txt.
    if os.path.islink(path):
        link_target = os.readlink(path)
        if not os.path.isabs(link_target):
            link_target = os.path.normpath(
                os.path.join(os.path.dirname(path), link_target)
            )
        archived = os.path.join(
            os.path.dirname(link_target), "old", os.path.basename(link_target)
        )
        if os.path.exists(archived):
            return os.path.abspath(archived), "archived_original_for_dangling_symlink"

    raise FileNotFoundError(path)


def read_lens(path):
    rows = {}
    with open(path, "r") as handle:
        for line_number, raw in enumerate(handle, 1):
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 3:
                raise ValueError(
                    "Malformed lens row {}:{}: {!r}".format(path, line_number, raw)
                )
            chromosome = str(parts[0])
            rows[chromosome] = {
                "chromosome_length_bp": int(float(parts[1])),
                "chromosome_gene_count": int(float(parts[2])),
            }
    return rows


def read_gff(path):
    genes = {}
    chromosome_max_order = {}
    with open(path, "r") as handle:
        for line_number, raw in enumerate(handle, 1):
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 6:
                raise ValueError(
                    "Malformed WGDI gff row {}:{}: {!r}".format(
                        path, line_number, raw
                    )
                )
            chromosome = str(parts[0])
            gene_id = parts[1]
            order = int(float(parts[5]))
            key = (chromosome, order)
            if key in genes:
                raise ValueError(
                    "Duplicate chromosome/order in {}: {} {}".format(
                        path, chromosome, order
                    )
                )
            genes[key] = {
                "gene_id": gene_id,
                "start_bp": int(float(parts[2])),
                "end_bp": int(float(parts[3])),
                "strand": parts[4],
            }
            chromosome_max_order[chromosome] = max(
                order, chromosome_max_order.get(chromosome, 0)
            )
    return genes, chromosome_max_order


def read_km_result(path):
    rows = []
    with open(path, "r") as handle:
        for line_number, raw in enumerate(handle, 1):
            parts = raw.rstrip("\n").split("\t")
            if len(parts) != 5:
                raise ValueError(
                    "Expected five columns in {}:{}: {!r}".format(
                        path, line_number, raw
                    )
                )
            rows.append(
                {
                    "chromosome": str(parts[0]),
                    "start_gene_index": int(float(parts[1])),
                    "end_gene_index": int(float(parts[2])),
                    "wgdi_color": parts[3],
                    "classification": int(float(parts[4])),
                }
            )
    rows.sort(
        key=lambda row: (
            chromosome_key(row["chromosome"]),
            row["start_gene_index"],
            row["end_gene_index"],
        )
    )
    return rows


def format_fraction(value):
    return "{:.8f}".format(value)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-root", default=ANALYSIS_ROOT_DEFAULT)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    analysis_root = os.path.abspath(args.analysis_root)
    output_root = os.path.abspath(args.output_dir)
    os.makedirs(output_root, exist_ok=True)

    all_blocks = []
    metadata_rows = []
    chromosome_qc_rows = []
    species_qc_rows = []
    manifest_rows = []
    palette_source = None
    fatal_errors = []

    def add_manifest(species, role, path, resolution="configured_path"):
        absolute = os.path.abspath(path)
        manifest_rows.append(
            {
                "species": species,
                "role": role,
                "path_resolution": resolution,
                "path": absolute,
                "realpath": os.path.realpath(absolute),
                "bytes": os.path.getsize(absolute),
                "sha256": sha256_file(absolute),
            }
        )

    wgdi_version = "unknown"
    wgdi_karyotype_mapping_source = ""
    try:
        import pkg_resources
        import wgdi.karyotype_mapping as wgdi_karyotype_mapping

        wgdi_version = pkg_resources.get_distribution("wgdi").version
        wgdi_karyotype_mapping_source = os.path.abspath(
            wgdi_karyotype_mapping.__file__
        )
        add_manifest(
            "__software__",
            "WGDI karyotype_mapping implementation",
            wgdi_karyotype_mapping_source,
            "installed_environment",
        )
    except Exception as exc:
        print("WARNING: could not inspect installed WGDI package: {}".format(exc))

    for species_order, (
        species,
        tree_tip,
        expected_chr_count,
        analysis_subdir,
    ) in enumerate(SPECIES, 1):
        output_dir = os.path.join(analysis_root, "output_Kary", analysis_subdir)
        config_path = os.path.join(output_dir, "total.conf")
        if not os.path.isdir(output_dir):
            fatal_errors.append("{}: missing {}".format(species, output_dir))
            continue
        if not os.path.isfile(config_path):
            fatal_errors.append("{}: missing {}".format(species, config_path))
            continue

        config = parse_section(config_path, "karyotype_mapping")
        required_keys = [
            "gff1",
            "the_other_lens",
            "the_other_ancestor_file",
            "ancestor_top",
        ]
        missing_keys = [key for key in required_keys if key not in config]
        if missing_keys:
            fatal_errors.append(
                "{}: total.conf missing {}".format(species, ",".join(missing_keys))
            )
            continue

        try:
            gff_path, gff_resolution = resolve_input(output_dir, config["gff1"])
            lens_path, lens_resolution = resolve_input(
                output_dir, config["the_other_lens"]
            )
            km_path, km_resolution = resolve_input(
                output_dir, config["the_other_ancestor_file"]
            )
            ancestor_path, ancestor_resolution = resolve_input(
                output_dir, config["ancestor_top"]
            )
            blockinfo_path, blockinfo_resolution = resolve_input(
                output_dir, config["blockinfo"]
            )
        except (FileNotFoundError, KeyError) as exc:
            fatal_errors.append("{}: unresolved input {}".format(species, exc))
            continue

        for role, path, resolution in [
            ("WGDI karyotype mapping output", km_path, km_resolution),
            ("WGDI chromosome lens", lens_path, lens_resolution),
            ("WGDI ordered gene table", gff_path, gff_resolution),
            ("WGDI run configuration", config_path, "configured_path"),
            (
                "WGDI ancestral karyotype definition",
                ancestor_path,
                ancestor_resolution,
            ),
            ("WGDI block-information input", blockinfo_path, blockinfo_resolution),
        ]:
            add_manifest(species, role, path, resolution)

        if palette_source is None:
            palette_source = ancestor_path

        lens = read_lens(lens_path)
        genes, gff_max_order = read_gff(gff_path)
        km_rows = read_km_result(km_path)
        chromosomes = sorted(
            set(row["chromosome"] for row in km_rows), key=chromosome_key
        )
        species_max_gene_count = max(
            lens[chromosome]["chromosome_gene_count"] for chromosome in chromosomes
        )

        species_ok = True
        species_messages = []
        if len(chromosomes) != expected_chr_count:
            species_ok = False
            species_messages.append(
                "observed {} chromosomes, expected {}".format(
                    len(chromosomes), expected_chr_count
                )
            )

        segment_number = 0
        for chromosome in chromosomes:
            segment_rows = [
                row for row in km_rows if row["chromosome"] == chromosome
            ]
            segment_rows.sort(key=lambda row: row["start_gene_index"])
            lens_row = lens.get(chromosome)
            if lens_row is None:
                fatal_errors.append(
                    "{}: chromosome {} absent from {}".format(
                        species, chromosome, lens_path
                    )
                )
                continue

            chromosome_gene_count = lens_row["chromosome_gene_count"]
            chromosome_length_bp = lens_row["chromosome_length_bp"]
            starts_at_one = bool(segment_rows) and (
                segment_rows[0]["start_gene_index"] == 1
            )
            ends_at_lens = bool(segment_rows) and (
                segment_rows[-1]["end_gene_index"] == chromosome_gene_count
            )
            contiguous = all(
                segment_rows[index]["start_gene_index"]
                == segment_rows[index - 1]["end_gene_index"] + 1
                for index in range(1, len(segment_rows))
            )
            gff_matches_lens = (
                gff_max_order.get(chromosome) == chromosome_gene_count
            )
            boundary_genes_complete = True
            unknown_colors = []

            for segment_index, row in enumerate(segment_rows, 1):
                segment_number += 1
                color_key = row["wgdi_color"].strip().lower()
                if color_key not in WGDI_COLOR_TO_AMK:
                    unknown_colors.append(row["wgdi_color"])
                    continue
                amk, source_hex, display_hex = WGDI_COLOR_TO_AMK[color_key]
                start_gene = genes.get((chromosome, row["start_gene_index"]))
                end_gene = genes.get((chromosome, row["end_gene_index"]))
                if start_gene is None or end_gene is None:
                    boundary_genes_complete = False
                    continue
                segment_gene_count = (
                    row["end_gene_index"] - row["start_gene_index"] + 1
                )
                all_blocks.append(
                    {
                        "species_order": species_order,
                        "species": species,
                        "tree_tip": tree_tip,
                        "basic_chromosome_number": expected_chr_count,
                        "analysis_directory": output_dir,
                        "analysis_subdirectory": analysis_subdir,
                        "chromosome": chromosome,
                        "segment_index": segment_index,
                        "AMK": amk,
                        "wgdi_color": row["wgdi_color"],
                        "source_hex_color": source_hex,
                        "display_hex_color": display_hex,
                        "wgdi_classification": row["classification"],
                        "start_gene_index": row["start_gene_index"],
                        "end_gene_index": row["end_gene_index"],
                        "segment_gene_count": segment_gene_count,
                        "chromosome_gene_count": chromosome_gene_count,
                        "chromosome_length_bp": chromosome_length_bp,
                        "species_max_gene_count": species_max_gene_count,
                        "start_fraction_of_species_max": format_fraction(
                            (row["start_gene_index"] - 1)
                            / float(species_max_gene_count)
                        ),
                        "end_fraction_of_species_max": format_fraction(
                            row["end_gene_index"] / float(species_max_gene_count)
                        ),
                        "segment_fraction_of_species_max": format_fraction(
                            segment_gene_count / float(species_max_gene_count)
                        ),
                        "start_fraction_of_chromosome": format_fraction(
                            (row["start_gene_index"] - 1)
                            / float(chromosome_gene_count)
                        ),
                        "end_fraction_of_chromosome": format_fraction(
                            row["end_gene_index"] / float(chromosome_gene_count)
                        ),
                        "start_boundary_gene_id": start_gene["gene_id"],
                        "start_boundary_gene_start_bp": start_gene["start_bp"],
                        "start_boundary_gene_end_bp": start_gene["end_bp"],
                        "end_boundary_gene_id": end_gene["gene_id"],
                        "end_boundary_gene_start_bp": end_gene["start_bp"],
                        "end_boundary_gene_end_bp": end_gene["end_bp"],
                        "source_km_result": km_path,
                        "source_lens": lens_path,
                        "source_gff": gff_path,
                        "source_total_conf": config_path,
                    }
                )

            chromosome_ok = (
                starts_at_one
                and ends_at_lens
                and contiguous
                and gff_matches_lens
                and boundary_genes_complete
                and not unknown_colors
            )
            if not chromosome_ok:
                species_ok = False
                species_messages.append(
                    "chromosome {} failed continuity/input checks".format(chromosome)
                )
            chromosome_qc_rows.append(
                {
                    "species": species,
                    "chromosome": chromosome,
                    "segment_count": len(segment_rows),
                    "starts_at_gene_1": starts_at_one,
                    "segments_are_contiguous": contiguous,
                    "ends_at_lens_gene_count": ends_at_lens,
                    "gff_max_order_matches_lens": gff_matches_lens,
                    "boundary_genes_complete": boundary_genes_complete,
                    "unknown_colors": ",".join(sorted(set(unknown_colors))),
                    "chromosome_gene_count": chromosome_gene_count,
                    "chromosome_length_bp": chromosome_length_bp,
                    "status": "PASS" if chromosome_ok else "FAIL",
                }
            )

        metadata_rows.append(
            {
                "species_order": species_order,
                "species": species,
                "tree_tip": tree_tip,
                "basic_chromosome_number": expected_chr_count,
                "analysis_subdirectory": analysis_subdir,
                "analysis_directory": output_dir,
                "observed_chromosome_count": len(chromosomes),
                "extracted_segment_count": segment_number,
                "species_max_gene_count": species_max_gene_count,
                "total_gene_count": sum(
                    lens[chromosome]["chromosome_gene_count"]
                    for chromosome in chromosomes
                ),
                "source_data_basis": (
                    "WGDI karyotype_mapping km_result.txt; coordinates are "
                    "ordered-gene indices normalized by the species maximum "
                    "chromosome gene count"
                ),
                "source_km_result": km_path,
                "source_lens": lens_path,
                "source_gff": gff_path,
                "source_total_conf": config_path,
            }
        )
        species_qc_rows.append(
            {
                "species": species,
                "expected_chromosome_count": expected_chr_count,
                "observed_chromosome_count": len(chromosomes),
                "segment_count": segment_number,
                "status": "PASS" if species_ok else "FAIL",
                "notes": "; ".join(species_messages),
            }
        )

    if fatal_errors:
        for message in fatal_errors:
            print("ERROR:", message, file=sys.stderr)
        return 2

    # Read the AMK definition directly from the actual WGDI ancestor input.
    palette_rows = []
    with open(palette_source, "r") as handle:
        for raw in handle:
            parts = raw.rstrip("\n").split("\t")
            if len(parts) != 5:
                raise ValueError("Malformed ancestor definition: {!r}".format(raw))
            color_key = parts[3].strip().lower()
            if color_key not in WGDI_COLOR_TO_AMK:
                raise ValueError("Unknown ancestral colour {}".format(parts[3]))
            amk, source_hex, display_hex = WGDI_COLOR_TO_AMK[color_key]
            palette_rows.append(
                {
                    "AMK": amk,
                    "wgdi_color": parts[3],
                    "source_hex_color": source_hex,
                    "display_hex_color": display_hex,
                    "ancestral_chromosome": parts[0],
                    "ancestral_start_gene_index": parts[1],
                    "ancestral_end_gene_index": parts[2],
                    "wgdi_classification": parts[4],
                    "source_ancestor_file": palette_source,
                }
            )
    palette_rows.sort(key=lambda row: int(row["AMK"].replace("AMK", "")))

    all_blocks.sort(
        key=lambda row: (
            int(row["species_order"]),
            chromosome_key(row["chromosome"]),
            int(row["segment_index"]),
        )
    )
    manifest_rows.sort(key=lambda row: (row["species"], row["role"], row["path"]))

    block_fields = [
        "species_order",
        "species",
        "tree_tip",
        "basic_chromosome_number",
        "analysis_directory",
        "analysis_subdirectory",
        "chromosome",
        "segment_index",
        "AMK",
        "wgdi_color",
        "source_hex_color",
        "display_hex_color",
        "wgdi_classification",
        "start_gene_index",
        "end_gene_index",
        "segment_gene_count",
        "chromosome_gene_count",
        "chromosome_length_bp",
        "species_max_gene_count",
        "start_fraction_of_species_max",
        "end_fraction_of_species_max",
        "segment_fraction_of_species_max",
        "start_fraction_of_chromosome",
        "end_fraction_of_chromosome",
        "start_boundary_gene_id",
        "start_boundary_gene_start_bp",
        "start_boundary_gene_end_bp",
        "end_boundary_gene_id",
        "end_boundary_gene_start_bp",
        "end_boundary_gene_end_bp",
        "source_km_result",
        "source_lens",
        "source_gff",
        "source_total_conf",
    ]
    write_tsv(
        os.path.join(output_root, "panel_a_ancestral_karyotype_blocks.tsv"),
        block_fields,
        all_blocks,
    )
    write_tsv(
        os.path.join(output_root, "panel_a_species_metadata.tsv"),
        list(metadata_rows[0].keys()),
        metadata_rows,
    )
    write_tsv(
        os.path.join(output_root, "panel_a_amk_palette.tsv"),
        list(palette_rows[0].keys()),
        palette_rows,
    )
    write_tsv(
        os.path.join(output_root, "panel_a_chromosome_QC.tsv"),
        list(chromosome_qc_rows[0].keys()),
        chromosome_qc_rows,
    )
    write_tsv(
        os.path.join(output_root, "panel_a_species_QC.tsv"),
        list(species_qc_rows[0].keys()),
        species_qc_rows,
    )
    write_tsv(
        os.path.join(output_root, "panel_a_input_manifest.tsv"),
        list(manifest_rows[0].keys()),
        manifest_rows,
    )

    failed_species = [row for row in species_qc_rows if row["status"] != "PASS"]
    failed_chromosomes = [
        row for row in chromosome_qc_rows if row["status"] != "PASS"
    ]
    summary_rows = [
        {"metric": "species_count", "value": len(metadata_rows)},
        {
            "metric": "chromosome_count",
            "value": len(chromosome_qc_rows),
        },
        {"metric": "segment_count", "value": len(all_blocks)},
        {
            "metric": "failed_species_count",
            "value": len(failed_species),
        },
        {
            "metric": "failed_chromosome_count",
            "value": len(failed_chromosomes),
        },
        {
            "metric": "Mpol_analysis_subdirectory",
            "value": "genome_Mpo_rerun_20260618",
        },
        {"metric": "WGDI_version", "value": wgdi_version},
        {
            "metric": "WGDI_karyotype_mapping_source",
            "value": wgdi_karyotype_mapping_source,
        },
        {
            "metric": "coordinate_definition",
            "value": (
                "WGDI ordered-gene indices; panel drawing uses fractions of "
                "the maximum chromosome gene count within each species"
            ),
        },
        {
            "metric": "image_or_PDF_measurement_used",
            "value": "FALSE",
        },
        {
            "metric": "overall_status",
            "value": (
                "PASS"
                if not failed_species and not failed_chromosomes
                else "FAIL"
            ),
        },
    ]
    write_tsv(
        os.path.join(output_root, "panel_a_extraction_summary.tsv"),
        ["metric", "value"],
        summary_rows,
    )

    print(
        "Extracted {} species, {} chromosomes and {} WGDI segments".format(
            len(metadata_rows), len(chromosome_qc_rows), len(all_blocks)
        )
    )
    if failed_species or failed_chromosomes:
        print("QC FAILED", file=sys.stderr)
        return 3
    print("QC PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
