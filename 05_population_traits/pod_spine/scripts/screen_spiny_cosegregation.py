#!/usr/bin/env python3
import argparse
import csv
import os
from collections import Counter, defaultdict


def read_tsv(path):
    with open(path, encoding="utf-8", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_inputs(path):
    paths = {}
    for row in read_tsv(path):
        paths[row["key"]] = row["path"]
    return paths


def load_phenotypes(path):
    rows = read_tsv(path)
    return {row["genome_id"]: row for row in rows}


def model_groups(phenotypes, model):
    spiny = []
    spineless = []
    for genome, row in phenotypes.items():
        use = True
        value = "NA"
        if model == "main_high_quality":
            use = row["use_main_model"] == "1"
            value = row["spiny_primary"]
        elif model == "cladeI_internal_high_quality":
            use = row["use_main_model"] == "1" and row["clade"] == "Clade-I"
            value = row["spiny_primary"]
        elif model == "include_low_quality":
            use = row["spiny_primary"] != "NA"
            value = row["spiny_primary"]
        elif model == "include_weak_spiny":
            use = row["spiny_inclusive"] != "NA"
            value = row["spiny_inclusive"]
        else:
            raise SystemExit(f"unknown model: {model}")

        if not use or value == "NA":
            continue
        if value == "1":
            spiny.append(genome)
        elif value == "0":
            spineless.append(genome)
    return sorted(spiny), sorted(spineless)


def infer_sample_columns(header, phenotypes):
    samples = [name for name in header if name in phenotypes]
    if not samples:
        raise SystemExit("no sample columns found in PAV matrix")
    first = min(header.index(sample) for sample in samples)
    return header[:first], header[first:]


def state_counts(row, samples, sample_index):
    values = {}
    for sample in samples:
        i = sample_index.get(sample)
        if i is None or i >= len(row):
            values[sample] = "NA"
        else:
            value = row[i]
            values[sample] = value if value in {"0", "1"} else "NA"
    return values


def summarize_group(values, samples):
    count = Counter(values[sample] for sample in samples)
    n = len(samples)
    called = count["0"] + count["1"]
    present = count["1"]
    absent = count["0"]
    missing = count["NA"]
    present_rate = float(present) / called if called else 0.0
    absent_rate = float(absent) / called if called else 0.0
    missing_rate = float(missing) / n if n else 0.0
    return {
        "n": n,
        "called": called,
        "present": present,
        "absent": absent,
        "missing": missing,
        "present_rate": present_rate,
        "absent_rate": absent_rate,
        "missing_rate": missing_rate,
    }


def is_candidate(spiny_summary, spineless_summary, args):
    min_called_spiny = max(args.min_called_spiny, int(round(spiny_summary["n"] * args.min_called_fraction_spiny + 0.499)))
    min_called_spineless = max(args.min_called_spineless, int(round(spineless_summary["n"] * args.min_called_fraction_spineless + 0.499)))
    if spiny_summary["called"] < min_called_spiny or spineless_summary["called"] < min_called_spineless:
        return None
    if spiny_summary["missing_rate"] > args.max_group_missing_rate or spineless_summary["missing_rate"] > args.max_group_missing_rate:
        return None
    if (
        spiny_summary["present_rate"] >= args.min_spiny_rate
        and spineless_summary["present_rate"] <= args.max_spineless_rate
    ):
        return "spiny_present"
    if (
        spiny_summary["present_rate"] <= args.max_spineless_rate
        and spineless_summary["present_rate"] >= args.min_spiny_rate
    ):
        return "spiny_absent"
    return None


def clade_rates(values, phenotypes):
    clade_i = [sample for sample, row in phenotypes.items() if row["clade"] == "Clade-I"]
    non_clade_i = [sample for sample, row in phenotypes.items() if row["clade"] != "Clade-I"]
    return summarize_group(values, clade_i), summarize_group(values, non_clade_i)


def format_samples(samples, values, state):
    return ",".join(sample for sample in samples if values.get(sample, "NA") == state)


def screen_matrix(dataset, matrix_path, phenotypes, models, out_dir, args):
    os.makedirs(out_dir, exist_ok=True)
    support = defaultdict(lambda: {"models": set(), "directions": set(), "datasets": set()})
    summary_rows = []
    with open(matrix_path, encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        meta_cols, sample_cols = infer_sample_columns(header, phenotypes)
        sample_index = {name: header.index(name) for name in sample_cols}

        model_files = {}
        writers = {}
        for model in models:
            model_dir = os.path.join(out_dir, model)
            os.makedirs(model_dir, exist_ok=True)
            out_path = os.path.join(model_dir, "candidates.tsv")
            out_handle = open(out_path, "w", encoding="utf-8", newline="")
            model_files[model] = out_handle
            spiny, spineless = model_groups(phenotypes, model)
            fieldnames = (
                ["dataset", "model", "association_direction"]
                + meta_cols
                + [
                    "spiny_n",
                    "spiny_called",
                    "spiny_present",
                    "spiny_absent",
                    "spiny_missing",
                    "spiny_present_rate",
                    "spineless_n",
                    "spineless_called",
                    "spineless_present",
                    "spineless_absent",
                    "spineless_missing",
                    "spineless_present_rate",
                    "cladeI_present_rate",
                    "noncladeI_present_rate",
                    "spiny_present_samples",
                    "spineless_present_samples",
                    "spiny_absent_samples",
                    "spineless_absent_samples",
                ]
            )
            writer = csv.DictWriter(out_handle, delimiter="\t", fieldnames=fieldnames)
            writer.writeheader()
            writers[model] = (writer, spiny, spineless, fieldnames)

        total_rows = 0
        candidate_counts = Counter()
        for line in handle:
            if not line.strip():
                continue
            total_rows += 1
            row = line.rstrip("\n").split("\t")
            values = state_counts(row, sample_cols, sample_index)
            clade_i_summary, non_clade_i_summary = clade_rates(values, phenotypes)
            meta = {name: row[i] if i < len(row) else "." for i, name in enumerate(meta_cols)}
            sv_id = meta.get("MetaSV_ID") or meta.get("SV_ID") or meta.get("ID") or row[0]

            for model, (writer, spiny, spineless, fieldnames) in writers.items():
                spiny_summary = summarize_group(values, spiny)
                spineless_summary = summarize_group(values, spineless)
                direction = is_candidate(spiny_summary, spineless_summary, args)
                if not direction:
                    continue
                out = {
                    "dataset": dataset,
                    "model": model,
                    "association_direction": direction,
                    "spiny_n": spiny_summary["n"],
                    "spiny_called": spiny_summary["called"],
                    "spiny_present": spiny_summary["present"],
                    "spiny_absent": spiny_summary["absent"],
                    "spiny_missing": spiny_summary["missing"],
                    "spiny_present_rate": "%.6f" % spiny_summary["present_rate"],
                    "spineless_n": spineless_summary["n"],
                    "spineless_called": spineless_summary["called"],
                    "spineless_present": spineless_summary["present"],
                    "spineless_absent": spineless_summary["absent"],
                    "spineless_missing": spineless_summary["missing"],
                    "spineless_present_rate": "%.6f" % spineless_summary["present_rate"],
                    "cladeI_present_rate": "%.6f" % clade_i_summary["present_rate"],
                    "noncladeI_present_rate": "%.6f" % non_clade_i_summary["present_rate"],
                    "spiny_present_samples": format_samples(spiny, values, "1"),
                    "spineless_present_samples": format_samples(spineless, values, "1"),
                    "spiny_absent_samples": format_samples(spiny, values, "0"),
                    "spineless_absent_samples": format_samples(spineless, values, "0"),
                }
                out.update(meta)
                writer.writerow({key: out.get(key, ".") for key in fieldnames})
                candidate_counts[(model, direction)] += 1
                support[sv_id]["models"].add(model)
                support[sv_id]["directions"].add(direction)
                support[sv_id]["datasets"].add(dataset)

        for model_handle in model_files.values():
            model_handle.close()
        for (model, direction), count in sorted(candidate_counts.items()):
            summary_rows.append(
                {
                    "dataset": dataset,
                    "model": model,
                    "association_direction": direction,
                    "candidate_count": count,
                    "total_tested": total_rows,
                }
            )
        summary_path = os.path.join(out_dir, "summary_counts.tsv")
        with open(summary_path, "w", encoding="utf-8", newline="") as out:
            writer = csv.DictWriter(
                out,
                delimiter="\t",
                fieldnames=["dataset", "model", "association_direction", "candidate_count", "total_tested"],
            )
            writer.writeheader()
            writer.writerows(summary_rows)
    return support


def write_model_groups(path, phenotypes, models):
    with open(path, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["model", "group", "n", "samples"])
        for model in models:
            spiny, spineless = model_groups(phenotypes, model)
            writer.writerow([model, "spiny", len(spiny), ",".join(spiny)])
            writer.writerow([model, "spineless", len(spineless), ",".join(spineless)])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-paths", required=True)
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--min-spiny-rate", type=float, default=0.75)
    parser.add_argument("--max-spineless-rate", type=float, default=0.10)
    parser.add_argument("--max-group-missing-rate", type=float, default=0.35)
    parser.add_argument("--min-called-spiny", type=int, default=3)
    parser.add_argument("--min-called-spineless", type=int, default=5)
    parser.add_argument("--min-called-fraction-spiny", type=float, default=0.75)
    parser.add_argument("--min-called-fraction-spineless", type=float, default=0.60)
    args = parser.parse_args()

    paths = load_inputs(args.input_paths)
    phenotypes = load_phenotypes(args.phenotype)
    models = [
        "main_high_quality",
        "cladeI_internal_high_quality",
        "include_low_quality",
        "include_weak_spiny",
    ]
    os.makedirs(args.out_dir, exist_ok=True)
    write_model_groups(os.path.join(args.out_dir, "model_groups.tsv"), phenotypes, models)

    datasets = [
        ("meta_publication", paths["meta_publication_pav"]),
        ("RefA", paths["refa_pav"]),
        ("RefB", paths["refb_pav"]),
    ]
    global_support = defaultdict(lambda: {"models": set(), "directions": set(), "datasets": set()})
    for dataset, matrix in datasets:
        support = screen_matrix(dataset, matrix, phenotypes, models, os.path.join(args.out_dir, dataset), args)
        for sv_id, item in support.items():
            global_support[(dataset, sv_id)]["models"].update(item["models"])
            global_support[(dataset, sv_id)]["directions"].update(item["directions"])
            global_support[(dataset, sv_id)]["datasets"].update(item["datasets"])

    with open(os.path.join(args.out_dir, "candidate_model_support.tsv"), "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["dataset", "sv_id", "models", "directions", "n_models", "main_support", "cladeI_internal_support"])
        for (dataset, sv_id), item in sorted(global_support.items(), key=lambda x: (x[0][0], x[0][1])):
            models_s = sorted(item["models"])
            directions_s = sorted(item["directions"])
            writer.writerow(
                [
                    dataset,
                    sv_id,
                    ",".join(models_s),
                    ",".join(directions_s),
                    len(models_s),
                    "1" if "main_high_quality" in item["models"] else "0",
                    "1" if "cladeI_internal_high_quality" in item["models"] else "0",
                ]
            )


if __name__ == "__main__":
    main()
