#!/usr/bin/env python3
import argparse
import re
from collections import Counter
from pathlib import Path

import pandas as pd


COPIA_ORDER = [
    "Ale", "Alesia", "Angela", "Bianca", "Ikeros", "Ivana",
    "Lyco", "SIRE", "TAR", "Tork", "Copia_unclassified",
]
GYPSY_ORDER = [
    "Athila", "CRM", "Galadriel", "Ogre", "Reina", "Retand",
    "Tekay", "Gypsy_unclassified",
]


def read_bed(path):
    rows = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            if re.match(r"^Chr[0-9]+$", f[0]):
                rows.append((f[0], int(f[1]), int(f[2])))
    return rows


def read_tesorter_cls(path):
    cls = {}
    if not Path(path).exists():
        return cls
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            raw, order, superfamily, clade, complete, strand, domains = f[:7]
            key = raw.split("#", 1)[0]
            cls[key] = {
                "order": order,
                "superfamily": superfamily,
                "clade": clade,
                "complete": complete,
                "strand": strand,
                "domains": domains,
                "has_rt": bool(re.search(r"(^|[|,; ])RT($|[|,; ])", domains)),
            }
    return cls


def classify(row, cls):
    rec = cls.get(row["seq_id"])
    if rec and rec["order"] == "LTR" and rec["superfamily"] in {"Copia", "Gypsy"}:
        clade = rec["clade"]
        if clade and clade not in {"-", "unknown", "unclassified"}:
            return rec["superfamily"], clade, rec["has_rt"], True
    sf = str(row["pass_superfamily"])
    if sf == "Copia":
        return "Copia", "Copia_unclassified", bool(rec and rec["has_rt"]), rec is not None
    if sf == "Gypsy":
        return "Gypsy", "Gypsy_unclassified", bool(rec and rec["has_rt"]), rec is not None
    return "Unknown intact LTR", "Unknown intact LTR", bool(rec and rec["has_rt"]), rec is not None


def overlap_any(row, regions):
    for chrom, start, end in regions:
        if row["chr"] == chrom and int(row["start"]) < end and int(row["end"]) > start:
            return True
    return False


def summarize(df, regions, genome_size_bp, centromere_size_bp):
    cent = df[df.apply(lambda r: overlap_any(r, regions), axis=1)].copy()
    out = {}
    for scope, sub, denom in [
        ("genome", df, genome_size_bp),
        ("centromere", cent, centromere_size_bp),
    ]:
        counts = Counter()
        for _, row in sub.iterrows():
            counts[(row["superfamily"], row["clade"])] += 1
        bp = float(sub["length_bp"].sum()) if len(sub) else 0
        out[scope] = {
            "counts": counts,
            "total_count": int(len(sub)),
            "total_size_mb": bp / 1_000_000,
            "ratio_pct": bp / denom * 100 if denom else 0,
            "rt_domain_count": int(sub["has_rt"].sum()) if len(sub) else 0,
            "tesorter_matched_count": int(sub["tesorter_matched"].sum()) if len(sub) else 0,
        }
    return out, cent


def ordered_rows(all_clades):
    rows = []
    for clade in COPIA_ORDER:
        rows.append(("Copia", clade))
    rows.extend(sorted((sf, c) for sf, c in all_clades if sf == "Copia" and c not in COPIA_ORDER))
    for clade in GYPSY_ORDER:
        rows.append(("Gypsy", clade))
    rows.extend(sorted((sf, c) for sf, c in all_clades if sf == "Gypsy" and c not in GYPSY_ORDER))
    rows.append(("Unknown intact LTR", "Unknown intact LTR"))
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    config = pd.read_csv(args.config, sep="\t")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    summaries = {}
    all_clades = set()
    details = []
    metrics = []

    for _, cfg in config.iterrows():
        label = cfg["label"]
        meta = pd.read_csv(cfg["pass_meta"], sep="\t")
        cls = read_tesorter_cls(cfg["tesorter_cls"])
        classified = meta.apply(lambda row: classify(row, cls), axis=1)
        meta["superfamily"] = [x[0] for x in classified]
        meta["clade"] = [x[1] for x in classified]
        meta["has_rt"] = [x[2] for x in classified]
        meta["tesorter_matched"] = [x[3] for x in classified]
        regions = read_bed(cfg["centromere_bed"])
        genome_size = int(cfg["genome_size_bp_chr_only"])
        centromere_size = sum(end - start for _, start, end in regions)
        summary, cent = summarize(meta, regions, genome_size, centromere_size)
        summaries[label] = summary
        all_clades.update(summary["genome"]["counts"])
        all_clades.update(summary["centromere"]["counts"])
        meta["genome"] = cfg["genome"]
        meta["label"] = label
        meta["in_functional_centromere"] = meta.apply(lambda r: overlap_any(r, regions), axis=1)
        details.append(meta)
        metrics.append({
            "genome": cfg["genome"],
            "label": label,
            "pass_meta": cfg["pass_meta"],
            "tesorter_cls": cfg["tesorter_cls"],
            "genome_size_bp_chr_only": genome_size,
            "functional_centromere_size_bp_chr_only": centromere_size,
            "intact_LTRRT_genome_count": summary["genome"]["total_count"],
            "intact_LTRRT_centromere_count": summary["centromere"]["total_count"],
            "intact_LTRRT_genome_size_mb": summary["genome"]["total_size_mb"],
            "intact_LTRRT_centromere_size_mb": summary["centromere"]["total_size_mb"],
            "intact_LTRRT_genome_ratio_pct": summary["genome"]["ratio_pct"],
            "intact_LTRRT_centromere_ratio_pct": summary["centromere"]["ratio_pct"],
            "intact_LTRRT_with_RT_domain_genome_count": summary["genome"]["rt_domain_count"],
            "intact_LTRRT_with_RT_domain_centromere_count": summary["centromere"]["rt_domain_count"],
            "tesorter_matched_genome_count": summary["genome"]["tesorter_matched_count"],
            "tesorter_matched_centromere_count": summary["centromere"]["tesorter_matched_count"],
        })

    table_rows = []
    for group, clade in ordered_rows(all_clades):
        record = {"Group": group, "Clade": clade}
        for label in config["label"].tolist():
            record[f"{label} centromere"] = summaries[label]["centromere"]["counts"].get((group, clade), 0)
            record[f"{label} genome"] = summaries[label]["genome"]["counts"].get((group, clade), 0)
        table_rows.append(record)

    for metric_name, key, fmt in [
        ("Total", "total_count", "int"),
        ("Total Size (Mb)", "total_size_mb", "float"),
        ("Ratio (%)", "ratio_pct", "float"),
        ("No. of intact LTR with RT domain", "rt_domain_count", "int"),
    ]:
        record = {"Group": "", "Clade": metric_name}
        for label in config["label"].tolist():
            for scope in ["centromere", "genome"]:
                value = summaries[label][scope][key]
                record[f"{label} {scope}"] = int(value) if fmt == "int" else round(float(value), 2)
        table_rows.append(record)

    table = pd.DataFrame(table_rows)
    table.to_csv(outdir / "intact_LTRRT_classification_summary.TableS8_style.flat.tsv", sep="\t", index=False)
    pd.concat(details, ignore_index=True).to_csv(outdir / "intact_LTRRT_classification_detail.tsv", sep="\t", index=False)
    pd.DataFrame(metrics).to_csv(outdir / "intact_LTRRT_summary_metrics.tsv", sep="\t", index=False, float_format="%.4f")

    count_rows = {"Total", "No. of intact LTR with RT domain"}
    with open(outdir / "intact_LTRRT_classification_summary.TableS8_style.tsv", "w") as handle:
        first = ["", ""]
        second = ["Group", "Clade"]
        for label in config["label"].tolist():
            first.extend([label, ""])
            second.extend(["centromere", "genome"])
        handle.write("\t".join(first) + "\n")
        handle.write("\t".join(second) + "\n")
        for _, row in table.iterrows():
            values = [str(row["Group"]), str(row["Clade"])]
            for label in config["label"].tolist():
                for scope in ["centromere", "genome"]:
                    value = row[f"{label} {scope}"]
                    if row["Group"] or row["Clade"] in count_rows:
                        values.append(str(int(round(float(value)))))
                    else:
                        values.append(str(value))
            handle.write("\t".join(values) + "\n")

    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
