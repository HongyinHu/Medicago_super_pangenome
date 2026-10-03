#!/usr/bin/env bash
set -eo pipefail

RUN_DIR="path/to/project/N_4.pod_spiny/26_Msa_SPL_family_phylogeny_20260710"
REF_PEP="path/to/project/N_4.pod_spiny/00_data/7.orthology_geneName/orthology_plant.pep"
MSA_PEP="path/to/project/16.T2T_ref_function_anno/output/blast_tair/T2T_tait/ref_Msa.T2T_ctg.pep"
CONDA_SH="path/to/home/anaconda3/etc/profile.d/conda.sh"

mkdir -p "${RUN_DIR}"/{data,scripts,logs,results,trees,figures,summary,tmp}
cd "${RUN_DIR}"

source "${CONDA_SH}"
conda activate biosofeware
set -u

echo "[INFO] run_dir=${RUN_DIR}"
echo "[INFO] started $(date)"

python3 - <<'PY'
from pathlib import Path
import re

run = Path("path/to/project/N_4.pod_spiny/26_Msa_SPL_family_phylogeny_20260710")
ref_pep = Path("path/to/project/N_4.pod_spiny/00_data/7.orthology_geneName/orthology_plant.pep")
msa_pep = Path("path/to/project/16.T2T_ref_function_anno/output/blast_tair/T2T_tait/ref_Msa.T2T_ctg.pep")

def read_fasta(path):
    records = []
    name = None
    header = None
    seq = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, header, "".join(seq)))
                header = line[1:]
                name = header.split()[0]
                seq = []
            else:
                seq.append(re.sub(r"\s+", "", line))
        if name is not None:
            records.append((name, header, "".join(seq)))
    return records

def write_fasta(records, path, width=80):
    with open(path, "w", encoding="utf-8") as out:
        for name, desc, seq in records:
            out.write(f">{name}")
            if desc:
                out.write(f" {desc}")
            out.write("\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i+width] + "\n")

ref_records = []
for rid, header, seq in read_fasta(ref_pep):
    ref_records.append((rid, "source=AtOs_reference_SPL", seq.replace("*", "")))
write_fasta(ref_records, run / "data" / "AtOs_SPL.reference.clean.pep")

msa_raw = read_fasta(msa_pep)
gene_counts = {}
parsed = []
for tid, header, seq in msa_raw:
    m = re.search(r"(?:^|\s)gene=([^\s]+)", header)
    gene = m.group(1) if m else tid
    gene_counts[gene] = gene_counts.get(gene, 0) + 1
    parsed.append((tid, gene, header, seq.replace("*", "")))

clean_records = []
map_rows = []
for tid, gene, header, seq in parsed:
    clean_id = gene if gene_counts[gene] == 1 else f"{gene}__{tid}"
    clean_records.append((clean_id, f"transcript={tid} original_header={header}", seq))
    map_rows.append((clean_id, gene, tid, len(seq), header))
write_fasta(clean_records, run / "data" / "Msa.T2T.proteins.geneid.pep")

with open(run / "data" / "Msa.T2T.proteins.geneid.map.tsv", "w", encoding="utf-8") as out:
    out.write("protein_id\tgene_id\ttranscript_id\tprotein_length\toriginal_header\n")
    for row in map_rows:
        out.write("\t".join(map(str, row)) + "\n")

with open(run / "logs" / "prepare_fasta.log", "w", encoding="utf-8") as out:
    out.write(f"AtOs reference records\t{len(ref_records)}\n")
    out.write(f"Msa raw protein records\t{len(msa_raw)}\n")
    out.write(f"Msa unique gene ids\t{len(set(x[1] for x in parsed))}\n")
    out.write(f"Chr23997 present\t{any(x[1] == 'Chr23997' for x in parsed)}\n")
PY

echo "[INFO] building SPL HMM from At/Os references"
mafft --auto --thread 8 data/AtOs_SPL.reference.clean.pep > data/AtOs_SPL.reference.clean.mafft.fa 2> logs/mafft_reference.log
hmmbuild tmp/AtOs_SPL.reference.hmm data/AtOs_SPL.reference.clean.mafft.fa > logs/hmmbuild_AtOs_SPL.log
hmmsearch --cpu 16 --tblout results/Msa_vs_AtOs_SPL.hmm.tbl --domtblout results/Msa_vs_AtOs_SPL.hmm.domtbl tmp/AtOs_SPL.reference.hmm data/Msa.T2T.proteins.geneid.pep > logs/hmmsearch_Msa_vs_AtOs_SPL.log

echo "[INFO] BLASTP At/Os SPL references against genome_Msa proteins"
makeblastdb -in data/Msa.T2T.proteins.geneid.pep -dbtype prot -out tmp/Msa.T2T.proteins.geneid >/dev/null
blastp \
  -query data/AtOs_SPL.reference.clean.pep \
  -db tmp/Msa.T2T.proteins.geneid \
  -evalue 1e-5 \
  -max_target_seqs 1000 \
  -num_threads 16 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp" \
  -out results/AtOs_SPL_vs_Msa.blastp.tsv

echo "[INFO] selecting Msa SPL candidates"
python3 - <<'PY'
from pathlib import Path
import math

run = Path("path/to/project/N_4.pod_spiny/26_Msa_SPL_family_phylogeny_20260710")

def read_fasta(path):
    records = {}
    name = None
    desc = ""
    seq = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = (desc, "".join(seq))
                parts = line[1:].split(maxsplit=1)
                name = parts[0]
                desc = parts[1] if len(parts) > 1 else ""
                seq = []
            else:
                seq.append(line.strip())
        if name is not None:
            records[name] = (desc, "".join(seq))
    return records

def write_fasta(ids, records, path, prefix="Msa_"):
    with open(path, "w", encoding="utf-8") as out:
        for x in ids:
            desc, seq = records[x]
            out.write(f">{prefix}{x} source=genome_Msa {desc}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")

msa_records = read_fasta(run / "data" / "Msa.T2T.proteins.geneid.pep")
candidate = {}

# hmmsearch --tblout columns:
# target name, accession, query name, accession, full E-value, full score, ...
with open(run / "results" / "Msa_vs_AtOs_SPL.hmm.tbl", "r", encoding="utf-8", errors="replace") as handle:
    for line in handle:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split()
        target = parts[0]
        evalue = float(parts[4])
        score = float(parts[5])
        if evalue <= 1e-5 and score >= 20:
            rec = candidate.setdefault(target, {"hmm_evalue": evalue, "hmm_score": score, "blast_hits": []})
            if evalue < rec["hmm_evalue"]:
                rec["hmm_evalue"] = evalue
                rec["hmm_score"] = score

with open(run / "results" / "AtOs_SPL_vs_Msa.blastp.tsv", "r", encoding="utf-8", errors="replace") as handle:
    for line in handle:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 15:
            continue
        qseqid, sseqid = parts[0], parts[1]
        pident = float(parts[2])
        aln_len = int(parts[3])
        evalue = float(parts[10])
        bitscore = float(parts[11])
        qlen = int(parts[12])
        slen = int(parts[13])
        qcov = float(parts[14])
        if evalue <= 1e-5 and bitscore >= 40 and aln_len >= 40:
            rec = candidate.setdefault(sseqid, {"hmm_evalue": math.nan, "hmm_score": math.nan, "blast_hits": []})
            rec["blast_hits"].append((qseqid, evalue, bitscore, pident, aln_len, qcov, qlen, slen))

rows = []
for gid, rec in candidate.items():
    hits = sorted(rec["blast_hits"], key=lambda x: (-x[2], x[1], -x[5]))
    best = hits[0] if hits else ("NA", math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan)
    rows.append({
        "gene_id": gid,
        "protein_id": gid,
        "protein_length": len(msa_records.get(gid, ("", ""))[1]),
        "hmm_evalue": rec.get("hmm_evalue", math.nan),
        "hmm_score": rec.get("hmm_score", math.nan),
        "blast_hit_count": len(hits),
        "best_reference_SPL": best[0],
        "best_blast_evalue": best[1],
        "best_blast_bitscore": best[2],
        "best_blast_pident": best[3],
        "best_blast_aln_len": best[4],
        "best_blast_qcov": best[5],
        "is_Chr23997": "yes" if gid == "Chr23997" else "no",
    })

rows.sort(key=lambda r: (
    0 if r["gene_id"] == "Chr23997" else 1,
    str(r["best_reference_SPL"]),
    -float(r["best_blast_bitscore"]) if str(r["best_blast_bitscore"]) != "nan" else 0,
    r["gene_id"],
))
header = list(rows[0].keys()) if rows else []
with open(run / "results" / "Msa_SPL_candidate_table.tsv", "w", encoding="utf-8") as out:
    out.write("\t".join(header) + "\n")
    for r in rows:
        out.write("\t".join(str(r[h]) for h in header) + "\n")

ids = [r["gene_id"] for r in rows]
missing = [x for x in ids if x not in msa_records]
if missing:
    raise SystemExit(f"missing records: {missing[:5]}")
write_fasta(ids, msa_records, run / "data" / "Msa_SPL_candidates.pep", prefix="Msa_")
with open(run / "summary" / "candidate_counts.txt", "w", encoding="utf-8") as out:
    out.write(f"Msa_SPL_candidates\t{len(rows)}\n")
    out.write(f"Chr23997_detected\t{any(r['gene_id'] == 'Chr23997' for r in rows)}\n")
    out.write("selection\tHMM evalue<=1e-5 score>=20 OR BLASTP evalue<=1e-5 bitscore>=40 aln_len>=40\n")
PY

cat data/AtOs_SPL.reference.clean.pep data/Msa_SPL_candidates.pep > data/AtOs_Msa_SPL.combined.pep

echo "[INFO] aligning combined At/Os/Msa SPL proteins"
mafft --auto --thread 16 data/AtOs_Msa_SPL.combined.pep > trees/AtOs_Msa_SPL.combined.mafft.fa 2> logs/mafft_combined.log

echo "[INFO] running IQ-TREE2 ML tree"
iqtree2 -s trees/AtOs_Msa_SPL.combined.mafft.fa -m MFP -bb 1000 -alrt 1000 -T 16 -pre trees/AtOs_Msa_SPL.ML > logs/iqtree2_AtOs_Msa_SPL.log 2>&1

echo "[INFO] plotting tree and summarizing Chr23997 placement"
Rscript - <<'RS'
suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
})
run <- "path/to/project/N_4.pod_spiny/26_Msa_SPL_family_phylogeny_20260710"
tree_file <- file.path(run, "trees", "AtOs_Msa_SPL.ML.treefile")
tree <- read.tree(tree_file)
tree <- midpoint(tree)
tip <- tree$tip.label
is_msa <- grepl("^Msa_", tip)
is_at <- grepl("^AtSPL", tip)
is_os <- grepl("^OsSPL", tip)
tip_col <- ifelse(tip == "Msa_Chr23997", "#D55E00", ifelse(is_msa, "#C00000", ifelse(is_at, "#1B9E77", "#386CB0")))
tip_cex <- ifelse(tip == "Msa_Chr23997", 0.82, ifelse(is_msa, 0.58, 0.55))

pdf(file.path(run, "figures", "AtOs_Msa_SPL_ML_tree_rectangular.pdf"), width=12, height=16, onefile=FALSE)
par(mar=c(1,1,1,1))
plot(tree, type="phylogram", cex=tip_cex, font=1, label.offset=0.006,
     tip.color=tip_col, edge.width=0.8, no.margin=TRUE)
if (!is.null(tree$node.label)) {
  labs <- tree$node.label
  labs[!grepl("/", labs)] <- ""
  nodelabels(labs, frame="none", cex=0.32, adj=c(1.1,-0.2), col="#444444")
}
legend("topleft", legend=c("genome_Msa SPL candidate", "Chr23997", "Arabidopsis SPL", "Rice SPL"),
       col=c("#C00000", "#D55E00", "#1B9E77", "#386CB0"), pch=19, bty="n", cex=0.75)
dev.off()

png(file.path(run, "figures", "AtOs_Msa_SPL_ML_tree_rectangular.png"), width=4200, height=5600, res=350)
par(mar=c(1,1,1,1))
plot(tree, type="phylogram", cex=tip_cex, font=1, label.offset=0.006,
     tip.color=tip_col, edge.width=0.8, no.margin=TRUE)
if (!is.null(tree$node.label)) {
  labs <- tree$node.label
  labs[!grepl("/", labs)] <- ""
  nodelabels(labs, frame="none", cex=0.32, adj=c(1.1,-0.2), col="#444444")
}
legend("topleft", legend=c("genome_Msa SPL candidate", "Chr23997", "Arabidopsis SPL", "Rice SPL"),
       col=c("#C00000", "#D55E00", "#1B9E77", "#386CB0"), pch=19, bty="n", cex=0.75)
dev.off()

pdf(file.path(run, "figures", "AtOs_Msa_SPL_ML_tree_fan.pdf"), width=14, height=14, onefile=FALSE)
par(mar=c(0,0,0,0))
plot(tree, type="fan", cex=0.48, font=1, tip.color=tip_col, edge.width=0.8, no.margin=TRUE)
legend("bottomleft", legend=c("genome_Msa SPL candidate", "Chr23997", "Arabidopsis SPL", "Rice SPL"),
       col=c("#C00000", "#D55E00", "#1B9E77", "#386CB0"), pch=19, bty="n", cex=0.8)
dev.off()

png(file.path(run, "figures", "AtOs_Msa_SPL_ML_tree_fan.png"), width=4800, height=4800, res=350)
par(mar=c(0,0,0,0))
plot(tree, type="fan", cex=0.48, font=1, tip.color=tip_col, edge.width=0.8, no.margin=TRUE)
legend("bottomleft", legend=c("genome_Msa SPL candidate", "Chr23997", "Arabidopsis SPL", "Rice SPL"),
       col=c("#C00000", "#D55E00", "#1B9E77", "#386CB0"), pch=19, bty="n", cex=0.8)
dev.off()

D <- cophenetic.phylo(tree)
target <- "Msa_Chr23997"
refs <- tip[grepl("^(At|Os)SPL", tip)]
msa_tips <- tip[grepl("^Msa_", tip)]
nearest_refs <- sort(D[target, refs])[1:min(12, length(refs))]
nearest_msa <- sort(D[target, setdiff(msa_tips, target)])[1:min(12, length(msa_tips)-1)]

clade_rows <- data.frame()
for (r in refs) {
  node <- getMRCA(tree, c(target, r))
  members <- extract.clade(tree, node)$tip.label
  clade_rows <- rbind(clade_rows, data.frame(
    reference=r,
    mrca_tip_count=length(members),
    mrca_msa_count=sum(grepl("^Msa_", members)),
    mrca_ref_count=sum(grepl("^(At|Os)SPL", members)),
    distance=D[target, r],
    stringsAsFactors=FALSE
  ))
}
clade_rows <- clade_rows[order(clade_rows$mrca_tip_count, clade_rows$distance),]

write.table(data.frame(reference=names(nearest_refs), patristic_distance=as.numeric(nearest_refs)),
            file=file.path(run, "summary", "Chr23997_nearest_reference_SPLs.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(data.frame(msa_gene=sub("^Msa_", "", names(nearest_msa)), patristic_distance=as.numeric(nearest_msa)),
            file=file.path(run, "summary", "Chr23997_nearest_Msa_SPL_candidates.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(clade_rows,
            file=file.path(run, "summary", "Chr23997_reference_MRCA_clades.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

cat("Chr23997 nearest reference SPLs\n", file=file.path(run, "summary", "Chr23997_naming_summary.txt"))
capture.output(print(head(data.frame(reference=names(nearest_refs), distance=as.numeric(nearest_refs)), 8)),
               file=file.path(run, "summary", "Chr23997_naming_summary.txt"), append=TRUE)
cat("\nSmallest MRCA clades with references\n", file=file.path(run, "summary", "Chr23997_naming_summary.txt"), append=TRUE)
capture.output(print(head(clade_rows, 12)),
               file=file.path(run, "summary", "Chr23997_naming_summary.txt"), append=TRUE)
RS

echo "[INFO] validating outputs"
python3 - <<'PY'
from pathlib import Path
run = Path("path/to/project/N_4.pod_spiny/26_Msa_SPL_family_phylogeny_20260710")
required = [
    run / "results" / "Msa_SPL_candidate_table.tsv",
    run / "trees" / "AtOs_Msa_SPL.ML.treefile",
    run / "figures" / "AtOs_Msa_SPL_ML_tree_rectangular.pdf",
    run / "figures" / "AtOs_Msa_SPL_ML_tree_rectangular.png",
    run / "figures" / "AtOs_Msa_SPL_ML_tree_fan.pdf",
    run / "figures" / "AtOs_Msa_SPL_ML_tree_fan.png",
    run / "summary" / "Chr23997_naming_summary.txt",
]
bad = [str(p) for p in required if not p.exists() or p.stat().st_size == 0]
if bad:
    raise SystemExit("Missing/empty outputs:\n" + "\n".join(bad))
with open(run / "summary" / "pipeline.done", "w") as out:
    out.write("done\n")
print("validated", len(required), "outputs")
PY

echo "[INFO] finished $(date)"
