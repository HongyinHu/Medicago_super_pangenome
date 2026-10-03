#!/usr/bin/env bash
set -eo pipefail

RUN_DIR="path/to/project/N_4.pod_spiny/27_R108_SPL_family_phylogeny_20260710"
REF_PEP="path/to/project/N_4.pod_spiny/00_data/7.orthology_geneName/orthology_plant.pep"
R108_PEP="path/to/project/9.T2T_gene_anno_new/output/4.evm_combind/genome_R108/test2_finally/R108_genome.anno.pep"
MT_SOURCE="path/to/project/N_4.pod_spiny/26_Msa_SPL_family_phylogeny_20260710/data/MtSPL.Wang2019.published.pep"
CONDA_SH="path/to/home/anaconda3/etc/profile.d/conda.sh"
TARGET="R10821044"

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

run = Path("path/to/project/N_4.pod_spiny/27_R108_SPL_family_phylogeny_20260710")
ref_pep = Path("path/to/project/N_4.pod_spiny/00_data/7.orthology_geneName/orthology_plant.pep")
r108_pep = Path("path/to/project/9.T2T_gene_anno_new/output/4.evm_combind/genome_R108/test2_finally/R108_genome.anno.pep")
mt_source = Path("path/to/project/N_4.pod_spiny/26_Msa_SPL_family_phylogeny_20260710/data/MtSPL.Wang2019.published.pep")


def read_fasta(path):
    records = []
    name = None
    header = None
    seq = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
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
                out.write(seq[i:i + width] + "\n")


ref_records = []
for rid, _header, seq in read_fasta(ref_pep):
    if not (rid.startswith("AtSPL") or rid.startswith("OsSPL")):
        raise SystemExit(f"Unexpected non-At/Os reference sequence: {rid}")
    ref_records.append((rid, "source=AtOs_reference_SPL", seq.replace("*", "").replace(".", "X")))
write_fasta(ref_records, run / "data" / "AtOs_SPL.reference.clean.pep")

mt_records = []
for rid, _header, seq in read_fasta(mt_source):
    if not rid.startswith("MtSPL"):
        raise SystemExit(f"Unexpected published MtSPL sequence: {rid}")
    mt_records.append((rid, "source=Wang2019_published_MtSPL", seq.replace("*", "").replace(".", "X")))
write_fasta(mt_records, run / "data" / "MtSPL.Wang2019.published.pep")

r108_raw = read_fasta(r108_pep)
gene_records = {}
map_rows = []
for tid, header, seq in r108_raw:
    match = re.search(r"(?:^|\s)gene=([^\s]+)", header)
    gene = match.group(1) if match else re.sub(r"\.\d+$", "", tid)
    clean_seq = seq.replace("*", "").replace(".", "X")
    if gene not in gene_records or len(clean_seq) > len(gene_records[gene][2]):
        gene_records[gene] = (tid, header, clean_seq)

clean_records = []
for gene in sorted(gene_records):
    tid, header, seq = gene_records[gene]
    clean_records.append((gene, f"transcript={tid} original_header={header}", seq))
    map_rows.append((gene, tid, len(seq), header))
write_fasta(clean_records, run / "data" / "R108.proteins.geneid.pep")

with open(run / "data" / "R108.proteins.geneid.map.tsv", "w", encoding="utf-8") as out:
    out.write("gene_id\ttranscript_id\tprotein_length\toriginal_header\n")
    for row in map_rows:
        out.write("\t".join(map(str, row)) + "\n")

with open(run / "logs" / "prepare_fasta.log", "w", encoding="utf-8") as out:
    out.write(f"AtOs reference records\t{len(ref_records)}\n")
    out.write(f"published MtSPL records\t{len(mt_records)}\n")
    out.write(f"R108 raw protein records\t{len(r108_raw)}\n")
    out.write(f"R108 unique gene records\t{len(clean_records)}\n")
    out.write(f"R10821044 present\t{'R10821044' in gene_records}\n")

if "R10821044" not in gene_records:
    raise SystemExit("Target R10821044 is missing from R108 protein annotation")
PY

echo "[INFO] building SPL HMM from At/Os references"
mafft --auto --thread 8 \
  data/AtOs_SPL.reference.clean.pep \
  > data/AtOs_SPL.reference.clean.mafft.fa \
  2> logs/mafft_reference.log
hmmbuild \
  tmp/AtOs_SPL.reference.hmm \
  data/AtOs_SPL.reference.clean.mafft.fa \
  > logs/hmmbuild_AtOs_SPL.log
hmmsearch \
  --cpu 16 \
  --tblout results/R108_vs_AtOs_SPL.hmm.tbl \
  --domtblout results/R108_vs_AtOs_SPL.hmm.domtbl \
  tmp/AtOs_SPL.reference.hmm \
  data/R108.proteins.geneid.pep \
  > logs/hmmsearch_R108_vs_AtOs_SPL.log

echo "[INFO] BLASTP At/Os SPL references against genome_R108 proteins"
makeblastdb \
  -in data/R108.proteins.geneid.pep \
  -dbtype prot \
  -out tmp/R108.proteins.geneid \
  > logs/makeblastdb_R108.log
blastp \
  -query data/AtOs_SPL.reference.clean.pep \
  -db tmp/R108.proteins.geneid \
  -evalue 1e-5 \
  -max_target_seqs 1000 \
  -num_threads 16 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp" \
  -out results/AtOs_SPL_vs_R108.blastp.tsv

echo "[INFO] selecting genome_R108 SPL candidates"
python3 - <<'PY'
from pathlib import Path
import math

run = Path("path/to/project/N_4.pod_spiny/27_R108_SPL_family_phylogeny_20260710")
target_gene = "R10821044"


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
                fields = line[1:].split(maxsplit=1)
                name = fields[0]
                desc = fields[1] if len(fields) > 1 else ""
                seq = []
            else:
                seq.append(line.strip())
        if name is not None:
            records[name] = (desc, "".join(seq))
    return records


def write_fasta(ids, records, path):
    with open(path, "w", encoding="utf-8") as out:
        for gene in ids:
            desc, seq = records[gene]
            out.write(f">{gene} source=genome_R108 {desc}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + "\n")


r108_records = read_fasta(run / "data" / "R108.proteins.geneid.pep")
candidate = {}

with open(run / "results" / "R108_vs_AtOs_SPL.hmm.tbl", "r", encoding="utf-8", errors="replace") as handle:
    for line in handle:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split()
        gene = parts[0]
        evalue = float(parts[4])
        score = float(parts[5])
        if evalue <= 1e-5 and score >= 20:
            candidate.setdefault(gene, {"hmm_evalue": evalue, "hmm_score": score, "blast_hits": []})

with open(run / "results" / "AtOs_SPL_vs_R108.blastp.tsv", "r", encoding="utf-8", errors="replace") as handle:
    for line in handle:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 15:
            continue
        qseqid, gene = parts[0], parts[1]
        pident = float(parts[2])
        aln_len = int(parts[3])
        evalue = float(parts[10])
        bitscore = float(parts[11])
        qlen = int(parts[12])
        slen = int(parts[13])
        qcov = float(parts[14])
        if evalue <= 1e-5 and bitscore >= 40 and aln_len >= 40:
            rec = candidate.setdefault(gene, {"hmm_evalue": math.nan, "hmm_score": math.nan, "blast_hits": []})
            rec["blast_hits"].append((qseqid, evalue, bitscore, pident, aln_len, qcov, qlen, slen))

rows = []
for gene, rec in candidate.items():
    hits = sorted(rec["blast_hits"], key=lambda x: (-x[2], x[1], -x[5]))
    best = hits[0] if hits else ("NA", math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan)
    rows.append({
        "gene_id": gene,
        "protein_length": len(r108_records.get(gene, ("", ""))[1]),
        "hmm_evalue": rec.get("hmm_evalue", math.nan),
        "hmm_score": rec.get("hmm_score", math.nan),
        "blast_hit_count": len(hits),
        "best_reference_SPL": best[0],
        "best_blast_evalue": best[1],
        "best_blast_bitscore": best[2],
        "best_blast_pident": best[3],
        "best_blast_aln_len": best[4],
        "best_blast_qcov": best[5],
        "is_target_R10821044": "yes" if gene == target_gene else "no",
    })

rows.sort(key=lambda row: (
    0 if row["gene_id"] == target_gene else 1,
    str(row["best_reference_SPL"]),
    -float(row["best_blast_bitscore"]) if str(row["best_blast_bitscore"]) != "nan" else 0,
    row["gene_id"],
))
if not rows:
    raise SystemExit("No R108 SPL candidates selected")
if target_gene not in {row["gene_id"] for row in rows}:
    raise SystemExit("Target R10821044 was not selected as an SPL candidate")

header = list(rows[0])
with open(run / "results" / "R108_SPL_candidate_table.tsv", "w", encoding="utf-8") as out:
    out.write("\t".join(header) + "\n")
    for row in rows:
        out.write("\t".join(str(row[field]) for field in header) + "\n")

ids = [row["gene_id"] for row in rows]
missing = [gene for gene in ids if gene not in r108_records]
if missing:
    raise SystemExit(f"Missing R108 protein records: {missing[:5]}")
write_fasta(ids, r108_records, run / "data" / "R108_SPL_candidates.pep")

with open(run / "summary" / "candidate_counts.txt", "w", encoding="utf-8") as out:
    out.write(f"R108_SPL_candidates\t{len(rows)}\n")
    out.write(f"R10821044_detected\t{target_gene in ids}\n")
    out.write("selection\tHMM evalue<=1e-5 score>=20 OR BLASTP evalue<=1e-5 bitscore>=40 aln_len>=40\n")
PY

cat \
  data/AtOs_SPL.reference.clean.pep \
  data/MtSPL.Wang2019.published.pep \
  data/R108_SPL_candidates.pep \
  > data/AtOs_Mt_R108_SPL.combined.pep

if grep -q '^>Msa_' data/AtOs_Mt_R108_SPL.combined.pep; then
  echo "[ERROR] Msa sequence detected in R108-only combined input" >&2
  exit 1
fi

echo "[INFO] aligning At/Os/published-Mt/genome_R108 SPL proteins"
mafft \
  --auto \
  --thread 16 \
  data/AtOs_Mt_R108_SPL.combined.pep \
  > trees/AtOs_Mt_R108_SPL.mafft.fa \
  2> logs/mafft_AtOs_Mt_R108_SPL.log

echo "[INFO] running genome_R108 IQ-TREE2 ML analysis"
iqtree2 \
  -s trees/AtOs_Mt_R108_SPL.mafft.fa \
  -m MFP \
  -bb 1000 \
  -alrt 1000 \
  -T 16 \
  -pre trees/AtOs_Mt_R108_SPL.ML \
  > logs/iqtree2_AtOs_Mt_R108_SPL.log 2>&1

echo "[INFO] plotting genome_R108 trees and summarizing R10821044 placement"
Rscript - <<'RS'
suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
})

run <- "path/to/project/N_4.pod_spiny/27_R108_SPL_family_phylogeny_20260710"
tree_file <- file.path(run, "trees", "AtOs_Mt_R108_SPL.ML.treefile")
tree <- midpoint(read.tree(tree_file))
target <- "R10821044"
tip <- tree$tip.label

if (any(grepl("^Msa_", tip))) {
  stop("Msa sequence detected in final R108-only tree")
}
if (!(target %in% tip)) {
  stop("Target R10821044 is absent from final tree")
}

is_r108 <- grepl("^R108", tip)
is_mt <- grepl("^MtSPL", tip)
is_at <- grepl("^AtSPL", tip)
tip_col <- ifelse(
  tip == target, "#D55E00",
  ifelse(is_r108, "#C00000", ifelse(is_mt, "#7A3E9D", ifelse(is_at, "#1B9E77", "#386CB0")))
)
tip_cex <- ifelse(tip == target, 0.78, ifelse(is_r108, 0.64, 0.50))

plot_rect <- function() {
  par(mar=c(1, 1, 1, 1))
  plot(tree, type="phylogram", cex=tip_cex, font=1, label.offset=0.005,
       tip.color=tip_col, edge.width=0.75, no.margin=TRUE)
  if (!is.null(tree$node.label)) {
    labs <- tree$node.label
    labs[!grepl("/", labs)] <- ""
    nodelabels(labs, frame="none", cex=0.25, adj=c(1.05, -0.15), col="#444444")
  }
  legend(
    "topleft",
    legend=c("genome_R108", "R10821044", "published MtSPL", "Arabidopsis SPL", "Rice SPL"),
    col=c("#C00000", "#D55E00", "#7A3E9D", "#1B9E77", "#386CB0"),
    pch=19, bty="n", cex=0.72
  )
}

pdf(file.path(run, "figures", "AtOs_Mt_R108_SPL_ML_tree_rectangular.pdf"), width=13, height=22, onefile=FALSE)
plot_rect()
dev.off()
png(file.path(run, "figures", "AtOs_Mt_R108_SPL_ML_tree_rectangular.png"), width=4600, height=7200, res=350)
plot_rect()
dev.off()

plot_fan <- function() {
  par(mar=c(0, 0, 0, 0))
  plot(tree, type="fan", cex=0.42, font=1, tip.color=tip_col, edge.width=0.7, no.margin=TRUE)
  legend(
    "bottomleft",
    legend=c("genome_R108", "R10821044", "published MtSPL", "Arabidopsis SPL", "Rice SPL"),
    col=c("#C00000", "#D55E00", "#7A3E9D", "#1B9E77", "#386CB0"),
    pch=19, bty="n", cex=0.72
  )
}

pdf(file.path(run, "figures", "AtOs_Mt_R108_SPL_ML_tree_fan.pdf"), width=16, height=16, onefile=FALSE)
plot_fan()
dev.off()
png(file.path(run, "figures", "AtOs_Mt_R108_SPL_ML_tree_fan.png"), width=5600, height=5600, res=350)
plot_fan()
dev.off()

D <- cophenetic.phylo(tree)
mt_refs <- tip[grepl("^MtSPL", tip)]
at_refs <- tip[grepl("^AtSPL", tip)]
r108_refs <- setdiff(tip[grepl("^R108", tip)], target)

nearest <- function(refs, n=12) {
  values <- sort(D[target, refs])
  values[seq_len(min(n, length(values)))]
}

nearest_mt <- nearest(mt_refs)
nearest_at <- nearest(at_refs)
nearest_r108 <- nearest(r108_refs)

mrca_rows <- data.frame()
for (reference in mt_refs) {
  node <- getMRCA(tree, c(target, reference))
  members <- extract.clade(tree, node)$tip.label
  support <- tree$node.label[node - Ntip(tree)]
  mrca_rows <- rbind(mrca_rows, data.frame(
    reference=reference,
    mrca_tip_count=length(members),
    mrca_r108_count=sum(grepl("^R108", members)),
    mrca_mt_count=sum(grepl("^MtSPL", members)),
    support=support,
    distance=D[target, reference],
    members=paste(members, collapse=","),
    stringsAsFactors=FALSE
  ))
}
mrca_rows <- mrca_rows[order(mrca_rows$mrca_tip_count, mrca_rows$distance), ]

write.table(
  data.frame(reference=names(nearest_mt), patristic_distance=as.numeric(nearest_mt)),
  file=file.path(run, "summary", "R10821044_nearest_published_MtSPLs.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)
write.table(
  data.frame(reference=names(nearest_at), patristic_distance=as.numeric(nearest_at)),
  file=file.path(run, "summary", "R10821044_nearest_AtSPLs.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)
write.table(
  data.frame(r108_gene=names(nearest_r108), patristic_distance=as.numeric(nearest_r108)),
  file=file.path(run, "summary", "R10821044_nearest_R108_SPLs.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)
write.table(
  mrca_rows,
  file=file.path(run, "summary", "R10821044_published_MtSPL_MRCA_clades.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

summary_file <- file.path(run, "summary", "R10821044_naming_summary.txt")
cat("R108-only corrected tree; no genome_Msa proteins included\n\n", file=summary_file)
cat("Nearest published MtSPL proteins\n", file=summary_file, append=TRUE)
capture.output(print(data.frame(reference=names(nearest_mt), distance=as.numeric(nearest_mt))),
               file=summary_file, append=TRUE)
cat("\nNearest Arabidopsis SPL proteins\n", file=summary_file, append=TRUE)
capture.output(print(data.frame(reference=names(nearest_at), distance=as.numeric(nearest_at))),
               file=summary_file, append=TRUE)
cat("\nNearest genome_R108 SPL candidates\n", file=summary_file, append=TRUE)
capture.output(print(data.frame(reference=names(nearest_r108), distance=as.numeric(nearest_r108))),
               file=summary_file, append=TRUE)
cat("\nSmallest MRCA clades containing R10821044 and published MtSPL proteins\n",
    file=summary_file, append=TRUE)
capture.output(print(head(mrca_rows, 12)), file=summary_file, append=TRUE)
RS

echo "[INFO] validating genome_R108-only outputs"
python3 - <<'PY'
from pathlib import Path

run = Path("path/to/project/N_4.pod_spiny/27_R108_SPL_family_phylogeny_20260710")
required = [
    run / "results" / "R108_SPL_candidate_table.tsv",
    run / "data" / "R108_SPL_candidates.pep",
    run / "data" / "AtOs_Mt_R108_SPL.combined.pep",
    run / "trees" / "AtOs_Mt_R108_SPL.ML.treefile",
    run / "trees" / "AtOs_Mt_R108_SPL.ML.contree",
    run / "figures" / "AtOs_Mt_R108_SPL_ML_tree_rectangular.pdf",
    run / "figures" / "AtOs_Mt_R108_SPL_ML_tree_rectangular.png",
    run / "figures" / "AtOs_Mt_R108_SPL_ML_tree_fan.pdf",
    run / "figures" / "AtOs_Mt_R108_SPL_ML_tree_fan.png",
    run / "summary" / "R10821044_naming_summary.txt",
]
bad = [str(path) for path in required if not path.exists() or path.stat().st_size == 0]
if bad:
    raise SystemExit("Missing or empty outputs:\n" + "\n".join(bad))

for path in [run / "data" / "AtOs_Mt_R108_SPL.combined.pep", run / "trees" / "AtOs_Mt_R108_SPL.ML.treefile"]:
    text = path.read_text(encoding="utf-8", errors="replace")
    if "Msa_" in text or "genome_Msa" in text:
        raise SystemExit(f"Forbidden genome_Msa sequence or label detected in {path}")

with open(run / "summary" / "no_genome_Msa_validation.txt", "w", encoding="utf-8") as out:
    out.write("combined_fasta_has_genome_Msa\tfalse\n")
    out.write("final_tree_has_genome_Msa\tfalse\n")

with open(run / "summary" / "pipeline.done", "w", encoding="utf-8") as out:
    out.write("done\n")
print("validated", len(required), "genome_R108-only outputs")
PY

echo "[INFO] finished $(date)"
