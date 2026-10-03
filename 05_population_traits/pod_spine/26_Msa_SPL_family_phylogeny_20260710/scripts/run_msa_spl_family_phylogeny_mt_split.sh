#!/usr/bin/env bash
set -eo pipefail

RUN_DIR="path/to/project/N_4.pod_spiny/26_Msa_SPL_family_phylogeny_20260710"
CONDA_SH="path/to/home/anaconda3/etc/profile.d/conda.sh"
MT_SUPP_URL="https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-019-5937-1/MediaObjects/12864_2019_5937_MOESM1_ESM.docx"
MT_SUPP_DOCX="${RUN_DIR}path/to/data.docx"

mkdir -p "${RUN_DIR}"/{data,scripts,logs,results,trees,figures,summary,tmp}
cd "${RUN_DIR}"

source "${CONDA_SH}"
conda activate biosofeware
set -u

echo "[INFO] started $(date)"
echo "[INFO] downloading the published Wang et al. 2019 MtSPL supplementary sequences"
curl -L -sS --retry 3 --max-time 180 "${MT_SUPP_URL}" -o "${MT_SUPP_DOCX}"
unzip -t "${MT_SUPP_DOCX}" > logs/Wang2019_MtSPL_Additional_file1.unzip_test.log

echo "[INFO] extracting published MtSPL proteins and splitting Chr23998"
python3 - <<'PY'
from collections import OrderedDict
from pathlib import Path
import html
import re
import zipfile
import xml.etree.ElementTree as ET

run = Path("path/to/project/N_4.pod_spiny/26_Msa_SPL_family_phylogeny_20260710")
docx = run / "data" / "Wang2019_MtSPL_Additional_file1.docx"


def read_fasta(path):
    records = OrderedDict()
    name = None
    desc = ""
    seq = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
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
                seq.append(re.sub(r"\s+", "", line))
        if name is not None:
            records[name] = (desc, "".join(seq))
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


with zipfile.ZipFile(docx) as archive:
    xml = archive.read("word/document.xml")
root = ET.fromstring(xml)

lines = []
for paragraph in root.iter("{http://schemas.openxmlformats.org/wordprocessingml/2006/main}p"):
    chunks = []
    for node in paragraph.iter():
        if node.tag.endswith("}t") and node.text:
            chunks.append(node.text)
        elif node.tag.endswith("}tab"):
            chunks.append("\t")
        elif node.tag.endswith("}br"):
            chunks.append("\n")
    text = html.unescape("".join(chunks)).strip()
    if text:
        lines.extend(x.strip() for x in text.splitlines() if x.strip())

published = OrderedDict()
current = None
for line in lines:
    if line.startswith(">"):
        current = line[1:].split()[0]
        if current.startswith("MtSPL"):
            published[current] = []
        continue
    if current in published:
        published[current].append(re.sub(r"[^A-Za-z*]", "", line).upper())

published = OrderedDict((name, "".join(parts).replace("*", "")) for name, parts in published.items())
if len(published) != 23:
    raise SystemExit(f"Expected 23 published MtSPL proteins, found {len(published)}: {list(published)}")
if any(not seq for seq in published.values()):
    raise SystemExit("One or more published MtSPL sequences are empty")

write_fasta(
    [(name, "source=Wang_et_al_2019_Additional_file1", seq) for name, seq in published.items()],
    run / "data" / "MtSPL.Wang2019.published.pep",
)
with open(run / "results" / "MtSPL.Wang2019.sequence_manifest.tsv", "w", encoding="utf-8") as out:
    out.write("published_name\tprotein_length\tsource\n")
    for name, seq in published.items():
        out.write(f"{name}\t{len(seq)}\tWang2019_Additional_file1\n")

msa = read_fasta(run / "data" / "Msa_SPL_candidates.pep")
target = "Msa_Chr23998"
if target not in msa:
    raise SystemExit(f"Missing {target} in Msa_SPL_candidates.pep")
desc, fused = msa.pop(target)
if len(fused) != 746:
    raise SystemExit(f"Expected Chr23998 protein length 746, observed {len(fused)}")
if fused[340] != "M":
    raise SystemExit(f"Expected second protein to start with M at aa341, observed {fused[340]}")

split_a = fused[:340]
split_b = fused[340:]
corrected = []
for name, (record_desc, seq) in msa.items():
    corrected.append((name, record_desc, seq))
corrected.extend([
    (
        "Msa_Chr23998a",
        "source=genome_Msa split_from=Chr23998 aa=1-340 genomic_block=Chr4:88989328-88991530(-) best_R108=R10821046",
        split_a,
    ),
    (
        "Msa_Chr23998b",
        "source=genome_Msa split_from=Chr23998 aa=341-746 genomic_block=Chr4:88984906-88987161(-) best_R108=R10821045",
        split_b,
    ),
])
write_fasta(corrected, run / "data" / "Msa_SPL_candidates.Chr23998_split.pep")

with open(run / "results" / "Chr23998_split_manifest.tsv", "w", encoding="utf-8") as out:
    out.write("new_id\toriginal_id\tprotein_interval\tprotein_length\tgenomic_block\tbest_R108\tevidence\n")
    out.write("Msa_Chr23998a\tMsa_Chr23998\t1-340\t340\tChr4:88989328-88991530(-)\tR10821046\tfirst_HMM_domain_and_BLAST_block\n")
    out.write("Msa_Chr23998b\tMsa_Chr23998\t341-746\t406\tChr4:88984906-88987161(-)\tR10821045\tsecond_HMM_domain_starts_with_M_and_BLAST_block\n")

with open(run / "summary" / "Chr23998_split_validation.txt", "w", encoding="utf-8") as out:
    out.write(f"original_length\t{len(fused)}\n")
    out.write(f"split_a_length\t{len(split_a)}\n")
    out.write(f"split_b_length\t{len(split_b)}\n")
    out.write(f"split_a_last20\t{split_a[-20:]}\n")
    out.write(f"split_b_first20\t{split_b[:20]}\n")
    out.write("HMM_domain_1\t5-303; envelope 1-338\n")
    out.write("HMM_domain_2\t341-705; envelope 341-739\n")
    out.write("R108_correspondence\tChr23998a->R10821046; Chr23998b->R10821045\n")
PY

cat \
  data/AtOs_SPL.reference.clean.pep \
  data/MtSPL.Wang2019.published.pep \
  data/Msa_SPL_candidates.Chr23998_split.pep \
  > data/AtOs_Mt_Msa_SPL.Chr23998_split.combined.pep

echo "[INFO] aligning At/Os/published-Mt/corrected-Msa SPL proteins"
mafft --auto --thread 16 \
  data/AtOs_Mt_Msa_SPL.Chr23998_split.combined.pep \
  > trees/AtOs_Mt_Msa_SPL.Chr23998_split.mafft.fa \
  2> logs/mafft_AtOs_Mt_Msa_SPL_Chr23998_split.log

echo "[INFO] running corrected IQ-TREE2 ML analysis"
iqtree2 \
  -s trees/AtOs_Mt_Msa_SPL.Chr23998_split.mafft.fa \
  -m MFP \
  -bb 1000 \
  -alrt 1000 \
  -T 16 \
  -pre trees/AtOs_Mt_Msa_SPL.Chr23998_split.ML \
  > logs/iqtree2_AtOs_Mt_Msa_SPL_Chr23998_split.log 2>&1

echo "[INFO] plotting corrected trees and summarizing Chr23997 placement"
Rscript - <<'RS'
suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
})

run <- "path/to/project/N_4.pod_spiny/26_Msa_SPL_family_phylogeny_20260710"
tree_file <- file.path(run, "trees", "AtOs_Mt_Msa_SPL.Chr23998_split.ML.treefile")
tree <- midpoint(read.tree(tree_file))
tip <- tree$tip.label

is_msa <- grepl("^Msa_", tip)
is_mt <- grepl("^MtSPL", tip)
is_at <- grepl("^AtSPL", tip)
is_os <- grepl("^OsSPL", tip)
tip_col <- ifelse(
  tip == "Msa_Chr23997", "#D55E00",
  ifelse(is_msa, "#C00000", ifelse(is_mt, "#7A3E9D", ifelse(is_at, "#1B9E77", "#386CB0")))
)
tip_cex <- ifelse(tip == "Msa_Chr23997", 0.78, ifelse(grepl("^Msa_Chr23998", tip), 0.66, 0.50))

plot_rect <- function() {
  par(mar=c(1,1,1,1))
  plot(tree, type="phylogram", cex=tip_cex, font=1, label.offset=0.005,
       tip.color=tip_col, edge.width=0.75, no.margin=TRUE)
  if (!is.null(tree$node.label)) {
    labs <- tree$node.label
    labs[!grepl("/", labs)] <- ""
    nodelabels(labs, frame="none", cex=0.25, adj=c(1.05,-0.15), col="#444444")
  }
  legend(
    "topleft",
    legend=c("genome_Msa", "Chr23997", "published MtSPL", "Arabidopsis SPL", "Rice SPL"),
    col=c("#C00000", "#D55E00", "#7A3E9D", "#1B9E77", "#386CB0"),
    pch=19, bty="n", cex=0.72
  )
}

pdf(file.path(run, "figures", "AtOs_Mt_Msa_SPL_Chr23998_split_ML_tree_rectangular.pdf"), width=13, height=22, onefile=FALSE)
plot_rect()
dev.off()
png(file.path(run, "figures", "AtOs_Mt_Msa_SPL_Chr23998_split_ML_tree_rectangular.png"), width=4600, height=7200, res=350)
plot_rect()
dev.off()

plot_fan <- function() {
  par(mar=c(0,0,0,0))
  plot(tree, type="fan", cex=0.42, font=1, tip.color=tip_col, edge.width=0.7, no.margin=TRUE)
  legend(
    "bottomleft",
    legend=c("genome_Msa", "Chr23997", "published MtSPL", "Arabidopsis SPL", "Rice SPL"),
    col=c("#C00000", "#D55E00", "#7A3E9D", "#1B9E77", "#386CB0"),
    pch=19, bty="n", cex=0.72
  )
}

pdf(file.path(run, "figures", "AtOs_Mt_Msa_SPL_Chr23998_split_ML_tree_fan.pdf"), width=16, height=16, onefile=FALSE)
plot_fan()
dev.off()
png(file.path(run, "figures", "AtOs_Mt_Msa_SPL_Chr23998_split_ML_tree_fan.png"), width=5600, height=5600, res=350)
plot_fan()
dev.off()

D <- cophenetic.phylo(tree)
target <- "Msa_Chr23997"
mt_refs <- tip[grepl("^MtSPL", tip)]
at_refs <- tip[grepl("^AtSPL", tip)]
msa_refs <- setdiff(tip[grepl("^Msa_", tip)], target)

nearest <- function(refs, n=12) {
  x <- sort(D[target, refs])
  x[seq_len(min(n, length(x)))]
}

nearest_mt <- nearest(mt_refs)
nearest_at <- nearest(at_refs)
nearest_msa <- nearest(msa_refs)

write.table(
  data.frame(reference=names(nearest_mt), patristic_distance=as.numeric(nearest_mt)),
  file=file.path(run, "summary", "Chr23997_nearest_published_MtSPLs.after_split.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)
write.table(
  data.frame(reference=names(nearest_at), patristic_distance=as.numeric(nearest_at)),
  file=file.path(run, "summary", "Chr23997_nearest_AtSPLs.after_split.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)
write.table(
  data.frame(msa_gene=names(nearest_msa), patristic_distance=as.numeric(nearest_msa)),
  file=file.path(run, "summary", "Chr23997_nearest_Msa_SPLs.after_split.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

clade_rows <- data.frame()
for (r in mt_refs) {
  node <- getMRCA(tree, c(target, r))
  members <- extract.clade(tree, node)$tip.label
  clade_rows <- rbind(clade_rows, data.frame(
    reference=r,
    mrca_tip_count=length(members),
    mrca_msa_count=sum(grepl("^Msa_", members)),
    mrca_mt_count=sum(grepl("^MtSPL", members)),
    distance=D[target, r],
    members=paste(members, collapse=","),
    stringsAsFactors=FALSE
  ))
}
clade_rows <- clade_rows[order(clade_rows$mrca_tip_count, clade_rows$distance),]
write.table(
  clade_rows,
  file=file.path(run, "summary", "Chr23997_published_MtSPL_MRCA_clades.after_split.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

summary_file <- file.path(run, "summary", "Chr23997_naming_summary.with_published_MtSPL_and_split_Chr23998.txt")
cat("Corrected tree: published MtSPL proteins added; Chr23998 split at aa340/341\n\n", file=summary_file)
cat("Nearest published MtSPL proteins\n", file=summary_file, append=TRUE)
capture.output(print(head(data.frame(reference=names(nearest_mt), distance=as.numeric(nearest_mt)), 10)), file=summary_file, append=TRUE)
cat("\nNearest Arabidopsis SPL proteins\n", file=summary_file, append=TRUE)
capture.output(print(head(data.frame(reference=names(nearest_at), distance=as.numeric(nearest_at)), 8)), file=summary_file, append=TRUE)
cat("\nNearest Msa SPL candidates\n", file=summary_file, append=TRUE)
capture.output(print(head(data.frame(reference=names(nearest_msa), distance=as.numeric(nearest_msa)), 8)), file=summary_file, append=TRUE)
cat("\nSmallest MRCA clades containing Chr23997 and published MtSPL proteins\n", file=summary_file, append=TRUE)
capture.output(print(head(clade_rows[, c("reference", "mrca_tip_count", "mrca_msa_count", "mrca_mt_count", "distance")], 12)), file=summary_file, append=TRUE)
RS

echo "[INFO] validating corrected outputs"
python3 - <<'PY'
from pathlib import Path

run = Path("path/to/project/N_4.pod_spiny/26_Msa_SPL_family_phylogeny_20260710")
required = [
    run / "data" / "MtSPL.Wang2019.published.pep",
    run / "data" / "Msa_SPL_candidates.Chr23998_split.pep",
    run / "results" / "MtSPL.Wang2019.sequence_manifest.tsv",
    run / "results" / "Chr23998_split_manifest.tsv",
    run / "trees" / "AtOs_Mt_Msa_SPL.Chr23998_split.ML.treefile",
    run / "figures" / "AtOs_Mt_Msa_SPL_Chr23998_split_ML_tree_rectangular.pdf",
    run / "figures" / "AtOs_Mt_Msa_SPL_Chr23998_split_ML_tree_rectangular.png",
    run / "figures" / "AtOs_Mt_Msa_SPL_Chr23998_split_ML_tree_fan.pdf",
    run / "figures" / "AtOs_Mt_Msa_SPL_Chr23998_split_ML_tree_fan.png",
    run / "summary" / "Chr23997_naming_summary.with_published_MtSPL_and_split_Chr23998.txt",
]
bad = [str(path) for path in required if not path.exists() or path.stat().st_size == 0]
if bad:
    raise SystemExit("Missing or empty corrected outputs:\n" + "\n".join(bad))

manifest = (run / "results" / "MtSPL.Wang2019.sequence_manifest.tsv").read_text().strip().splitlines()
if len(manifest) != 24:
    raise SystemExit(f"Expected header plus 23 MtSPL records, observed {len(manifest)} lines")

(run / "summary" / "pipeline_mt_split.done").write_text("done\n")
print("validated", len(required), "corrected outputs")
PY

echo "[INFO] finished $(date)"
