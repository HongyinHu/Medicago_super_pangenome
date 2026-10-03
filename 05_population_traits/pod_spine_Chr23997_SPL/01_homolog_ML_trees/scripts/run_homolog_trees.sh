#!/usr/bin/env bash
set -eo pipefail

ROOT="path/to/project/41.Chr23997_realName"
QUERY="$ROOT/00_data/7.orthology_geneName/SPL_orthology.pep.fa"
MSA="$ROOT/00_data/genome_Msa.pep"
R108="$ROOT/00_data/genome_R108.pep"
OUT="$ROOT/01_homolog_ML_trees"
THREADS=32
TOP_N=3

source path/to/home/anaconda3/etc/profile.d/conda.sh
conda activate biosofeware
set -u
mkdir -p "$OUT"/{01_blast/db,02_fasta,03_alignment,04_iqtree,logs,scripts}

{
  echo -e "item\tvalue"
  echo -e "query_fasta\t$QUERY"
  echo -e "msa_proteome\t$(readlink -f "$MSA")"
  echo -e "r108_proteome\t$(readlink -f "$R108")"
  echo -e "query_records\t$(grep -c '^>' "$QUERY")"
  echo -e "msa_records\t$(grep -c '^>' "$MSA")"
  echo -e "r108_records\t$(grep -c '^>' "$R108")"
  echo -e "blastp\t$(blastp -version | head -1)"
  echo -e "mafft\t$(mafft --version 2>&1 | head -1)"
  echo -e "iqtree2\t$(iqtree2 --version | head -1)"
  echo -e "homolog_filter\tevalue<=1e-10; qcovs>=50; pident>=25; top_${TOP_N}_per_query_per_genome"
} > "$OUT/logs/input_summary.tsv"

makeblastdb -in "$MSA" -dbtype prot -parse_seqids -out "$OUT/01_blast/db/genome_Msa" >"$OUT/logs/makeblastdb_Msa.log" 2>&1
makeblastdb -in "$R108" -dbtype prot -parse_seqids -out "$OUT/01_blast/db/genome_R108" >"$OUT/logs/makeblastdb_R108.log" 2>&1

FMT='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs'
blastp -query "$QUERY" -db "$OUT/01_blast/db/genome_Msa" -evalue 1e-10 -max_target_seqs 50 -num_threads "$THREADS" -outfmt "$FMT" > "$OUT/01_blast/SPL_vs_Msa.tsv"
blastp -query "$QUERY" -db "$OUT/01_blast/db/genome_R108" -evalue 1e-10 -max_target_seqs 50 -num_threads "$THREADS" -outfmt "$FMT" > "$OUT/01_blast/SPL_vs_R108.tsv"

OUT="$OUT" QUERY="$QUERY" MSA="$MSA" R108="$R108" TOP_N="$TOP_N" python - <<'PY'
import os
from pathlib import Path
from collections import defaultdict

out = Path(os.environ["OUT"])
query_path = Path(os.environ["QUERY"])
msa_path = Path(os.environ["MSA"])
r108_path = Path(os.environ["R108"])
top_n = int(os.environ["TOP_N"])
columns = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qlen","slen","qcovs"]

def parse_fasta(path):
    records = {}
    ident = None
    header = None
    seq = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if ident is not None:
                    records[ident] = (header, "".join(seq).replace("*", ""))
                header = line[1:]
                ident = header.split()[0]
                if ident in records:
                    raise RuntimeError(f"Duplicate FASTA identifier in {path}: {ident}")
                seq = []
            else:
                seq.append(line.replace(" ", ""))
        if ident is not None:
            records[ident] = (header, "".join(seq).replace("*", ""))
    if not records:
        raise RuntimeError(f"No records in {path}")
    return records

queries = parse_fasta(query_path)
db_records = {"Msa": parse_fasta(msa_path), "R108": parse_fasta(r108_path)}
selected = {}
all_selected_rows = []
summary = {q: {"Msa": [], "R108": []} for q in queries}

for source, table in [("Msa", out/"01_blast/SPL_vs_Msa.tsv"), ("R108", out/"01_blast/SPL_vs_R108.tsv")]:
    grouped = defaultdict(list)
    with open(table) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(columns):
                continue
            row = dict(zip(columns, fields))
            grouped[row["qseqid"]].append(row)
    selected[source] = {}
    for qid in queries:
        rows = sorted(grouped.get(qid, []), key=lambda x: (-float(x["bitscore"]), float(x["evalue"]), -float(x["qcovs"])))
        passing = [r for r in rows if float(r["evalue"]) <= 1e-10 and float(r["qcovs"]) >= 50 and float(r["pident"]) >= 25]
        keep = passing[:top_n] if passing else rows[:1]
        selected[source][qid] = keep
        summary[qid][source] = keep
        for rank, row in enumerate(keep, 1):
            row = dict(row)
            row["genome"] = source
            row["rank"] = str(rank)
            row["selection"] = "threshold_pass" if row in passing else "best_available_fallback"
            all_selected_rows.append(row)

with open(out/"01_blast/homolog_hits.tsv", "w") as fh:
    fh.write("\t".join(["genome","rank","selection"] + columns) + "\n")
    for row in all_selected_rows:
        fh.write("\t".join(row[x] for x in ["genome","rank","selection"] + columns) + "\n")

def write_tree_fasta(name, sources):
    seen = set()
    members = []
    path = out/"02_fasta"/f"{name}.faa"
    with open(path, "w") as fh:
        for qid, (header, seq) in queries.items():
            label = f"TARGET|{qid}"
            fh.write(f">{label} source=target original_header={header}\n{seq}\n")
            seen.add(label)
            members.append((name, "target", label, qid))
        for source in sources:
            for qid in queries:
                for row in selected[source][qid]:
                    sid_raw = row["sseqid"]
                    candidates = [sid_raw] + [x for x in sid_raw.split("|") if x]
                    sid = next((x for x in candidates if x in db_records[source]), None)
                    if sid is None:
                        raise RuntimeError(f"{source} hit not found in proteome: {sid_raw}")
                    label = f"{source}|{sid}"
                    if label in seen:
                        continue
                    header, seq = db_records[source][sid]
                    fh.write(f">{label} source={source} matched_query={qid} original_header={header}\n{seq}\n")
                    seen.add(label)
                    members.append((name, source, label, sid))
    return members

members = []
members += write_tree_fasta("Msa_background", ["Msa"])
members += write_tree_fasta("R108_background", ["R108"])
members += write_tree_fasta("combined_Msa_R108", ["Msa", "R108"])

with open(out/"02_fasta/members.tsv", "w") as fh:
    fh.write("tree_input\tsource\tsequence_id\toriginal_id\n")
    for r in members:
        fh.write("\t".join(r) + "\n")

with open(out/"summary.tsv", "w") as fh:
    fh.write("query_id\tMsa_homologs_selected\tR108_homologs_selected\n")
    for qid in queries:
        fh.write(f"{qid}\t{len(summary[qid]['Msa'])}\t{len(summary[qid]['R108'])}\n")
PY

for fasta in "$OUT"/02_fasta/*.faa; do
  name=$(basename "$fasta" .faa)
  mafft --auto --thread "$THREADS" "$fasta" > "$OUT/03_alignment/$name.aln.fa" 2> "$OUT/logs/$name.mafft.log"
done

: > "$OUT/logs/tree_status.tsv"
pids=()
for aln in "$OUT"/03_alignment/*.aln.fa; do
  name=$(basename "$aln" .aln.fa)
  nseq=$(grep -c '^>' "$aln")
  if [[ "$nseq" -lt 4 ]]; then
    echo -e "$name\tSKIPPED\t$nseq sequences (<4)" >> "$OUT/logs/tree_status.tsv"
    continue
  fi
  (
    iqtree2 -s "$aln" -m MFP -B 1000 -alrt 1000 -nt "$THREADS" -pre "$OUT/04_iqtree/$name" -redo
    echo -e "$name\tCOMPLETED\t$nseq sequences" >> "$OUT/logs/tree_status.tsv"
  ) > "$OUT/logs/$name.iqtree.stdout.log" 2>&1 &
  pids+=("$!")
done
status=0
for pid in "${pids[@]}"; do
  wait "$pid" || status=1
done
if [[ "$status" -ne 0 ]]; then
  echo "ERROR: one or more IQ-TREE runs failed; inspect $OUT/logs/*.iqtree.stdout.log" >&2
  exit 1
fi

{
  echo "# SPL homolog and ML tree analysis"
  echo
  echo "Queries: 00_data/7.orthology_geneName/SPL_orthology.pep.fa."
  echo "For each query and each genome proteome, up to the three highest-scoring BLASTP hits passing evalue <= 1e-10, query coverage >= 50%, and identity >= 25% were retained."
  echo "Tree inputs contain all 58 query proteins plus the non-redundant selected homologs for Msa, R108, or both."
  echo "MAFFT --auto was used for alignment; IQ-TREE 2 used ModelFinder Plus, 1,000 ultrafast bootstraps, and 1,000 SH-aLRT replicates."
  echo
  echo "Final consensus trees: 04_iqtree/*.contree."
} > "$OUT/README.md"

touch "$OUT/COMPLETED"
