#!/usr/bin/env bash
set -euo pipefail

BASE=path/to/project/38.medicago_resequence
OUT=${BASE}/2.call_SNP_new
DATA_DIR=${BASE}/1.resequenc_data/ncbi_medicago_submit
SAMPLE_DIR=${OUT}/01_sample_list
LOG_DIR=${SAMPLE_DIR}/logs

mkdir -p "${SAMPLE_DIR}" "${LOG_DIR}"

python3 - "${DATA_DIR}" "${SAMPLE_DIR}" > "${LOG_DIR}/make_sample_list.log" 2>&1 <<'PY'
import collections
import pathlib
import re
import sys

data_dir = pathlib.Path(sys.argv[1])
out_dir = pathlib.Path(sys.argv[2])
pattern = re.compile(r"(.+)_([12])\.f(?:ast)?q\.gz$")

pairs = collections.defaultdict(dict)
ignored = []
for path in sorted(data_dir.glob("*.fq.gz")):
    match = pattern.match(path.name)
    if not match:
        ignored.append(str(path))
        continue
    sample, mate = match.group(1), match.group(2)
    pairs[sample][mate] = str(path.resolve())

sample_info = out_dir / "sample_info.tsv"
sample_ids = out_dir / "sample_ids.txt"
missing = out_dir / "missing_pairs.tsv"

bad = []
with sample_info.open("w", encoding="utf-8") as handle, sample_ids.open("w", encoding="utf-8") as ids:
    handle.write("sample\tfq1\tfq2\n")
    for sample in sorted(pairs):
        mates = pairs[sample]
        if "1" not in mates or "2" not in mates:
            bad.append((sample, mates.get("1", ""), mates.get("2", "")))
            continue
        handle.write(f"{sample}\t{mates['1']}\t{mates['2']}\n")
        ids.write(f"{sample}\n")

with missing.open("w", encoding="utf-8") as handle:
    handle.write("sample\tfq1\tfq2\n")
    for sample, fq1, fq2 in bad:
        handle.write(f"{sample}\t{fq1}\t{fq2}\n")

print(f"fastq_files={sum(len(v) for v in pairs.values())}")
print(f"complete_samples={sum(1 for v in pairs.values() if '1' in v and '2' in v)}")
print(f"incomplete_samples={len(bad)}")
print(f"ignored_files={len(ignored)}")
if ignored:
    print("ignored:")
    for item in ignored:
        print(item)
if bad:
    raise SystemExit("Incomplete FASTQ pairs found; see missing_pairs.tsv")
PY
