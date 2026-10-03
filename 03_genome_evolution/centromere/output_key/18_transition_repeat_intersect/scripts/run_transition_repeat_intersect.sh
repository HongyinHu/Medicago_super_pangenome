#!/usr/bin/env bash
set -euo pipefail

ROOT="path/to/project/10.centromere_analysis/output_key"
RUN="${ROOT}/18_transition_repeat_intersect"
SCRIPT="${RUN}/scripts/analyze_transition_repeat_intersect.py"
FA="${ROOT}path/to/data"
TARGET_FA="${RUN}path/to/data.fa"
TRF_DAT="${RUN}path/to/data.fa.2.7.7.80.10.50.2000.dat"

mkdir -p "${RUN}"/{scripts,data,results,figures,logs}

echo "[INFO] Extract transition intervals"
samtools faidx "${FA}" Chr3:22119603-22565328 > "${TARGET_FA}.tmp"
samtools faidx "${FA}" Chr5:25489882-26760313 >> "${TARGET_FA}.tmp"
python - "${TARGET_FA}.tmp" "${TARGET_FA}" <<'PY'
import sys
inp, out = sys.argv[1], sys.argv[2]
name_map = {
    "Chr3:22119603-22565328": "Mpo_Chr3_CEN5_to_CEN6_transition",
    "Chr5:25489882-26760313": "Mpo_Chr5_CEN5_to_CEN6_transition",
}
with open(inp) as ih, open(out, "w") as oh:
    for line in ih:
        if line.startswith(">"):
            raw = line[1:].strip().split()[0]
            oh.write(">" + name_map.get(raw, raw) + "\n")
        else:
            oh.write(line)
PY
rm -f "${TARGET_FA}.tmp"

echo "[INFO] Run targeted TRF if needed"
if [[ ! -s "${TRF_DAT}" ]]; then
  (
    cd "${RUN}/data"
    timeout 900 conda run -n genome_repeat trf "$(basename "${TARGET_FA}")" 2 7 7 80 10 50 2000 -d -h > "${RUN}/logs/trf.stdout" 2> "${RUN}/logs/trf.stderr" || true
  )
fi

echo "[INFO] Summarize repeat intersections"
conda run -n biosofeware python "${SCRIPT}" \
  --run-root "${ROOT}" \
  --outdir "${RUN}/results" \
  --trf-dat "${TRF_DAT}" \
  --random-n 1000

cp "${RUN}/results"/Mpo_Chr5_to_Chr6_transition_repeat_content.* "${RUN}/figures"/
echo "[INFO] Done"
