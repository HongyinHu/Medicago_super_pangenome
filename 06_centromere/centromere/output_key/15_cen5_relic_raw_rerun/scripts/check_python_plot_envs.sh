#!/usr/bin/env bash
set -euo pipefail

CONDA=path/to/home/anaconda3/bin/conda

for env_name in biosofeware chromap_env pygenometracks HIC_AB_TAD CENH3_env base; do
  echo "ENV:${env_name}"
  if "${CONDA}" run -n "${env_name}" python -c "import matplotlib, numpy; import pandas; print('matplotlib', matplotlib.__version__); print('numpy', numpy.__version__); print('pandas', pandas.__version__)"
  then
    exit 0
  fi
done

exit 1
