cp genome_R108.40000.cool genome_R108.40000.raw.cool

python - <<'PY'
import h5py

cool = "genome_R108.40000.raw.cool"

with h5py.File(cool, "a") as f:
    if "bins/weight" in f:
        del f["bins/weight"]
        print("Removed bins/weight")
    else:
        print("No bins/weight found")
PY
