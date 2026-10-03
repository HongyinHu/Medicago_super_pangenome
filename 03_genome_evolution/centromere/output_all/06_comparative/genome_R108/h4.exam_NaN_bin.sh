python - <<'PY'
import cooler
import pandas as pd

cool = "genome_R108.100000.cool"
chrom = "Chr1"

c = cooler.Cooler(cool)
bins = c.bins().fetch(chrom)

bad = bins[bins["weight"].isna()]
print("total bins:", len(bins))
print("bad bins:", len(bad))
print(bad[["chrom", "start", "end", "weight"]].head(50))
PY
