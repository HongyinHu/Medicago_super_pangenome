python - <<'PY'
import cooler
import numpy as np

cool = "genome_R108.100000.cool"
chrom = "Chr1"

c = cooler.Cooler(cool)
mat = c.matrix(balance=False).fetch(chrom)

vals = mat[np.isfinite(mat) & (mat > 0)]

print("nonzero pixels:", vals.size)
for q in [50, 80, 90, 95, 98, 99, 99.5]:
    v = np.percentile(vals, q)
    print(q, v, "log1p =", np.log1p(v))
PY
