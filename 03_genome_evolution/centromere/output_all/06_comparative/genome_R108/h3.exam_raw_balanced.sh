python - <<'PY'
import cooler
import numpy as np

cool = "genome_R108.100000.cool"
chrom = "Chr1"   # 改成你实际画图的染色体名，比如 A01 或 Chr1

c = cooler.Cooler(cool)

for balance in [False, True]:
    try:
        mat = c.matrix(balance=balance).fetch(chrom)
    except Exception as e:
        print("balance =", balance, "failed:", e)
        continue

    vals = mat[np.isfinite(mat) & (mat > 0)]
    print("\nBalance =", balance)
    print("nonzero pixels:", vals.size)
    for q in [50, 80, 90, 95, 98, 99, 99.5]:
        v = np.percentile(vals, q)
        print(q, v, "log1p =", np.log1p(v))
PY