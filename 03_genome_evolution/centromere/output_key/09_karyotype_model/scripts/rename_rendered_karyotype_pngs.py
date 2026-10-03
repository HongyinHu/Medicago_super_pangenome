#!/usr/bin/env python3
from pathlib import Path
import shutil

src = Path("path/to/project/10.centromere_analysis/output_key/09_karyotype_model/rendered")
dst = Path("path/to/project/10.centromere_analysis/output_key/09_karyotype_model/rendered_ascii")
dst.mkdir(exist_ok=True)
for i, p in enumerate(sorted(src.glob("*.png")), 1):
    out = dst / f"karyotype_rendered_{i}.png"
    shutil.copy2(p, out)
    print(f"{p.name}\t{out.name}")
