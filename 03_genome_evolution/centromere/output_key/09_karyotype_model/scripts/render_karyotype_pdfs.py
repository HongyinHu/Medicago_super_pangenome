#!/usr/bin/env python3
"""Render local karyotype/chromosome-painting PDFs to PNG previews."""

from __future__ import annotations

import subprocess
from pathlib import Path


def main() -> None:
    root = Path("path/to/project/10.centromere_analysis/output_key/09_karyotype_model")
    pdf_dir = root / "painting_pdf"
    out_dir = root / "rendered"
    out_dir.mkdir(parents=True, exist_ok=True)
    for pdf in sorted(pdf_dir.glob("*.pdf")):
        out_prefix = out_dir / pdf.stem
        subprocess.run(["pdftoppm", "-png", "-r", "220", "-singlefile", str(pdf), str(out_prefix)], check=True)
        print(out_prefix.with_suffix(".png"))


if __name__ == "__main__":
    main()
