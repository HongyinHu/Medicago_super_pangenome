import sys
import os
import pandas as pd
import numpy as np

windows_file = os.environ["windows_file"]
meth_file = os.environ["meth_file"]
out_file = os.environ["out_file"]

# windows: chr start end
windows = pd.read_csv(
    windows_file, sep="\t", header=None,
    names=["chr", "start", "end"]
)

meth = pd.read_csv(
    meth_file, sep="\t", header=None,
    names=["chr", "start", "end", "score", "cov"]
)

meth["bin_start"] = (meth["start"] // 50000) * 50000
meth["weighted"] = meth["score"] * meth["cov"]

summary = (
    meth.groupby(["chr", "bin_start"], as_index=False)
        .agg(weighted_sum=("weighted", "sum"),
             cov_sum=("cov", "sum"))
)

summary["value"] = summary["weighted_sum"] / summary["cov_sum"]

summary = summary.rename(columns={"bin_start": "start"})
summary["end"] = summary["start"] + 50000

result = windows.merge(
    summary[["chr", "start", "value"]],
    on=["chr", "start"],
    how="left"
)

result["value"] = result["value"].fillna(0)

result[["chr", "start", "end", "value"]].to_csv(
    out_file, sep="\t", header=False, index=False,
    float_format="%.6f"
)

print("write file to %s"%out_file)
