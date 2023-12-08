import pandas as pd
import os, subprocess, sys

l = [
    pd.read_csv(e, compression="gzip", sep="\t")
    for j, e in enumerate(list(snakemake.input))
]
pd.concat(l).to_csv(
    snakemake.output[0], sep="\t", compression="gzip", mode="w", index=False
)
