import zarr
import pandas as pd
import numpy as np

fh = zarr.open_group(zarr.ZipStore(snakemake.input.zarr), "r")
chrom = snakemake.wildcards.chrom

samples = fh[chrom]["samples"][:].astype("<U8").tolist()

df = pd.DataFrame(columns=["ID1", "ID2", "missing"])
df["ID1"] = samples
df["ID2"] = samples
df["missing"] = np.zeros(len(samples), dtype="int")

df.loc[-1] = 0, 0, 0 
df = df.sort_index().reset_index(drop=True)

df.to_csv(
    snakemake.output.txt,
    sep="\t", 
    mode="w",
    index=False, 
    compression="gzip")
