import zarr
import pandas as pd

fh = zarr.open_group(zarr.ZipStore(snakemake.input.zarr), "r")

chrom = snakemake.wildcards.chrom
chunksize = 100000

genotypes = fh[chrom]["calldata/genotype"]

positions = fh[chrom]["variants/POS"]
reference = fh[chrom]["variants/REF"][:].astype("<U8")
alternate = fh[chrom]["variants/ALT"][:].astype("<U8")

nvar, nsam, ploidy = genotypes.shape

samples = [str(i) for i in range(nsam)]

assert nvar == positions.shape[0] == reference.shape[0]

mode = "w"
for start in range(0, nvar, chunksize):
    stop = start + chunksize
    df = pd.DataFrame(index=positions[start:stop], columns=["chrom", "id", "pos", "ref", "alt"])
    df["pos"] = positions[start:stop]
    df["ref"] = reference[start:stop]
    df["alt"] = alternate[start:stop]
    df["id"] = "."
    df["chrom"] = chrom

    for i, sid in enumerate(samples):
        df[sid+"A"] = genotypes[start:stop, i, 0]
        df[sid+"B"] = genotypes[start:stop, i, 1]

    df.to_csv(
        snakemake.output.txt,
        sep=" ", 
        mode=mode,
        index=False, 
        compression="gzip",
        header=False)

    mode = "a"
