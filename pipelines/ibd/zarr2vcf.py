## Nick Harding
## Snakemake python script.

import allel
import pandas as pd
import numpy as np

import os
import zarr

def write_vcf(path, chrom, fh, sample_idx=None, chunk_size=100000):
    
    # write header.
    VCF_FIXED_FIELDS = 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'
    nrows = fh[chrom]["variants/POS"].shape[0]
    samples = tuple(fh["samples"][:].astype("<U8").tolist())
    inc_samples = tuple(np.take(samples, sample_idx).tolist())
    
    assert len(samples) == fh[chrom]["calldata/genotype"].shape[1]
    assert nrows == fh[chrom]["calldata/genotype"].shape[0]
    
    with open(path, mode="w") as gz:
        
        f = '##fileformat=VCFv4.1\n'
        gz.write(f)
        f = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        gz.write(f)

        h = "#" + "\t".join(VCF_FIXED_FIELDS + inc_samples)
        gz.write(h + "\n")
        
        chunks = np.arange(0, nrows + chunk_size, chunk_size)
        assert chunks.max() > nrows

        for start, stop in zip(chunks[:-1], chunks[1:]):
            sl = slice(start, stop)
            positions = fh[chrom]['variants/POS'][sl]
            reference = fh[chrom]['variants/REF'][sl].astype("<U1")
            alternate = fh[chrom]['variants/ALT'][sl].astype("<U1")
            genotypes = fh[chrom]['calldata/genotype'][sl].astype("<U2").take(sample_idx, axis=1)
            
            multiple_alts = alternate.ndim > 1

            for pos, ref, alt, gt in zip(positions, reference, alternate, genotypes):
                filterstring = 'PASS'

                # alt may be an np array, with several entries.
                if multiple_alts:
                    alt = ",".join(x for x in alt if x != '')

                try:
                    gstr = ["/".join(x) for x in gt]
                    genotype_str = "\t".join(gstr)
                    genotype_str = genotype_str.replace("-1/-1", "./.")

                    line = "\t".join(
                        [chrom, str(pos)] + ['.', ref, alt, '0', '.', '.', 'GT'] + [genotype_str]) + "\n"

                    gz.write(line)

                except TypeError:
                    print(pos)
                    print(ref)
                    print(alt)
                    print(gt)
                    raise TypeError("Bad data type.")
                    
meta = pd.read_table(snakemake.input.meta)
meta.index.name = "idx"
ix = meta.reset_index().set_index("ox_code").loc[snakemake.params.samples].idx.values

zh = zarr.open_group(snakemake.input.zarr, "r")
write_vcf(snakemake.output.vcf, snakemake.wildcards.chrom, zh, sample_idx=ix)
