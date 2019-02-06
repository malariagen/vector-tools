__author__ = "Nicholas Harding"
# This script takes a single zarr zipstore, and computes summary stats
#  of the alignment  from the AD and GT fields of calldata/
import zarr
import allel
import pandas as pd
import sys
import argparse
import numpy as np


def log(*msg):
    print(*msg, file=sys.stderr)
    sys.stderr.flush()


def main():

    parser = argparse.ArgumentParser(
        description='Given a zipped Zarr store summarize information about the coverage/calls')

    parser.add_argument(
        '--input', required=True, help='path to Zarr store containing calldata for a single sample')

    parser.add_argument(
        '--seqid', required=False, help='name of chromosomes or contigs to process',
        nargs="+", default=["2R", "2L", "3R", "3L", "X"], type=str)

    parser.add_argument(
        '--output', required=True, help='path to output basename for stats table and coverage histogram')

    try:
        args = {
            "input": snakemake.input.input,
            "seqid": snakemake.params.seqid,
            "output": snakemake.params.stem
        }
        log("Args read via snakemake")
    except NameError:
        args = vars(parser.parse_args())
        log("Args read via command line")

    zfn = args['input']
    store = zarr.ZipStore(zfn, mode="r")
    callset = zarr.Group(store)

    csv_out = args['output'] + ".callstats.csv"
    npy_out = args['output'] + ".covhist.npz"

    # Holders for data
    df_cols = ["nSitesCalled", "nHomRef", "nHet", "nHetRef", "nHomAlt", "nNonRefAlleles"]
    ad_df = pd.DataFrame(
        index=args['seqid'],
        columns=df_cols)

    cov_hist = dict()

    sample = next(iter(callset))

    for chrom in args['seqid']:

        gt = allel.GenotypeArray(
            callset[sample][chrom]["calldata/GT"])

        is_called = np.squeeze(gt.is_called())

        # if no coverage at all
        if is_called.sum() == 0:
            ad_df.loc[chrom] = 0
            cov_hist[chrom] = np.array([is_called.shape[0]])
            continue

        all_allele_depths = callset[sample][chrom]["calldata/AD"][:, 0]

        # To calculate DP at each position/sample we need to filter missing genotpyes.
        allele_depths_with_cov = np.compress(is_called, all_allele_depths, axis=0)

        # Histogram of coverage
        sum_by_alt = allele_depths_with_cov.sum(axis=1)
        bc = np.bincount(sum_by_alt, minlength=251)
        bc[0] = (~is_called).sum()
        cov_hist[chrom] = bc

        nHomRef = gt.count_hom_ref()
        nHet = gt.count_het()
        nHetRef = gt.count_het(allele=0)
        nHomAlt = gt.count_hom_alt()

        nNonRefAlleles = (2 * nHomAlt) + nHetRef + (2 * (nHet - nHetRef))
        total_calls = is_called.sum()

        ad_df.loc[chrom] = [total_calls, nHomRef, nHet, nHetRef, nHomAlt, nNonRefAlleles]

    ad_df.to_csv(csv_out)
    np.savez(npy_out, **cov_hist)


if __name__ == '__main__':
    main()

