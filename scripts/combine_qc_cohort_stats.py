import pandas as pd
import numpy as np
import sys
import argparse


## functions to compute mean, median, mode from histogram counts
# ignore 0s and 1s.
# compute mean
def mean_binc(x):
    x = np.array(x)[2:]
    return np.average(np.arange(0, x.shape[0], dtype="int") + 2, weights=x)


# compute median
def median_binc(x):
    x = np.array(x)[2:]
    return np.searchsorted(np.cumsum(x), np.sum(x) // 2) + 2


# compute mode
def mode_binc(x):
    x = np.array(x)[2:]
    return np.argmax(x) + 2


def log(*msg):
    print(*msg, file=sys.stderr)
    sys.stderr.flush()


def main():

    parser = argparse.ArgumentParser(
        description='Given a set of manifest samples and expected filepaths merge all data and apply defined filters')

    parser.add_argument(
        '--manifest', required=True, help='path to csv containing sample IDs in the first column')

    parser.add_argument(
        '--input-path', required=False,
        help='Location of csv files including a {sample} key. The file endings ".contamination.csv", ".covstats.csv" '
             'and "covhist.npz" are assumed"', default="{sample}", type=str)

    parser.add_argument(
        '--output', required=True, help='path to output location for table')

    try:
        args = {
            "manifest": snakemake.input.manifest,
            "output": snakemake.output.csv,
            "input_path": snakemake.params.path
        }
        log("Args read via snakemake")
    except NameError:
        args = vars(parser.parse_args())
        log("Args read via command line")

    with open(args["manifest"], mode="r") as sin:
        sids = [x.strip() for x in sin.readlines()]

    contam_suffix = ".contamination.csv"
    align_suffix = ".callstats.csv"
    cov_suffix = ".covhist.npz"

    df = pd.concat(
        {sid: pd.read_csv(args["input_path"].format(sample=sid) + align_suffix, index_col=0) for sid in sids})

    coverage_dict = {}
    for sid in sids:
        v = np.load(args["input_path"].format(sample=sid) + cov_suffix)
        for kv, y in v.items():
            coverage_dict[sid, kv] = pd.Series(y, dtype=int)

    cov_frame = pd.concat(coverage_dict, axis=1).T

    # this groups by individual. (and sums over chromosomes).
    cov_grouped = cov_frame.groupby(level=0)
    coverage_all_sum = cov_grouped.agg(sum)

    qc_frame = pd.DataFrame(index=sids)
    qc_frame.index.name = "derived_sample_id"

    qc_frame["mean_cov"] = coverage_all_sum.apply(mean_binc, axis=1).round(2)
    qc_frame["median_cov"] = coverage_all_sum.apply(median_binc, axis=1)
    qc_frame["modal_cov"] = coverage_all_sum.apply(mode_binc, axis=1)

    # report per chromosome mean/medians.
    cov_grouped_by_chrom = cov_frame.groupby(level=1)
    for (_chrom, _gdf) in cov_grouped_by_chrom:
        _gdf.index = _gdf.index.droplevel(level=1)
        qc_frame["mean_cov_" + _chrom] = _gdf.apply(mean_binc, axis=1).round(2) 
        qc_frame["median_cov_" + _chrom] = _gdf.apply(median_binc, axis=1)
        qc_frame["mode_cov_" + _chrom] = _gdf.apply(mode_binc, axis=1)

    qc_frame['frac_gen_cov'] = (coverage_all_sum.T[1:].sum(axis=0) / coverage_all_sum.sum(axis=1)).round(3)

    allele_calls_summed = df.groupby(level=0).agg(sum)

    qc_frame['divergence'] = (allele_calls_summed["nNonRefAlleles"] / (2 * allele_calls_summed["nSitesCalled"])).round(5)

    contam_df = pd.concat(
        [pd.read_csv(args["input_path"].format(sample=sid) + contam_suffix, index_col=0) for sid in sids])

    output_frame = pd.concat(
        [qc_frame, contam_df[['pc_contam', 'LLR']].round(3)], axis=1)

    output_frame.to_csv(args["output"], sep="\t")


if __name__== "__main__":
    main()

