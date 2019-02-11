import pandas as pd
import numpy as np
import sys
import argparse


# mean_median, mode
def mean_binc(x):
    x = np.array(x)
    return np.average(np.arange(0, x.shape[0], dtype="int"), weights=x)


# mean_median, mode
def median_binc(x):
    x = np.array(x)
    return np.searchsorted(np.cumsum(x), np.sum(x) // 2)


# mean_median, mode
def mode_binc(x):
    x = np.array(x)
    return np.argmax(x[1:]) + 1


def make_sex_call(ratio, female_min=0.8, female_max=1.2, male_max=0.6, male_min=0.4):
    if np.isnan(ratio):
        return "NA"
    elif female_max >= ratio >= female_min:
        return "F"
    elif male_min <= ratio <= male_max:
        return "M"
    else:
        return "UNK"


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

    parser.add_argument('--min-frac-genome-covered', default=0.5, required=False, type=float)
    parser.add_argument('--min-median-coverage', default=10, required=False, type=int)
    parser.add_argument('--min-female-xratio', default=0.8, required=False, type=float)
    parser.add_argument('--max-female-xratio', default=1.1, required=False, type=float)
    parser.add_argument('--min-male-xratio', default=0.4, required=False, type=float)
    parser.add_argument('--max-male-xratio', default=0.6, required=False, type=float)
    parser.add_argument('--max-pc-contamination', default=4.5, required=False, type=float)

    try:
        args = {
            "manifest": snakemake.input.manifest,
            "output": snakemake.output.csv,
            "input_path": snakemake.params.path,
            "min_frac_genome_covered": snakemake.params.min_fgc,
            "min_median_coverage": snakemake.params.min_med_cov,
            "min_female_xratio":snakemake.params.min_female_xratio,
            "max_male_xratio":snakemake.params.max_male_xratio,
            "max_female_xratio":snakemake.params.max_female_xratio,
            "min_male_xratio":snakemake.params.min_male_xratio,
            "max_pc_contamination": snakemake.params.max_pc_contam
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

    coverage_all_sum = cov_frame.groupby(level=0).agg(sum)

    qc_frame = pd.DataFrame(index=sids)
    qc_frame.index.name = "derived_sample_id"

    qc_frame["mean_cov"] = coverage_all_sum.apply(mean_binc, axis=1).round(2)
    qc_frame['median_cov'] = coverage_all_sum.apply(median_binc, axis=1)

    # mode sum not in output
    coverage_all_sum.apply(mode_binc, axis=1)

    qc_frame['coverage_ratio_mean'] = (cov_frame.xs('X', level=1).apply(mean_binc, axis=1) / cov_frame.xs('3L', level=1).apply(mean_binc, axis=1)).round(3)
    qc_frame['coverage_ratio_mode'] = (cov_frame.xs('X', level=1).apply(mode_binc, axis=1) / cov_frame.xs('3L', level=1).apply(mode_binc, axis=1)).round(3)

    qc_frame['sex_call'] = qc_frame['coverage_ratio_mode'].apply(
        make_sex_call, args=(args['min_female_xratio'], args['max_female_xratio'] , args['min_male_xratio'], args['max_male_xratio']))

    qc_frame['frac_gen_cov'] = (coverage_all_sum.T[1:].sum(axis=0) / coverage_all_sum.sum(axis=1)).round(3)

    allele_calls_summed = df.groupby(level=0).agg(sum)

    qc_frame['divergence'] = (allele_calls_summed["nNonRefAlleles"] / (2 * allele_calls_summed["nSitesCalled"])).round(3)

    contam_df = pd.concat(
        [pd.read_csv(args["input_path"].format(sample=sid) + contam_suffix, index_col=0) for sid in sids])

    output_frame = pd.concat(
        [qc_frame, contam_df[['pc_contam', 'LLR']].round(3)], axis=1)

    output_frame["FILTER_frac_genome_cov"] = output_frame["frac_gen_cov"] >= args["min_frac_genome_covered"]
    output_frame["FILTER_median_cov"] = output_frame["median_cov"] >= args["min_median_coverage"]
    output_frame["FILTER_contamination"] = output_frame["pc_contam"] <= args["max_pc_contamination"]
    output_frame["FILTER_nosexcall"] = output_frame["sex_call"].apply(lambda y: y in ["M", "F"])

    output_frame.to_csv(args["output"])


if __name__== "__main__":
    main()

