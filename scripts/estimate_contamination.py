import zarr
import allel
import pandas as pd
import numpy as np
import sys
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import product
from scipy.optimize import minimize_scalar


# This defines the probabilities of observing alleles REF/ALT/ERR,
# when there is no sequence error in the read, and when there is error. (inner dim)
# for each genotype HOMREF, HET, HOMALT. (outer dim)
# Only biallelic alleles are supported.
ALLELE_OBSERVATION_PROBABILITIES = [
    (np.array((1, 0, 0)), np.array((0, 1/3, 2/3))),
    (np.array((0.5, 0.5, 0)), np.array((1/6, 1/6, 2/3))),
    (np.array((0, 1, 0)), np.array((1/3, 0, 2/3)))
]


# scikit-allele does not yet have a locate intersection that works on multi-indexes
# this is a bit of a hacky way to solve this problem, based on the known indices of the different seq ids
# for each seqid in turn find the intersection and replace in the master intersection data
def locate_intersection(positions_a, lengths_a, positions_b, lengths_b):

    log("Computing position overlap")
    loc_a = np.zeros(positions_a.shape, dtype=bool)
    loc_b = np.zeros(positions_b.shape, dtype=bool)

    ix_b, ix_a = 0, 0

    for va, vb in zip(lengths_a, lengths_b):

        positions_given_seq_a = allel.SortedIndex(positions_a[ix_a:(ix_a + va)])
        positions_given_seq_b = allel.SortedIndex(positions_b[ix_b:(ix_b + vb)])

        temp_loc_a, temp_loc_b = positions_given_seq_a.locate_intersection(positions_given_seq_b)
        loc_a[ix_a:(ix_a + va)] = temp_loc_a
        loc_b[ix_b:(ix_b + vb)] = temp_loc_b
        ix_a += va
        ix_b += vb

    return loc_a, loc_b


# single function: requires: input handle holding seq ids, the seq ids to concat, paths to arrays
# idea of this is to make it a general function, then sort the complexity with multi-indexes later.
def concatenate_arrays(callset, seqids, paths):

    # ensure seq ids is a list
    if isinstance(seqids, str):
        seqids = [seqids]
    else:
        assert isinstance(seqids, list), "`seq_ids` must be a string, or a list of strings"

    assert isinstance(paths, list) and len(paths) > 0, "Paths must be a list of at least one string"

    path_holder_dict = {p: np.concatenate([callset[seq][p] for seq in seqids]) for p in paths}

    shapes = []
    for seq in seqids:
        shape = [callset[seq][p].shape[0] for p in paths]
        assert len(set(shape)) == 1, "Inconsistent axis size: {0} for {1}".format(shape, seq)
        shapes.append(shape[0])

    return path_holder_dict, shapes


# There are 9 possible genotypes of the sample and putative contaminant
# Using fixed alpha and error rates we can calculate a 9 x 3 matrix of
# the probabilities of observing each of REF/ALT/ERR
def compute_probability_matrix(alpha, error_rate, definitions):
    assert isinstance(alpha, float)
    assert isinstance(error_rate, float)
    assert isinstance(definitions, list)

    prob_matrix = np.zeros((9, 3))

    for i, j in product(range(len(definitions)), repeat=2):
        no_error, error = definitions[i]
        gt1 = ((1 - error_rate) * no_error) + (error_rate * error)
        gt1 *= (1 - alpha)

        no_error, error = definitions[j]
        gt2 = ((1 - error_rate) * no_error) + (error_rate * error)
        gt2 *= alpha

        prob_matrix[(i * 3) + j] = gt1 + gt2

    return prob_matrix


# Given the 9 possible genotypes we have some prior expectation of the probabilities of
# those genotypes. This needs to be computed once, and given the genotypes of the sample
# and contaminant are independent, we use the product of the expected genotypes.
# Expected genotypes are calculated based on HW proportions
def determine_weights(afv):

    afv = np.array(afv, dtype=float)
    assert np.all(afv < 1.0) & np.all(afv > 0.0), \
        "frequencies must be withing 0.0 > < 1.0"

    genotype_freqs = np.stack(
        [afv ** 2, 2 * (afv * (1 - afv)), (1 - afv) ** 2],
        axis=1)

    # step 3 create the GT1/GT2 object.
    return np.matmul(
        genotype_freqs[:, :, np.newaxis],
        genotype_freqs[:, np.newaxis, :]).reshape((-1, 9)).T


# This is the main workhorse function that is optimized.
# Given an alpha, we compute the likelihood of the observed data.
def compute_likelihood(alpha, seq_error, allele_depths, gt_weights, logging=False, logstem=None):

    # compute the probability matrix given the function above
    prob_mat = compute_probability_matrix(alpha, seq_error, ALLELE_OBSERVATION_PROBABILITIES)

    assert np.all(gt_weights < 0.0), "Weights must be log space"

    # use some log maths and some numpy indexing to take the prob mat to the power of the allele depths
    log_master_matrix = np.log(prob_mat[:, None, :]) * allele_depths[None, :, :]

    # log_gts_over_sites takes the product across observed alleles, so sum in log space.
    log_gts_over_sites = log_master_matrix.sum(axis=2)

    # un-log so we can sum over possible GTs.
    # At this point we have marginalized out potential genotypes
    vector_over_sites = np.exp(gt_weights + log_gts_over_sites).sum(axis=0)

    # in some places, this will return a 0.0, which is fine.
    # problems only occur when *all* genotypes very unlikely.
    # This can happen if there is a 3rd allele at high frequency.
    # The only way the model has of accounting for this is
    # as a sequence error. Therefore we have the filter for 'probably biallelic'
    try:
        assert np.all(vector_over_sites < 1.0), "All values must be below 1.0"
        assert np.all(vector_over_sites > 0.0), "All values must be above 0.0"
    except AssertionError:
        log(vector_over_sites.min(), vector_over_sites.max())
        mini, maxi = vector_over_sites.argmin(), vector_over_sites.argmax()
        log(mini, allele_depths[mini])
        log(maxi, allele_depths[maxi])
        raise AssertionError("bad values")

    # This is where we make a log file when we want to make very detailed logs,
    # ie expaining why some alphas are better than others at explaining the data
    if logging:
        log_frame = pd.DataFrame(
            index=range(allele_depths.shape[0]),
            columns=["REFcount", "ALTcount", "ERRcount", "LL"])

        log_frame["REFcount"] = allele_depths[:, 0]
        log_frame["ALTcount"] = allele_depths[:, 1]
        log_frame["ERRcount"] = allele_depths[:, 2]
        log_frame["LL"] = np.log(vector_over_sites)

        fn = logstem.format(alpha=alpha)
        log_frame.to_csv(fn)

    # return the negative, as the function is _minimized_
    return - np.log(vector_over_sites).sum()


def log(*msg):
    print(*msg, file=sys.stderr)
    sys.stderr.flush()


def plot_allele_balance(genotypes, allele_depths, path, annotation):

    is_het = np.squeeze(genotypes.is_het())
    het_calls = genotypes.compress(is_het, axis=0)[:, 0]
    assert het_calls.shape[1] == 2
    het_ad = np.compress(is_het, allele_depths, axis=0)[:, 0]
    assert het_ad.shape[1] == 4

    altref = np.array([(hc[ref], hc[alt]) for hc, (ref, alt) in zip(het_ad, het_calls)], dtype=int)
    assert np.all(altref > 0), "ALL cells in het alt/ref should be above zero."

    fig, ax = plt.subplots(figsize=(8, 8))

    # compute site frequency spectrum
    mac1, mac2 = altref[:, 0], altref[:, 1]
    m = 250
    n = 250
    tmp = (mac1 * n + mac2).astype(int, copy=False)
    s = np.bincount(tmp)
    s.resize(m, n)

    # commented line below allows log of frequencies.
    # needs s += 1
    cax = ax.imshow(
        s[:80, :80],
        cmap="Reds")
    fig.colorbar(cax, ax=ax, orientation='horizontal', fraction=0.08, extend='max')
    sns.despine(ax=ax)

    ax.set_xlabel("REF count")
    ax.set_ylabel("ALT count")
    ax.set_xlim((0, 80))
    ax.set_ylim((0, 80))

    ax.text(65, 75, r"%$\alpha$ : {0:.2f}".format(annotation["pc_contam"]), color='k')
    ax.text(65, 70, r"$\Lambda$ : {0:.1f}".format(annotation["LLR"]), color='k')

    fig.savefig(path, bbox_inches="tight")


def main():

    parser = argparse.ArgumentParser(
        description='This script takes a single zarr zipstore, and estimates the contamination rate, providing a log'
        'likelihood ratio vs the null model')

    parser.add_argument(
        '--input', required=True,
        help='Path to zarr file containing genotypes and allele depths, zipped Zarr file with data for a single sample.'
             'This should follow the standard format of {sample}/{seqid}/calldata/GT and {sample}/{seqid}/calldata/AD.')

    parser.add_argument(
        '--sites', required=True,
        help='Path to zarr describing which sites in `input` were genotyped. This is used to match the `input` to the'
             'allele frequencies below. variants/POS is required.')

    parser.add_argument(
        '--allele-frequencies', required=True,
        help='path to zarr file describing allele frequencies. This has two purposes: 1) to select SNPs to downsample'
             'to, based on the `minimum_af` argument. 2) To provide a prior expectation on the frequency of genotypes.'
             'The first level of the zarr file should be groups for seqids, with each containing `POS` (position) and'
             'AF (allele frequencies). The shape of the AF array must be Nx4, where N is the size of the 1D POS array.'
             'The order of alleles *must* correspond to the coding in the input data. There is no requirement to have a'
             'similar shape to the input genotypes, although a minimum level of intersection is required!')

    parser.add_argument(
        '--seqid', required=True, nargs='+', help='name of chromosome(s) or contig(s) to process. ')

    parser.add_argument(
        '--output', required=True, help='path to output file stem')

    parser.add_argument(
        '--downsample', required=False, default=20000, help='number of sites to consider.', type=int)

    parser.add_argument(
        '--minimum-af', required=False, default=0.05,
        help='minimum minor allele frequency in reference population to consider. Sites with higher MAF are more '
             'powerful at detecting contamination', type=float)

    parser.add_argument(
        '--sequence-error-rate', required=False, default=1e-3, help='probability of observing a non REF/ALT base',
        type=float)

    parser.add_argument(
        '--minimum-coverage', required=False, default=10,
        help='minimum read depth to use. Low depths have low power to detect contamination', type=int)

    parser.add_argument('--plot', dest='plot', action='store_true')
    parser.add_argument('--no-plot', dest='plot', action='store_false')

    parser.add_argument('--log', dest='log', action='store_true')
    parser.add_argument('--no-log', dest='log', action='store_false')

    parser.set_defaults(plot=True, log=False)

    try:
        args = {
            "input": snakemake.input.input,
            "sites": snakemake.input.sites,
            "allele_frequencies": snakemake.input.allele_frequencies,
            "seqid": snakemake.params.seqid,
            "output": snakemake.params.stem,
            "minimum_af": snakemake.params.minimum_af,
            "minimum_coverage": snakemake.params.minimum_coverage,
            "sequence_error_rate": snakemake.params.seq_err_rate,
            "downsample": snakemake.params.downsample,
            "plot": snakemake.params.plot,
            "log": snakemake.params.log
        }
        log("Args read via snakemake")
    except NameError:
        args = vars(parser.parse_args())
        log("Args read via command line")

    seqids = args['seqid']
    sequence_error_rate = args['sequence_error_rate']
    downsample_n = args["downsample"]
    minimum_minor_af = args["minimum_af"]

    output_csv = args['output'] + ".contamination.csv"
    output_png = args['output'] + ".allele_balance.png"
    output_log = args["output"] + ".{alpha}.log"

    sample_store = zarr.ZipStore(args["input"], mode="r")
    sample_callset = zarr.Group(sample_store)

    sites = zarr.ZipStore(args["sites"], mode="r")
    variant_sites = zarr.Group(sites)
    sample = next(iter(sample_callset))

    concatenated_sample_callset, _ = concatenate_arrays(
        sample_callset[sample],
        seqids,
        paths=["calldata/GT", "calldata/AD"])

    gt = allel.GenotypeArray(concatenated_sample_callset["calldata/GT"])
    ad = concatenated_sample_callset["calldata/AD"]

    concatenated_sites, concatenated_site_shapes = concatenate_arrays(
        variant_sites,
        seqids,
        ["variants/POS"])
    pos = concatenated_sites["variants/POS"]
    assert pos.shape[0] == gt.shape[0] == ad.shape[0], "Shape inconsistency. {0}, {1}, {2}".format(
        pos.shape, gt.shape, ad.shape)

    # load allele frequencies required to compute weights
    allele_frequencies_z = zarr.open_group(args['allele_frequencies'], "r")
    concatenated_af_arrays, concatenated_af_shapes = concatenate_arrays(
        allele_frequencies_z,
        seqids,
        ["POS", "AF"])
    af_pos = concatenated_af_arrays["POS"]
    # This is a 2D array of the frequency of the ALT allele in some other dataset.
    af_val = concatenated_af_arrays["AF"]
    assert af_val.shape[1] == 4, "Allele frequencies must contain all 4 alleles, even if unobserved."

    # for the sample_gt: Keep if
    # a) in af, b) is_called and c) is_biallelic
    # step 1 find the intersection this works on multi indexes
    loc_gt, loc_af = locate_intersection(pos, concatenated_site_shapes, af_pos, concatenated_af_shapes)

    flt_af_val = np.compress(loc_af, af_val, axis=0)
    flt_gt = np.compress(loc_gt, gt, axis=0)
    flt_ad = np.compress(loc_gt, ad, axis=0)

    # now we need to filter both by is biallelic and is called.
    is_bial_ref_pop = np.count_nonzero(flt_af_val, axis=1) == 2
    is_called = flt_gt.is_called()[:, 0]

    # compress the intersection by the AND of these
    keep_loc = is_called & is_bial_ref_pop
    alt_frequency_pass = np.compress(keep_loc, flt_af_val, axis=0)
    allele_depth_pass = np.compress(keep_loc, flt_ad, axis=0)

    # recode the allele depth to 0/1.
    # find the "alt" column.
    log("Ordering alleles by frequency for REF/ALT/ERR")
    min_cov_reached = allele_depth_pass[:, 0].sum(axis=1) >= args['minimum_coverage']

    major_allele_indexes = np.argsort(alt_frequency_pass, axis=1)[:, ::-1]
    allele_depth_pass_reordered = np.array([np.take(q, ix) for q, ix in zip(allele_depth_pass, major_allele_indexes)])

    # Define allele counts: sum final 2 columns, representing ref/alt/error
    allele_depths = allele_depth_pass_reordered[:, :3]
    allele_depths[:, 2] = allele_depths[:, 2] + allele_depth_pass_reordered[:, 3]
    assert allele_depths.shape[1] == 3

    # issue with some samples having a third allele (ie not in phase 2) discovered at high frequency
    # Filter sites where more than 10% of reads look like errors.
    probably_biallelic = allele_depth_pass_reordered[:, 2] < (.1 * allele_depth_pass_reordered.sum(axis=1))

    # step 2 create the 0/1/2 from the allele frequencies.
    major_af = alt_frequency_pass.max(axis=1)

    # select the values with the highest MAF.
    log("Selecting variants on which to perform analysis")
    while True:
        eligible = probably_biallelic & min_cov_reached & ((1 - major_af) > minimum_minor_af)
        if eligible.sum() > downsample_n:
            break

        minimum_minor_af -= 0.01
        if minimum_minor_af < 0:
            log("Insufficient variants meet criteria to compute contamination. n={0}, min={1}".format(
                eligible.sum(), downsample_n))
            break

    res = pd.DataFrame(index=[sample], columns=["LLR", "LL", "pc_contam"])

    if eligible.sum() > downsample_n:
        log("Downsample from {0} to {1}".format(eligible.sum(), downsample_n))

        ix_ds = np.sort(np.random.choice(
            np.where(eligible)[0], size=downsample_n))

        major_af = np.take(major_af, ix_ds, axis=0)
        allele_depths = np.take(allele_depths, ix_ds, axis=0)

        genotype_weights = np.log(determine_weights(major_af))

        log("estimating contamination...")
        xv = minimize_scalar(
            compute_likelihood,
            args=(sequence_error_rate, allele_depths, genotype_weights, args["log"], output_log),
            bounds=(0, 0.5),
            method="Bounded",
            options={"xatol": 1e-6})

        # compute the likelihood at alpha = 0, to report likelihood ratio.
        null = compute_likelihood(0.0, sequence_error_rate, allele_depths, genotype_weights, args["log"], output_log)

        # return the llr / ll / estimate
        res.loc[sample] = -min(xv.fun - null, 0), -xv.fun, xv.x * 100

        if args['plot']:
            plot_allele_balance(flt_gt, flt_ad, output_png, res.iloc[0])

    res.to_csv(output_csv)



if __name__ == '__main__':
    main()
