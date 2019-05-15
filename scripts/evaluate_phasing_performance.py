
__author__ = 'Nicholas Harding'
__version__ = 'v1.3.0'

# This script takes a list of samples,
# two sets of phased data and evaluates one against the other

# to produce:
# - A summary table of number of errors in each window.

import argparse
import h5py
import numpy as np
import pandas as pd
import allel
import zarr
import os
import pyfasta
import re
from intervaltree import IntervalTree

def read_roh_h5(path, names, chrom):
    roh_fh = h5py.File(path, "r")
    roh = dict()
    for s in names:
        try:
            roh[s] = roh_fh[chrom][s][:]
        except KeyError:
            roh[s] = None
    return roh

def read_roh_zarr(path, names, chrom):
    roh_fh = zarr.open_group(path, "r")
    roh = dict()
    for s in names:
        try:
            roh[s] = np.column_stack(
                [roh_fh[chrom][s]["start"][:], roh_fh[chrom][s]["stop"][:]])
        except KeyError:
            roh[s] = None
    return roh

def haldane(r):
    return -np.log(1-(2*r))/2

def check_equivalency(a, b):
    assert a.ndim == b.ndim == 2
    assert np.array_equal(a.shape, b.shape)
    assert np.array_equal(a.sum(axis=1), b.sum(axis=1))

def compute_allele_agreement(hap, geno):

    a = hap == np.take(geno, 0, axis=1)
    b = hap == np.take(geno, 1, axis=1)
    # a = 0, b = 1, a & b == err, not a | b = 2.

    out = 4 - ((a + 1) * (b + 2))
    # 2 = unkn, 1 = 2nd allele, 0 = 1st allele, -2 = Fail- both!
    assert not np.any(-2 == out), "Error: some homozygotes must be present"
    return out

def determine_haplotype_switch(geno_a, geno_b):
    """
    Operates on two alternative haplotype assertions on two identical genotypes.
    Returns haplotype switch array, as integers to avoid confusion
    The information is encoded in runs of 0s and 1s. a 0 indicates that the
    0th allele from the first genotype array matches the 0th allele from the
    second. A 1 indicates that the 0th allele from the first matches the
    first from the second. As the order of the 3rd axis is arbitrary the
    distinction between 0/1s is unimportant, we may as well work from the 1st
    axis and invert the result. Will error if any inconsistencies.

    :return: np.array() of length N
    """

    # check for equivalence
    check_equivalency(geno_a, geno_b)
    return compute_allele_agreement(np.take(geno_a, 0, axis=1), geno_b)


def find_reference_gaps(genome_fa):

    # fix on a gapsize of 1kb, as ld breaks down over this distance and
    # no chance of read based phasing
    # do chunks of 10k at a time. look for consecutive Ns.
    size = 10000
    gap_size = 1000
    match = "N{{{n}}}?".format(n=gap_size)

    gaps = list()

    for i in range(0, len(genome_fa), size - (2 * gap_size)):
        text = genome_fa[i:i+size]

        for m in re.finditer(match, text, flags=re.IGNORECASE):
            gaps.append((i + m.start(), i + m.end()))

    tree = IntervalTree.from_tuples(gaps)
    tree.merge_overlaps()

    return np.array([[iv.begin, iv.end] for iv in sorted(tree)])


def calc_marker_dist(pos, homozygosity=None):

    if pos.size == 0:
        return 0
    else:
        start, stop = pos[0], pos[-1]

    if homozygosity is None:
        sum_roh = 0
    else:
        int_tree = IntervalTree(IntervalTree.from_tuples(homozygosity))
        int_tree.slice(start)
        int_tree.slice(stop)
        sum_roh = np.sum([ival.length() for ival in int_tree[start:stop]])

    return stop - start - sum_roh


def evaluate_markers(markers, error_positions):
    """
    function that retuns all marker distances, and whether they are a SE.
    start: is the genomic window start pos
    markers: an arr of genetic positions of the markers
    error_pos: an arr of error positions
    window_size: window size
    """

    if markers.size == 0:
        return np.empty(0, dtype="bool")

    # don't care if first marker is an error as window based ie pertains to
    # gap between marker and immediately before
    errors = np.in1d(markers[1:], error_positions)

    return errors

def derive_position_switch_array(alleles):
    # The position switch array describes the runs of allele ids.
    # Any non inherited alleles are ignored...just masked.
    # consistent_alleles = alleles != 2
    # a = np.compress(consistent_alleles, alleles)
    # if consistent_alleles.any():
    #     print "{0} non consistent site(s), either mutation or genotype " \
    #           "error. These are ignored for purposes of switches.".format(
    #           np.invert(consistent_alleles).sum())
    #
    m = alleles[:-1] != alleles[1:]
    co = np.where(np.concatenate(([False], m, [True])))[0]
    ans = np.diff(np.insert(co, 0, 0))
    return ans


def calculate_switch_distances(windows, switch_array, marker_pos, rohz, gaps):

    marker_pos = allel.SortedIndex(marker_pos)
    gap_mp = np.mean(gaps, axis=1)

    marker_count = np.zeros(windows.shape[0], dtype="int")
    marker_dist = np.zeros(windows.shape[0], dtype="float")
    error_count = np.zeros(windows.shape[0], dtype="int")

    pos_sw = derive_position_switch_array(switch_array)
    pos_errors = np.take(marker_pos, pos_sw[:-1].cumsum())

    for i, (start, stop) in enumerate(windows):

        # this is the code I need to change
        # A don't count error if immediately after GAP
        # B don't count towards distance
        try:
            ix = marker_pos.locate_range(start, stop)

        except KeyError:
            marker_dist[i] = 0.0
            marker_count[i] = 0
            error_count[i] = 0
            continue

        # how many separate gaps between first and last ix?
        gap_ix = np.searchsorted(marker_pos[ix], gap_mp)

        # interested in number of gaps
        gap_pos = np.unique(
            np.compress((gap_ix < marker_pos[ix].size) & (gap_ix > 0), gap_ix))

        # now insert 0 and pos size at beginning and end
        cuts = np.concatenate([[0], gap_pos, [marker_pos[ix].size]])
        assert cuts.size >= 2

        for p, q in zip(cuts[:-1], cuts[1:]):

            first, last = marker_pos[ix][p], marker_pos[ix][q-1]

            error_count[i] += np.sum(evaluate_markers(marker_pos[ix][p:q], pos_errors))

            marker_dist[i] += calc_marker_dist(marker_pos[ix][p:q], rohz)
            # just one marker is not informative.
            marker_count[i] += (q - p - 1)

    return np.vstack([marker_dist, marker_count, error_count])

def open_group(path):
    if path.endswith("h5"):
        return h5py.File(path, "r")
    elif path.endswith(("zarr2", "zarr")):
        return zarr.open_group(path, "r")
    else:
        raise ValueError("Bad filepath provided: {0}. Only hdf5/zarr supported.".format(path))

parser = argparse.ArgumentParser(
    description='Evaluation pipeline for phasing evaluation')

# data:
parser.add_argument('--test', '-T', help='input hdf5 file', action='store',
                    dest="test", default=None, type=str)

parser.add_argument('--eval', '-E', help='input hdf5 file to evaluate against',
                    action='store', dest="eval", default=None, type=str)

parser.add_argument('--output', '-O', help='output filename', action="store",
                    dest="output", default=None, type=str, required=True)

parser.add_argument('--fasta', '-F', action='store', default=None,
                    dest='fasta', help='FASTA file', type=str,
                    required=True)

parser.add_argument('--chr', '-C', default=None, required=True, action="store",
                    dest='chrom', help='Which contig to evaluate')

# non-required parameters
parser.add_argument('--accessibility', '-A', action='store', default=None,
                    dest='accessibility', help='Accessibility h5', type=str,
                    required=False)

parser.add_argument('--samples', '-S', action='store', nargs="+", default=None,
                    dest='samples', required=False,
                    help='Which samples to evaluate.')

parser.add_argument('--roh', '-R', action='store',
                    default=None, dest='roh', type=str,
                    help='Path to ROH file for calculating switch distance')

parser.add_argument('--overlap', '-wo', action='store',
                    default=250000, dest='overlap', type=int,
                    help='Overlap between windows.')

parser.add_argument('--windowsize', '-ws', action='store', default=500000,
                    type=int, dest='winsize',
                    help='Evenly spaced windows across genome')

args = parser.parse_args()
genome = pyfasta.Fasta(args.fasta)
contig_length = len(genome[args.chrom])
reference_gaps = find_reference_gaps(genome[args.chrom])

test_fh = open_group(args.test)
eval_fh = open_group(args.eval)

dtype = eval_fh[args.chrom]["samples"][:].dtype
eval_samples = eval_fh[args.chrom]["samples"][:].tolist()
test_samples = test_fh[args.chrom]['samples'][:].astype(dtype).tolist()

if args.samples is not None:
    test_idx = [test_samples.index(s.encode()) for s in args.samples]
    eval_idx = [eval_samples.index(s.encode()) for s in args.samples]
    sample_names = [s.encode() for s in args.samples]
else:
    intersection = np.intersect1d(eval_samples, test_samples)
    test_idx = [test_samples.index(s) for s in intersection]
    eval_idx = [eval_samples.index(s) for s in intersection]
    sample_names = intersection

if args.roh is not None:
    if args.roh.endswith(".h5"):
        roh = read_roh_h5(args.roh, sample_names, args.chrom)
    elif args.roh.endswith(".zarr2"):
        roh = read_roh_zarr(args.roh, sample_names, args.chrom)
    else:
        raise ValueError("Bad filepath provided: {0}. Only hdf5/zarr supported.".format(path))
else:
    roh = {s: None for s in sample_names}


test_genotypes = allel.GenotypeChunkedArray(test_fh[args.chrom]["calldata"]["genotype"])
eval_genotypes = allel.GenotypeChunkedArray(eval_fh[args.chrom]["calldata"]["genotype"])

# now lets line up postions...
e_pos = allel.SortedIndex(eval_fh[args.chrom]['variants']['POS'])
t_pos = allel.SortedIndex(test_fh[args.chrom]['variants']['POS'])
loc_e, loc_t = e_pos.locate_intersection(t_pos)

# Subset data
test_gt_subset = test_genotypes.subset(loc_t, test_idx)
eval_gt_subset = eval_genotypes.subset(loc_e, eval_idx)

positions = e_pos.compress(loc_e)

# conversely all eval MUST be in test set.
# Not true - as multiallelics are present too.
assert eval_gt_subset.shape == test_gt_subset.shape, \
    ("Not same shape:", eval_gt_subset.shape, test_gt_subset.shape)

is_missing = eval_gt_subset.is_missing()

scan_windows = allel.stats.window.position_windows(pos=None, start=1,
    stop=contig_length, size=args.winsize, step=args.overlap)

res = np.empty((len(sample_names), 3, scan_windows.shape[0]))

for idx, sid in enumerate(sample_names):

    # Markers are hets in the evaluation set.
    # There is a check for equivalency 
    hz = eval_gt_subset[:, idx].is_het()
    marker_positions = np.compress(hz, positions)
    sampleval_gt_subset = np.compress(hz, test_gt_subset[:, idx], axis=0)
    sample_gs = np.compress(hz, eval_gt_subset[:, idx], axis=0)

    switch = determine_haplotype_switch(sample_gs, sampleval_gt_subset)

    res[idx, :, :] = calculate_switch_distances(
        windows=scan_windows, 
        switch_array=switch, 
        marker_pos=marker_positions,
        rohz=roh[sid], 
        gaps=reference_gaps)

# now summarize across the samples dimension.
sum_r = res.sum(axis=0)

distance, n_markers, n_errors = sum_r
n_markers = n_markers.astype("int")
n_errors = n_errors.astype("int")

marker_err_rate = n_errors/(n_markers + n_errors)

mean_marker_d = distance/n_markers
mean_switch_d = mean_marker_d/marker_err_rate

df = pd.DataFrame.from_items((("start", scan_windows.T[0]),
                              ("stop", scan_windows.T[1]),
                              ("n_markers", n_markers),
                              ("n_errors", n_errors),
                              ("err_rate", marker_err_rate),
                              ("distance", distance),
                              ("mean_marker_dist", mean_marker_d),
                              ("mean_switch_dist", mean_switch_d)))

df.to_csv(args.output, sep="\t", index=False)
