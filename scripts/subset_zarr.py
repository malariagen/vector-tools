# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division
import sys
import zarr
import argparse
import allel
import numpy as np
import pandas as pd
import numcodecs


def log(*msg):
    print(*msg, file=sys.stderr)
    sys.stderr.flush()


def output_to_vcf(path, seq_id, positions, ref, alt, n_vars, gt, sample_id):

    log('Output VCF to {!r} ...'.format(path))

    # use pandas because supposed to be fast CSV output
    a1 = pd.Series(gt[:, 0, 0].astype('U1')).replace('-', '.')
    a2 = pd.Series(gt[:, 0, 1].astype('U1')).replace('-', '.')
    gt_text = (a1 + '/' + a2)
    df = pd.DataFrame.from_items([
        ('#CHROM', [seq_id] * n_vars),
        ('POS', positions),
        ('ID', ['.'] * n_vars),
        ('REF', ref.astype('U')),
        ('ALT', alt[:, 0].astype('U')),
        ('QUAL', ['.'] * n_vars),
        ('FILTER', ['.'] * n_vars),
        ('INFO', ['.'] * n_vars),
        ('FORMAT', ['GT'] * n_vars),
        (sample_id.replace('_', '-'), gt_text)
    ])
    with open(path, mode='w', encoding='ascii') as f:
        print('##fileformat=VCFv4.0', file=f)
        df.to_csv(f, index=False, sep='\t')


def output_to_zarr(path, seq_id, sample_id, arrays, cname, clevel, shuffle):

    log('Output zarr to {!r} ...'.format(path))

    store = zarr.ZipStore(path, mode='w')
    root = zarr.group(store=store)
    callset = root.create_group(sample_id)
    seq_group = callset.require_group(seq_id)
    calldata_group = seq_group.require_group('calldata')
    variants_group = seq_group.require_group('variants')

    compressor = numcodecs.Blosc(cname=cname, clevel=clevel, shuffle=shuffle)

    for key, value in arrays.items():
        calldata_group.create_dataset(key, data=value, compressor=compressor)
        log('Created output array: ' + repr(key))

    store.close()


def main():

    parser = argparse.ArgumentParser(description='Extract data from a Zarr store for a given set '
                                                 'of variants, recoding alleles, and output to '
                                                 'VCF.')
    parser.add_argument('--genotypes', required=True,
                        help='path to genotypes file, should be a zipped Zarr file with data for a '
                             'single sample')
    parser.add_argument('--sites', required=True,
                        help='path to Zarr store containing the sites the sample was '
                             'genotyped at')
    parser.add_argument('--target-sites', required=True,
                        help='path to Zarr store containing the sites to subset to')
    parser.add_argument('--seqid', required=True,
                        help='name of chromosome or contig to process')
    parser.add_argument('--output', required=True,
                        help='path to output location for VCF file')
    parser.add_argument('--cname', required=False, default='zstd',
                        help='name of compressor library, defaults to zstd')
    parser.add_argument('--clevel', required=False, default=1,
                        help='compression level, defaults to 1')
    parser.add_argument('--shuffle', required=False, default=0,
                        help='shuffle filter; 0 = no shuffle (default), 1 = byte shuffle, '
                             '2 = bit shuffle')
    args = parser.parse_args()

    log('Load target sites from {!r} ...'.format(args.target_sites))
    target_callset = zarr.open_group(args.target_sites, mode='r')
    seqid = args.seqid
    target_variants = target_callset[seqid]['variants']
    target_pos = allel.SortedIndex(target_variants['POS'])
    n_variants = target_pos.shape[0]
    log('Found {:,} target sites.'.format(n_variants))

    log('Load genotype sites from {!r} ...'.format(args.sites))
    if args.sites.endswith("zip"):
      sites_store = zarr.ZipStore(args.sites, mode='r')
      sites_callset = zarr.group(store=sites_store)
    else:
      sites_callset = zarr.open_group(args.sites, mode='r')

    sites_variants = sites_callset[seqid]['variants']
    pos = allel.SortedIndex(sites_variants['POS'])
    log('Found {:,} genotype sites.'.format(pos.shape[0]))

    log('Locate subset ...')
    loc_subset, _ = pos.locate_intersection(target_pos)
    assert np.count_nonzero(loc_subset) == n_variants, (np.count_nonzero(loc_subset), n_variants)

    log('Map alleles ...')
    ref = sites_variants['REF'].oindex[loc_subset]
    alt = sites_variants['ALT'].oindex[loc_subset]
    target_ref = target_variants['REF'][:]
    target_alt = target_variants['ALT'][:]
    target_alleles = np.column_stack([target_ref, target_alt])
    allele_mapping = allel.create_allele_mapping(ref, alt, target_alleles)

    log('Load genotypes from {!r} ...'.format(args.genotypes))
    if args.genotypes.endswith("zip"):
      genotypes_store = zarr.ZipStore(args.genotypes, mode='r')
      genotypes_callset = zarr.group(store=genotypes_store)
    else:
      genotypes_callset = zarr.open_group(args.genotypes, mode='r')

    # assume sample ID is only child of hierarchy root
    sample_id = next(iter(genotypes_callset))
    gt = genotypes_callset[sample_id][seqid]['calldata/GT'].oindex[loc_subset]
    gt = allel.GenotypeArray(gt)
    assert gt.shape[0] == n_variants
    assert gt.shape[1] == 1
    log('No. variant alleles: {:,}'.format(np.count_nonzero(gt > 0)))

    log('Recode genotypes ...')
    gt_recoded = gt.map_alleles(allele_mapping)
    # ensure there are no half-missing genotypes, convert to full missing
    loc_missing = np.any(gt_recoded < 0, axis=(1, 2))
    gt_recoded[loc_missing] = -1
    log('No. variant alleles after recoding: {:,}'.format(np.count_nonzero(gt_recoded > 0)))

    if args.output.endswith(".vcf"):
        output_to_vcf(args.output,
                      seqid,
                      target_pos,
                      target_ref,
                      target_alt,
                      n_variants,
                      gt_recoded,
                      sample_id)

    else:
        log('Load GQ from {!r} ...'.format(args.genotypes))
        gq = genotypes_callset[sample_id][seqid]['calldata/GQ'].oindex[loc_subset]
        log('Load AD from {!r} ...'.format(args.genotypes))
        ad = allel.AlleleCountsChunkedArray(
            genotypes_callset[sample_id][seqid]['calldata/AD'].oindex[loc_subset][:, 0])
        allele_mapping = allele_mapping.astype(ad.dtype)

        
        # option would be to change the target alleles?
        
        # most rows need a 2 and a 3. Some just need a 3.
        # first replace with 2s.
        # But one per row where you have two -1s. (This is most).
        # is -1 & is sum_neg == 2 & 
        is_neg = allele_mapping == -1
        cum_sum_neg = is_neg.cumsum(axis=1)
        allele_mapping[is_neg & (cum_sum_neg == 1)] = 3
        allele_mapping[(allele_mapping == -1)] = 2

        # now all should be 
        assert np.all(allele_mapping.sum(1) == 6), "all rows must sum to 6"

        ad_recoded = ad.map_alleles(allele_mapping)
        assert ad.min() >= -1, "AD min must be a minimum of -1 : {minv}".format(minv=ad.min())
        assert ad_recoded.min() >= -1, "AD recoded min must be a minimum of -1 : {minv}".format(minv=ad_recoded.min())
        ad_recoded = np.reshape(ad_recoded, (-1, 1, 4))

        assert gq.shape[0] == n_variants == ad_recoded.shape[0]
        assert gq.shape[1] == 1
        assert ad_recoded.shape[1] == 1, ad_recoded.shape
        assert ad_recoded.shape[2] == 4, ad_recoded.shape

        log('Load MQ from {!r} ...'.format(args.genotypes))
        mq = genotypes_callset[sample_id][seqid]['variants/MQ'].oindex[loc_subset][:, None]
        assert mq.shape[0] == n_variants
        assert mq.shape[1] == 1

        output_to_zarr(args.output,
                       seqid,
                       sample_id,
                       {"GT": gt_recoded, "AD": ad_recoded, "GQ": gq, "MQ": mq},
                       cname=args.cname,
                       clevel=args.clevel,
                       shuffle=args.shuffle)

    log('All done.')


if __name__ == '__main__':
    main()
