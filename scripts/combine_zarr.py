# -*- coding: utf-8 -*-
import argparse
import os
import sys
import zarr
import numcodecs
import dask
import dask.array as da
from dask.diagnostics import ProgressBar
import numpy as np


def log(*msg):
    print(*msg, file=sys.stderr)
    sys.stderr.flush()


def main():

    parser = argparse.ArgumentParser(description='Combine multiple per-sample Zarr files into a '
                                                 'single multi-sample Zarr store.')
    parser.add_argument('--samples', required=True,
                        help='path to tab-separated file containing manifest of samples to be processed, '
                             'file expected to have a header and sample names as a first column', )
    parser.add_argument('--input-pattern', required=True,
                        help='path to location of per-sample Zarr files, including "{sample}" '
                             'where sample ID should be substituted, e.g., "/path/to/{sample}.zip"')
    parser.add_argument('--output', required=True,
                        help='path to output Zarr directory')
    parser.add_argument('--seqid', required=True,
                        help='name of chromosome or contig to process')
    parser.add_argument('--field', required=True,
                        help='name of field to extract, e.g., "calldata/GT". If the root not given e.g. "GT" defaults to "calldata/GT".')
    parser.add_argument('--chunk-width', required=False, default=50,
                        help='chunk width, defaults to 50')
    parser.add_argument('--cname', required=False, default='zstd',
                        help='name of compressor library, defaults to zstd')
    parser.add_argument('--clevel', required=False, default=1,
                        help='compression level, defaults to 1')
    parser.add_argument('--shuffle', required=False, default=0,
                        help='shuffle filter; 0 = no shuffle (default), 1 = byte shuffle, '
                             '2 = bit shuffle')
    parser.add_argument('--num-workers', required=False, default=None, type=int,
                        help='number of parallel workers to use, defaults to number of cores')

    try:
        args = {
          "samples": snakemake.params.samples,
          "input_pattern": snakemake.params.input_pattern,
          "output": snakemake.params.output,
          "seqid": snakemake.params.seqid,
          "field": snakemake.params.field,
          "num_workers": snakemake.params.num_workers,
          "cname": snakemake.params.cname,
          "clevel": snakemake.params.clevel,
          "chunk_width": snakemake.params.chunk_width,
          "shuffle": snakemake.params.shuffle
        }
        log("Args read via snakemake")
    except NameError:
        args = vars(parser.parse_args())
        log("Args read via command line")

    log('Combine sample zarrs')
    log('====================')
    log('Begin, using zarr {}, numcodecs {}, dask {}.'
        .format(zarr.__version__, numcodecs.__version__, dask.__version__))

    samples = load_samples(path=args["samples"])
    check_genotypes_files(samples=samples, input_pattern=args["input_pattern"])

    # check field.
    if len(args["field"].split("/")) == 1:
        renamedfield = os.path.join("calldata", args["field"])
        log("Assuming given field {0} refers to {1}".format(args["field"], renamedfield))
        args["field"] = renamedfield

    arr = check_array_setup(samples=samples, input_pattern=args["input_pattern"], seqid=args["seqid"],
                            field=args["field"])
    output_arr = setup_output(output_path=args["output"], seqid=args["seqid"], field=args["field"],
                              example_arr=arr, samples=samples, cname=args["cname"],
                              clevel=args["clevel"], shuffle=args["shuffle"], chunk_width=args["chunk_width"])
    input_arr = setup_input(samples=samples, input_pattern=args["input_pattern"], seqid=args["seqid"],
                            field=args["field"])
    copy_data(input_arr=input_arr, output_arr=output_arr, num_workers=args["num_workers"])
    log('All done.')


def load_samples(path):
    log('Loading sample manifest from {!r} ...'.format(path))
    with open(path, mode='r') as f:
        samples = [s.strip().split('\t')[0] for s in f.readlines()[1:]]
    log('Found {:,} samples.'.format(len(samples)))
    return samples


def check_genotypes_files(samples, input_pattern):
    log('Checking genotypes files exist for all samples ...')
    for s in samples:
        input_path = input_pattern.format(sample=s)
        if not os.path.exists(input_path):
            raise Exception('genotypes file not found for sample {!r} at path {!r}'
                            .format(s, input_path))
    log('Found all genotypes files.')


def check_array_setup(samples, input_pattern, seqid, field):
    log('Determining number of variants ...')
    path = input_pattern.format(sample=samples[0])
    callset = zarr.group(zarr.ZipStore(path, mode='r'))
    # expect sample name in hierarchy
    array = callset[samples[0]][seqid][field]
    n_variants = array.shape[0]
    log('Found {:,} variants.'.format(n_variants))
    return array


def setup_output(output_path, seqid, field, example_arr, samples, cname, clevel, shuffle,
                 chunk_width):
    log('Setting up output at {!r} ...'.format(output_path))
    callset = zarr.open_group(output_path, mode='a')
    seq_group = callset.require_group(seqid)
    field_root, field_id = field.split("/")
    root_group = seq_group.require_group(field_root)
    output_shape = (example_arr.shape[0], len(samples)) + example_arr.shape[2:]

    c1 = 2 **26 // np.prod((chunk_width,) + example_arr.shape[2:])
    output_chunks = (c1, chunk_width) + example_arr.chunks[2:]

    compressor = numcodecs.Blosc(cname=cname, clevel=clevel, shuffle=shuffle)
    output_arr = root_group.empty_like(
      field_id, example_arr, shape=output_shape, 
      chunks=output_chunks, overwrite=True, compressor=compressor)
    log('Created output array: ' + repr(output_arr))
    return output_arr


def setup_input(samples, input_pattern, seqid, field):
    log('Setting up input array ...')
    input_paths = [input_pattern.format(sample=s) for s in samples]
    input_stores = [zarr.ZipStore(ip, mode='r') for ip in input_paths]
    input_roots = [zarr.group(store) for store in input_stores]
    input_arrays = [root[s][seqid][field] for root, s in zip(input_roots, samples)]
    input_arrays = [da.from_array(a, chunks=a.chunks) for a in input_arrays]

    # here we add a dim to allow the hstack to work. must share the shape (X, 1, )
    input_arrays = [a[:, None] if a.ndim == 1 else a for a in input_arrays]

    input_array = da.hstack(input_arrays)
    log('Input array:', input_array)
    return input_array


def copy_data(input_arr, output_arr, num_workers):
    log('Copying data ...')
    input_arr = input_arr.rechunk(output_arr.chunks)
    with ProgressBar(out=sys.stderr):
        input_arr.store(output_arr, lock=False, num_workers=num_workers)
    log('Copying done.')
    log(output_arr.info)


if __name__ == '__main__':
    main()
