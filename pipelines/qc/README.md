# QC pipeline guide

## Overview

This pipeline takes a sampleset run through the vo-agam pipeline, and runs scripts to estimate contamination and to calculate some alignment summary stats.
As part of this it applies filters which blacklist failing samples and makes a sex call based on the ratio of coverage of 3L:X.
The output is a `tsv` file that contains the pertinent statistics for each sample, as well as filters to be applied.

Filter definitions are specified in `config.yml`

## Running the pipeline

### 1. Ensure the data are available in the relevant environment. In this case the Oxford cluster.

The input to this pipeline are a set of zipped zarr files, in the standard `calldata/GT` etc. format.

In the case of *Anopheles gambiae* these are copied from the Cambridge cluster to a location on the Oxford cluster nightly. 

An rsync command may look something like this:

`rsync -options othercluster.ac.uk:/path/on/other/cluster/ /path/to/inputdir`

The command used for the AG1000G project (to foxtrot from malsrv2) is:

```
# first set up a port forward to local host to navigate two-step Sanger system
ssh -N -f -L localhost:2222:malsrv2:22 nh7@ssh.sanger.ac.uk

# remove the "n" flag to take off dry-run mode
rsync --progress -nrogptvLe "ssh -p 2222" \
  nh7@localhost:/lustre/scratch118/malaria/team112/pipelines/setups/vo_agam/symlink/vo_agam_vcf2zarr/convert_vcf_to_zarr_serial_contig/ \
  /kwiat/vector/mirror/observatory/vo_agam/genotype_zarrs/
```

### 2. Ensure your repository contains `vector-tools` as a submodule. 

The second thing to check is that the `vector-tools` repository has been added as a submodule of your repository. In this example case, the repository for the AG1000G project is `vector-ops`. All commands run here are done so from the root directory of the `vector-ops` repository.

### 3. Create the configfile from the example file provided.

In `vector-tools/pipelines/qc` the `config.yml` file is that used by the AG1000G project. Copy this file into an appropriate place in your repo (for `vector-ops` we used `pipelines/qc/config.yml`, and edit the filepaths and filter settings where necessary. 

### 4. Run the pipeline!

To run the pipeline (e.g. on the Uganda sampleset) invoke snakemake using the `SGE` system on the Oxford cluster (from `vector-ops`):

```
bash vector-tools/pipelines/qc/submit.sh qc_report \
  --config fofn=${pwd}/tracking/AG1000G-UG/agam.fofn.tsv \
  --configfile vector-tools/pipelines/qc/config.yml \
  -d /kwiat/vector/observatory/production/vo_agam/qc
```

The path to the `fofn` must be specified. This path is used to interpolate the name of the sampleset, which is the name of the directory containing the `fofn`. The sampleset name is used to define the output directory.

The `-d` argument specifies the output directory. The `-d` flag is used instead of a config entry to ensure that any additional snakemake files are created in the designated output directory.

This pipeline creates a file in `/kwiat/vector/observatory/production/vo_agam/qc/AG1000G-UG`, `qc_summary.tsv`.

### 5. Further analysis

Once the output file `qc_summary.tsv` has been created in the output directory, a new branch should be created in your respository, and the file committed. You may find it easier to push to the remote, and then pull the branch down to a local machine, where a jupyter notebook visualisation can be carried out. 

For example, from the cluster:
```
git checkout -b xxx-qc-uganda
cp /kwiat/vector/observatory/production/vo_agam/qc/AG1000G-UG/qc_summary.tsv tracking/AG1000G-UG/qc/
cp /kwiat/vector/observatory/production/vo_agam/qc/AG1000G-UG/qc_config.yml tracking/AG1000G-UG/qc/
git add tracking/AG1000G-UG/qc/
git commit -m "add qc summary and config file for Uganda sample set"
git push -u origin xxx-qc-uganda

```

Then from a local machine:
```
git fetch
git checkout xxx-qc-uganda
```
Then using the notebook in `vector-ops/docs/guides/qc-visualisation.ipynb` as a base, proceed to examine the QC stats and identify outliers etc. The new notebook should also be added and committed to the branch.


## Requirements

- up to date git repositories `vector-ops` with the `vector-tools` submodule
- access to the set of genotype files (in zipped zarr format) from the agam pipeline.
- An environment containing `snakemake v5.1` or higher. We anticipate using `vector-ops/binder`.

## Notes on the config

Some of the entries in the `config.yml` merit some explanation. 

```
## Contamination
# Number of sites to consider for the contamination analysis
downsample: 50000

# minimum allele frequency of sites to consider for contamination analysis. (higher MAF SNPs are more informative).
minimum_af: 0.01

# minimum coverage to include a site in the contamination analysis.
minimum_coverage: 10

# which chromosomes/contigs to perform contamination analysis on.
seqid_contam: ["3L", "3R"]

# sequence error rate, required for contamination estimation algorithm
ser: 0.001

# sites at which samples have been genotyped
genotyped_sites: /kwiat/vector/observatory/resources/agam/ag.allsites.nonN.zarr.zip

# allele frequencies of above sites for contamination
allele_frequencies: /kwiat/vector/observatory/resources/agam/phase2_biallelic_variantsites_af.zarr

## Alignment summary statistics
# which chromosomes to calculate alignment summary stats over.
seqid_stats: ["3L", "3R", "X", "2L", "2R"]

# input/zip zarrs (for both contamination and alignment summary)
input_location: /kwiat/vector/mirror/observatory/vo_agam/genotype_zarrs
```
