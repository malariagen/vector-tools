## QC pipeline guide

### Overview

This pipeline takes a sampleset run through the vo-agam pipeline, and runs scripts to estimate contamination and to calculate some alignment summary stats.
As part of this it applies filters which blacklist failing samples and makes a sex call based on the ratio of coverage of 3L:X.
The output is a `tsv` file that contains the pertinent statistics for each sample, as well as filters to be applied.

Filter definitions are specified in `config.yml`

### Running the pipeline

To run the pipeline (e.g. on the Uganda sampleset) invoke snakemake using the `SGE` system on the Oxford cluster (from `vector-tools/pipelines/qc`):

```
bash submit.sh qc_report -p -d /kwiat/vector/observatory/production/vo_agam/qc --config fofn={path/to/fofn} --configfile config.yml
```

The path to the `fofn` must be specified. This path is used to interpolate the name of the sampleset using a regex. The sampleset is used to define the output directory.

The `-d` argument specifies the output directory. The `-d` flag is used instead of a config entry to ensure that any additional snakemake files are created in the designated output directory.

This pipeline creates a file in `/kwiat/vector/observatory/production/vo_agam/qc/AG1000G-UG `, `qc_summary.tsv` which should be copied to `vector-ops/tracking/AG1000G-UG/qc` (or equivalent).

### Requirements

- up to date git repositories `vector-ops` and `vector-tools`
- access to the set of genotype files from the agam pipeline.
- An environment containing `snakemake v5.1` or higher. We anticipate using `vector-ops/binder`.

### Notes on the config

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

# QC analysis
# minimum fraction of genome covered (by at least 1x) to not blacklist a sample.
min_fgc: 0.5

# minimum median coverage below which to blacklist samples
min_med_cov: 10

# min/max values of X/3L modal coverage for sex calling
min_female_xratio: 0.8
max_female_xratio: 1.1
min_male_xratio: 0.4
max_male_xratio: 0.6

# max % contamination before blacklisted for contamination
max_pc_contam: 4.5

# input/zip zarrs (for both contamination and alignment summary)
input_location: /kwiat/vector/mirror/observatory/vo_agam/genotype_zarrs

# output location for the QC work
output_location: /kwiat/vector/observatory/production/vo_agam/qc/{sampleset}/
```

### Visualisation

In `docs/guides` there is a `qc_visualisation.ipynb`. This runs through the diagnostics of an example sample set. 

This forms the basis of the diagnostics for each sample set and copied to `vector-ops/tracking/{sampleset}/qc` when complete.

