#! /bin/bash
snakemake --snakefile `dirname $0`/Snakefile $@ \
    --jn "{rulename}.{jobid}.sh"
