#! /bin/bash
snakemake --snakefile `dirname $0`/Snakefile $@ \
    --cluster 'qsub -v PATH="binder/deps/conda/bin::$PATH" -j y -b n -l {params.req} -S /bin/bash -o $(pwd)' \
    --jobs 50 \
    --latency-wait 60 \
    --restart-times 0 \
    --jn "{rulename}.{jobid}.sh" \
    --use-conda
