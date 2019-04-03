#! /bin/bash
set -e
mkdir -p _log

snakemake -p $@ \
    --cluster 'qsub -v PATH="/home/njh/miniconda3/bin:$PATH" -j y -o $(pwd)/_log/ -b n -l {params.req} -S /bin/bash' \
    --jobs 99 \
    --latency-wait 600 \
    --restart-times 0 \
    --jn "{rulename}.{jobid}.sh" \
    --use-conda

