#! /bin/bash

mkdir -p _log
snakemake --snakefile `dirname $0`/Snakefile $@ \
    --cluster 'qsub -v PATH="/home/njh/miniconda3/bin:$PATH" -j y -o $(pwd)/_log/ -b n -l {params.req},h=!foxtrot.well.ox.ac.uk -S /bin/bash' \
    --jobs 50 \
    --latency-wait 300 \
    --restart-times 0 \
    --jn "{rulename}.{jobid}.sh" \
    --use-conda
