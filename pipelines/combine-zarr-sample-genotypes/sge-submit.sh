#! /bin/bash
snakemake --snakefile `dirname $0`/Snakefile $@ \
    --cluster 'qsub -v PATH="binder/deps/conda/bin::$PATH" -j y -b n -l {params.req},h=!foxtrot.well.ox.ac.uk -S /bin/bash' \
    --jobs 50 \
    --ri \
    --latency-wait 300 \
    --max-jobs-per-second 0.1 \
    --max-status-checks-per-second 0.1 \
    --restart-times 1 \
    --jn "{rulename}.{jobid}.sh" \
    --use-conda
