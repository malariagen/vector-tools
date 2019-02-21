#! /bin/bash
snakemake --snakefile `dirname $0`/Snakefile $@ \
    --cluster 'qsub -v PATH="binder/deps/conda/bin::$PATH" -j y -b n -l {params.req},h=!foxtrot.well.ox.ac.uk -S /bin/bash' \
    --jobs 50 \
    --latency-wait 300 \
    --restart-times 0 \
    --jn "{rulename}.{jobid}.sh" \
    --use-conda
