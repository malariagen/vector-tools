set -e

if [ -z ${SGE_log_dir+x} ]; then echo "SGE_log_dir is unset: default to home" && SGE_log_dir=$HOME; else echo "log dir is set to '$SGE_log_dir'"; fi

mkdir -p $SGE_log_dir

#! /bin/bash
snakemake $@ \
    --cluster 'qsub -v PATH="binder/deps/conda/bin::$PATH" -j y -o $SGE_log_dir -b n -l {params.req} -S /bin/bash' \
    --jobs 10 \
    --latency-wait 60 \
    --restart-times 0 \
    --jn "{rulename}.{jobid}.sh" \
    --use-conda
