set -e

if [[ -z "${SGE_LOG}" ]]; then
  SGE_LOG=$HOME
  echo "SGE_LOG not set. The log location defaults to $SGE_LOG. Remember this variable needs to be exported."
else
  echo "SGE_LOG set as $SGE_LOG"
fi

mkdir -p $SGE_LOG

#! /bin/bash
snakemake --snakefile `dirname $0`/Snakefile $@ \
    --cluster 'qsub -v PATH="binder/deps/conda/bin::$PATH" -j y -o $SGE_LOG -b n -l {params.req} -S /bin/bash' \
    --jobs 10 \
    --latency-wait 60 \
    --restart-times 0 \
    --jn "{rulename}.{jobid}.sh" \
    --use-conda
