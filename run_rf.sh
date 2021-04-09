#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=s.a.g.e.vandenbrand@uu.nl
RUN_NUMBER=$1


source ./vars.sh
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  exit 1
else
  ITERATION=$SLURM_ARRAY_TASK_ID
fi

echo "going to run $ITERATION for run $RUN_NUMBER"

Rscript -e 'source("./src/one_run_rf.R")' "$RUN_NUMBER" "$ITERATION"
