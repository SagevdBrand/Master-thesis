#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=s.a.g.e.vandenbrand@uu.nl
START_ITERATION=$1
END_ITERATION=$2

source ./vars.sh
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  exit 1
else
  RUN_NUMBER=$SLURM_ARRAY_TASK_ID
fi

echo "going to run $START_ITERATION until $END_ITERATION for run $RUN_NUMBER"

Rscript -e 'source("./src/one_run.R")' "$RUN_NUMBER" "$START_ITERATION" "$END_ITERATION"
