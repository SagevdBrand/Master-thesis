#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=s.a.g.e.vandenbrand@uu.nl

source ./vars.sh
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  exit 1
else
  RUN_NUMBER=$SLURM_ARRAY_TASK_ID
fi

Rscript -e 'source("./src/one_run.R")' "$RUN_NUMBER"
