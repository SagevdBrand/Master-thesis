#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=s.a.g.e.vandenbrand@uu.nl
source ./vars.sh
#for i in {1..5000}; do
#Rscript -e 'source("./src/one_run.R")' --run_number="$i"
Rscript -e 'source("./src/one_run.R")'
#done
