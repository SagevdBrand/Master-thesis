# Simulation scripts

## Folder structure
In this main folder one can find all the R-scripts (`.R`) and bash scripts
(`.sh`) that were used to run the simulation studies for the thesis.

In the folder `old` the older versions of the scripts and testing scripts are
stored as part of the process during the year I wrote the code.

## Bash scripts
There are several bash scripts present. The most important being `run.sh` and
`vars.sh`. The latter was used to inform the High Performance Computer which
version of R to use and where to find the packages needed to run the simulation.
In case of reproducing the results, please adjust this file accordingly.

The `run.sh` file is the basic script used to run the simulations. Here three
arguments are created based on the input of an array job. The run-id is taken
from the array input. For example when running the command `sbatch --array=1-60
run.sh 1-500`. All 60 scenarios will be run for 500 times in the array job.

However, please note that some scenarios require a long time for a single
iteration. Meaning that the computer will not be able to finish 500 iterations
within the given 24 hours. Therefore, one could either reduce the number of
iterations to run within the command above, or increase the computation time
allowed.

Running the largest sample size scenario for the random forest
model may take more than 36 hours to complete a single iteration. Therefore,
there are also bash files created specifically for the random forest scenario.
Where `run_rf_short.sh`, asks for 12 hours for a single iteration (used
  for the random forest with n/2 sample size),
`run_rf_middle.sh` (can be used for the sample size n) for 24 hours,
 and `run_rf.sh` for 48 hours (for the n*2 scenario).
An example command to run this script is `sbatch --array=1-500
run_rf_short.sh 58`, which will create 500 jobs to be run for scenario 58, where
each job has a time limit for 12 hours. Each job refers to a single iteration.

## R scripts
Most R scripts quite self explanatory in their names. However, to be crystal
clear here is a short description of each script:
- `beta optimization script.R`: A script used to optimize regression coefficients
for the data generating mechanism. Requires study scenarios to be specified
regarding event-fraction and AUC value.
- `check_results.R`: A script meant to calculate the performance measures from
all of the simulation estimands. Requires the simulation studies to be performed
 first, and having called `error handling.R` to count errors that occured during
 the simulations and deal with for example missing calibration slopes.
- `create validation data`: A script to generate validation datasets belonging
to each scenario. Requires the simulation scenarios to be specified, by running
`study scenarios.R`.
- `data generation functions.R`: Functions necessary for the optimization of the
regression coefficients and for the creation of both development datasets and
the validation datasets.
- `error handling.R`: After all results of the simulations are stored, this
script binds all files from the
[all estimands folder](https://github.com/SagevdBrand/Master-thesis/tree/master/results/output/estimands/all)
into a single .rds file.
Thereafter, the errors that occured are counted and in those cases where the
calibration sloped could not be estimated because of degenerate linear
predictors. The values are filled with the maximum value of the calibration
slope found in the same scenario. Eventually this script saves the processed
results as an .rds file within
[the main estimands folder](https://github.com/SagevdBrand/Master-thesis/tree/master/results/output/estimands).
- `estimand functions.R`: This script defines all the functions used in
`one_run.R ` and `one_run_rf.R` to compute the final estimands to be processed
in `error handling.R`.
- `libraries.R`: Simply calls all the necessary packages to be used for all
scripts.
- `one_run.R` and `one_run_rf.R`: Both are the main script to run the
simulations. They require the study scenarios to be specified within
`study scenarios` and the validation data to be generated from `create
validation data.R`.
- `setup_packages.R`: An initial script to be run when not all packages within
the `libraries.R` are installed yet.
-  `setup.R`: This script defines all the file paths called for in the other
scripts. When changing the folder structure, please make sure to change this
script accordingly.
- `study scenarios.R`: Within this script the study scenarios are defined. It
requires that the regression coefficients of the data generating mechanism have
already been defined through the `beta optimization script.R`.
- `visualization.R`: Within this script the visualization of the results are
created and then stored in the
[figures folder](https://github.com/SagevdBrand/Master-thesis/tree/master/results/figures). 
