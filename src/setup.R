################# INITIAL SETUP ##################
## Script to define paths to the simulated data ##
##################################################
source("./src/libraries.R")
#####################
### General paths ###
#####################

setting_path <- "./data/simulation settings/"
val_data_path <- "./data/validation data/"

### Results
## Visualization
full_estimand_plots <- "./results/figures/full estimands/"
full_performance_plots <- "./results/figures/full performance measures/"
estimand_plots <- "./results/figures/estimands thesis/"
performance_measures_plots <- "./results/figures/performance measures thesis/"

## Estimands
estimands_general_path <- "./results/output/estimands/"
estimands_path <- "./results/output/estimands/all/"
app_ext_path <- "./results/output/estimands/app_ext/"
cv_10_fold_path <- "./results/output/estimands/10_fold_cv/"
cv_5_fold_path <- "./results/output/estimands/5_fold_cv/"
cv_10x10_fold_path <- "./results/output/estimands/10x10_fold_cv/"
bootstrap_path <- "./results/output/estimands/bootstrap/"

## Performance measures
performance_general_path <- "./results/output/performance/"

## errors occurred durings sims
errors_path <- "./results/output/errors/"

#######################################
## Results names for estimand matrix ##
#######################################
results_estimands_names <-
  c(
    "seed",
    "iteration",
    "study",
    "scenario",
    "dim",
    "prev",
    "model",
    "pred_selection",
    "n",
    "expected events",
    "observed events",
    "approach",
    "auc",
    "auc_se",
    "calib_int",
    "calib_int_se",
    "calib_slope",
    "calib_slope_se",
    "Tjur",
    "Tjur_se",
    "R2_CS",
    "R2_CS_se",
    "eci",
    "eci_se",
    "mape",
    "mape_se",
    "rmspe",
    "rmspe_se",
    "error_info"
  )

apparent_col_names <-
  c(
    "approach",
    "auc",
    "calib_int",
    "calib_slope",
    "Tjur",
    "R2_CS",
    "eci",
    "mape",
    "rmspe",
    "error_info"
    )

estimands_names <- 
  c(
    "auc",
    "calib_int",
    "calib_slope",
    "Tjur",
    "R2_CS",
    "eci",
    "mape",
    "rmspe"
  )

estimands_se_names <- 
  c(
    "auc_se",
    "calib_int_se",
    "calib_slope_se",
    "Tjur_se",
    "R2_CS_se",
    "eci_se",
    "mape_se",
    "rmspe_se"
  )

iv_colnames <- c(
  "approach",
  "auc",
  "auc_se",
  "calib_int",
  "calib_int_se",
  "calib_slope",
  "calib_slope_se",
  "Tjur",
  "Tjur_se",
  "R2_CS",
  "R2_CS_se",
  "eci",
  "eci_se",
  "mape",
  "mape_se",
  "rmspe",
  "rmspe_se",
  "error_info"
)

study_info <- c("dim",
    "prev",
    "model",
    "pred_selection",
    "n")

#####################
###### study 1 ######
#####################

s1_performance <- "./results/output/performance/study 1/"
s1_figures <- "./results/figures/study 1/"

#####################
###### study 2 ######
#####################

s2_performance <- "./results/output/performance/study 2/"
s2_figures <- "./results/figures/study 2/"

#####################
###### study 3 ######
#####################

s3_performance <- "./results/output/performance/study 3/"
s3_figures <- "./results/figures/study 3/"


### RUN THE FOLLOWING ONCE ON THE HPC
# sessionInfo()
