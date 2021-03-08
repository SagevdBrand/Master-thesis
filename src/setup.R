################# INITIAL SETUP ##################
## Script to define paths to the simulated data ##
##################################################
source("./src/libraries.R")

#####################
### General paths ###
#####################

setting_path <- "./data/simulation settings/"
data_path <- "./data/simulation data/"
val_data_path <- "./data/validation data/"

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
    "shrinkage",
    "prev",
    "model",
    "pred_selection",
    "n",
    "expected events",
    "approach",
    "auc",
    "auc_se",
    "auc_ci_lower",
    "auc_ci_upper",
    "calib_int",
    "calib_int_se",
    "calib_slope",
    "calib_slope_se",
    "R2_CS",
    "R2_CS_se",
    "eci",
    "eci_se",
    "mape",
    "mape_se",
    "error_info"
  )

apparent_col_names <-
  c("approach",
    "auc",
    "calib_int",
    "calib_slope",
    "R2_CS",
    "eci",
    "mape",
    "error_info"
    )

iv_colnames <- 
  results_estimands_names[12:27]

study_info <- c("dim",
    "shrinkage",
    "prev",
    "model",
    "pred_selection",
    "n")

#####################
###### study 1 ######
#####################

study_1_settings <- "./data/simulation settings/s1.Rds"
study_1_data <- "./data/simulation data/study 1/"
study_1_val_data <- "./data/validation data/study 1/"
s1_estimands <- "./results/output/estimands/study 1/"
s1_performance <- "./results/output/performance/study 1/"
s1_figures <- "./results/figures/study 1/"

#####################
###### study 2 ######
#####################

#####################
###### study 3 ######
#####################

#####################
###### study 4 ######
#####################


# for (i in 1:4){
# study_specific <- paste0("study ", i, "/")
# dir.create(file.path(paste0("Output/", study_specific)), recursive = TRUE)
#  }

# datasets <- list.dirs(path = "data", recursive = F, full.names = F)
# data_files <- glue(datasets, ".csv")
# data_dir <- glue("data/", datasets)


