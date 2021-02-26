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


