##################################
##################################
##################################
## get estimands 1 scenario:
set.seed(123)

############ Load necessary stuff ############
source("scripts/setup.R")
source("scripts/libraries.R")
source("scripts/estimand functions.R")

#### Load data #####
data_files <- list.files(path = scenario_1_data, recursive = T, full.names = F)
df <- lapply(paste0(scenario_1_data,data_files),readRDS,.GlobalEnv)
names(df) <- data_files
s1 <- read_rds(paste0(scenario_1_settings,"s1.Rds"))

## Obtain apparent performance results ##
system.time(results_app <- get_app_results(scenario = s1, df = df))

## Obtain cross validation results ##
system.time(results_10_cv <- get_cv_results(scenario = s1, df = df, V = 10))
system.time(results_5_cv <- get_cv_results(scenario = s1, df = df, V = 5))










