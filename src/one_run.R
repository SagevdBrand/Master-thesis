#####################################################
#####################################################
#####################################################
##### CODE FOR N ITERATIONS OF 1 SCENARIO

args = commandArgs(trailingOnly = TRUE)
scenario_id <- as.numeric(args[1])
start_iteration <- as.numeric(args[2])
end_iteration <- as.numeric(args[3])

############################################################################
############ Load libraries, validation data and study settings ############
############################################################################
## Libraries, file paths and functions
source("./src/setup.R")
source("./src/estimand functions.R")
source("./src/data generation functions.R")

## Load scenario settings
scenarios <- readRDS(paste0(setting_path, "studies.RDS"))
scenario <- scenarios[scenario_id,]

## Load validation data
val_data <- readRDS(paste0(val_data_path, "val_data_", scenario_id, ".Rds"))

## How many iterations?
num_iterations <- 1000
set.seed(123)
seed_state <- sample(1:500000, num_iterations)

#########################################################
############## LET THE SIMULATION COMMENCE ##############
#########################################################

for (j in start_iteration:end_iteration){
errors_during_sim <- ErrorsWarnings({
  print(j)
  
  set.seed(seed_state[j]) # set a new seed for each iteration
  
  ## Create and load simulation data
  sim_data <- generate_data(scenario, validation = FALSE)
  
  ## Print run time of a single run:
  runtime <- system.time({
  ############################################
  ## Obtain apparent and external estimands ##
  ############################################
  results_app_ext <- get_app_ext_results(study = scenario, df = sim_data, df_val = val_data, studyname = as.character(scenario_id))
  
  # Save just to be sure
  #saveRDS(results_app_ext, paste0(app_ext_path, "app_ext_estimands_", scenario_id, "_seed_", seed_state[j], "_iteration_", j, ".Rds")) 
  print("app_ext_done")
  ##############################
  ## 10 fold cross-validation ##
  ##############################
  results_10_cv <- get_cv_results(study = scenario, df = sim_data, V = 10, studyname = as.character(scenario_id))
  
  # Save just to be sure
  #saveRDS(results_10_cv, paste0(cv_10_fold_path, "10fcv_estimands_", scenario_id, "_seed_", seed_state[j], "_iteration_", j, ".Rds"))
  print("10_cv_done")
  #############################
  ## 5 fold cross-validation ##
  #############################
  results_5_cv <- get_cv_results(study = scenario, df = sim_data, V = 5, studyname = as.character(scenario_id))
  
  # Save just to be sure
  #saveRDS(results_5_cv, paste0(cv_5_fold_path, "5fcv_estimands_", scenario_id, "_seed_", seed_state[j], "_iteration_", j, ".Rds")) 
  print("5_cv_done")
  #################################
  ## 10X10 fold cross-validation ##
  #################################
  results_10x10_cv <- get_10x10_results(study = scenario, df = sim_data, V = 10, studyname = as.character(scenario_id))
  
  # Save just to be sure
  #saveRDS(results_10x10_cv, paste0(cv_10x10_fold_path, "10x10fcv_estimands_", scenario_id, "_seed_", seed_state[j], "_iteration_", j, ".Rds"))
  print("10x10_cv_done")
  ###############
  ## Bootstrap ##
  ###############
  results_bootstrap <- get_bootstrap_results(study = scenario, df = sim_data, nboot = 500, studyname = as.character(scenario_id))
  ## Save just to be sure
  #saveRDS(results_bootstrap, paste0(bootstrap_path, "bootstrap_estimands_", scenario_id, "_seed_", seed_state[j], "_iteration_", j, ".Rds"))
  print("bootstrap_done")
  })
  print(paste("Single iteration run time:", runtime))
  
  #################################################
  ########## Wrangling into nice format ###########
  #################################################
  
  ## Bind all results together 
  results_estimands <-
    rbind(results_app_ext,
          results_10_cv,
          results_5_cv,
          results_10x10_cv,
          results_bootstrap
          ) 
  
  ## Filling in missing details:
  results_estimands$iteration <- j
  results_estimands$seed <- seed_state[j]
  results_estimands <- results_estimands %>% mutate(`expected events` = n * prev)
  
  # Saving estimands
  saveRDS(results_estimands, file = paste0(estimands_path, "estimands_", scenario$study, "_", scenario$scenario, "_seed_", seed_state[j], ".Rds"))
  })# close Error warnings

} # close for loop





