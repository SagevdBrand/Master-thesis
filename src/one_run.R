#####################################################
#####################################################
#####################################################
##### CODE FOR A SINGLE ITERATION
##### SEED IS DEPENDENT ON RUN-ID


############################################################################
############ Load libraries, validation data and study settings ############
############################################################################
## Libraries, file paths and functions
source("./src/setup.R")
source("./src/estimand functions.R")
source("./src/data generation functions.R")

## Load scenario settings
s1 <- readRDS(study_1_settings)
s2 <- readRDS(study_2_settings)
s3 <- readRDS(study_3_settings)

## Load validation data
s1_val_data <- readRDS(study_1_val_data)
s2_val_data <- readRDS(study_2_val_data)
s3_val_data <- readRDS(study_3_val_data)

#########################################################
############## LET THE SIMULATION COMMENCE ##############
#########################################################
args = commandArgs(trailingOnly = TRUE)

run_id <- as.numeric(args[1])

errors_during_sim <- ErrorsWarnings({
  
  set.seed(run_id) # for each run the job_id will be used as seed
  
  ## Create and load simulation data
  s1_data <- generate_data(s1, validation = FALSE)
  s2_data <- generate_data(s2, validation = FALSE)
  s3_data <- generate_data(s3, validation = FALSE)
  
  ## Obtain apparent and external estimands ##
  runtime <- system.time({
  results_app_ext_study_1 <- get_app_ext_results(study = s1, df = s1_data, df_val = s1_val_data, studyname = "study_1")
  results_app_ext_study_2 <- get_app_ext_results(study = s2, df = s2_data, df_val = s2_val_data, studyname = "study_2")
  results_app_ext_study_3 <- get_app_ext_results(study = s3, df = s3_data, df_val = s3_val_data, studyname = "study_3")
  })
  print(paste("app_ext", runtime))
  
  ## Save just to be sure
  saveRDS(results_app_ext_study_1, paste0(app_ext_path, "app_ext_estimands_study_1_seed_", run_id, ".Rds")) 
  saveRDS(results_app_ext_study_2, paste0(app_ext_path, "app_ext_estimands_study_2_seed_", run_id, ".Rds")) 
  saveRDS(results_app_ext_study_3, paste0(app_ext_path, "app_ext_estimands_study_3_seed_", run_id, ".Rds")) 
  
  ## Obtain internal validation estimands ##
  # 10 fold cross-validation
  runtime <- system.time({
  results_10_cv_s1 <- get_cv_results(study = s1, df = s1_data, V = 10, studyname = "study_1")
  results_10_cv_s2 <- get_cv_results(study = s2, df = s2_data, V = 10, studyname = "study_2")
  results_10_cv_s3 <- get_cv_results(study = s3, df = s3_data, V = 10, studyname = "study_3")
  })
  print(paste("10_cv", runtime))
  
  ## Save just to be sure
  saveRDS(results_10_cv_s1, paste0(cv_10_fold_path, "10fcv_estimands_study_1_seed_", run_id, ".Rds")) 
  saveRDS(results_10_cv_s2, paste0(cv_10_fold_path, "10fcv_estimands_study_2_seed_", run_id, ".Rds")) 
  saveRDS(results_10_cv_s3, paste0(cv_10_fold_path, "10fcv_estimands_study_3_seed_", run_id, ".Rds")) 
  
  
  # 5 fold cross-validation
  runtime <- system.time({
  results_5_cv_s1 <- get_cv_results(study = s1, df = s1_data, V = 5, studyname = "study_1")
  results_5_cv_s2 <- get_cv_results(study = s2, df = s2_data, V = 5, studyname = "study_2")
  results_5_cv_s3 <- get_cv_results(study = s3, df = s3_data, V = 5, studyname = "study_3")
  })
  print(paste("5_cv", runtime))
  
  ## Save just to be sure
  saveRDS(results_5_cv_s1, paste0(cv_5_fold_path, "5fcv_estimands_study_1_seed_", run_id, ".Rds")) 
  saveRDS(results_5_cv_s2, paste0(cv_5_fold_path, "5fcv_estimands_study_2_seed_", run_id, ".Rds")) 
  saveRDS(results_5_cv_s3, paste0(cv_5_fold_path, "5fcv_estimands_study_3_seed_", run_id, ".Rds")) 
  
  # 10X10 fold cross-validation 
  runtime <- system.time({
  results_10x10_cv_s1 <- get_10x10_results(study = s1, df = s1_data, V = 10, studyname = "study_1")
  results_10x10_cv_s2 <- get_10x10_results(study = s2, df = s2_data, V = 10, studyname = "study_2")
  results_10x10_cv_s3 <- get_10x10_results(study = s3, df = s3_data, V = 10, studyname = "study_3")
  })
  print(paste("10x10cv", runtime))
  
  ## Save just to be sure
  saveRDS(results_10x10_cv_s1 , paste0(cv_10x10_fold_path, "10x10fcv_estimands_study_1_seed_", run_id, ".Rds")) 
  saveRDS(results_10x10_cv_s2 , paste0(cv_10x10_fold_path, "10x10fcv_estimands_study_2_seed_", run_id, ".Rds")) 
  saveRDS(results_10x10_cv_s3 , paste0(cv_10x10_fold_path, "10x10fcv_estimands_study_3_seed_", run_id, ".Rds"))
  
  # Bootstrap
  runtime <- system.time({
  results_bootstrap_s1 <- get_bootstrap_results(study = s1, df = s1_data, nboot = 500, studyname = "study_1")
  results_bootstrap_s2 <- get_bootstrap_results(study = s2, df = s2_data, nboot = 500, studyname = "study_2")
  results_bootstrap_s3 <- get_bootstrap_results(study = s3, df = s3_data, nboot = 500, studyname = "study_3")
  })
  print(paste("res_bootstrap", runtime))
  
  ## Save just to be sure
  saveRDS(results_bootstrap_s1, paste0(bootstrap_path, "bootstrap_estimands_study_1_seed_", run_id, ".Rds")) 
  saveRDS(results_bootstrap_s2, paste0(bootstrap_path, "bootstrap_estimands_study_2_seed_", run_id, ".Rds")) 
  saveRDS(results_bootstrap_s3, paste0(bootstrap_path, "bootstrap_estimands_study_3_seed_", run_id, ".Rds"))
  
  
  #################################################
  ########## Wrangling into nice format ###########
  #################################################
  
  ## Bind all results together 
  results_estimands <-
    rbind(results_app_ext_study_1,
          results_app_ext_study_2,
          results_app_ext_study_3,
          results_10_cv_s1,
          results_10_cv_s2,
          results_10_cv_s3,
          results_5_cv_s1,
          results_5_cv_s2,
          results_5_cv_s3,
          results_10x10_cv_s1,
          results_10x10_cv_s2,
          results_10x10_cv_s3,
          results_bootstrap_s1,
          results_bootstrap_s2,
          results_bootstrap_s3
          ) 
  
  ## Filling in missing details:
  results_estimands$iteration <- run_id
  results_estimands$seed <- run_id
  results_estimands <- results_estimands %>% mutate(`expected events` = n * prev)
  
  
  # Saving estimands
  saveRDS(results_estimands, file = paste0(estimands_path, "estimands_seed_", run_id, ".Rds"))
  })# close Error warnings







