##################################
##################################
##################################
## GET ESTIMANDS SCENARIO 1

############ TO FIX/DO #################
## [x] IMPLEMENT STUDY 1 FOR ESTIMANDS
## [ ] IMPLEMENT STUDY 2 FOR ESTIMANDS
## [ ] IMPLEMENT STUDY 3 FOR ESTIMANDS

## [ ] PERFORMANCE MEASURES IN DATAFRAME
##      [ ] ADD SEED
##      [ ] ADD STUDY SCENARIO
##      [ ] ADD COLUMN FOR IV/APP.
##      [ ] ADD PERFORMANCE MEASURE (COMPARED TO EXT.)

## [ ] BOOTSTRAP ESTIMAND FUNCTION

## [ ] ADD OPTIONS FOR OTHER MODELS
##      [ ] MACHINE LEARNING -> CODE IS BEING DEVELOPED

## [ ] BUILD IN ERROR HANDLING AS SPECIFIED IN PROTOCOL!
##      [X] CHECK FOR VAR(LP) == 0 in BE:
##            [ ] RETURN HIGHEST VALUE FOR CALIBRATION SLOPE WITHIN THAT SCENARIO
##      [X] CHECK FOR VAR(LP) == 0 in LASSO:
##            [ ] RETURN HIGHEST VALUE FOR CALIBRATION SLOPE WITHIN THAT SCENARIO


## DONE:
## [X] CREATE DIFFERENT DGM-PAR IN DIFFERENT STUDIES
## [X] FIX SPAN ISSUES WITH LOESS -> MAKE SPAN WIDER
## [X] OBTAIN DATA FOR STUDY 2
## [X] OBTAIN DATA FOR STUDY 3
## [x] ADD MAPE 
## [x] ADD RMSPE
## [ ] BUILD IN ERROR HANDLING AS SPECIFIED IN PROTOCOL!
##      [X] IF ERROR OCCURS, MAKE SURE IT CONTINUES AND JUST RETURNS AN ERROR WITHIN THE RESULTS VECTOR
##      [x] cHECK FOR SEPARATION ISSUES
##            [x] FOR ML
##      [X] CHECK FOR CONVERGENCE ISSUES 
## [X] RIDGE -> CODE IS READY, ONLY NEEDS IMPLEMENTATION
## [X] LASSO -> CODE IS READY ONLY NEEDS IMPLEMENTATION

## [X] CHECK FOR VAR(Y) == 0 |SUM(Y) < 8 | N - SUM(Y) < 8  FOR LASSO AND RIDGE REGRESSION
## [X] ADD OBSERVED NUMBER OF EVENTS
## [X] ADD TJUR (MAKE SURE RESULTS ARE STILL IN RIGHT COLUMNS)
## [X] CHECK FOR VAR(LP) == 0
## [X] RESULTS IN DATAFRAME:
##      [X] ADD SEED
##      [X] ADD STUDY SCENARIO
##      [X] ADD COLUMN FOR EACH IV/APP/EXT.
##      [X] ADD ESTIMANDS
##      [X] ADD ERROR MESSAGE COLUMN
##      [X] GET RID OF CI AUC
##      [X] CHECK FOR EVENTS IN RESAMPLING
## [X] CHECK FOR EVENTS VAR(Y) == 0 IN DEVELOPMENT DATASET
## [X] RETURN VECTOR OF NA/ERROR AS RESULT
## [x] 10 & 5 FOLD CV
## [x] 10X10 FOLD CV
## [x] PIPELINE FOR EASY GENERATION OF DATA
## [x] PIPELINE FOR STORING ESTIMAND RESULTS
## [X] OUTPUT SEED
## [X] INTEGRATE CREATE DATA FUNCTION SO WE DONT HAVE THE SAME DATA ALL THE TIME :')
## [X] LOOCV FUNCTION or not?
## [X] MAKE SURE IT WORKS WITH MORE OR LESS PREDICTORS
## [X] FUNCTION FOR EXTERNAL VALIDATION OF ALL MODELS USED IN SCENARIO

########################################


############################################################################
############ Load libraries, validation data and study settings ############
############################################################################

## Libraries, file paths and functions
source("./src/study scenarios.R")
source("./src/estimand functions.R")
source("./src/data generation functions.R")
source("./src/setup.R")

#########################################################
############## LET THE SIMULATION COMMENCE ##############
#########################################################

# Reproduction seed
set.seed(123)

# Store seed values
n_sim <- 1 # how many iterations?
seed_state <- sample(1:50000, n_sim)

system.time({for(j in 1:n_sim){
  set.seed(seed_state[j]) # for each run the next value in the state vector will be chosen (and saved!)
  
  ## Create and load simulation data
  s1_data <- generate_data(s1, validation = FALSE)
  s2_data <- generate_data(s2, validation = FALSE)
  s3_data <- generate_data(s3, validation = FALSE)
  ## Create and load validation data
  ## Create and load simulation data
  s1_val_data <- generate_data(s1, validation = TRUE)
  s2_val_data <- generate_data(s2, validation = TRUE)
  s3_val_data <- generate_data(s3, validation = TRUE)
  
  ## Obtain apparent and external estimands ##
  results_app_ext_s1 <- get_app_ext_results(study = s1[1:3,], df = s1_data[1:3], df_val = s1_val_data[1:3], studyname = "Study 1")
  results_app_ext_s2 <- get_app_ext_results(study = s2[1:3,], df = s2_data[1:3], df_val = s2_val_data[1:3], studyname = "Study 2")
  results_app_ext_s3 <- get_app_ext_results(study = s3[1:3,], df = s3_data[1:3], df_val = s3_val_data[1:3], studyname = "Study 3")
  ## Obtain internal validation estimands ##
  # 10 fold cross-validation
  results_10_cv_s1 <- get_cv_results(study = s1[1:3,], df = s1_data[1:3], V = 10, studyname = "Study 1")
  results_10_cv_s2 <- get_cv_results(study = s2[1:3,], df = s2_data[1:3], V = 10, studyname = "Study 2")
  results_10_cv_s3 <- get_cv_results(study = s3[1:3,], df = s3_data[1:3], V = 10, studyname = "Study 3")
  
  # 5 fold cross-validation
  results_5_cv_s1 <- get_cv_results(study = s1[1:3,], df = s1_data[1:3], V = 5, studyname = "Study 1")
  results_5_cv_s2 <- get_cv_results(study = s2[1:3,], df = s2_data[1:3], V = 5, studyname = "Study 2")
  results_5_cv_s3 <- get_cv_results(study = s3[1:3,], df = s3_data[1:3], V = 5, studyname = "Study 3")
  
  # 10X10 fold cross-validation 
  results_10x10_cv_s1 <- get_10x10_results(study = s1[1:3,], df = s1_data[1:3], V = 10, studyname = "Study 1")
  results_10x10_cv_s2 <- get_10x10_results(study = s2[1:3,], df = s2_data[1:3], V = 10, studyname = "Study 2")
  results_10x10_cv_s3 <- get_10x10_results(study = s3[1:3,], df = s3_data[1:3], V = 10, studyname = "Study 3")
  
  # Bootstrap 3 varieties in one go

  
  #################################################
  ########## Wrangling into nice format ###########
  #################################################
  
  ## Bind all results together 
  results_estimands_s1 <-
    rbind(results_app_ext, results_10_cv, results_5_cv, results_10x10_cv) ## ADD OTHER RESULTS FROM VALIDATION APPROACHES
  
  ## Filling in missing details:
  results_estimands_s1$iteration <- j
  results_estimands_s1$seed <- seed_state[j]
  results_estimands_s1 <- results_estimands_s1 %>% mutate(`expected events` = n * prev)
  
  
  # Saving estimands
  saveRDS(results_estimands_s1, file = paste0(s1_estimands, "s1_estimands_seed_", seed_state[j], ".Rds"))
  
  } # close simulation for loop
  
  }) # Close timing function



## Obtain external validation estimands based on models fitted


#################################
## Obtain performance measures ##
#################################











