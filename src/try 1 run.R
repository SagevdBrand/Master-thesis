##################################
##################################
##################################
## GET ESTIMANDS SCENARIO 1

############ TO FIX/DO #################
## [ ] PERFORMANCE MEASURES IN DATAFRAME
##      [ ] ADD SEED
##      [ ] ADD STUDY SCENARIO
##      [ ] ADD COLUMN FOR IV/APP.
##      [ ] ADD PERFORMANCE MEASURE (COMPARED TO EXT.)

## [ ] BUILD IN ERROR HANDLING AS SPECIFIED IN PROTOCOL!
##      [X] CHECK FOR VAR(LP) == 0 in BE:
##            [ ] RETURN HIGHEST VALUE FOR CALIBRATION SLOPE WITHIN THAT SCENARIO
##      [X] CHECK FOR VAR(LP) == 0 in LASSO:
##            [ ] RETURN HIGHEST VALUE FOR CALIBRATION SLOPE WITHIN THAT SCENARIO

## [ ] BUILD ERROR HANDLING FOR SIMULATION RUNS!!

## DONE:
## [x] BOOTSTRAP ESTIMAND FUNCTION
## [x] IMPLEMENT STUDY 1 FOR ESTIMANDS
## [X] IMPLEMENT STUDY 2 FOR ESTIMANDS
## [X] IMPLEMENT STUDY 3 FOR ESTIMANDS
## [X] ADD OPTIONS FOR OTHER MODELS
##      [X] MACHINE LEARNING -> CODE IS BEING DEVELOPED
##          [X] RF
##          [X] SVM
##          [X] CART
##          [X] ANN
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
source("./src/setup.R")
source("./src/estimand functions.R")
source("./src/data generation functions.R")

## Obtain scenario settings
s1 <- readRDS(study_1_settings)
s2 <- readRDS(study_2_settings)
s3 <- readRDS(study_3_settings)

## Obtain validation data
s1_val_data <- readRDS(study_1_val_data)
s2_val_data <- readRDS(study_2_val_data)
s3_val_data <- readRDS(study_3_val_data)

#########################################################
############## LET THE SIMULATION COMMENCE ##############
#########################################################

# Reproduction seed
set.seed(123)

# Store seed values
n_sim <- 1 # how many iterations?
seed_state <- sample(1:50000, n_sim)

system.time({
set.seed(seed_state[1]) # for each run the next value in the state vector will be chosen (and saved!)
  
  ## Create and load simulation data
  s1_data <- generate_data(s1, validation = FALSE)
  s2_data <- generate_data(s2, validation = FALSE)
  s3_data <- generate_data(s3, validation = FALSE)
  
  progression<-winProgressBar(title = "Progress bar", min = 0,max = 100 , width = 300)
  l <- seq(1, 100, length.out = 15)
  
  
  ## Obtain apparent and external estimands ##
  results_app_ext_s1 <- get_app_ext_results(study = s1, df = s1_data, df_val = s1_val_data, studyname = "study_1")
  setWinProgressBar(progression, l[1], title=paste(l[1],"% done"))
  
  results_app_ext_s2 <- get_app_ext_results(study = s2, df = s2_data, df_val = s2_val_data, studyname = "study_2")
  setWinProgressBar(progression, l[2], title=paste(l[2],"% done"))
  
  results_app_ext_s3 <- get_app_ext_results(study = s3, df = s3_data, df_val = s3_val_data, studyname = "study_3")
  setWinProgressBar(progression, l[3], title=paste(l[3],"% done"))
  
  saveRDS(results_app_ext_s1, paste0(estimands_path, "app_estimands_for_bootstrap_study1.RDS")) #if bootstrap is not done yet, we can retest
  saveRDS(results_app_ext_s2, paste0(estimands_path, "app_estimands_for_bootstrap_study2.RDS")) #if bootstrap is not done yet, we can retest
  saveRDS(results_app_ext_s3, paste0(estimands_path, "app_estimands_for_bootstrap_study3.RDS")) #if bootstrap is not done yet, we can retest
  
  saveRDS(p_app_study_1, paste0(estimands_path, "app_preds_for_bootstrap_study1.RDS")) #if bootstrap is not done yet, we can retest
  saveRDS(p_app_study_2, paste0(estimands_path, "app_preds_for_bootstrap_study2.RDS")) #if bootstrap is not done yet, we can retest
  saveRDS(p_app_study_3, paste0(estimands_path, "app_preds_for_bootstrap_study3.RDS")) #if bootstrap is not done yet, we can retest
  
  
  ## Obtain internal validation estimands ##
  # 10 fold cross-validation
  results_10_cv_s1 <- get_cv_results(study = s1, df = s1_data, V = 10, studyname = "study_1")
  setWinProgressBar(progression, l[4], title=paste(l[4],"% done"))
  
  results_10_cv_s2 <- get_cv_results(study = s2, df = s2_data, V = 10, studyname = "study_2")
  setWinProgressBar(progression, l[5], title=paste(l[5],"% done"))
  
  results_10_cv_s3 <- get_cv_results(study = s3, df = s3_data, V = 10, studyname = "study_3")
  setWinProgressBar(progression, l[6], title=paste(l[6],"% done"))
  
  # 5 fold cross-validation
  results_5_cv_s1 <- get_cv_results(study = s1, df = s1_data, V = 5, studyname = "study_1")
  setWinProgressBar(progression, l[7], title=paste(l[7],"% done"))
  
  results_5_cv_s2 <- get_cv_results(study = s2, df = s2_data, V = 5, studyname = "study_2")
  setWinProgressBar(progression, l[8], title=paste(l[8],"% done"))
  
  results_5_cv_s3 <- get_cv_results(study = s3, df = s3_data, V = 5, studyname = "study_3")
  setWinProgressBar(progression, l[9], title=paste(l[9],"% done"))
  
  # 10X10 fold cross-validation 
  results_10x10_cv_s1 <- get_10x10_results(study = s1, df = s1_data, V = 10, studyname = "study_1")
  setWinProgressBar(progression, l[10], title=paste(l[10],"% done"))
  
  results_10x10_cv_s2 <- get_10x10_results(study = s2, df = s2_data, V = 10, studyname = "study_2")
  setWinProgressBar(progression, l[11], title=paste(l[11],"% done"))
  
  results_10x10_cv_s3 <- get_10x10_results(study = s3, df = s3_data, V = 10, studyname = "study_3")
  setWinProgressBar(progression, l[12], title=paste(l[12],"% done"))
  
  # Bootstrap 3 varieties in one go
  p_app_study_1 <- readRDS(paste0(estimands_path, "app_preds_for_bootstrap_study1.RDS")) 
  p_app_study_2 <- readRDS(paste0(estimands_path, "app_preds_for_bootstrap_study2.RDS"))
  p_app_study_3 <- readRDS(paste0(estimands_path, "app_preds_for_bootstrap_study3.RDS"))
  
  results_app_ext_s1 <- readRDS(paste0(estimands_path, "app_estimands_for_bootstrap_study1.RDS"))
  results_app_ext_s2 <- readRDS(paste0(estimands_path, "app_estimands_for_bootstrap_study2.RDS"))
  results_app_ext_s3 <- readRDS(paste0(estimands_path, "app_estimands_for_bootstrap_study3.RDS"))
  
  results_bootstrap_s1 <- get_bootstrap_results(study = s1, df = s1_data, nboot = 500, studyname = "study_1")
  setWinProgressBar(progression, l[13], title=paste(l[13],"% done"))
  
  results_bootstrap_s2 <- get_bootstrap_results(study = s1, df = s2_data, nboot = 500, studyname = "study_2")
  setWinProgressBar(progression, l[14], title=paste(l[14],"% done"))
  
  results_bootstrap_s3 <- get_bootstrap_results(study = s1, df = s3_data, nboot = 500, studyname = "study_3")
  setWinProgressBar(progression, l[15], title=paste(l[15],"% done, now 'only' save it!"))
  
  saveRDS(results_bootstrap_s1, paste0(estimands_path, "bootstrap_results_study1.RDS"))
  
  #################################################
  ########## Wrangling into nice format ###########
  #################################################
  
  ## Bind all results together 
  results_estimands <-
    rbind(results_app_ext_s1,
          results_app_ext_s2,
          results_app_ext_s3,
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
  results_estimands$iteration <- j
  results_estimands$seed <- seed_state[j]
  results_estimands <- results_estimands %>% mutate(`expected events` = n * prev)
  
  
  # Saving estimands
  saveRDS(results_estimands, file = paste0(estimands_path, "estimands_seed_", seed_state[j], ".Rds"))

}) # Close timing function

#################################
## Obtain performance measures ##
#################################











