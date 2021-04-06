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
##      [X] CHECK FOR CART WHETHER ANY SPLITS WERE MADE
##            [ ] RETURN HIGHEST VALUE FOR CALIBRATION SLOPE WITHIN THAT SCENARIO 

## [X] BUILD ERROR HANDLING FOR SIMULATION RUNS!!

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
## [X] ADD MAPE 
## [X] ADD RMSPE
## [X] BUILD IN ERROR HANDLING AS SPECIFIED IN PROTOCOL!
##      [X] IF ERROR OCCURS, MAKE SURE IT CONTINUES AND JUST RETURNS AN ERROR WITHIN THE RESULTS VECTOR
##      [X] CHECK FOR SEPARATION ISSUES
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

run_id <- 999999
set.seed(999)
## Create and load simulation data
s1_data <- generate_data(s1, validation = FALSE)
s2_data <- generate_data(s2, validation = FALSE)
s3_data <- generate_data(s3, validation = FALSE)


## create one test set for all studies, shortes scenarios:
df <- s1_data[1]#c(s1_data[c(7)], s2_data[c(10)], s3_data[c(16)])
df_val <- s1_val_data[1]#, s2_val_data[c(10)], s3_val_data[c(16)])
study <- s1[1,]#rbind(s1[c(7),], s2[c(10),], s3[c(16),])
study$scenario <- "Scenario 1" #c("Scenario 1", "Scenario 2", "Scenario 3")
model <- "Firth"
pred_selection <- "none"
dgm_par <-
  c(study[1, ]$par1, 
    rep(study[1, ]$par2 * 3, round(0.3 * study[1, ]$dim)), 
    rep(study[1, ]$par2, round(0.5 *  study[1, ]$dim)), 
    rep(study[1, ]$par2 * 0, round(0.2 * study[1, ]$dim)))

app_estimands_for_bootstrap_study1 <- readRDS("~/GitHub/Master-thesis/results/output/estimands/old/app_estimands_for_bootstrap_study1.RDS")
app_preds_for_bootstrap_study1 <- readRDS("~/GitHub/Master-thesis/results/output/estimands/old/app_preds_for_bootstrap_study1.RDS")
p_app_study1 <- app_preds_for_bootstrap_study1
results_app_ext_study1 <- app_estimands_for_bootstrap_study1
nboot <- 10

# Test apparent & external:
system.time(results_app_ext_test <- get_app_ext_results(study = study, df = df, df_val = df_val, studyname = "test"))
results_app_ext_test$scenario <- c(paste("Scenario", seq(1:3)))
saveRDS(results_app_ext_test, paste0(app_ext_path, "app_ext_estimands_test_seed_", run_id, ".Rds")) 
## Time: 670.2
## Time: 643
## Time: 542.41
## Time (num.trees@250): 529.77
## Time (num.trees@250): 498.73 

system.time(test_10cv_short <- get_cv_results(study = study, df = df, V = 3, studyname = "test"))
saveRDS(test_10cv_short, paste0(cv_10_fold_path, "10fcv_estimands_test_seed_", run_id, ".Rds")) 
saveRDS(test_10cv_short, paste0(cv_5_fold_path, "5fcv_estimands_test_seed_", run_id, ".Rds")) 
## Time (num.trees@250): 102.55

system.time(test_boot_short <- get_bootstrap_results(study = study, df = df, nboot = 50, studyname = "study1"))
saveRDS(test_boot_short, paste0(bootstrap_path, "bootstrap_estimands_test_seed_", run_id, ".Rds")) 
## Time (num.trees@250): 134.11 | 134.95


system.time(test_10x10_short <- get_10x10_results(study = study, df = df, V = 3, studyname = "test"))
saveRDS(test_10x10_short, paste0(cv_10x10_fold_path, "10x10fcv_estimands_test_seed_", run_id, ".Rds")) 
## Time (num.trees@250): 1100/12
