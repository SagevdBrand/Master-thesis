##################################
##################################
##################################
## GET ESTIMANDS SCENARIO 1

############ TO FIX/DO #################

## [ ] FUNCTION FOR EXTERNAL VALIDATION OF ALL MODELS USED IN SCENARIO - FOR PERFORMANCE MEASURES
## [ ] BOOTSTRAP ESTIMAND FUNCTION
## [ ] ADD OPTIONS FOR OTHER MODELS
##      [ ] RIDGE -> CODE IS READY, ONLY NEEDS IMPLEMENTATION
##      [ ] LASSO -> CODE IS READY ONLY NEEDS IMPLEMENTATION
##      [ ] MACHINE LEARNING -> CODE IS BEING DEVELOPED, BUT WAIT FOR FEEDBACK
## [ ] CREATE DIFFERENT DGM-PAR IN DIFFERENT STUDIES

## [x] IMPLEMENT STUDY 1 FOR ESTIMANDS
## [ ] IMPLEMENT STUDY 2 FOR ESTIMANDS
## [ ] IMPLEMENT STUDY 3 FOR ESTIMANDS
## [ ] IMPLEMENT STUDY 4 FOR ESTIMANDS

## [ ] IF ERROR OCCURS, MAKE SURE IT CONTIUES AND JUST RETURNS AN ERROR WITHIN THE RESULTS VECTOR
## [ ] BUILD IN ERROR HANDLING AS SPECIFIED IN PROTOCOL!

##      [ ] CHECK FOR EVENTS IN RESAMPLING
##      [ ] CHECK FOR VAR(Y) == 0 |SUM(Y) < 8 | N - SUM(Y) < 8  FOR LASSO AND REGRESSION
##      [ ] CHECK FOR VAR(LP) == 0
##      [ ] RETURN HIGHEST VALUE FOR CALIBRATION SLOPE WITHIN THAT SCENARIO
##      [ ] 


## DONE:
## [X] CHECK FOR EVENTS IN DEVELOPMENT DATASET
## [X] RETURN VECTOR OF NA/ERROR AS RESULT
## [x] 10 & 5 FOLD CV
## [x] 10X10 FOLD CV
## [x] PIPELINE FOR EASY GENERATION OF DATA
## [x] PIPELINE FOR STORING ESTIMAND RESULTS
## [X] OUTPUT SEED
## [X] INTEGRATE CREATE DATA FUNCTION SO WE DONT HAVE THE SAME DATA ALL THE TIME :')
## [X] LOOCV FUNCTION or not?
## [X] MAKE SURE IT WORKS WITH MORE OR LESS PREDICTORS

########################################


############################################################################
############ Load libraries, validation data and study settings ############
############################################################################

## Libraries, file paths and functions
source("./src/setup.R")
source("./src/estimand functions.R")

## Load scenario settings
s1 <- read_rds(study_1_settings)

## Load validation data
source("./src/validation data generation study 1.R") # Have it run at least once, so that there are files in the folder.
val_data_files <- list.files(path = study_1_val_data, recursive = T, full.names = F)

source("./src/validation data generation study 1.R") # Have it run at least once, so that there are files in the folder.
data_files <- list.files(path = study_1_data, recursive = T, full.names = F) # get the data names

#########################################################
############## LET THE SIMULATION COMMENCE ##############
#########################################################

# Reproduction seed
set.seed(123)

# Store seed values
n_sim <- 1 # how many iterations?
state <- c(1:5000) # Create a vector of seed states

system.time({for(j in 1:n_sim){
  set.seed(state[j]) # for each run the next value in the state vector will be chosen (and saved!)
  
  ## Create and load simulation data
  source("./src/data generation study 1.R") # Generate data, and save temporarily
  names(s1_data) <- data_files # Change the names of each element in the list, to be sure it corresponds to the right scenario
  
  ## Create and load validation data
  source("./src/validation data generation study 1.R") # Generate data, and save temporarily
  names(s1_val_data) <- val_data_files
  
  ## Obtain apparent estimands ##
  results_app <- get_app_results(study = s1, df = s1_data)
  
  ## Obtain internal validation estimands ##
  # 10 fold cross-validation
  results_10_cv <- get_cv_results(study = s1, df = s1_data, V = 10)
  
  # 5 fold cross-validation
  results_5_cv <- get_cv_results(study = s1, df = s1_data, V = 5)
  
  # 10X10 fold cross-validation 
  results_10x10_cv <- get_10x10_results(study = s1, df = s1_data, V = 10)
  
  # Bootstrap 3 varieties in one go
  
  
  ## Obtain external validation estimands
  
  
  # External validation
  

  ## Make a list of all results 
  results_lists <- list("results_app" = results_app,
                        "results_10_cv" = results_10_cv,
                        "results_5_cv"  = results_5_cv,
                        "results_10x10_cv"= results_10x10_cv) ## ADD OTHER RESULTS FROM VALIDATION APPROACHES
  
  # SAVING THE ESTIMANDS
  for(i in 1:nrow(s1)){ # For however many scenarios there are within the study 
  
  # Save an RDS object in the shape of a vector
  saveRDS(
    # create an object with name "s1_results_i
    object = assign(paste0("s1_estimands", i), 
                    # Which contains the seed state, 
                    c("seed_state" = state[j], 
                      # and results belonging to the scenario
                      unlist(
                        lapply(results_lists, "[[", i) 
                      ))), # Closes assign()
    # Save the vector as an Rds file in the correct folder
    file = paste0(s1_estimands, "s1_scen", i, "_seed", state[j], ".Rds")
  )  # close saveRDS
    
  } # close saving for loop
  
} # close simulation for loop
  
  }) # Close timing function



## Obtain external validation estimands based on models fitted


#################################
## Obtain performance measures ##
#################################











