##################################
##################################
##################################
## GET ESTIMANDS SCENARIO 1

############ TO FIX/DO #################

## [ ] FUNCTION FOR EXTERNAL VALIDATION OF ALL MODELS USED IN SCENARIO - FOR PERFORMANCE MEASURES
## [ ] BOOTSTRAP ESTIMAND FUNCTION
## [ ] 10X10 CV FUNCTION 
## [X] LOOCV FUNCTION or not?
## [ ] ADD OPTIONS FOR OTHER MODELS
## [ ] MAKE SURE IT WORKS WITH MORE OR LESS PREDICTORS
## [ ] BUILD IN CHECKS AS SPECIFIED IN PROTOCOL!
## [ ] IF ERROR OCCURS, MAKE SURE IT CONTIUES AND JUST RETURNS AN ERROR WITHIN THE RESULTS VECTOR
## [X] OUTPUT SEED!!!!
## [X] INTEGRATE CREATE DATA FUNCTION SO WE DONT HAVE THE SAME DATA ALL THE TIME :')
## [ ] CREATE DIFFERENT DGM-PAR IN DIFFERENT STUDIES

########################################


###############################################################################
############ Load libraries, validation data and scenario settings ############
###############################################################################

## Libraries, file paths and functions
source("scripts/setup.R")
source("scripts/estimand functions.R")

## Load scenario settings
s1 <- read_rds(paste0(scenario_1_settings,"s1.Rds"))

## Load validation data
val_data_files <- list.files(path = scenario_1_val_data, recursive = T, full.names = F)
val_df <- lapply(paste0(scenario_1_val_data,val_data_files),readRDS,.GlobalEnv)
names(val_df) <- val_data_files

#########################################################
############## LET THE SIMULATION COMMENCE ##############
#########################################################

# Reproduction seed
set.seed(123)

# Store seed values
n_sim <- 2 # how many iterations?
state <- floor(runif(n_sim, 0, 10000)) # Create a vector of seed states

system.time(for(j in 1:n_sim){
  set.seed(state[j]) # for each run the next value in the state vector will be chosen (and saved!)
  
  ## Create and load simulation data
  source("scripts/data generation scenario 1.R") # Generate data, and save temporarily
  data_files <- list.files(path = scenario_1_data, recursive = T, full.names = F) # get the data path 
  df <- lapply(paste0(scenario_1_data,data_files),readRDS,.GlobalEnv) # Load the data
  names(df) <- data_files # Change the names of each element in the list, to be sure it corresponds to the right scenario
  
  ## Obtain apparent estimands ##
  results_app <- get_app_results(scenario = s1, df = df)
  
  ## Obtain internal validation estimands ##
  # 10 fold cross-validation
  results_10_cv <- get_cv_results(scenario = s1, df = df, V = 10)
  
  # 5 fold cross-validation
  results_5_cv <- get_cv_results(scenario = s1, df = df, V = 5)
  
  # 10X10 fold cross-validation 
  results_10x10_cv <- get_10x10_results(scenario = s1, df = df, V = 10)
  
  # Bootstrap 3 varieties in one go
  
  
  ## Obtain external validation estimands
  # Exteral validation
  

  ## Make a vector of all results + state
  for(i in 1:length(results_app)){ # For however many scenarios there are within the scenario 
  results_lists <- list(results_app, results_10_cv, results_5_cv, results_10x10_cv) ## ADD OTHER RESULTS FROM VALIDATION APPROACHES
  saveRDS(
    object = assign(paste0("s1_results", i), # create an object with name "s1_results_i
                    c("seed_state" = state[j], #Fill it with the seed state, 
                      unlist(
                        lapply(results_lists, "[[", i) # and results belonging to the subscenario
                      ))),
    file = paste0(s1_output, "s1_", i, "_", state[j], ".Rds")
  )  # close saveRDS
  } # close saving for loop
} # close simulation for loop
)

######################
## Obtain estimands ##
######################



## Obtain external validation estimands based on models fitted


#################################
## Obtain performance measures ##
#################################











