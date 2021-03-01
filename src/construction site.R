##########################################################
################## CONSTRUCTION SITE #####################
##########################################################

## Get libraries, paths and a dataset:

## Libraries, file paths and functions
source("./src/setup.R")
source("./src/estimand functions.R")

## Load scenario settings
s1 <- read_rds(study_1_settings)
study <- s1

### Load data as df ###
# Have it run at least once, so that there are files in the folder.
source("./src/validation data generation study 1.R") 

# get the data names
data_files <- list.files(path = study_1_data, recursive = T, full.names = F)

# Generate data, and save temporarily
source("./src/data generation study 1.R")

# Change the names of each element in the list
# to be sure it corresponds to the right scenario
names(s1_data) <- data_files
df <- s1_data[[1]]

#### TESTING SITE FOR ERROR HANDLING
# Set i and V
i <- 1
V <- 5

#Preset which model, predictor selection methods and dgm
model <- study[i, ]$model
pred_selection <- study[i, ]$pred_selection
dgm_par <- c(study[i, ]$par1, 
             rep(study[i, ]$par2 * 3, round(0.3 * study[i, ]$dim)),  # strong
             rep(study[i, ]$par2,     round(0.5 * study[i, ]$dim)),  # weak
             rep(study[i, ]$par2 * 0, round(0.2 * study[i, ]$dim)))  # noise

get_cv_estimands(df = df[[i]], 
                 model = model,
                 dgm_par = dgm_par,
                 pred_selection = pred_selection,
                 V = V,
                 x10 = FALSE)

## BY MAARTEN:
# tryCatch.W.E <- function(expr){
#   W <- NULL
#   w.handler <- function(w){
#     W <<- w
#     invokeRestart("muffleWarning")}
#   list(value = withCallingHandlers(tryCatch(expr, error = function(e) e), warning = w.handler),warning = W)
# }


## FIRST CHECK ##
## NO EVENTS SAMPLED IN GENERATED SIMULATION DATA
## ESPECIALLY PROBLEM FOR LASSO AND RIDGE?
######## Setting up #########
## set.seed(123)
source("./src/data generation functions.R")

s1 <- read_rds(study_1_settings)
system.time(s1_data <- generate_data(s1, validation = FALSE))

check_events <- lapply(lapply(s1_data,"[[", "y"), var) 
check_events[[3]] <- 0 

lapply(check_events, function(x) {if(x == 0) warning("No events sampled")})



## For whatever many scenarios are present, generate a dataset
for(i in 1:length(s1_data)) {
  # Instead of saving a number, save by letters, to keep the correct order of files
  ind <- letters[1:length(s1_data)]
  
  # Save an rds object for each generated dataset
  # Names of the objects are "s1_data_a"
  saveRDS(object = assign(paste0("s1_", i), s1_data[[i]]),
          file = paste0(study_1_data, "s1_data_", ind[i],".Rds")) 
  
  # Remove from environment to avoid clutter
  rm(list = ls(pattern = paste0("s1_",i)))
  
}



  