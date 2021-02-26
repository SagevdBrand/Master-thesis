################################
################################
################################
### Try simulate data scenario 1 

######## Setting up #########
## set.seed(123)
source("./src/setup.R")
source("./src/data generation functions.R")

s1 <- read_rds(study_1_settings)
system.time(s1_data <- generate_data(s1, validation = FALSE))
lapply(lapply(s1_data,'[[', 11), mean) # check whether the prevalence is somewhat okay 
# - Get rid after fixing error check

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


