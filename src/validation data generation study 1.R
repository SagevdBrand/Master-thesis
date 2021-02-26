#######################################
#######################################
#######################################
### Simulate validation data study 1 

######## Setting up #########
set.seed(123)
source("./src/setup.R")
source("./src/data generation functions.R")

s1 <- read_rds(study_1_settings)
system.time(s1_val_data <- generate_data(s1, validation = TRUE)) 
lapply(lapply(s1_val_data,'[[', 11), mean)

for(i in 1:length(s1_val_data)) {
  ind <- letters[1:length(s1_val_data)]
  saveRDS(assign(paste0("s1_", i), s1_val_data[[i]]), file = paste0(study_1_val_data, "s1_val_data_", ind[i],".Rds")) # add name of file to path
  rm(list = ls(pattern = paste0("s1_",i)))
  }

