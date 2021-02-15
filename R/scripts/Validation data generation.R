#######################################
#######################################
#######################################
### Simulate validation data scenario 1 

######## Setting up #########
set.seed(123)
source("scripts/Setup.R")
source("scripts/Data generation functions.R")

s1 <- read_rds(paste0(scenario_1_settings,"s1.Rds"))
system.time(s1_val_data <- generate_data(s1, validation = TRUE)) 
lapply(lapply(s1_val_data,'[[', 11), mean)

for(i in 1:length(s1_val_data)) {
  ind <- letters[1:length(s1_val_data)]
  saveRDS(assign(paste0("s1_", i), s1_val_data[[i]]), file = paste0(scenario_1_validation_data, "s1_val_data_", ind[i],".Rds")) # add name of file to path
}

