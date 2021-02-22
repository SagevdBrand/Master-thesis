################################
################################
################################
### Try simulate data scenario 1 

######## Setting up #########
## set.seed(123)
source("scripts/Setup.R")
source("scripts/Data generation functions.R")

s1 <- read_rds(paste0(scenario_1_settings,"s1.Rds"))
system.time(s1_data <- generate_data(s1, validation = FALSE))
lapply(lapply(s1_data,'[[', 11), mean) # check whether the prevalence is somewhat okay


## Maybe not save it like this though?
for(i in 1:length(s1_data)) {
  ind <- letters[1:length(s1_data)]
  saveRDS(assign(paste0("s1_", i), s1_data[[i]]), file = paste0(scenario_1_data, "s1_data_", ind[i],".Rds")) # add name of file to path
  rm(list = ls(pattern = paste0("s1_",i)))
  }


