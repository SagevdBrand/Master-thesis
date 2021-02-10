################################
################################
################################
### Try simulate data scenario 1 

######## Setting up #########
library(tidyverse)
library(MASS)

set.seed(123)

source("scripts/Scenarios.R")
source("scripts/Data generation functions.R")

system.time(s1_data <- generate_data(s1))

lapply(lapply(s1_data,'[[', 11), mean) # check whether the prevalence is somewhat okay

for(i in 1:length(s1_data)) {
  s1_path <- paste0("Data/simulation data/scenario 1/", i) # create right path
  saveRDS(assign(paste0("s1_", i), s1_data[[i]]), file = paste0(s1_path,"/s1_", i)) # add name of file to path
}

