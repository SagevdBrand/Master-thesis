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
