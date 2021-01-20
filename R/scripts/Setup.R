################# INITIAL SETUP ##################
## Script to define paths to the simulated data ##
##################################################
library(tidyverse)
library(glue)

datasets <- list.dirs(path = "data", recursive = F, full.names = F)
data_files <- glue(datasets, ".csv")
data_dir <- glue("data/", datasets)


################
## Vocabulary ##
################

# Within each of the scenarios, the specific settings used are abbreviated.
# These abbreviations mean the following:

########################
###### Scenario 1 ######
########################

########################
###### Scenario 2 ######
########################

########################
###### Scenario 3 ######
########################

########################
###### Scenario 4 ######
########################



