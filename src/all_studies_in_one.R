## All studies in one frame
############################################################################
############ Load libraries, validation data and study settings ############
############################################################################
## Libraries, file paths and functions
source("./src/setup.R")
source("./src/estimand functions.R")
source("./src/data generation functions.R")

## Load scenario settings
s1 <- readRDS(study_1_settings)
s2 <- readRDS(study_2_settings)
s3 <- readRDS(study_3_settings)

s1$study <- "Study_1"
s2$study <- "Study_2"
s3$study <- "Study_3"

## Load validation data
s1_val_data <- readRDS(study_1_val_data)
s2_val_data <- readRDS(study_2_val_data)
s3_val_data <- readRDS(study_3_val_data)

studies <- rbind(s1,s2,s3)
saveRDS(studies, paste0(setting_path, "studies.RDS"))


validation_data <- list(s1_val_data, s2_val_data, s3_val_data)
validation_data <- flatten(validation_data)
str(validation_data)
