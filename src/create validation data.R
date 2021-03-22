source("./src/study scenarios.R")
source("./src/data generation functions.R")
source("./src/setup.R")

set.seed(123)

## Create and load validation data
s1_val_data <- generate_data(s1, validation = TRUE)
s2_val_data <- generate_data(s2, validation = TRUE)
s3_val_data <- generate_data(s3, validation = TRUE)

## Save results
saveRDS(s1_val_data, study_1_val_data)
saveRDS(s2_val_data, study_2_val_data)
saveRDS(s3_val_data, study_3_val_data)


######################################
######################################
##### END SCRIPT