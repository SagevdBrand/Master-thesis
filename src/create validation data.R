source("./src/data generation functions.R")
source("./src/setup.R")

set.seed(123)

studies <- readRDS(paste0(setting_path, "studies.RDS"))
## Create and load validation data
for (i in 1:nrow(studies)){
val_data <- generate_data(studies[i,], validation = TRUE)
saveRDS(val_data, 
        paste0(
        val_data_path, paste("val_data", i, sep = "_"), ".Rds"))
}


######################################
######################################
##### END SCRIPT