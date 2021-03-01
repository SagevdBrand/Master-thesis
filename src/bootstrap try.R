#########################
## BOOTSTRAP FUNCTIONS ##
#########################

## Create and load simulation data
source("./src/setup.R")
source("./src/estimand functions.R")
source("./src/data generation study 1.R") # Generate data, and save temporarily
names(s1_data) <- data_files # Change the names of each element in the list, to be sure it corresponds to the right scenario

df <- s1_data[[i]]
model <- s1[i, ]$model
pred_selection <- s1[i, ]$pred_selection

dgm_par <- c(s1[i, ]$par1, 
               rep(s1[i, ]$par2 * 3, round(0.3 * s1[i, ]$dim)),  # strong
               rep(s1[i, ]$par2,     round(0.5 * s1[i, ]$dim)),  # medium
               rep(s1[i, ]$par2 * 0, round(0.2 * s1[i, ]$dim)))  # noise
  

