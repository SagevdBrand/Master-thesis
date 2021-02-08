#############################################################
### Optimizing coefficients for C-statistic and Prevalence ##
#############################################################

############# To do:
## Depending on decisions regarding to predictors, effects of predictors etc.
## Expand script to get betas for all situation.


##################
### Setting up ###
##################

library(MASS)
library(tidyverse)

set.seed(123)

source("scripts/Scenarios.R") 
source("scripts/Data generation functions.R")


#######################
### data generation ###
#######################

n <- 30000 # setting n to develop betas on.
n_pred <- 10 # For the first 2 scenario's at least
sigma <- matrix(0.2, ncol = n_pred, nrow = n_pred) # create covariance matrix to be used as input
diag(sigma) <- 1 # set the diagonal to 1
mu <- c(rep(0, n_pred)) # provide a vector of values for mu -> standard normal

X <- mvrnorm(n = n, mu = mu, Sigma = sigma) # create 3 predictor columns
dm <- cbind(1, X) # Putting the above in a data matrix, including intercept

# There is only in scenario 1 a change in prevalence. 
# Therefore, only 2 runs are needed. 

system.time(beta_prev_0.05 <- optim_beta(prev = 0.05, R2 = R2[1]))
system.time(beta_prev_0.2 <- optim_beta(prev = 0.2, R2 = R2[2]))

#####################################################################
######### Validate results on independent validation set ############
#####################################################################

## Create validation dataset ##
n_val <- 100000 # setting n
X_val <- mvrnorm(n = n_val, mu = mu, Sigma = sigma) # create 10 predictor columns
dm_val <- cbind(1, X_val) # Putting the above in a data matrix, including intercept

## Use the checking function for the validation set:
checking_val(par = beta_prev_0.05)
checking_val(par = beta_prev_0.2)

## Save beta's in an object!
saveRDS(beta_prev_0.05, file = "Data/Simulation settings/Scenario 1/Betas/Betas_prev_0.05_halfstrong.Rds")
saveRDS(beta_prev_0.2, file = "Data/Simulation settings/Scenario 1/Betas/Betas_prev_0.2_halfstrong.Rds")

####################################
####################################
####################################
### END SCRIPT
