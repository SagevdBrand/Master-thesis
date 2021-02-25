############################################
############################################
### Generate data for optimization

set.seed(123)

############################################
## Generate data to develop parameters on ##
############################################


# setting n to develop parameters on.
n <- 30000 

# How many predictors?
n_pred <- 10 

# create covariance matrix to be used as input
sigma <- matrix(0.2, ncol = n_pred, nrow = n_pred) 
diag(sigma) <- 1 # set the diagonal to 1

# provide a vector of values for mu -> standard normal
mu <- c(rep(0, n_pred)) 

# create n_pred predictor columns
X <- mvrnorm(n = n, mu = mu, Sigma = sigma)

# Putting the above in a data matrix, including intercept
dm <- cbind(1, X) 

###############################
## Create validation dataset ##
###############################

# setting n
n_val <- 100000 

# create n_pred predictor columns
X_val <- mvrnorm(n = n_val, mu = mu, Sigma = sigma)

# Putting the above in a data matrix, including intercept
dm_val <- cbind(1, X_val)

##########################
#### Save as raw data ####
##########################

saveRDS(object = dm, file = "./data/raw/dm_development.RDS")
saveRDS(object = dm_val, file = "./data/raw/dm_validation.RDS")

#######################################################
#######################################################
#### END SCRIPT
