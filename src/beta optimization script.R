####################################################
####################################################
####### TESTING OPTIMIZATION REGRESSION COEFFICIENTS
#### Libraries #####
library(tidyverse)
library(MASS)

source("./src/data generation functions.R")

###############
## Define R2 ##
############### 

set.seed(123)
# These values were determined in a different part of my study.
# They were numerically approximated given an AUC of 0.75 and event-rate (either 0.05, 0.2 or 0.5).
R2 <- c(approximate_R2(auc = AUC1, prev = 0.05)$R2.coxsnell,
        approximate_R2(auc = AUC1, prev = 0.2)$R2.coxsnell, 
        approximate_R2(auc = AUC1, prev = 0.5)$R2.coxsnell)


#############
## STUDY 1 ##
#############

############################################
## Generate data to develop parameters on ##
############################################

# setting n to develop parameters on.
n <- 30000 
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


############# OPTIMIZATION PROCESS ############

# Optimizing the dist_R2_prev to determine the beta coefficients for which the squared distance 
# between the observed and preferred values of R2 and prevalence are minimal. 
# par_i are just initial values for the coefficients. 

## Through trial and error, these were found to give the best results:
## Initial coefficients for event rate at 0.05, AUC = 0.75 (n_pred = 10)
## c(-3.37, 0.13)
## Initial coefficients for event rate at 0.2, AUC  = 0.75 (n_pred = 10)
## c(-1.66, 0.13)
## Initial coefficients for event rate at 0.5, AUC  = 0.75 (n_pred = 10)  
## c(0.01, 0.13))

##################
### PREV = .05 ###
##################
set.seed(123)
# Define initial parameters
par1_i <- c(-3.37, 0.13)

# Run optimization script n times, to get different values for regression parameters
system.time(results1 <- replicate(n = 20, 
                                 optim(par1_i, 
                                       dist_R2_prev,
                                       pref_R2 = R2[1], 
                                       pref_prev = 0.05,
                                       noise_contrib = "default",
                                       n_pred = n_pred,
                                       control = list(maxit = 1000)
                                 ),
                                 simplify = F))

# Use the checking function to get the AUC and prevalence, when using the combination of params
# for each repetition of the above.
results_check1 <- apply(sapply(results1, "[[", 1), 2, checking, noise_contrib = "default", n_pred = n_pred)

# Give a summary of the AUC values (Cstat) and prevalence (prev)
apply(results_check1, 1, summary)

# If it looks good, take the median of the params from results1
par1 <- apply(sapply(results1,  "[[", 1), 1, median) 

# Double check whether these values really come down to the preferred AUC value and prevalence
# AUC should be about 0.75 and prevalence 0.05
checking_val(par = par1, noise_contrib = "default", n_pred = n_pred)

##################
### PREV = 0.2 ###
##################

# Define initial parameters
par2_i <- c(-1.66, 0.13) # These have been found by trial and error

# Run optimization script n times, to get different values for regression parameters
system.time(results2 <- replicate(n = 20, 
                                 optim(par2_i, 
                                       dist_R2_prev,
                                       pref_R2 = R2[2],
                                       pref_prev = 0.2,
                                       noise_contrib = "default",
                                       n_pred = n_pred,
                                       control = list(maxit = 1000)
                                 ),
                                 simplify = F))

# Use the checking function to get the AUC and prevalence, when using the combination of params
# for each repetition of the above.
results_check2 <- apply(sapply(results2,  "[[", 1), 2, checking, noise_contrib = "default", n_pred = n_pred)

# Give a summary of the AUC values (Cstat) and prevalence (prev)
apply(results_check2, 1, summary)

# If it looks good, take the median of the params from results2
par2 <- apply(sapply(results2,  "[[", 1), 1, median) 

# Double check whether these values really come down to the preferred AUC value and prevalence
# AUC should be about 0.75 and prevalence 0.2
checking_val(par = par2, noise_contrib = "default", n_pred = n_pred)


##################
### PREV = 0.5 ###
##################

# Define initial parameters
par3_i <- c(0.01, 0.13)# These have been found by trial and error

# Run optimization script n times, to get different values for regression parameters
system.time(results3 <- replicate(n = 20, 
                                 optim(par3_i, 
                                       dist_R2_prev,
                                       pref_R2 = R2[3],
                                       pref_prev = 0.5,
                                       noise_contrib = "default",
                                       n_pred = n_pred,
                                       control = list(maxit = 1000)
                                 ),
                                 simplify = F))

# Use the checking function to get the AUC and prevalence, when using the combination of params
# for each repetition of the above.
results_check3 <- apply(sapply(results3, "[[", 1), 2, checking, noise_contrib = "default", n_pred = n_pred)

# Give a summary of the AUC values (Cstat) and prevalence (prev)
apply(results_check3, 1, summary)

# If it looks good, take the median of the params from results3
par3 <- apply(sapply(results3, "[[", 1), 1, median) 

# Double check whether these values really come down to the preferred AUC value and prevalence
# AUC should be about 0.75 and prevalence 0.5
checking_val(par = par3, noise_contrib = "default", n_pred = n_pred)



##############################################
################### STUDY 2 ##################
##############################################

## what are the relevant factors
## incluencing the parameters values?

dim2 <- c(6, 30, 60)
noise_2 <- c("none", "half")

# Put them in a matrix
s2_pred_set <- expand.grid(dim = dim2,
                           noise = noise_2
)

# To be filled with each iteration
s2_pred_set$par1 <- NA
s2_pred_set$par2 <- NA

# Which R2 belongs to the study:
s2_R2 <- R2[2]


####################
## GET PARAMETERS ##
####################

for (i in 1:nrow(s2_pred_set)){
n_pred <- s2_pred_set$dim[i]
noise_contrib <- s2_pred_set$noise[i]


# create covariance matrix to be used as input
sigma <- matrix(0.2, ncol = n_pred, nrow = n_pred) 
diag(sigma) <- 1 # set the diagonal to 1

# provide a vector of values for mu -> standard normal
mu <- c(rep(0, n_pred)) 

# create n_pred predictor columns
X <- mvrnorm(n = n, mu = mu, Sigma = sigma)

# Putting the above in a data matrix, including intercept
dm <- cbind(1, X) 


####################
## Validation set ##
####################


# setting n
n_val <- 100000 

# create n_pred predictor columns
X_val <- mvrnorm(n = n_val, mu = mu, Sigma = sigma)

# Putting the above in a data matrix, including intercept
dm_val <- cbind(1, X_val)


#####################
## OPTIMIZATION S2 ##
#####################

par_i <- c(-1.61, 0.2)

results <- replicate(n = 20, optim(par_i, 
                                        dist_R2_prev,
                                        pref_R2 = R2_s2, 
                                        pref_prev = 0.2,
                                        n_pred = n_pred,
                                        noise_contrib = noise_contrib,
                                        n_pred = n_pred,
                                        control = list(maxit = 1000)
                                  ),
                                  simplify = F)


results_check <- apply(sapply(results, "[[", 1), 2, checking, noise_contrib = noise_contrib, n_pred = n_pred)

# Give a summary of the AUC values (Cstat) and prevalence (prev)
apply(results_check, 1, summary)

# If it looks good, take the median of the params from results
par_results <- apply(sapply(results, "[[", 1), 1, median) 

# Double check whether these values really come down to the preferred AUC value and prevalence
# AUC should be about 0.75 and prevalence 0.2
checking_val(par = par_results, noise_contrib = noise_contrib, n_pred = n_pred)


s2_pred_set$par1[i] <- par_results[1] 
s2_pred_set$par2[i] <- par_results[2]

}

#############
## RESULTS ##
#############
#dim noise      par1        par2
#1   6  none -1.670588  0.20894297
#2  30  none -1.667019 -0.05344094
#3  60  none -1.652066  0.02734855
#4   6  half -1.671603  0.27602491
#5  30  half -1.623487  0.07781071
#6  60  half -1.657480 -0.04277436

#############
## STUDY 3 ##
#############

############################################
## Generate data to develop parameters on ##
############################################
set.seed(123)
# setting n to develop parameters on.
n <- 30000 

# How many predictors?
n_pred <- 20 

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

##################
## Optimization ##
##################
# Use some basis for the function below
# By trial and error:
par_i <- c(-1.64187732, -0.07185555)

results_s3 <- replicate(n = 20, 
                      optim(par_i, 
                            dist_R2_prev,
                            pref_R2 = R2[2],
                            pref_prev = 0.2, 
                            noise_contrib = "default",
                            n_pred = n_pred,
                            control = list(maxit = 1000)
                      ),
                      simplify = F)


# Use the checking function to get the AUC and prevalence, when using the combination of params
# for each repetition of the above.
results_check_s3 <- apply(sapply(results_s3, "[[", 1), 2, checking, noise_contrib = "default", n_pred = n_pred)

# Give a summary of the AUC values (Cstat) and prevalence (prev)
apply(results_check_s3, 1, summary)

# If it looks good, take the median of the params from results3
par_s3 <- apply(sapply(results_s3, "[[", 1), 1, median) 

# Double check whether these values really come down to the preferred AUC value and prevalence
# AUC should be about 0.75 and prevalence 0.5
checking_val(par = par_s3, noise_contrib = "default", n_pred = n_pred)

#############
## RESULTS ##
#############
## [1]  -1.64975835 -0.07082455


########################
########################
######## END SCRIPT