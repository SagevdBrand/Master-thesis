####################################################
####################################################
####### TESTING OPTIMIZATION REGRESSION COEFFICIENTS
set.seed(123)

#### Libraries #####
library(tidyverse)
library(MASS)

#### Defining functions ####

###  Determines the squared distance between the preferred and observed R2 and prevalence.
dist_R2_prev <- function(par, pref_R2, pref_prev){
  # par is a vector with initial guesses for both
  # the intercept and the regression coefficients.
  # However since regression have been restricted
  # only one value is necessary.
  
  # Providing the intercept and regression coefficients as specified in the par object. 
  dgm_par <-
    c(par[1], 
      rep(par[2] * 3, round(0.3 * n_pred)),  # strong
      rep(par[2], round(0.5 *  n_pred)),  # weaker
      rep(par[2] * 0, round(0.2 * n_pred))) # noise
  
  # Obtain values for y based on Bernoulli distribution, with input p
  p <- 1/(1+exp(-dm %*% dgm_par))
  y <- as.numeric(p > runif(length(p))) 
  
  # Obtain observed values of c-statistic and 
  # average predicted probability of an event
  
  obs_R2 <- pseudo_Rsqrs(p = p, y = y) # obtain R2cs based on p and y
  obs_prev <- mean(y) # "Observed" prevalence
  
  # Sum of squared differences of both values:
  (obs_R2-pref_R2)^2 + (obs_prev-pref_prev)^2
  
}

#### Checks whether the optimized parameters look like we want them to:
checking <- function(par){
  ## What do the observed prevalence and c-statistic look like?
  # Providing the  and regression coefficients as specified in the par object. 
  
  # Defining candidate predictors
  dgm_par <-
    c(par[1], 
      rep(par[2] * 3, round(0.3 * n_pred)), # Strong
      rep(par[2], round(0.5 *  n_pred)),  # Weaker
      rep(par[2] * 0, round(0.2 * n_pred))) # Noise
  
  # Obtain values for y based on uniform distribtuin, with input p
  p <- 1/(1+exp(-dm %*% dgm_par))
  y <- as.numeric(p > runif(length(p))) 
  
  # Obtain observed values
  obs_cstat <- fastAUC(p = p, y = y)
  obs_prev <- mean(y) 
  c("cstat" = obs_cstat, "prev" = obs_prev)
}

### Same idea as the checking function, however, using the validation data:
checking_val <- function(par){
  
  # Providing the intercept and regression coefficients as specified in the par object. 
  dgm_par_val <-
    c(par[1], 
      rep(par[2] * 3, round(0.3 * n_pred)), 
      rep(par[2], round(0.5 *  n_pred)), 
      rep(par[2] * 0, round(0.2 * n_pred)))
  
  
  # Obtain values for y based on uniform distribution, with input p
  p_val <- 1/(1+exp(-dm_val %*% dgm_par_val))
  y_val <- as.numeric(p_val > runif(length(p_val)))
  
  # Obtain observed values
  #obs_cstat <- c_stat2(preds = p_val, outcome = y_val) # obtain c-statistic based on p and y
  obs_cstat <- fastAUC(p = p_val, y = y_val)
  obs_prev <- mean(y_val)
  
  # return results
  c("cstat" = obs_cstat, "prev" = obs_prev)
}

## R^2 Cox-Snell
pseudo_Rsqrs <- function(p, y){ 
  
  .LL <- function(p, y){
    sum(y*log(p)+(1-y)*log(1-p))
  }
  
  LL_fit  <- .LL(p=p, y=y) 
  LL_null <- .LL(p=mean(y), y=y)
  
  cox <- 1-exp(-(LL_fit-LL_null)*2/length(y)) 
  cox_max <- 1 - exp(2 * length(y) ^ (-1) * LL_null)
  c("cox"=cox)
  
}

## AUC (C-statistic)
fastAUC <- function(p, y) {
  x1 = p[y==1]; n1 = length(x1); 
  x2 = p[y==0]; n2 = length(x2);
  r = rank(c(x1,x2))  
  auc = (sum(r[1:n1]) - n1*(n1+1)/2) / n1 / n2
  return(auc)
}

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

###############
## Define R2 ##
############### 
# These values were determined in a different part of my study.
# They were numerically approximated given an AUC of 0.75 and event-rate (either 0.05, 0.2 or 0.5).
R2 <- c(0.04131983, 0.12658143, 0.18407039) 

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

# Define initial parameters
par1_i <- c(-3.37, 0.13)

# Run optimization script n times, to get different values for regression parameters
system.time(results1 <- replicate(n = 5, 
                                 optim(par1_i, 
                                       dist_R2_prev,
                                       pref_R2 = R2[1], 
                                       pref_prev = 0.05, 
                                       control = list(maxit = 1000)
                                 ),
                                 simplify = F))

# Use the checking function to get the AUC and prevalence, when using the combination of params
# for each repetition of the above.
results_check1 <- apply(sapply(results1, '[[', 1), 2, checking)

# Give a summary of the AUC values (Cstat) and prevalence (prev)
apply(results_check1, 1, summary)

# If it looks good, take the median of the params from results1
par1 <- apply(sapply(results1, '[[', 1), 1, median) 

# Double check whether these values really come down to the preferred AUC value and prevalence
# AUC should be about 0.75 and prevalence 0.05
checking_val(par = par1)

##################
### PREV = 0.2 ###
##################

# Define initial parameters
par2_i <- c(-1.66, 0.13) # These have been found by trial and error

# Run optimization script n times, to get different values for regression parameters
system.time(results2 <- replicate(n = 5, 
                                 optim(par2_i, 
                                       dist_R2_prev,
                                       pref_R2 = R2[2],
                                       pref_prev = 0.2, 
                                       control = list(maxit = 1000)
                                 ),
                                 simplify = F))

# Use the checking function to get the AUC and prevalence, when using the combination of params
# for each repetition of the above.
results_check2 <- apply(sapply(results2, '[[', 1), 2, checking)

# Give a summary of the AUC values (Cstat) and prevalence (prev)
apply(results_check2, 1, summary)

# If it looks good, take the median of the params from results2
par2 <- apply(sapply(results2, '[[', 1), 1, median) 

# Double check whether these values really come down to the preferred AUC value and prevalence
# AUC should be about 0.75 and prevalence 0.2
checking_val(par = par2)


##################
### PREV = 0.5 ###
##################

# Define initial parameters
par3_i <- c(0.01, 0.13)# These have been found by trial and error

# Run optimization script n times, to get different values for regression parameters
system.time(results3 <- replicate(n = 5, 
                                 optim(par3_i, 
                                       dist_R2_prev,
                                       pref_R2 = R2[3],
                                       pref_prev = 0.5, 
                                       control = list(maxit = 1000)
                                 ),
                                 simplify = F))

# Use the checking function to get the AUC and prevalence, when using the combination of params
# for each repetition of the above.
results_check3 <- apply(sapply(results3, '[[', 1), 2, checking)

# Give a summary of the AUC values (Cstat) and prevalence (prev)
apply(results_check3, 1, summary)

# If it looks good, take the median of the params from results3
par3 <- apply(sapply(results3, '[[', 1), 1, median) 

# Double check whether these values really come down to the preferred AUC value and prevalence
# AUC should be about 0.75 and prevalence 0.5
checking_val(par = par3)



######################################################################################################
######################################################################################################
######## END SCRIPT

  


