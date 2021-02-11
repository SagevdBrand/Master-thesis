############################################################
######################### TESTING ##########################
############################################################
library(MASS)
library(tidyverse)

set.seed(123)

source("scripts/Scenarios.R") 
source("scripts/Data generation functions.R")

dist_R2_prev <- function(par,pref_R2, pref_prev, effects = c("half strong", "equal")){
  # par is a vector with initial guesses for both
  # the intercept and the beta1-3 coefficients.
  # However since beta1-3 have been restricted (half is strong, other half weak)
  # only one value is necessary.
  
  # Providing the beta0 and beta1-3 as specified in the par object. 
    if (effects == "half strong"){
    # Half of the predictors is strong
    dgm_par <- c(par[1], rep(par[2], round(n_pred/2)), rep((par[2]*3), round(n_pred/2))) 
    
  } else if (effects == "equal"){
    # All predictors are equally strong
    dgm_par <- c(par[1], rep(par[2], n_pred)) 
  } else {
    dgm_par <- c(par[1], rep(par[2], n_pred)) 
    print("No effects supplied, so all are set to be equally strong")
  }
  
  
  # Obtain values for y based on Bernoulli distribution, with input p
  #system.time(p <- exp(dm %*% dgm_par)/(1+exp(dm %*% dgm_par)))
  p <- 1/(1+exp(-dm %*% dgm_par))# SAME THING, JUST SLIGHTLY LESS COMPUTATION
  
  #system.time(y <- rbinom(length(p),1,p))
  y <-as.numeric(p>runif(length(p))) # SAME THING, JUST SLIGHTLY FASTER
  
  # Obtain observed values of c-statistic and 
  # average predicted probability of an event
  
  obs_R2 <- pseudo_Rsqrs(p = p, y = y) # obtain R2cs based on p and y
  obs_prev <- mean(y) # "Observed" prevalence
  
  # Sum of absolute differences of both values:
  #abs(obs_cstat-pref_cstat) + abs(obs_prev-pref_prev)
  (obs_R2-pref_R2)^2 + (obs_prev-pref_prev)^2# alternative, not sure which one is better
  
}


checking <- function(par, effects = c("half strong", "equal")){
  ## What do the observed prevalence and c-statistic look like?
  # Providing the beta0 and beta1-3 as specified in the par object. 
  if (effects == "half strong"){
    # Half of the predictors is strong
    dgm_par <- c(par[1], rep(par[2], round(n_pred/2)), rep((par[2]*3), round(n_pred/2))) 
    
  } else if (effects == "equal"){
    # All predictors are equally strong
    dgm_par <- c(par[1], rep(par[2], n_pred)) 
  } else {
    dgm_par <- c(par[1], rep(par[2], n_pred)) 
    print("No effects supplied, so all are set to be equally strong")
  }
  
  # Obtain values for y based on Bernoulli distribution, with input p
  p <- plogis(dm %*% dgm_par)
  y <- rbinom(length(p),1,p)
  
  # Obtain observed values
  #obs_cstat <- c_stat2(preds = p, outcome = y) # obtain c-statistic based on p and y
  obs_cstat <- fastAUC(p = p, y = y)
  obs_prev <- mean(y) # THE OBSERVED PREVALENCE IS NOT A FUNCTION OF JUST THE INTERCEPT
  c("cstat" = obs_cstat, "prev" = obs_prev)
}



n <- 30000 # setting n to develop betas on.
n_pred <- 10 # For the first and third scenario's at least
sigma <- matrix(0.2, ncol = n_pred, nrow = n_pred) # create covariance matrix to be used as input
diag(sigma) <- 1 # set the diagonal to 1
mu <- c(rep(0, n_pred)) # provide a vector of values for mu -> standard normal

X <- mvrnorm(n = n, mu = mu, Sigma = sigma) # create 3 predictor columns
dm <- cbind(1, X) # Putting the above in a data matrix, including intercept


  # Optimizing the dist_R2_prev to determine the beta coefficients for which the squared distance 
  # between the observed and preferred values of R2 and prevalence are minimal. 
  # par_i are just initial values for the coefficients. 
  
#### 30 repetitions #####
  system.time(results <- replicate(n = 10, 
                       optim(c(-2, 0.2), ## depends a lot on the input :(
                             dist_R2_prev,
                             pref_R2 = R2[1],
                             pref_prev = 0.05, 
                             effects = "equal",
                             control = list(maxit = 1000)),
                       simplify = F))
  # Computation time is about 1 minute
  
  # Getting the estimated optimal coefficients from each list.
  # Thereafter, use the checking function to determine the c-statistic
  # and prevalence belonging to the estimated betas
  results_check <- apply(sapply(results, '[[', 1), 2, checking, effects = "equal")
  
  # What do the results look like?
  apply(results_check, 1, summary)
  
  # Taking the median of the coefficients,
  # to reduce influence of potential outliers 
  par <- apply(sapply(results, '[[', 1), 1, median) 

  
## Create validation dataset ##
  n_val <- 100000 # setting n
  X_val <- mvrnorm(n = n_val, mu = mu, Sigma = sigma) # create 10 predictor columns
  dm_val <- cbind(1, X_val) # Putting the above in a data matrix, including intercept

  
# Same idea as the checking function, however, using the validation data:
checking_val <- function(par, effects = c("half strong", "equal")){
    
    # Providing the beta0 and beta1-3 as specified in the par object. 
    if (effects == "half strong"){
      # Half of the predictors is strong
      dgm_par_val <- c(par[1], rep(par[2], round(n_pred/2)), rep((par[2]*3), round(n_pred/2))) 
      
    } else if (effects == "equal"){
      # All predictors are equally strong
      dgm_par_val <- c(par[1], rep(par[2], n_pred)) 
    } else {
      dgm_par_val <- c(par[1], rep(par[2], n_pred)) 
      print("No effects supplied, so all are set to be equally strong")
    }
    
    p_val <- plogis(dm_val %*% dgm_par_val)
    y_val <- rbinom(length(p_val),1,p_val)
    
    # Obtain observed values
    #obs_cstat <- c_stat2(preds = p_val, outcome = y_val) # obtain c-statistic based on p and y
    obs_cstat <- fastAUC(p = p_val, y = y_val)
    obs_prev <- mean(y_val) # THE OBSERVED PREVALENCE IS NOT A FUNCTION OF JUST THE INTERCEPT
    c("cstat" = obs_cstat, "prev" = obs_prev)
  }
  
## Use the checking function for the validation set:
checking_val(par = par, effects = "equal")

