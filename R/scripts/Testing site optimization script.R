############################################################
######################### TESTING ##########################
############################################################
set.seed(123)
source("scripts/libraries.R")
source("scripts/Scenarios.R") 
source("scripts/Data generation functions.R")

dist_AUC_prev <- function(par, pref_auc, pref_prev){
  # par is a vector with initial guesses for both
  # the intercept and the beta1-3 coefficients.
  # However since beta1-3 have been restricted (half is strong, other half weak)
  # only one value is necessary.
  
  # Providing the beta0 and beta1-3 as specified in the par object. 
  dgm_par <-
    c(par[1], 
      rep(par[2] * 3, round(0.4 * n_pred)),  # strong
      rep(par[2], round(0.4 *  n_pred)),  # medium
      rep(par[2] * 0, round(0.2 * n_pred))) # noise
  
  # Obtain values for y based on Bernoulli distribution, with input p
  p <- 1/(1+exp(-dm %*% dgm_par))
  y <- as.numeric(p > runif(length(p))) 
  
  # Obtain observed values of c-statistic and 
  # average predicted probability of an event
  
  obs_auc <- fastAUC(p = p, y = y) # obtain R2cs based on p and y
  obs_prev <- mean(y) # "Observed" prevalence
  
  # Sum of absolute differences of both values:
  #abs(obs_cstat-pref_cstat) + abs(obs_prev-pref_prev)
  (obs_auc-pref_auc)^2 + (obs_prev-pref_prev)^2# alternative, not sure which one is better
  
}

dist_R2_prev <- function(par, pref_R2, pref_prev){
  # par is a vector with initial guesses for both
  # the intercept and the beta1-3 coefficients.
  # However since beta1-3 have been restricted (half is strong, other half weak)
  # only one value is necessary.
  
  # Providing the beta0 and beta1-3 as specified in the par object. 
  dgm_par <-
    c(par[1], 
      rep(par[2] * 3, round(0.4 * n_pred)),  # strong
      rep(par[2], round(0.4 *  n_pred)),  # medium
      rep(par[2] * 0, round(0.2 * n_pred))) # noise
  
  # Obtain values for y based on Bernoulli distribution, with input p
  p <- 1/(1+exp(-dm %*% dgm_par))
  y <- as.numeric(p > runif(length(p))) 
  
  # Obtain observed values of c-statistic and 
  # average predicted probability of an event
  
  obs_R2 <- pseudo_Rsqrs(p = p, y = y) # obtain R2cs based on p and y
  obs_prev <- mean(y) # "Observed" prevalence
  
  # Sum of absolute differences of both values:
  #abs(obs_cstat-pref_cstat) + abs(obs_prev-pref_prev)
  (obs_R2-pref_R2)^2 + (obs_prev-pref_prev)^2# alternative, not sure which one is better
  
}


checking <- function(par){
  ## What do the observed prevalence and c-statistic look like?
  # Providing the beta0 and beta1-3 as specified in the par object. 
  
  # 40% of candidate predictors is strong (*3), 40% is medium, and 20% is noise.
  dgm_par <-
    c(par[1], 
      rep(par[2] * 3, round(0.4 * n_pred)), 
      rep(par[2], round(0.4 *  n_pred)), 
      rep(par[2] * 0, round(0.2 * n_pred)))
  
  # Obtain values for y based on uniform distribtuin, with input p
  p <- 1/(1+exp(-dm %*% dgm_par))
  y <- as.numeric(p > runif(length(p))) 
  
  # Obtain observed values
  #obs_cstat <- c_stat2(preds = p, outcome = y) # obtain c-statistic based on p and y
  obs_cstat <- fastAUC(p = p, y = y)
  obs_prev <- mean(y) # THE OBSERVED PREVALENCE IS NOT A FUNCTION OF JUST THE INTERCEPT
  c("cstat" = obs_cstat, "prev" = obs_prev)
}

# Same idea as the checking function, however, using the validation data:
checking_val <- function(par){
  
  # Providing the beta0 and beta1-3 as specified in the par object. 
  dgm_par_val <-
    c(par[1], 
      rep(par[2] * 3, round(0.4 * n_pred)), 
      rep(par[2], round(0.4 *  n_pred)), 
      rep(par[2] * 0, round(0.2 * n_pred)))
  
  
  # Obtain values for y based on uniform distribution, with input p
  p_val <- 1/(1+exp(-dm_val %*% dgm_par_val))
  y_val <- as.numeric(p_val > runif(length(p_val)))
  
  # Obtain observed values
  #obs_cstat <- c_stat2(preds = p_val, outcome = y_val) # obtain c-statistic based on p and y
  obs_cstat <- fastAUC(p = p_val, y = y_val)
  obs_prev <- mean(y_val) # THE OBSERVED PREVALENCE IS NOT A FUNCTION OF JUST THE INTERCEPT
  c("cstat" = obs_cstat, "prev" = obs_prev)
}




## Development dataset ##
n <- 30000 # setting n to develop betas on.
n_pred <- 10 # For the first and third scenario's at least
sigma <- matrix(0.2, ncol = n_pred, nrow = n_pred) # create covariance matrix to be used as input
diag(sigma) <- 1 # set the diagonal to 1
mu <- c(rep(0, n_pred)) # provide a vector of values for mu -> standard normal

X <- mvrnorm(n = n, mu = mu, Sigma = sigma) # create 3 predictor columns
dm <- cbind(1, X) # Putting the above in a data matrix, including intercept


## Create validation dataset ##
n_val <- 100000 # setting n
X_val <- mvrnorm(n = n_val, mu = mu, Sigma = sigma) # create 10 predictor columns
dm_val <- cbind(1, X_val) # Putting the above in a data matrix, including intercept

  # Optimizing the dist_R2_prev to determine the beta coefficients for which the squared distance 
  # between the observed and preferred values of R2 and prevalence are minimal. 
  # par_i are just initial values for the coefficients. 
  
#### 30 repetitions #####
set.seed(123)
system.time(results <- replicate(n = 30, 
                       optim(c(0.01, 0.1), 
                             dist_R2_prev,
                             pref_R2 = R2[3],
                             pref_prev = 0.5, 
                             control = list(maxit = 1000)
                             ),
                       simplify = F))

results_check <- apply(sapply(results, '[[', 1), 2, checking)
apply(results_check, 1, summary)
par <- apply(sapply(results, '[[', 1), 1, median) 
checking_val(par = par)

set.seed(123)
system.time(results <- replicate(n = 30, 
                                 optim(c(-1.66, 0.11), 
                                       dist_AUC_prev,
                                       pref_auc = 0.75,
                                       pref_prev = 0.5, 
                                       control = list(maxit = 1000)
                                 ),
                                 simplify = F))

results_check <- apply(sapply(results, '[[', 1), 2, checking)
apply(results_check, 1, summary)
par <- apply(sapply(results, '[[', 1), 1, median) 
checking_val(par = par)

## Initial coefficients for event rate at 0.05, AUC = 0.75 (n_pred = 10)
   ## c(-3.36, 0.1)
## Initial coefficients for event rate at 0.2, AUC  = 0.75 (n_pred = 10)
   ## c(-1.66, 0.11)
## Initial coefficients for event rate at 0.5, AUC  = 0.75 (n_pred = 10)  
   ## c(0.01, 0.1))
  


