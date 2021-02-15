########################################
########################################
########################################
### Simulation functions
### Data generation

######### To do:
## Add options for when difference in predictor effects
## Add options for different number of predictors (scenario 3)
## Adapt initial coefficients for different number of predictors (scenario 3)
## How to determine initial coefficients for optimization?

## DATA GEN:
## Add a check for how many y-values are drawn as an event!

#####################################
### Approximation of R2 Cox-Snell ###
#####################################

##based on predefined AUC and prevalence
## Using code from Riley et al., see title below!
source("scripts/libraries.R")

## Approximate R^2 from a desired AUC of 0.75 and prev of 0.2 & 0.05
# “A note on estimating the Cox-Snell R2 from a reported C-statistic 
# (AUROC) to inform sample size calculations for developing 
# a prediction model with a binary outcome” 
approximate_R2 <- function(auc, prev, n = 1000000){
  # define mu as a function of the C statistic
  mu <- sqrt(2) * qnorm(auc)
  # simulate large sample linear prediction based on two normals
  # for non-eventsN(0, 1), events and N(mu, 1)
  LP <- c(rnorm(prev*n, mean=0, sd=1), rnorm((1-prev)*n, mean=mu, sd=1))
  y <- c(rep(0, prev*n), rep(1, (1-prev)*n))
  # Fit a logistic regression with LP as covariate;
  # this is essentially a calibration model, and the intercept and
  # slope estimate will ensure the outcome proportion is accounted
  # for, without changing C statistic
  fit <- lrm(y~LP)
  
  max_R2 <- function(prev){
    1-(prev^prev*(1-prev)^(1-prev))^2
  }
  
  return(list(R2.nagelkerke = as.numeric(fit$stats['R2']),
              R2.coxsnell = as.numeric(fit$stats['R2']) * max_R2(prev)))
}

##########################
## Performance measures ##
##########################

## AUC
## Adapted from: 1.	https://gist.github.com/traversc/1446ebe1dcc2d84bccdca781d4f1fa2a
fastAUC <- function(p, y) {
  x1 = p[y==1]; n1 = length(x1); 
  x2 = p[y==0]; n2 = length(x2);
  r = rank(c(x1,x2))  
  auc = (sum(r[1:n1]) - n1*(n1+1)/2) / n1 / n2
  return(auc)
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

###########################
#### Distance function ####
###########################

## Defining a function to get the sum of 
## squared differences between the preferred and observed 
## values of both the R^2CS and prevalence.
dist_R2_prev <- function(par,pref_R2, pref_prev){
  # par is a vector with initial guesses for both
  # the intercept and the beta1-3 coefficients.
  # However since beta1-3 have been restricted (half is strong, other half weak)
  # only one value is necessary.
  
  # Providing the beta0 and beta1-3 as specified in the par object. 
  # Half of the predictors is strong
  dgm_par <-
    c(par[1], 
      rep(par[2] * 3, round(0.4 * n_pred)),  # strong
      rep(par[2], round(0.4 *  n_pred)),  # medium
      rep(par[2] * 0, round(0.2 * n_pred))) # noise

  # Obtain values for y based on Bernoulli distribution, with input p
  p <- 1/(1+exp(-dm %*% dgm_par))
  y <- as.numeric(p > runif(length(p))) 

  # Obtain observed values of R2 and 
  # average predicted probability of an event
  
  obs_R2 <- pseudo_Rsqrs(p = p, y = y) # obtain R2cs based on p and y
  obs_prev <- mean(y) # "Observed" prevalence
  
  # Sum of absolute differences of both values:
  #abs(obs_cstat-pref_cstat) + abs(obs_prev-pref_prev)
  (obs_R2-pref_R2)^2 + (obs_prev-pref_prev)^2# alternative, not sure which one is better
  
}

###############################
#### Optimization function ####
###############################

# function to optimize beta coefficients
optim_beta <- function(prev_scenario, R2_scenario){
  
  # depending on predictor effects and the prevalence, start with different initial coefficients
  # With trial and error, it was discovered that the optimization would get
  # stuck at a local minimum.
  if (prev_scenario == 0.05){ 
    par_i <- c(-3.36, 0.1) 
  } else if (prev_scenario == 0.2) { 
    par_i <- c(-1.66, 0.11)
  } else {
    par_i <- c(0.01, 0.1)
    }

  # Optimizing the dist_R2_prev to determine the beta coefficients for which the squared distance 
  # between the observed and preferred values of R2 and prevalence are minimal. 
  # par_i are just initial values for the coefficients. 
  
  #### 20 repetitions #####
  results <- replicate(n = 20, 
                       optim(par_i, 
                             dist_R2_prev,
                             pref_R2 = R2_scenario,
                             pref_prev = prev_scenario),
                       simplify = F)
  # Computation time is about 1 minute
  
  # Taking the median of the coefficients,
  # to reduce influence of potential outliers 
  par <- apply(sapply(results, '[[', 1), 1, median) 
  
  return(par)
}

#################################
####  Check results function ####
#################################

## To check whether the coefficients as given by optim()
## return the desired C-stat and Prevalence
checking <- function(par){

  dgm_par <-
    c(par[1], 
      rep(par[2] * 3, round(0.4 * n_pred)), 
      rep(par[2], round(0.4 *  n_pred)), 
      rep(par[2] * 0, round(0.2 * n_pred)))
  
  # Obtain values for y based on uniform distribution, with input p
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
  
  # Half strong
  dgm_par_val <- c(par[1], 
                   rep(par[2] * 3, round(0.4 * n_pred)), 
                   rep(par[2], round(0.4 *  n_pred)), 
                   rep(par[2] * 0, round(0.2 * n_pred)))
  
  # Obtain values for y based on uniform distribution, with input p
  p_val <- 1/(1+exp(-dm_val %*% dgm_par_val))
  y_val <- as.numeric(p_val > runif(length(p_val)))
  
  # Obtain observed values
  obs_cstat <- fastAUC(p = p_val, y = y_val)
  obs_prev <- mean(y_val) 
  c("cstat" = obs_cstat, "prev" = obs_prev)
}


#######################################################
### Generate data script once betas have been found ###
#######################################################


generate_data <- function(scenario, validation = c(TRUE, FALSE)){
  
  s_list <- split(scenario, seq(nrow(scenario))) 
  
  .per_element <- function(s_list){
    
    sigma <- matrix(0.2, ncol = s_list$dim, nrow = s_list$dim) # create covariance matrix to be used as input
    diag(sigma) <- 1 # set the diagonal to 1
    mu <- c(rep(0, s_list$dim))
    
    # create candidate predictors
    # If it is the validation set, n = always 100000. Otherwise, take n as specified in each scenario
    if (validation == TRUE) {
      X <- mvrnorm(n = 100000, mu = mu, Sigma = sigma)
    }
    else {
      X <-
        mvrnorm(n = s_list$n, mu = mu, Sigma = sigma)
    } 
    dm <- cbind(1, X) # Putting the above in a data matrix, including intercept

    dgm_par <- c(s_list$par1, 
                 rep(s_list$par2 * 3, round(0.4 * s_list$dim)),
                 rep(s_list$par2, round(0.4 * s_list$dim)),
                 rep(0, round(0.2 * s_list$dim))
                 ) 
    
    
    # Obtain values for y based on Bernoulli distribution, with input p
    p <- 1/(1+exp(-dm %*% dgm_par))
    y <- as.numeric(p > runif(length(p)))
    
    #######################################################
    ## BUILT CHECK FOR HOW MANY EVENTS HAVE BEEN SAMPLED ##
    #######################################################
    
    return(as.data.frame(cbind(X, y)))
  }
  
  data_in_lists <- lapply(s_list, .per_element) #Apply function above
  return(data_in_lists)
  
}

####################################
####################################
####################################
### END SCRIPT