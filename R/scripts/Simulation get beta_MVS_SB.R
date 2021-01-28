#############################
### try-out simulate data ###
#############################

################### Setting up #############################
library(MASS)
library(tidyverse)

set.seed(123)

#### Functions #####

c_stat2 <- function(preds, outcome){
  preds <- as.matrix(preds)
  cats <- sort(unique(outcome))
  n_cat <- length(cats)
  n0   <- sum(outcome == cats[2])
  n1   <- length(outcome) - n0
  r <- rank(preds[,1])
  S0 <- sum(as.numeric(r[outcome == cats[2]]))
  (S0 - n0 * (n0 + 1)/2)/(as.numeric(n0) * as.numeric(n1))
}

#########
## R^2 ##
#########

pseudo_Rsqrs <- function(p, SIMxy){ 
  
  .LL <- function(p, SIMxy){
    sum(SIMxy[,"y"]*log(p)+(1-SIMxy[,"y"])*log(1-p))
  }
  
  LL_fit  <- .LL(p=p, SIMxy=SIMxy) 
  LL_null <- .LL(p=mean(SIMxy[,"y"]), SIMxy=SIMxy)
  
  cox <- 1-exp(-(LL_fit-LL_null)*2/nrow(SIMxy)) 
  cox_max <- 1 - exp(2 * nrow(SIMxy) ^ (-1) * LL_null)
  c("cox"=cox)
  
}

#######################
### data generation ###
#######################

sigma <- matrix(0.2, ncol = 3, nrow = 3) # create covariance matrix to be used as input
diag(sigma) <- 1 # set the diagonal to 1
mu <- c(0,0, 0) # provide a vector of values for mu
n <- 10000 #setting n
X <- mvrnorm(n = n, mu = mu, Sigma = sigma) # create 3 predictor columns
dm <- cbind(1, X) # Putting the above in a data matrix, including intercept


##########################################
##########################################
##########################################
## Defining a function to get the sum of absolute 
## differences between the preferred and observed 
## values of both the C-statistic and prevalence

par <- c(0.5, 0.5)

func <- function(par,pref_cstat,pref_prev){
# par is a vector with initial guesses for both
# the intercept and the beta1-3 coefficients.
# However since beta1-3 have been restricted to be equal
# only one value is necessary.

# Providing the beta0 and beta1-3 as specified in the par object. 
# Again, beta1-3 are restricted to be equal.
dgm_par <- c(par[1], rep(par[2], 3)) 

# Obtain values for y based on Bernoulli distribution, with input p
#system.time(p <- exp(dm %*% dgm_par)/(1+exp(dm %*% dgm_par)))
system.time(p <- 1/(1+exp(-dm %*% dgm_par))) # SAME THING, JUST SLIGHTLY LESS COMPUTATION

#system.time(y <- rbinom(length(p),1,p))
system.time(y <-as.numeric(p>runif(length(p)))) # SAME THING, JUST SLIGHTLY FASTER

# Obtain observed values of c-statistic and 
# average predicted probability of an event

obs_cstat <- c_stat2(preds = p, outcome = y) # obtain c-statistic based on p and y
obs_prev <- mean(y) # KEEP IT SIMPLE ;)

# Sum of absolute differences of both values:
#abs(obs_cstat-pref_cstat) + abs(obs_prev-pref_prev)
(obs_cstat-pref_cstat)^2 + (obs_prev-pref_prev)^2 # alternative, not sure which one is better

}

## Try to optimize using initial guesses 0.5 and 0.5 
system.time(results <- optim(c(0.5, 0.5), func, pref_cstat = 0.75, pref_prev = 0.2))

###################
## check results ##
###################

## What do the observed prevalence and c-statistic look like?
dgm_par <- c(results$par[1], rep(results$par[2], 3)) 

# Obtain values for y based on Bernoulli distribution, with input p
p <- plogis(dm %*% dgm_par)
y <- rbinom(length(p),1,p)

# Obtain observed values
(obs_cstat <- c_stat2(preds = p, outcome = y)) # obtain c-statistic based on p and y
(obs_prev <- mean(y)) # THE OBSERVED PREVALENCE IS NOT A FUNCTION OF JUST THE INTERCEPT

