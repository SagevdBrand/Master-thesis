#############################
### try-out simulate data ###
#############################

################### Setting up #############################

library(MASS)
library(tidyverse)

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

#####################################
### Initial guesses and framework ###
#####################################

sigma <- matrix(0.2, ncol = 3, nrow = 3) # create covariance matrix to be used as input
diag(sigma) <- 1 # set the diagonal to 1
mu <- rep(1, 3) # provide a vector of mu
n <- 1000
X <- mvrnorm(n = n, mu = mu, Sigma = sigma) # create 3 predictor columns
dm <- cbind(1, X)

# random guess for the coefficients: they're all equal
# add a guess for the intercept

dgm_par_i <- c(0.2, rep(0.3, 3))

# Obtain initial values for the y
p <- exp(dm %*% dgm_par_i)/(1+exp(dm %*% dgm_par_i))
y <- rbinom(length(p),1,p)

glm_fit_i <- glm(y ~ X, family = "binomial")
obs_p <- predict(fit_i, type =  "response")


obs_cstat <- c_stat2(preds = obs_p, outcome = y)
obs_prev <- exp(glm_fit$coefficients[1]) / (1 + exp(glm_fit$coefficients[1]))

# Sum of absolute differences of both values:
sum(abs(obs_cstat-pref_cstat), abs(obs_prev-pref_prev))



##### Putting the above together:
func <- function(par){

# What are the values we are looking for?
pref_cstat <- 0.75
pref_prev <- 0.2


# Providing the beta0 and beta1-3 as specified in the par object. 
dgm_par <- c(par[1], rep(par[2], 3))

# Obtain initial values for the y
p <- exp(dm %*% dgm_par)/(1+exp(dm %*% dgm_par))
y <- rbinom(length(p),1,p)

# fit a model to retrieve predicted probability (for c-statistic)
glm_fit <- glm(y ~ X, family = "binomial")
obs_p <- predict(glm_fit, type =  "response")

# Obtain observed values
obs_cstat <- c_stat2(preds = obs_p, outcome = y)
obs_prev <- exp(glm_fit$coefficients[1]) / (1 + exp(glm_fit$coefficients[1]))

# Sum of absolute differences of both values:
sum(abs(obs_cstat-pref_cstat), abs(obs_prev-pref_prev))

}

optim(c(0.2, 0.4), func)
func(par = results)



exp(results[1]) / (1 + exp(results[1]))
