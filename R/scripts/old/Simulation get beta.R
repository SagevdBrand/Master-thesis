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


#######################
### data generation ###
#######################

sigma <- matrix(0.2, ncol = 3, nrow = 3) # create covariance matrix to be used as input
diag(sigma) <- 1 # set the diagonal to 1
mu <- c(rep(3,3)) # provide a vector of values for mu
n <- 10000 # setting n
X <- mvrnorm(n = n, mu = mu, Sigma = sigma) # create 3 predictor columns
dm <- cbind(1, X) # Putting the above in a data matrix, including intercept


##########################################
##########################################
##########################################
## Defining a function to get the sum of absolute 
## differences between the preferred and observed 
## values of both the C-statistic and prevalence

func <- function(par){
# par is a vector with initial guesses for both
# the intercept and the beta1-3 coefficients.
# However since beta1-3 have been restricted to be equal
# only one value is necessary.
  
# What are the preferred values of the c-statistic and prevalence?
pref_cstat <- 0.75
pref_prev <- 0.2


# Providing the beta0 and beta1-3 as specified in the par object. 
# Again, beta1-3 are restricted to be equal.
dgm_par <- c(par[1], rep(par[2], 3)) 

# Obtain values for y based on Bernoulli distribution, with input p
p <- exp(dm %*% dgm_par)/(1+exp(dm %*% dgm_par))
y <- rbinom(length(p),1,p)

# fit a model to retrieve predicted probability
glm_fit <- glm(y ~ X, family = "binomial")

#obs_p <- predict(glm_fit, type =  "response") 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Not sure for the c-statistic whether I need to use obs_p or p? #
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Obtain observed values of c-statistic and 
# average predicted probability of an event

obs_cstat <- c_stat2(preds = p, outcome = y) # obtain c-statistic based on p and y
obs_prev <- 1 / (1 + exp(-glm_fit$coefficients[1])) #transform the log-odds to probability

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Maybe a dumb question, but is this the way to get the average predicted probability for an event? #
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Sum of absolute differences of both values:
sum(abs(obs_cstat-pref_cstat), abs(obs_prev-pref_prev))

}

## Try to optimize using initial guesses 0.5 and 0.5 
(results <- optim(c(0.1, 0.3), func, method = "BFGS"))

###################
## check results ##
###################

## What do the observed prevalence and c-statistic look like?
dgm_par <- c(results$par[1], rep(results$par[2], 3)) 

# Obtain values for y based on Bernoulli distribution, with input p
p <- exp(dm %*% dgm_par)/(1+exp(dm %*% dgm_par))
y <- rbinom(length(p),1,p)

# fit a model to retrieve predicted probability (for c-statistic)
glm_fit <- glm(y ~ X, family = "binomial")
#obs_p <- predict(glm_fit, type =  "response")

# Obtain observed values
(obs_cstat <- c_stat2(preds = p, outcome = y)) # obtain c-statistic based on p and y

(obs_prev <- 1 / (1 + exp(-glm_fit$coefficients[1]))) #transform the log-odds to probability

sum(y)/length(y)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# It appears that the cstat does what we want it to do, but the prevalence does not..                #
# Maybe need to add some weight to the absolute difference of the preferred and observed prevalence? #
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
