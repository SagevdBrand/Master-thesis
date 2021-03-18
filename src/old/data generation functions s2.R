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

# Based on predefined AUC and prevalence
# Using code from Riley et al.:
# Approximate R^2 from a desired AUC of 0.75 and prev of 0.2 & 0.05
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


####################################
####################################
####################################
### END SCRIPT