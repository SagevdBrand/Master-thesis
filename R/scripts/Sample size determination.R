######## Sample size #######
## libraries:
library(pmsampsize)
library(rms)

set.seed(123)

## Approximate R^2 from a desired AUC of 0.75 and prev of 0.2 & 0.05

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

pref_cstat <- 0.75
pref_prev <- .2

Rcs_prev_.2 <- approximate_R2(auc = pref_cstat, prev = 0.2)$R2.coxsnell
Rcs_prev_.05 <- approximate_R2(auc = pref_cstat, prev = 0.05)$R2.coxsnell

n_estimate <- pmsampsize(type = "b", parameters = 3, prevalence = 0.2, rsquared = Rcs_prev_.2)$sample_size