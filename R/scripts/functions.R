########################################
########################################
########################################
### Simulation functions
### Data generation

## AUC
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

#### distance function ####

## Defining a function to get the sum of 
## squared differences between the preferred and observed 
## values of both the R^2CS and prevalence.
d_R2_prev <- function(par,pref_R2, pref_prev){
  # par is a vector with initial guesses for both
  # the intercept and the beta1-3 coefficients.
  # However since beta1-3 are restricted to be equal
  # only one value is necessary.
  
  # Providing the beta0 and beta1-3 as specified in the par object. 
  # Again, beta1-3 are restricted to be equal.
  dgm_par <- c(par[1], rep(par[2], n_pred)) 
  
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

####  Check results function ####

## To check whether the coefficients as given by optim()
## return the desired C-stat and Prevalence
checking <- function(par){
  ## What do the observed prevalence and c-statistic look like?
  dgm_par <- c(par[1], rep(par[2], n_pred)) 
  
  # Obtain values for y based on Bernoulli distribution, with input p
  p <- plogis(dm %*% dgm_par)
  y <- rbinom(length(p),1,p)
  
  # Obtain observed values
  #obs_cstat <- c_stat2(preds = p, outcome = y) # obtain c-statistic based on p and y
  obs_cstat <- fastAUC(p = p, y = y)
  obs_prev <- mean(y) # THE OBSERVED PREVALENCE IS NOT A FUNCTION OF JUST THE INTERCEPT
  c("cstat" = obs_cstat, "prev" = obs_prev)
}

# Same idea as the checking function, however, using the validation data:
checking_val <- function(par){
  dgm_par_val <- c(par[1], rep(par[2], n_pred)) 
  
  p_val <- plogis(dm_val %*% dgm_par_val)
  y_val <- rbinom(length(p_val),1,p_val)
  
  # Obtain observed values
  #obs_cstat <- c_stat2(preds = p_val, outcome = y_val) # obtain c-statistic based on p and y
  obs_cstat <- fastAUC(p = p_val, y = y_val)
  obs_prev <- mean(y_val) # THE OBSERVED PREVALENCE IS NOT A FUNCTION OF JUST THE INTERCEPT
  c("cstat" = obs_cstat, "prev" = obs_prev)
}
