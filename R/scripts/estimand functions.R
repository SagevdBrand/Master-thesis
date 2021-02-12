########################################
########################################
########################################
### Simulation functions
### Estimands

###### Apparent performance ######
## Depends on which model is used.
get_app_estimands <- function(df, model) {
  # Fit model depending on scenario
  # And get predicted probabilities
  if (model == "OLS") {
    fit_app <- glm(y ~ ., family = "binomial", data = df) #%>%
    #step(direction = "backward", trace = F) # This should be changeable depending on which thing you're using
    p_app <- predict(fit_app, type = "response")
    app_matrix <- model.matrix(object = fit_app$formula, data = df)
    
  } else { # If model = Firth (or svm for now)
    fit_app <- logistf(y ~ ., data = df, flic = T)
    app_matrix <- model.matrix(object = fit_app$formula, data = df)
    p_app <- 1 / (1 + exp(-app_matrix %*% fit_app$coefficients))
  }
  
  # Performance measures
  auc_app <- fastAUC(p = p_app, y = df$y)
  R2_app <- pseudo_Rsqrs(p = p_app, y = df$y)
  
  calib_app <-
    calib(
      modelmatrix = app_matrix,
      data = df,
      coefs = fit_app$coefficients
    )
  
  eci_app <-
    eci_bvc(
      data = df,
      modelmatrix = app_matrix,
      coefs = fit_app$coefficients,
      preds = p_app
    )
  
  # Save results
  results <-
    c(
      "AUC_app" = auc_app,
      "calib_int_app" = calib_app['intercept'],
      "calib_slope_app" = calib_app['slope'],
      "R2_cox_snell_app" = R2_app,
      "ECI_app" = eci_app
    )
}

get_app_results <- function(scenario, df) {
  results_app <- list() #object to store results
  
  for (i in 1:length(df)) {
    model <- s1[i, ]$model
    results_app[[i]] <- get_app_estimands(df = df[[i]], model = model)
    
    }
  names(results_app) <- c(1:length(df))
  return(results_app)
}


###### Cross-validation  approaches ######
## Functions called within:

# obtain matrix for ciAUC
auc_mat <- function(ciAUC){
  ciAUC <- as.matrix(unlist(ciAUC))
  ciAUC <- as.matrix(ciAUC[-5,])
  rownames(ciAUC) <- c("cvmean", "se", "ci_lower", "ci_upper")
  colnames(ciAUC) <- ("cstat_bias_c")
  return(ciAUC)
}

## Store other results in a nice table
results_table <- function(results_folds) {
  # create matrix to store results
  results <- matrix(NA, nrow = 4, ncol = ncol(results_folds), dimnames = list(c("cvmean", "se", "ci_lower", "ci_upper"), c(colnames(results_folds))))
  
  for (i in 1:ncol(results_folds)){
    cvmean <- mean(results_folds[,i])
    results[1,i] <- cvmean
    se <- sd(results_folds[,i]) / sqrt(length(results_folds[,i]))
    results[2,i] <- se
    
    if (length(results_folds) < 30) {
      results[3,i] <- cvmean - qt(0.975, df = (length(results_folds[,i]) -1)) * se
      results[4,i] <- cvmean + qt(0.975, df = (length(results_folds[,i]) -1)) * se
      
    } else {
      results[3,i] <- cvmean - qnorm(0.975) * se
      results[4,i] <- cvmean + qnorm(0.975) * se
    }
    
  }
  return(results)
}

#########
## R^2 ##
#########

Rcoxsnell_for_folds <- function(folds, p, df){
  
  .pseudo_Rsqrs <- function(p, y){ 
    
    .LL <- function(p, y){
      sum(y*log(p)+(1-y)*log(1-p))
    }
    
    LL_fit  <- .LL(p=p, y=y) 
    LL_null <- .LL(p=mean(y), y=y)
    
    cox <- 1-exp(-(LL_fit-LL_null)*2/length(y)) 
    cox_max <- 1 - exp(2 * length(y) ^ (-1) * LL_null)
    c("cox"=cox)
    
  }
  n.folds <- length(folds)
  
  Rcoxsnell_folds <- matrix(NA, nrow = n.folds, ncol = 1, dimnames = list(c(), c("R2cox-snell_bias_c")))
  
  for (i in 1:n.folds){
    Rcoxsnell_folds[i,] <- .pseudo_Rsqrs(p = p[[i]], y = df[unlist(folds[[i]]),]$y)
  }
  
  results<- results_table(Rcoxsnell_folds)
  names(results) <- colnames(Rcoxsnell_folds)
  return(results)
  
}

#################
## Calibration ##
#################





## 10 fold cv
# write a function to obtain estimates per df, so first try with subscenario 1a

get_10cv_estimands <- function(df, model){
V <- 10
####################
## Splitting data ##
####################
  cvFoldsB <- function(Y, V){  #Create Balanced CV folds (stratify by outcome)
    Y0 <- split(sample(which(Y=="0")), rep(1:V, length=length(which(Y==0))))
    Y1 <- split(sample(which(Y=="1")), rep(1:V, length=length(which(Y==1))))
    folds <- vector("list", length=V)
    for (v in seq(V)) {folds[[v]] <- c(Y0[[v]], Y1[[v]])}		
    return(folds)
  }
  
  foldsb <- cvFoldsB(Y = df$y, V = V)

#################################################
## Getting predictions depending on model used ##
#################################################
  doFit <- function(V, folds, model){  #Train/test glm for each fold
      
    if (model == "OLS") {
      fit <- glm(y~., data=df[-foldsb[[V]],], family=binomial) #%>%
      #step(direction = "backward", trace = F) # This should be changeable depending on which thing you're using
      iv_matrix <- model.matrix(object = fit$formula, data = df[foldsb[[V]],])
      p <- predict(fit, newdata=df[foldsb[[V]],], type = "response")
      
    } else { # If model = Firth (or svm for now)
      fit <- logistf(y ~ ., data = df[-foldsb[[V]],], flic = T)
      iv_matrix <- model.matrix(object = fit$formula, data = df[foldsb[[V]],]) #iv matrix of test fold
      p <- 1 / (1 + exp(-iv_matrix %*% fit$coefficients))
    }
    
    coefs <- coef(fit)
    results <- list(p, coefs, iv_matrix)
  
  }
  
  results <- (lapply(seq(V), doFit, folds = foldsb, model = model))  
  p <- unlist(sapply(results, "[[", 1)) # Get out pred values as a single thing
  p[unlist(foldsb)] <- p #Re-order pred values
  
  p_per_fold <- sapply(results, "[[", 1) # getting the predictions per fold
  coefs <- (sapply(results, "[[", 2)) # Get out model coefficients
  iv_matrix <- (sapply(results, "[[", 3)) # get out the test matrices
  
  ## Optain estimands ##
  ciAUC <- ci.cvAUC(predictions=p, labels=df$y, folds=foldsb, confidence=0.95)
  auc <- auc_mat(ciAUC = ciAUC)
 
  R2 <- Rcoxsnell_for_folds(folds = foldsb, p = p_per_fold, df = df)
  calibration_results <- calibration_for_folds(folds = foldsb, modelmatrix = mm_full, data = df, coefs = coefs)
  
  return(results)
  }

test <- get_10cv_estimands(df = df, model = model)

######################################################################################
##########################################OLD#########################################
######################################################################################


fastAUC <- function(p, y) {
  x1 = p[y==1]; n1 = length(x1); 
  x2 = p[y==0]; n2 = length(x2);
  r = rank(c(x1,x2))  
  auc = (sum(r[1:n1]) - n1*(n1+1)/2) / n1 / n2
  return(auc)
}

# obtain matrix for ciAUC
auc_mat <- function(ciAUC){
  ciAUC <- as.matrix(unlist(ciAUC))
  ciAUC <- as.matrix(ciAUC[-5,])
  rownames(ciAUC) <- c("cvmean", "se", "ci_lower", "ci_upper")
  colnames(ciAUC) <- ("cstat_bias_c")
  return(ciAUC)
}

############################################# 
## obtain nice matrix for each performance ##
#############################################

#################
## Calibration ##
#################

calibration <- function(lp,data){
  coef(glm(data$y~lp, family = "binomial"))
}


calib <- function(modelmatrix, data, coefs) {
  phat <- 1/(1+exp(-(iv_matrix[[1]]%*%coefs[,1]))) # risk score of individual patients
  slope <- c(coef(glm(df$y ~ log(phat/(1-phat)),family="binomial"))[2])
  intercept <- coef(glm(data$y ~ offset(log(phat/(1-phat))),family="binomial"))
  results <- c(intercept, slope)
  names(results) <- c("intercept", "slope")
  return(results)
}


calibration_for_folds <- function(folds, modelmatrix, data, coefs){
  n.folds <- length(folds)
  cal_folds <- matrix(NA, nrow = n.folds, ncol = 2, dimnames = list(c(), c("calib_int_bias_c", "calib_slope_bias_c")))
  
  for (i in 1:n.folds){
    cal_folds[i,] <- calib(modelmatrix = modelmatrix, data = data, coefs = coefs[,i])
  }
  
  results<- results_table(cal_folds)
  names(results) <- colnames(cal_folds)
  return(results)
}


#########
## ECI ##
#########
## In contrast, the ECI is the average squared difference of the predicted
# probabilities ^pn with the estimated observed probabilities ^on instead
# of the actual outcomes.


#### Code by Ben: #####
# ECI own function
eci_bvc <- function(data, modelmatrix, coefs, preds){
  
  .obtain_calout <- function(data, modelmatrix, coefs){
    phat <- 1/(1+exp(-(modelmatrix%*%coefs)))
    calout <- loess(data$y ~ log(phat/(1-phat)), data = data)
    return(calout)
  }
  
  calout <- .obtain_calout(data = data, modelmatrix = modelmatrix, coefs = coefs)
  eci <- (mean((preds-fitted(calout))*(preds-fitted(calout))))*(100)
  
  return(eci)
}
# ECI own function, compared with random model
# eci_rel <- function(calout,preds,k,outc){
#   prevm=matrix((table(outc)/length(outc))[1:k],nrow=dim(preds)[1],ncol=k,byrow=T)
#   ecir=mean((preds-prevm)*(preds-prevm))
#   ecim=mean((preds-fitted(calout))*(preds-fitted(calout)))
#   return(ecim/ecir)
# }


