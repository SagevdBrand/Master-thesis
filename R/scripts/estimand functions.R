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
    results_app[[i]] <- get_app_estimands(df = as.data.frame(df[[i]]), model = model)
  }
  
  return(results_app)
}


###### Cross-validation  approaches ######
get_cv_estimands(data, scenario){
  
  
  
  
}





######################################################################################
##########################################OLD#########################################
######################################################################################


##########################
### Splitting the data ###
##########################

###################
## CV approaches ##
###################

## either balanced ##
cvFoldsB <- function(Y, V){  #Create Balanced CV folds (stratify by outcome)
  Y0 <- split(sample(which(Y=="0")), rep(1:V, length=length(which(Y==0))))
  Y1 <- split(sample(which(Y=="1")), rep(1:V, length=length(which(Y==1))))
  folds <- vector("list", length=V)
  for (v in seq(V)) {folds[[v]] <- c(Y0[[v]], Y1[[v]])}		
  return(folds)
}

## or unbalanced ##
cvFoldsUB <- function(Y, V) {
  folds <- as.vector(split(sample(row_number(Y)), rep(1:V, length = length(Y))))
  return(folds) 
}


########################
## applying the model ##
########################

## Using the same model across all the training sets and only getting the predictions out ##
predicting <- function(folds, data, formula) {
  V <- length(folds)
  
  doFit <- function(v, folds, data){  #Train/test glm for each fold
    fit <-  glm(modelform, data=data[-folds[[v]],], family=binomial)
    coefs <- coef(fit)
    pred <- predict(fit, newdata=data[folds[[v]],], type="response")
    lp <- fit$linear.predictors
    results <- list(pred, coefs, lp)
    return(results)
  }
  
  results <- (lapply(seq(V), doFit, folds=foldsb, data=xydata)) #CV train/predict
  predictions <- unlist(sapply(results, "[[", 1)) # Get out pred values as a single thing
  predictions[unlist(folds)] <- predictions #Re-order pred values
  
  p_per_fold <- sapply(results, "[[", 1) #getting the predictions per fold
  coefs <- (sapply(results, "[[", 2)) # Get out model coefficients
  lp <- (sapply(results, "[[", 3)) # get out the linear predictors
  
  results <-  list("predictions" = predictions, "coefficients" = coefs, "linear predictors" = lp, "p_per_fold" = p_per_fold)
  return(results)
}


####################################
## obtaining performance measures ##
####################################
c.stat2 <- function(preds, outcome){
  preds <- as.matrix(preds)
  cats <- sort(unique(outcome))
  n_cat <- length(cats)
  n0   <- sum(outcome == cats[2])
  n1   <- length(outcome) - n0
  r <- rank(preds[,1])
  S0 <- sum(as.numeric(r[outcome == cats[2]]))
  (S0 - n0 * (n0 + 1)/2)/(as.numeric(n0) * as.numeric(n1))
}

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

Rcoxsnell_for_folds <- function(folds, p, SIMxy){
  n.folds <- length(folds)
  Rcoxsnell_folds <- matrix(NA, nrow = n.folds, ncol = 1, dimnames = list(c(), c("R2cox-snell_bias_c")))
  
  for (i in 1:n.folds){
    Rcoxsnell_folds[i,] <- pseudo_Rsqrs(p = p[[i]], SIMxy = SIMxy[unlist(folds[[i]]),])
  }
  
  results<- results_table(Rcoxsnell_folds)
  names(results) <- colnames(Rcoxsnell_folds)
  return(results)
  
}

#################
## Calibration ##
#################

calibration <- function(lp,data){
  coef(glm(data$y~lp, family = "binomial"))
}


calib <- function(modelmatrix, data, coefs) {
  phat <- 1/(1+exp(-(modelmatrix%*%coefs))) # risk score of individual patients
  slope <- c(coef(glm(data$y ~ log(phat/(1-phat)),family="binomial"))[2])
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


