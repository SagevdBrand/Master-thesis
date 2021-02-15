########################################
########################################
########################################
### Simulation functions
### Estimands

######################################
###### Basic estimand functions ######
######################################

#########
## AUC ##
#########

fastAUC <- function(p, y) {
  x1 = p[y==1]; n1 = length(x1); 
  x2 = p[y==0]; n2 = length(x2);
  r = rank(c(x1,x2))  
  auc = (sum(r[1:n1]) - n1*(n1+1)/2) / n1 / n2
  return(auc)
}

#################
## Calibration ##
#################

calib <- function(modelmatrix, data, coefs) {
  phat <- 1/(1+exp(-(modelmatrix%*%coefs))) # risk score of individual patients
  slope <- c(coef(glm(data$y ~ log(phat/(1-phat)),family="binomial"))[2])
  intercept <- coef(glm(data$y ~ offset(log(phat/(1-phat))),family="binomial"))
  results <- c(intercept, slope)
  names(results) <- c("intercept", "slope")
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

#########
## ECI ##
#########

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

##########################################
###### Cross-validation  approaches ######
##########################################

## Function called within:



## 10 fold cv
# write a function to obtain estimates per df, so first try with subscenario 1a

get_cv_estimands <- function(df, model, V){

####################
## Splitting data ##
####################
  .cvFoldsB <- function(Y, V) {  #Create Balanced CV folds (stratify by outcome)
    Y0 <- split(sample(which(Y=="0")), rep(1:V, length=length(which(Y==0))))
    Y1 <- split(sample(which(Y=="1")), rep(1:V, length=length(which(Y==1))))
    folds <- vector("list", length=V)
    for (v in seq(V)) {folds[[v]] <- c(Y0[[v]], Y1[[v]])}		
    return(folds)
  }
  
  foldsb <- .cvFoldsB(Y = df$y, V = V)

#################################################
## Getting predictions depending on model used ##
#################################################
  .doFit <- function(V, folds, model){  #Train/test glm for each fold
      
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
  
  results <- (lapply(seq(V), .doFit, folds = foldsb, model = model))  
  p <- unlist(sapply(results, "[[", 1)) # Get out pred values as a single thing
  p[unlist(foldsb)] <- p #Re-order pred values
  
  p_per_fold <- sapply(results, "[[", 1) # getting the predictions per fold
  coefs <- (sapply(results, "[[", 2)) # Get out model coefficients
  iv_matrix <- (sapply(results, "[[", 3)) # get out the test matrices
  
  ## Obtain estimands ##
  ## for ECI
    .obtain_calout <- function(data, modelmatrix, coefs){
      phat <- 1/(1+exp(-(modelmatrix%*%coefs)))
      calout <- loess(data$y ~ log(phat/(1-phat)), data = data)
      return(calout)
    }
    
    ## Empty objects to store results
    cal_folds <- matrix(NA, nrow = V, ncol = 2)
    R2_folds <- c()
    eci_folds <- c()
    
    # For each fold, calculate calibration intercept & slope, R2 and ECI
    for (i in 1:V){
      cal_folds[i,] <- calib(modelmatrix = iv_matrix[[i]], data = df[unlist(foldsb[[i]]),], coefs = coefs[,i])
      R2_folds[i] <- pseudo_Rsqrs(p = p_per_fold[[i]], y = df[unlist(foldsb[[i]]),]$y)
      calout <- .obtain_calout(data = df[unlist(foldsb[[i]]),], modelmatrix = iv_matrix[[i]], coefs = coefs[,i])
      eci_folds[i] <- (mean((p_per_fold[[i]]-fitted(calout))*(p_per_fold[[i]]-fitted(calout))))*(100)
    }
    
    ## AUC
    auc_results <- as.vector(unlist(ci.cvAUC(predictions=p, labels=df$y, folds=foldsb, confidence=0.95)))[-5]
    names(auc_results) <- c(paste0("AUC_mean_", V, "fcv" ), paste0("AUC_se_", V, "fcv"), paste0("AUC_ci_lower_", V, "fcv"), paste0("AUC_ci_upper_", V, "fcv"))
    
    ## Get mean and standard error over all other results and store in a single vector:
    ## Calibration
    means <- apply(cal_folds, 2, mean)
    se <- apply(cal_folds, 2, function(x) (sd(x)/(sqrt(V) - 1)))
    
    calibration_results <-  c(means[1], se[1], means[2], se[2])
    names(calibration_results) <- c(paste0("calib_int_mean_", V, "fcv" ),
                                    paste0("calib_int_se_", V, "fcv"),
                                    paste0("calib_slope_mean_", V, "fcv" ),
                                    paste0("calib_slope_se_", V, "fcv"))
    
    ## R2 cox snell
    R2 <- c(mean(R2_folds), (sd(R2_folds)/(sqrt(V) - 1)))
    names(R2) <- c(paste0("R2_mean_", V, "fcv" ), paste0("R2_se_", V, "fcv"))
    
    ## ECI
    eci <- c(mean(eci_folds), (sd(eci_folds)/(sqrt(V) - 1)))
    names(eci) <- c(paste0("eci_mean_", V, "fcv" ), paste0("eci_se_", V, "fcv"))
    
    results <- c(auc_results, calibration_results, R2, eci)
  
  return(results)
  }


get_cv_results <- function(scenario, df, V) {
  results_cv <- list() #object to store results
  
  for (i in 1:length(df)) {
    model <- s1[i, ]$model
    results_cv[[i]] <- get_cv_estimands(df = df[[i]], model = model, V = V)
    
  }
  names(results_cv) <- c(1:length(df))
  return(results_cv)
}



