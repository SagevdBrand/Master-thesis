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

##########
## MAPE ##
##########

MAPE <- function(p,iv_matrix, dgm_par){
  p_true <- 1 / (1 + exp(-iv_matrix %*% dgm_par))
  mean(abs(p_true-p))
}

##################################
###### Apparent performance ######
##################################

get_app_estimands <- function(df, model, dgm_par, pred_selection) {
  
  # Fit model depending on scenario
  # And get predicted probabilities and modelmatrix
  if (model == "Firth" & pred_selection == "none") {
    
    fit_app <- logistf(y ~ ., data = df, flic = T)
    app_matrix <- model.matrix(object = fit_app$formula, data = df)
    p_app <- 1 / (1 + exp(-app_matrix %*% fit_app$coefficients))
    
  } else if (model == "Firth" & pred_selection == "<0.05") {
    assign("df", as.data.frame(df), envir = .GlobalEnv)
    
    
    model_form <- as.formula(paste0("y ~", paste(colnames(df)[!colnames(df)%in%"y"], collapse = "+" )))
    fit_app <- logistf(formula = model_form, data = df, flic = T)
    fit_app <- backward(fit_app, trace = FALSE)
    
    app_matrix <- model.matrix(object = fit_app$formula, data = df)
    p_app <- 1 / (1 + exp(-app_matrix %*% fit_app$coefficients))
    
    # Get the elements of the design generating mechanism that are
    # belonging to the model after backwards elimination
    # Always get the first element as this is the intercept
    ind <- na.omit(c(1, (as.numeric(str_extract_all(colnames(app_matrix), "(?<=V).*"))
                         # Add 1, because the indices of columns exclude the intercept
                         + 1)))
    
    dgm_par <- dgm_par[ind]
    
    
  } else {
    
    fit_app <- glm(y ~ ., family = "binomial", data = df) 
    p_app <- predict(fit_app, type = "response")
    app_matrix <- model.matrix(object = fit_app$formula, data = df)
    
  }
  
  # Performance measures
  auc_app <- fastAUC(p = p_app, y = df$y)
  R2_app <- pseudo_Rsqrs(p = p_app, y = df$y)
  MAPE_app <- MAPE(p = p_app, dgm_par = dgm_par, iv_matrix = app_matrix) 
  
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
      "ECI_app" = eci_app,
      "MAPE_app" = MAPE_app
    )
}

################################
## Function to obtain results ##
################################

# Given the study, obtain apparent estimands
# for each scenario

get_app_results <- function(study, df) {
  
  # object to store results
  results_app <- list()
  
  # For each scenario within the study
  for (i in 1:nrow(study)) {
    
    # determine which model & pre_selection is specified
    model <- study[i, ]$model
    pred_selection <- study[i, ]$pred_selection
    
    # What do the model parameters look like in the
    # Data generating mechanism?
    dgm_par <- c(study[i, ]$par1, 
                 rep(study[i, ]$par2 * 3, round(0.3 * s1[i, ]$dim)),  # strong
                 rep(study[i, ]$par2,     round(0.5 * s1[i, ]$dim)),  # weaker
                 rep(study[i, ]$par2 * 0, round(0.2 * s1[i, ]$dim)))  # noise
    
    results_app[[i]] <- get_app_estimands(df = df[[i]], model = model, dgm_par = dgm_par, pred_selection = pred_selection)
    
  }
  
  names(results_app) <- c(1:nrow(study))
  return(results_app)
  
}

##########################################
###### Cross-validation  approaches ######
##########################################

###########################
## 10, 5 & 10x10 fold cv ##
###########################

get_cv_estimands <- function(df, model, dgm_par, pred_selection, V, x10 = c(FALSE, TRUE)){

  # Create Balanced CV folds (stratify by outcome)
  .cvFoldsB <- function(Y, V) {  
    
    Y0 <- split(sample(which(Y == "0")), 
                rep(1:V, length = length(which(Y == 0))))
    
    Y1 <- split(sample(which(Y == "1")), 
                rep(1:V, length = length(which(Y == 1))))
    
    folds <- vector("list", length=V)
    
    for (v in seq(V)) { 
      folds[[v]] <- c(Y0[[v]], Y1[[v]])
      }
    
    return(folds)
  
    }
  
  folds <- .cvFoldsB(Y = df$y, V = V)
  
  # Getting predictions depending on model used
  .doFit <- function(V, folds, model, pred_selection){ 
    
    # Train/test glm for each fold
    # Fit model depending on scenario
    # And get predicted probabilities and modelmatrix
    if (model == "Firth" & pred_selection == "none") { # If model = Firth , no predictor selection
      
      fit <- logistf(y ~ ., data = df[-folds[[V]],], flic = T)
      iv_matrix <- model.matrix(object = fit$formula, data = df[folds[[V]],])
      p <- 1 / (1 + exp(-iv_matrix %*% fit$coefficients))
      dgm_par_folds <- dgm_par
      
      } else if (model == "Firth" & pred_selection == "<0.05") { # If model = Firth , backwards predictor selection
        
        df_train <- as.data.frame(df[-folds[[V]],])
        assign("df_train", df_train, envir = .GlobalEnv)
       
        model_form <- as.formula(paste0("y ~", paste(colnames(df)[!colnames(df)%in%"y"], collapse = "+" )))
        fit <- logistf(formula = model_form, data = df_train, flic = T)
        fit <- backward(fit, trace = FALSE)
        
        iv_matrix <- model.matrix(object = fit$formula, data = df[folds[[V]],])
        p <- 1 / (1 + exp(-iv_matrix %*% fit$coefficients))
        
        # Get the elements of the design generating mechanism that are
        # belonging to the model after backwards elimination
        # Always get the first element as this is the intercept
        ind <- na.omit(c(1, (as.numeric(str_extract_all(colnames(iv_matrix), "(?<=V).*"))
                                   # Add 1, because the indices of columns exclude the intercept
                                    + 1)))
       
        dgm_par_folds <- dgm_par[ind]
       
        
      } else {
        
      fit <- glm(y ~ ., family = "binomial", data = df[-folds[[V]],])
      p <- predict(fit, newdata = df[folds[[V]],] , type = "response")
      iv_matrix <- model.matrix(object = fit$formula, data = df[folds[[V]],])
      dgm_par_folds <- dgm_par
      
      }
    
    results <- list(p, iv_matrix, dgm_par_folds)
  }
  
  results <- (lapply(seq(V), .doFit, folds = folds, model = model, pred_selection = pred_selection)) # obtain model results for all folds
  p <- c(unlist(sapply(results, "[[", 1))) # Get out pred values as a single vector
  p[unlist(folds)] <- p #Re-order pred values, so they are in line with y
  
  p_per_fold <- lapply(results, "[[", 1) # getting the predictions per fold
  iv_matrix <- lapply(results, "[[", 2) # iv_matrix per fold
  dgm_par_folds <- lapply(results, "[[", 3) # Dgm_par per fold
  
  ######################  
  ## Obtain estimands ##
  ######################

    ## Empty objects to store results
    slope_folds <- c()
    intercept_folds <- c()
    R2_folds <- c()
    eci_folds <- c()
    MAPE_folds <- c()
    
    # For each fold, calculate calibration intercept & slope, R2, ECI and MAPE
    for (v in 1:V){
      # create object to shorten future code:
      data <- df[unlist(folds[[v]]),]
      ppf <- p_per_fold[[v]]
      dgm_par_ppf <- dgm_par_folds[[v]]
      
      # Calibration slope for folds
      slope_folds[v] <- c(coef(glm(data$y ~ log(ppf/(1-ppf)), family="binomial"))[2])
      # Calibration intercept for folds
      intercept_folds[v] <- coef(glm(data$y ~ offset(log(ppf/(1-ppf))), family="binomial"))
      # R2 for folds
      R2_folds[v] <- pseudo_Rsqrs(p = ppf, y = data$y)
      # ECI for folds
      calout <- loess(y ~ log(ppf/(1-ppf)), data = data)
      eci_folds[v] <- (mean((ppf-fitted(calout))*(ppf-fitted(calout))))*(100)
      # MAPE for folds
      MAPE_folds[v] <- MAPE(p = ppf, iv_matrix = iv_matrix[[v]], dgm_par = dgm_par_ppf)
    }
    
    ## AUC, function by ledell already gives bias corrected se and ci as well
    auc_results <- as.vector(unlist(ci.cvAUC(predictions=p, labels=df$y, folds=folds, confidence=0.95)))[-5]
    names(auc_results) <- c(paste0("AUC_mean_", V, "fcv" ), 
                            paste0("AUC_se_", V, "fcv"), 
                            paste0("AUC_ci_lower_", V, "fcv"),
                            paste0("AUC_ci_upper_", V, "fcv"))
    
    ## Get mean and standard error over all other results and store in a single vector:
    ## Calibration
    intercept <- c(mean(intercept_folds), (sd(intercept_folds)/(sqrt(V) - 1)))
    names(intercept) <- c(paste0("calib_int_mean_", V, "fcv" ), paste0("calib_int_se_", V, "fcv"))
                                    
    slope <-  c(mean(slope_folds), (sd(slope_folds)/(sqrt(V) - 1)))
    names(slope) <- c(paste0("calib_slope_mean_", V, "fcv" ), paste0("calib_slope_se_", V, "fcv"))
    
    ## R2 cox snell
    R2 <- c(mean(R2_folds), (sd(R2_folds)/(sqrt(V) - 1)))
    names(R2) <- c(paste0("R2_mean_", V, "fcv" ), paste0("R2_se_", V, "fcv"))
    
    ## ECI
    eci <- c(mean(eci_folds), (sd(eci_folds)/(sqrt(V) - 1)))
    names(eci) <- c(paste0("eci_mean_", V, "fcv" ), paste0("eci_se_", V, "fcv"))
    
    ## MAPE
    MAPE_results <- c(mean(MAPE_folds), (sd(MAPE_folds)/(sqrt(V) - 1)))
    names(MAPE_results) <- c(paste0("MAPE_mean_", V, "fcv" ), paste0("MAPE_se_", V, "fcv"))
    
    # If the it is anything other than 10x10 cv:
    if (x10 == FALSE){
      
      results <- c(auc_results, intercept, slope, R2, eci, MAPE_results)
    
      } else {
        
      # If it is 10x10 cv, only keep the raw results of each fold and name them:  
      names(intercept_folds) <- c(rep("intercept", length(intercept_folds)))
      names(slope_folds) <- c(rep("slope", length(slope_folds)))
      names(R2_folds) <- c(rep("R2", length(R2_folds)))
      names(eci_folds) <- c(rep("eci", length(eci_folds)))
      names(MAPE_folds) <- c(rep("MAPE", length(MAPE_folds)))
      results <- list("auc" = auc_results, 
                      "intercept" = intercept_folds, 
                      "slope" = slope_folds, 
                      "R2" = R2_folds,
                      "eci" = eci_folds,
                      "MAPE" = MAPE_folds)
      }
 
  return(results)
    
  }

## Function to apply get_cv_estimands to all datasets in a study:
get_cv_results <- function(study, df, V) {
  results_cv <- list() #object to store results
  
  for (i in 1:nrow(study)) {
    print(i)
    model <- study[i, ]$model
    pred_selection <- study[i, ]$pred_selection
    dgm_par <- c(study[i, ]$par1, 
                 rep(study[i, ]$par2 * 3, round(0.3 * study[i, ]$dim)),  # strong
                 rep(study[i, ]$par2,     round(0.5 * study[i, ]$dim)),  # weak
                 rep(study[i, ]$par2 * 0, round(0.2 * study[i, ]$dim)))  # noise
    
    results_cv[[i]] <- get_cv_estimands(df = df[[i]], 
                                        model = model,
                                        dgm_par = dgm_par,
                                        pred_selection = pred_selection,
                                        V = V,
                                        x10 = FALSE)
  }
  names(results_cv) <- c(1:length(df))
  return(results_cv)
}

##############
## 10x10 cv ##
##############

get_10x10_results <- function(study, df, V){
  results_cv <- list()
  
  for (i in 1:length(df)) {
    print(i)
    # Settings for get cv estimands function:
    model <- study[i, ]$model
    pred_selection <- study[i, ]$pred_selection
    data <- df[[i]]
    dgm_par <- c(study[i, ]$par1, 
                 rep(study[i, ]$par2 * 3, round(0.3 * study[i, ]$dim)),  # strong
                 rep(study[i, ]$par2,     round(0.5 * study[i, ]$dim)),  # medium
                 rep(study[i, ]$par2 * 0, round(0.2 * study[i, ]$dim)))  # noise
    
    # obtain the results of each fold separately 10 times and stored for each replication
    results <- replicate(n = 10, # represents the 10x10 part
                         expr = get_cv_estimands(df = data, model = model, 
                                                 dgm_par = dgm_par,
                                                 pred_selection = pred_selection,
                                                 V = V,
                                                 x10 = TRUE), 
                         simplify = F)
    
    # obtain mean results for auc:
    auc_results <- rowMeans(sapply(results, "[[", 1))
    names(auc_results) <- c(paste0("AUC_mean_", V, "x10fcv" ), 
                            paste0("AUC_se_", V, "x10fcv"), 
                            paste0("AUC_ci_lower_", V, "x10fcv"),
                            paste0("AUC_ci_upper_", V, "x10fcv"))
    
    # Get the results as vectors for each estimand:
    intercept_results <- c(sapply(results, "[[", 2))
    slope_results <- c(sapply(results, "[[", 3))
    R2_results <- c(sapply(results, "[[", 4))
    eci_results <- c(sapply(results, "[[", 5))
    MAPE_results <- c(sapply(results, "[[", 6))
    
    ## Get mean and standard error over all other results and store in a single vector:
    ## Calibration
    intercept <- c(mean(intercept_results), (sd(intercept_results)/(sqrt(V) - 1)))
    names(intercept) <- c(paste0("calib_int_mean_", V, "x10fcv" ), paste0("calib_int_se_", V, "fcv"))
    
    slope <-  c(mean(slope_results), (sd(slope_results)/(sqrt(V) - 1)))
    names(slope) <- c(paste0("calib_slope_mean_", V, "x10fcv" ), paste0("calib_slope_se_", V, "fcv"))
    
    ## R2 cox snell
    R2 <- c(mean(R2_results), (sd(R2_results)/(sqrt(V) - 1)))
    names(R2) <- c(paste0("R2_mean_", V, "x10fcv" ), paste0("R2_se_", V, "fcv"))
    
    ## ECI
    eci <- c(mean(eci_results), (sd(eci_results)/(sqrt(V) - 1)))
    names(eci) <- c(paste0("eci_mean_", V, "x10fcv" ), paste0("eci_se_", V, "fcv"))
    
    ## MAPE
    MAPE_results <- c(mean(MAPE_results), (sd(MAPE_results)/(sqrt(V) - 1)))
    names(MAPE_results) <- c(paste0("MAPE_mean_", V, "x10fcv" ), paste0("MAPE_se_", V, "fcv"))
    
    results_cv[[i]] <- c(auc_results, intercept, slope, R2, eci, MAPE_results)
    
  }
  
  names(results_cv) <- c(1:length(df))
  return(results_cv)
}


#########################
## BOOTSTRAP FUNCTIONS ##
#########################





