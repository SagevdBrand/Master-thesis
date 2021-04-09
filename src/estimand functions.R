########################################
########################################
########################################
### Simulation functions
### Estimands

####################
## Error handling ##
####################
## from https://stackoverflow.com/questions/3903157/how-can-i-check-whether-a-function-call-results-in-a-warning
## and https://github.com/MvanSmeden/Beyond-EPV/blob/master/sim_LRMs.R
ErrorsWarnings <- function(expr) {
  
  W <- NULL
  
  w_handler <- function(w){
    W <<- w
    invokeRestart("muffleWarning")
  }
  
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e), warning = w_handler), warning = W)
} 

#########################
#### Model functions ####
#########################
# Lasso and Ridge as defined 
# by Van Calster et al., 2020
# alpha = 1  LASSO
# alpha = 0 = Ridge
Pen_reg_VC <- function(df, alpha = c(0,1), nfolds = c(10, nrow(df))){
  
lambda <- c(exp(seq(log(64), log(0.00001), length.out = 250)), 0) # DOUBLE CHECK THIS
fit <- cv.glmnet(y ~., 
                 family = "binomial",
                 lambda = lambda,
                 alpha = alpha,
                 data = df,
                 nfolds = nfolds
)

return(fit)

}

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

## Two lose lines in the code!

###############
## Tjur's R2 ##
###############
tjur <- function(p, y) {

categories <- unique(y)
m1 <- mean(p[which(y == categories[1])], na.rm = TRUE)
m2 <- mean(p[which(y == categories[2])], na.rm = TRUE)

tjur_d <- abs(m2 - m1)
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

## Two lose lines in the code!

##################
## MAPE & rMSPE ##
##################

MAPE_rMSPE <- function(p,iv_matrix, dgm_par){
  p_true <- 1 / (1 + exp(-iv_matrix %*% dgm_par))
  
  mape <- mean(abs(p_true-p))
  rmspe <- sqrt(mean((p_true-p)^2))
  
  results <-  c("MAPE" = mape, "rMSPE" = rmspe)
  return(results)
}

##################################
###### Apparent performance ######
##################################
get_app_ext_estimands <- function(df, df_val, model, dgm_par, pred_selection){
  
  error_info <- NA
  
  errors_warnings <- ErrorsWarnings({
    
    # obtaining model matrices
    app_matrix <- as.matrix(df[,-which(colnames(df)=="y")])
    app_matrix <- as.matrix(cbind("(Intercept)" = 1, app_matrix))
    
    ext_matrix <- as.matrix(df_val[, -which(colnames(df)=="y")])
    ext_matrix <- as.matrix(cbind("(Intercept)" = 1, ext_matrix))
    # Fit model depending on scenario
    # And get predicted probabilities
    
    if (model == "Firth" | model == "ML"){ 
      # Specifiy model formula (necessary for logistf)
      model_form <- as.formula(paste0("y ~", paste(colnames(df)[!colnames(df)%in%"y"], collapse = "+" )))
      
      ## If model is Firth without backwards elimination:
      if (model == "Firth" & pred_selection == "none") {
        
        fit_app <- logistf(formula = model_form, data = df, flic = T, firth = T)
        
      #   ## If model is Firth with backwards elimination:  
      # } else if (model == "Firth" & pred_selection == "<0.157") {
      #   df_app <- df
      #   assign("df_app", as.data.frame(df_app), envir = .GlobalEnv)
      #   
      #   fit_app <- logistf(formula = model_form, data = df_app, flic = T, firth = T)
      #   fit_app <- backward(fit_app, trace = FALSE)
      #   
      #   # Check whether any predictors have been selected.
      #   if (var(fit_app$linear.predictors) == 0) {
      #     error_info <- paste(toString(error_info), paste0("No predictors selected -> no calibration slope"), sep = " + ")
        }
        
      ## If model is using ML :  
      if (model == "ML"){
        
        if (pred_selection == "none") {
        
        fit_app <- glm(y ~ ., family = "binomial", data = df) 
        
      } else if (model == "ML" & pred_selection == "<0.157"){
        
        fit_app <- glm(y ~ ., family = "binomial", data = df) 
        step(fit_app, trace = FALSE)
        
        # Check whether any predictors have been selected.
           if (var(fit_app$linear.predictors) == 0) {
             error_info <- paste(toString(error_info), paste0("No predictors selected -> no calibration slope"), sep = " + ")
           } # close check for predictors
        
      } # close Which specific type of ML 
      
      # Check for separation
      if(any(sqrt(diag(summary(fit_app)$cov.unscaled)*summary(fit_app)$dispersion) > 70)){
        error_info <- paste(toString(error_info), paste0("Data separation might have occured"), sep = " + ")
      } # Close separation detection
      } # Close ML
      # Predictions:
      # Apparent
      fit_app_matrix <- model.matrix(object = fit_app$formula, data = df)
      p_app <- 1 / (1 + exp(-fit_app_matrix %*% fit_app$coefficients))
      
      # External
      fit_ext_matrix <- model.matrix(object = fit_app$formula, data = df_val)
      p_ext <- 1 / (1 + exp(-fit_ext_matrix %*% fit_app$coefficients))
      
      # If the model uses penalized regression  
    } else if (model == "Lasso" | model == "Ridge") {
      
      alpha <- ifelse(model == "Lasso", 1, 0) # alpha = 1 for Lasso
      
      # Make sure that there are at least 8 events or no-events:
      if (sum(df$y) < 8 | sum(1-df$y) <8 ){
        # Fit the model depending on alpha with LOOCV
        fit_app <-  Pen_reg_VC(df = df, alpha = alpha, nfolds = nrow(df))
        error_info <- paste(toString(error_info), 
                            paste0("Too few (non-)events for tuning -> LOOCV"),
                            sep = " + ")
      } else {
        # Otherwise just use 10-f-cv for tuning of lambda
        fit_app <-  Pen_reg_VC(df = df, alpha = alpha, nfolds = 10)
      }  # close error handling too few events
      
      ## Check whether predictors have been selected
      # Create linear predictor
      lp <- predict(fit_app, newdata = df, s = "lambda.min") 
      
      if (var(lp) == 0) {
        error_info <- paste(toString(error_info), paste0("No predictors selected -> no calibration slope"), sep = " + ")
      }
      # Apparent
      p_app <- predict(fit_app, newdata = df, s = "lambda.min", type = "response")
      # External
      p_ext <- predict(fit_app, newdata = df_val, s = "lambda.min", type = "response")
      
    } else  if (model == "RF"| model == "CART") {
      
      if (model == "RF"){
        # pre-defining a tune grid
        rf_grid <-
          expand.grid(.mtry = seq(1, ncol(df) - 1, by = 2),
                      .splitrule = "gini",
                      .min.node.size = 1)
        
        df_factor <- df
        df_factor$y <- ifelse(df$y == 1, "yes", "no")
        
        df_val_factor <- df_val
        df_val_factor$y <- ifelse(df_val$y == 1, "yes", "no")
        
        # Check number of events
        if (sum(df$y) < 8 | sum(1-df$y) <8 ){
          # Fit the model with LOOCV
          fit_app <- caret::train(as.factor(y) ~.,
                              data = df_factor, 
                              method = 'ranger',
                              tuneGrid = rf_grid,
                              num.trees = 250,
                              trControl = trainControl(method = "LOOCV", classProbs = T))
          # Save the error:
          error_info <- paste(toString(error_info), 
                              paste0("Too few (non-)events for tuning -> LOOCV"),
                              sep = " + ")
          
        } else {
          
        fit_app <- caret::train(as.factor(y) ~.,
                                data = df_factor, 
                                method = 'ranger',
                                tuneGrid = rf_grid,
                                num.trees = 250,
                                trControl = trainControl(method = "cv", classProbs = T)
        )
        }
        
        p_app <- predict(fit_app$finalModel, df_factor, type = "response")$predictions[,2]
        p_ext <- predict(fit_app$finalModel, df_val_factor, type = "response")$predictions[,2]
        
        # Otherwise fit CART
        
      } else {
        # Pre-defining a tuning grid
        rpart_grid <- expand.grid(cp = seq(0, 0.01, by = 0.00034))
        
        # Fitting the model
        # Check number of events
        if (sum(df$y) < 8 | sum(1-df$y) <8 ){
        fit_app <- caret::train(as.factor(y) ~.,
                                data = df,
                                method = 'rpart',
                                minsplit = 1,
                                minbucket = 1,
                                tuneGrid = rpart_grid,
                                trControl = trainControl(method = "LOOCV"))
        
        # Save the error:
        error_info <- paste(toString(error_info), 
                            paste0("Too few (non-)events for tuning -> LOOCV"),
                            sep = " + ")
        
        # Check whether there were more than 1 splits
        if (nrow(fit_app$finalModel$cptable) == 1) {
          error_info <- paste(toString(error_info), paste0("No predictors selected -> no calibration slope"), sep = " + ")
          } # close handling if there is only one split in the tree
        
        } else {
          
          fit_app <- caret::train(as.factor(y) ~.,
                                  data = df,
                                  method = 'rpart',
                                  minsplit = 1,
                                  minbucket = 1,
                                  tuneGrid = rpart_grid,
                                  trControl = trainControl(method = "cv"))
          
        
        # Check whether there were more than 1 splits
        if (nrow(fit_app$finalModel$cptable) == 1) {
          error_info <- paste(toString(error_info), paste0("No predictors selected -> no calibration slope"), sep = " + ")
          } # Close handling if there is only one split in the tree
        
        } # Close handling for checking events.
        
        # Apparent
        p_app <- predict(fit_app$finalModel, newdata = df, type = "prob")[,2]
        # External 
        p_ext <- predict(fit_app$finalModel, newdata = df_val, type = "prob")[,2]
      } # close which tree

      
    # } else if (model == "SVM"){
    #   
    #   # Pre-defining a tuning grid
    #   svm_grid <- expand.grid(.C = seq(0.001, 10, length.out = 4),
    #                           .sigma = seq(0.001, 1, length.out = 4))
    #   # Fitting the model, making sure that no output
    #   # is printed!
    #   invisible(capture.output(fit_app <- caret::train(as.factor(y) ~.,
    #                                                    data = df,
    #                                                    method = 'svmRadial',
    #                                                    tuneGrid = svm_grid,
    #                                                    trControl = trainControl(method = "cv"),
    #                                                    prob.model = TRUE
    #   )))
    #   
    #   # Apparent
    #   p_app <- predict(fit_app$finalModel, newdata =  df[, which(!colnames(df) %in% "y")], type = "prob")[,2]
    #   # External 
    #   p_ext <- predict(fit_app$finalModel, newdata = df_val[, which(!colnames(df_val) %in% "y")], type = "prob")[,2]
    #   
    # } else if (model == "ANN"){
    #   # Pre-defining a tuning grid
    #   nnet_grid <- expand.grid(size = seq(from = 1, to = 10, by = 2),
    #                            decay = seq(from = 0.1, to = 2, by = 0.4))
    #   # Fitting the model
    #   fit_app <- caret::train(as.factor(y) ~.,
    #                           data = df, 
    #                           method = 'nnet', 
    #                           tuneGrid = nnet_grid, 
    #                           trace = F, 
    #                           trControl = trainControl(method = "cv"),
    #                           preProcess = c("center", "scale"))
    #   # Checking for convergence
    #   if (fit_app$finalModel$convergence == 1){
    #     error_info <- paste(toString(error_info), paste0("nnet: Maximum number of iterations was reached"), sep = " + ")
    #   } # close check for convergence
    #   # Apparent
    #   p_app <- predict(fit_app$finalModel, newdata = df, type='raw')
    #   # External 
    #   p_ext <- predict(fit_app$finalModel, newdata = df_val, type = 'raw')
      
    } # Close model if else statements
    
    # Check whether any of the predictions are 0
    if (any(p_app == 1) | any(p_app == 0) | any(p_ext == 1) | any(p_ext == 0)) {
      error_info <- paste(toString(error_info), paste0("probabilities of 0 or 1 occured"), sep = " + ")
    }
    
    p_app <-  ifelse(p_app == 0, 0.000001, p_app)
    p_app <- ifelse(p_app == 1,  0.999999, p_app)
    
    p_ext <-  ifelse(p_ext == 0, 0.000001, p_ext)
    p_ext <- ifelse(p_ext == 1,  0.999999, p_ext)
    
    # Obtain apparent results
    auc_app <- fastAUC(p = p_app, y = df$y)
    R2_app <- pseudo_Rsqrs(p = p_app, y = df$y)
    MAPE_rMSPE_app <- MAPE_rMSPE(p = p_app, dgm_par = dgm_par, iv_matrix = app_matrix) 
    tjur_app <- tjur(p = p_app, y = df$y)
    slope_app <- c(coef(glm(df$y ~ log(p_app/(1-p_app)), family="binomial"))[2])
    intercept_app <- coef(glm(df$y ~ offset(log(p_app/(1-p_app))), family="binomial"))
    # ECI calculations
    calout <- loess(y ~ log(p_app/(1-p_app)), data = df, span = 10)
    ecir <- mean((p_app - mean(df$y))*(p_app - mean(df$y)))
    ecim <- (mean((p_app-fitted(calout))*(p_app-fitted(calout))))
    eci_app <- ecim/ecir
    
    # obtain external results
    auc_ext <- fastAUC(p = p_ext, y = df_val$y)
    R2_ext <- pseudo_Rsqrs(p = p_ext, y = df_val$y)
    MAPE_rMSPE_ext <- MAPE_rMSPE(p = p_ext, dgm_par = dgm_par, iv_matrix = ext_matrix) 
    tjur_ext <- tjur(p = p_ext, y = df_val$y)
    slope_ext <- c(coef(glm(df_val$y ~ log(p_ext/(1-p_ext)), family="binomial"))[2])
    intercept_ext <- coef(glm(df_val$y ~ offset(log(p_ext/(1-p_ext))), family="binomial"))
    # ECI calculations
    calout <- loess(y ~ log(p_ext/(1-p_ext)), data = df_val, span = 10)
    ecir <- mean((p_ext - mean(df_val$y))*(p_ext - mean(df_val$y)))
    ecim <- (mean((p_ext-fitted(calout))*(p_ext-fitted(calout))))
    eci_ext <- ecim/ecir
    
  # all results together
  results <- list(
    c(
      "Apparent",
      auc_app,
      intercept_app,
      slope_app,
      tjur_app,
      R2_app,
      eci_app,
      MAPE_rMSPE_app['MAPE'],
      MAPE_rMSPE_app['rMSPE'],
      "error_info" = error_info
    ),
    c(
      "External",
      auc_ext,
      intercept_ext,
      slope_ext,
      tjur_ext,
      R2_ext,
      eci_ext,
      MAPE_rMSPE_ext['MAPE'],
      MAPE_rMSPE_ext['rMSPE'],
      "error_info" = error_info
    ),
    p_app # Save the prediction to be used in bootstrap!
  )
  
  }) #Close warnings and error check
  
  # If there are warnings, paste those in the error_info
  if (!is.null(errors_warnings$warning)) {
    results[[1]]["error_info"] <- paste(toString(error_info), toString(errors_warnings$warning), sep = " + ")
    results[[2]]["error_info"] <- paste(toString(error_info), toString(errors_warnings$warning), sep = " + ")
  }  
  # If there are error messages print this in the error_info
  if ("message" %in% errors_warnings$value) {
    results[[1]]["error_info"] <- paste(toString(error_info), toString(errors_warnings$value$message), sep = " + ")
    results[[2]]["error_info"] <- paste(toString(error_info), toString(errors_warnings$value$message), sep = " + ")
  }
  
  return(results)
} # End function

################################
## Function to obtain results ##
################################

# Given the study, obtain apparent estimands
# for each scenario

get_app_ext_results <- function(study, df, df_val, studyname) {
  
  # objects to store results
  results_app <- as.data.frame(matrix(NA, nrow = nrow(study), 
                                      ncol = length(results_estimands_names), 
                                      dimnames = list(c(), results_estimands_names)))
  results_ext <- as.data.frame(matrix(NA, nrow = nrow(study), 
                                      ncol = length(results_estimands_names), 
                                      dimnames = list(c(), results_estimands_names)))
  p_app_preds <- list()
  # Add the information of the study to the respective columns
  results_app[, which(colnames(results_app) %in% study_info)] <- study[, which(colnames(study) %in% study_info)]
  results_ext[, which(colnames(results_ext) %in% study_info)] <- study[, which(colnames(study) %in% study_info)]
  # The n for external validation will be overwritten below
  
  # For each scenario within the study
  for (i in 1:nrow(study)) {
    # Fill in some information for each scenario:
    results_app[i, 'scenario'] <- study$scenario[i]
    results_ext[i, 'scenario'] <- study$scenario[i]
    # Overwrite n of the development data
    results_ext[i, 'n'] <- nrow(df_val[[i]])
    # How many events observed in validation set?
    results_ext[i, 'observed events'] <- sum(df_val[[i]]$y)
    
    # Check for each scenario whether there are events sampled
    if (any(str_detect(names(df[[i]]),"Error: No events sampled"))){
      
      # If no events were sampled, then the following will be the result:
      results_app[i, 'error_info'] <- c("Error: No events sampled")
      results_app[i, 'approach'] <- c("Apparent")
      # For the external results, there will be the following messages
      results_ext[i, "error_info"] <- c("Error: No events sampled in development set")
      results_ext[i, "approach"] <- c("External")
      # and all the other columns of the estimands remain NA
      
    } else {
      
      # Otherwise do this:
      results_app[i, 'observed events'] <- sum(df[[i]]$y)
      # determine which model & pre_selection is specified
      model <- study[i, ]$model
      pred_selection <- study[i, ]$pred_selection
      
      # Obtain dgm, depending on the dimensionality,
      # noise contribution, and therefore dgm settings:
      if (study[i, ]$noise == "default"){
        dgm_par <-
          c(study[i, ]$par1, 
            rep(study[i, ]$par2 * 3, round(0.3 * study[i, ]$dim)), 
            rep(study[i, ]$par2, round(0.5 *  study[i, ]$dim)), 
            rep(study[i, ]$par2 * 0, round(0.2 * study[i, ]$dim)))
      }
      if (study[i, ]$noise == "half"){
        dgm_par <-
          c(study[i, ]$par1, 
            rep(study[i, ]$par2 * 3, round(1/6 * study[i, ]$dim)), # strong
            rep(study[i, ]$par2, round(2/6 *  study[i, ]$dim)),    # weak
            rep(study[i, ]$par2 * 0, round(3/6 * study[i, ]$dim))) # noise
      } 
      if (study[i, ]$noise == "none"){
        dgm_par <-
          c(study[i, ]$par1, 
            rep(study[i, ]$par2 * 3, round(1/6 * study[i, ]$dim)), # strong
            rep(study[i, ]$par2, round(5/6 *  study[i, ]$dim))     # weak
          )
      }
      
      # Fill the columns with apparent results, no SE!
      results <- get_app_ext_estimands(
        df = df[[i]],
        df_val = df_val[[i]],
        model = model,
        dgm_par = dgm_par,
        pred_selection = pred_selection
      )
      
      # Save the results of each list into the correct matrix
      results_app[i, which(colnames(results_app) %in% apparent_col_names)] <- results[[1]]
      results_ext[i, which(colnames(results_ext) %in% apparent_col_names)] <- results[[2]]
      # Bind the matrices together
      results_app_ext <- rbind(results_app, results_ext)
      # Fill in details of the study
      results_app_ext$study <- study$study
      
      # Save the p_apps to the environment
      p_app_preds[[i]] <- c(results[[3]])
      
    } # end else statement
  } # end for loop
  assign(paste0("p_app_", studyname), p_app_preds, envir = .GlobalEnv)
  return(results_app_ext)
} # end function

##########################################
###### Cross-validation  approaches ######
##########################################

###########################
## 10, 5 & 10x10 fold cv ##
###########################

get_cv_estimands <- function(df, model, dgm_par, pred_selection, V, x10 = c(FALSE, TRUE)){
  
  # If nothing 'bad' happens, the error info remains as NA
  error_info <- NA
  original_V <- V
  
  errors_warnings <- ErrorsWarnings({
    # Create Balanced CV folds (stratify by outcome)
    # Code by Erin Ledell:
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
    
    #  Check for var(y) == 0 within each fold:
    if (any(lapply(seq(V), function(V) var(df[folds[[V]],]$y)) == 0)) {
      # Find out which fold has no events:
      bad_folds_ind <- which(lapply(seq(V), function(V) var(df[folds[[V]],]$y)) == 0)
      bad_folds <- folds[bad_folds_ind]
      # Remove the bad folds from the folds:
      folds <- folds[-bad_folds_ind]
      V <- V - length(bad_folds)
      # Turn the data to NA
      df[unlist(bad_folds),] <- NA
      # And exclude it from the data, keeping row indices intact
      df <- na.exclude(df)
      # Return informative error message, 
      # adding to any that might have been there before
      error_info <- paste(toString(error_info), paste0("No events sampled in folds ",
                                                       paste(bad_folds_ind, collapse = " & ")),
                          sep = " + ")
    } # Close error handling folds without events 
    # (should not happen because of stratified folds, but just to be sure)
    
    # Getting predictions depending on model used
    .doFit <- function(V, folds, model, pred_selection){ 
      
      # Create the training and test folds
      df_train <- as.data.frame(df[-(which(rownames(df) %in% folds[[V]])),])
      df_test <-  as.data.frame(df[which(rownames(df) %in% folds[[V]]),])
      assign("df_train", df_train, envir = .GlobalEnv)
      assign("df_test", df_test, envir = .GlobalEnv)
      
      iv_matrix <- as.matrix(df_test[,-which(colnames(df_test)=="y")])
      iv_matrix <- as.matrix(cbind("(Intercept)" = 1, iv_matrix))
      
      # Train/test glm for each fold
      # Fit model depending on scenario
      # And get predicted probabilities and modelmatrix
      if (model == "Firth" | model == "ML"){ 
        
        # Specifiy model formula (necessary for logistf)
        model_form <- as.formula(paste0("y ~", paste(colnames(df)[!colnames(df)%in%"y"], collapse = "+" )))
        
        if (model == "Firth") { # If model = Firth 
          fit <- logistf(formula = model_form, data = df_train, flic = T, firth = T)
        } # Close if for firth and no predictor selection
        
        # else if (model == "Firth" & pred_selection == "<0.157"){
        #   # Fit Firth model with Backward elimination
        #   fit <- logistf(formula = model_form, data = df_train, flic = T, firth = T)
        #   fit <- backward(fit, trace = FALSE)
        #   
        #   # Check whether any predictors have been selected at all. 
        #   if (var(fit$linear.predictors) == 0) {
        #     error_info <- paste(toString(error_info), paste0("No predictors selected -> no calibration slope"), sep = " + ")
        #     
        #   } # close error handling predictor selection
      
      
      if (model == "ML"){ 
          
          if(pred_selection == "none") {
            
            fit <- glm(y ~ ., family = "binomial", data = df_train) 
            
          } else if (pred_selection == "<0.157"){
            
            fit <- glm(y ~ ., family = "binomial", data = df_train) 
            step(fit, trace = FALSE)
            if (var(fit$linear.predictors) == 0) {
              error_info <-paste(toString(error_info), paste0("No predictors selected -> no calibration slope"), sep = " + ")
            } # close error handling predictor selction
          } # Close for specific models
          
          # Check for separation
          if(any(sqrt(diag(summary(fit)$cov.unscaled)*summary(fit)$dispersion) > 70)){
            error_info <- paste(toString(error_info), paste0("Data separation might have occured"), sep = " + ")
          } # Close separation detection
       
      }
        # Obtain model matrix
        fit_matrix <- model.matrix(object = fit$formula, data = df_test)
        # Obtain predictions
        p <- 1 / (1 + exp(-fit_matrix %*% fit$coefficients))
        
        # If the model is penalized regression  
      } else if (model == "Lasso" | model == "Ridge") {
        
        alpha <- ifelse(model == "Lasso", 1, 0) # alpha = 1 for Lasso
        
        # Make sure that there are at least 8 events or no-events:
        if (sum(df_train$y) < 8 | sum(1-df_train$y) <8 ){
          # Fit the model depending on alpha with LOOCV
          fit <-  Pen_reg_VC(df = df_train, alpha = alpha, nfolds = nrow(df_train))
          error_info <- paste(toString(error_info), 
                              paste0("Too few (non-)events for tuning -> LOOCV"),
                              sep = " + ")
        } else {
          # Otherwise just use 10-f-cv for tuning of lambda
          fit <-  Pen_reg_VC(df = df_train, alpha = alpha, nfolds = 10)
        }  # close error handling too few events
        
        # Check whether predictors have been selected
        # Create linear predictor
        lp <- predict(fit, newdata = df_train, s = "lambda.min") 
        if (var(lp) == 0) {
          error_info <- paste(toString(error_info), paste0("No predictors selected -> no calibration slope"), sep = " + ")
        }
        
        # Obtain predictions:
        p <- predict(fit, newdata = df_test, s = "lambda.min", type = "response")
        
      } else if (model == "RF"| model == "CART") {
        
        if (model == "RF"){
          # pre-defining a tune grid
          rf_grid <-
            expand.grid(.mtry = seq(1, ncol(df_train) - 1, by = 2),
                        .splitrule = "gini",
                        .min.node.size = 1)
          
          df_train_factor <- df_train
          df_train_factor$y <- ifelse(df_train$y == 1, "yes", "no")
          df_test_factor <- df_test
          df_test_factor$y <- ifelse(df_test$y == 1, "yes", "no")
          
          # Make sure that there are at least 8 events or no-events:
          if (sum(df_train$y) < 8 | sum(1-df_train$y) <8 ){
            # Fit the model with LOOCV
            fit <- caret::train(as.factor(y) ~.,
                                data = df_train_factor, 
                                method = 'ranger',
                                tuneGrid = rf_grid,
                                num.trees = 250,
                                trControl = trainControl(method = "LOOCV", classProbs = T))
            # Save the error:
            error_info <- paste(toString(error_info), 
                                paste0("Too few (non-)events for tuning -> LOOCV"),
                                sep = " + ")
            
          } else {
          # Fitting the appropriate model:
          fit <- caret::train(as.character(as.factor(y)) ~.,
                              data = df_train_factor, 
                              method = 'ranger',
                              tuneGrid = rf_grid,
                              num.trees = 250,
                              trControl = trainControl(method = "cv", classProbs = T))
          
         
          } # Close if else for events check
          p <- predict(fit$finalModel, df_test_factor, type = "response")$predictions[,2]
          
        } else { # the model is CART
          rpart_grid <- expand.grid(cp = seq(0, 0.01, by = 0.00034))
          # Make sure that there are at least 8 events or no-events:
          if (sum(df_train$y) < 8 | sum(1-df_train$y) <8 ){
          # Fit the model with LOOCV
          fit <- caret::train(as.factor(y) ~.,
                              data = df_train,
                              method = 'rpart',
                              minsplit = 1,
                              minbucket = 1,
                              tuneGrid = rpart_grid,
                              trControl = trainControl(method = "LOOCV"))
          # Save the error:
          error_info <- paste(toString(error_info), 
                              paste0("Too few (non-)events for tuning -> LOOCV"),
                              sep = " + ")
          
         } else {
            fit <- caret::train(as.factor(y) ~.,
                                data = df_train,
                                method = 'rpart',
                                minsplit = 1,
                                minbucket = 1,
                                tuneGrid = rpart_grid,
                                trControl = trainControl(method = "cv"))
          } # Close events check
          
          # Check whether there were more than 1 splits
          if (nrow(fit$finalModel$cptable) == 1) {
            error_info <- paste(toString(error_info), paste0("No predictors selected -> no calibration slope"), sep = " + ")
          } # close check for only 1 splits
          p <- predict(fit$finalModel, newdata = df_test, type = "prob")[,2]
        } # Close fitting tree models 
        
      # } else if ( model == "SVM") {
      #   # pre-defining a tune grid
      #   svm_grid <- expand.grid(.C = seq(0.001, 10, length.out = 4),
      #                           .sigma = seq(0.001, 1, length.out = 4))
      #   
      #   # Fit svm with Radial Basis Function kernel
      #   invisible(capture.output(fit <- caret::train(as.factor(y) ~.,
      #                                                data = df_train,
      #                                                method = 'svmRadial',
      #                                                tuneGrid = svm_grid,
      #                                                trControl = trainControl(method = "cv"),
      #                                                prob.model = TRUE
      #   )))
      #   
      #   # Obtain predictions
      #   p <- predict(fit$finalModel, newdata = df_test[!colnames(df_test)%in%"y"], type = "prob")[,2]
      #   
      # } else if ( model == "ANN") {
      #   # pre-defining a tune grid
      #   nnet_grid <- expand.grid(size = seq(from = 1, to = 10, by = 2),
      #                            decay = seq(from = 0.1, to = 2, by = 0.4))
      #   
      #   # Fit Neural network
      #   fit <- caret::train(as.factor(y) ~.,
      #                       data = df_train, 
      #                       method = 'nnet', 
      #                       tuneGrid = nnet_grid, 
      #                       trace = F, 
      #                       trControl = trainControl(method = "cv"),
      #                       preProcess = c("center", "scale"))
      #   
      #   # Check convergence of model
      #   if (fit$finalModel$convergence == 1){
      #     error_info <- paste(toString(error_info), paste0("nnet: Maximum number of iterations was reached"), sep = " + ")
      #   }
      #   
      #   # Obtain predictions
      #   p <- predict(fit$finalModel, newdata = df_test)
        
      } # Close model if else statements
      results <- list(p, iv_matrix)
    } # Close function .doFit
    
    results <- (lapply(seq(V), .doFit, folds = folds, model = model, pred_selection = pred_selection)) # obtain model results for all folds
    p <- c(unlist(sapply(results, "[[", 1))) # Get out pred values as a single vector
    p[which(rownames(df) %in% unlist(folds))] <- p #Re-order pred values, so they are in line with y
    
    p_per_fold <- lapply(results, "[[", 1) # getting the predictions per fold
    iv_matrix <- lapply(results, "[[", 2) # iv_matrix per fold
    
    ######################  
    ## Obtain estimands ##
    ######################
    
    ## Empty objects to store results
    auc_folds <- c()
    slope_folds <- c()
    intercept_folds <- c()
    tjur_folds <- c()
    R2_folds <- c()
    eci_folds <- c()
    mape_folds <- c()
    rmspe_folds <- c()
    
    # For each fold, calculate calibration intercept & slope, R2, ECI and MAPE
    for (v in 1:V){
      # create objects to shorten future code:
      data <- df[unlist(which(rownames(df) %in% folds[[v]])),]
      ppf <- p_per_fold[[v]]
      ppf_true <- 1 / (1 + exp(-iv_matrix[[v]] %*% dgm_par))
      
      # Check whether any of the predictions are 0
      if (any(ppf == 1) | any(ppf == 1)) {
        error_info <- paste(toString(error_info), paste0("probabilities of 0 or 1 occured"), sep = " + ")
      }
      ppf <-  ifelse(ppf == 0, 0.000001, ppf)
      ppf <- ifelse(ppf == 1,  0.999999, ppf)
      
      # AUC for folds:
      auc_folds[v] <- fastAUC(p = ppf, y = data$y)
      # Calibration slope for folds
      slope_folds[v] <- c(coef(glm(data$y ~ log(ppf/(1-ppf)), family="binomial"))[2])
      # Calibration intercept for folds
      intercept_folds[v] <- coef(glm(data$y ~ offset(log(ppf/(1-ppf))), family="binomial"))
      # Tjur's R2 for folds
      tjur_folds[v] <- tjur(p = ppf, y = data$y)
      # R2 for folds
      R2_folds[v] <- pseudo_Rsqrs(p = ppf, y = data$y)
      # ECI for folds
      calout <- loess(y ~ log(ppf/(1-ppf)), data = data, span = 10)
      ecir <- mean((ppf - mean(data$y))*(ppf - mean(data$y)))
      ecim <- mean((ppf-fitted(calout))*(ppf-fitted(calout)))
      eci_folds[v] <- ecim/ecir
      # MAPE for folds
      mape_folds[v] <- mean(abs(ppf_true-ppf))
      # rMSPE for folds
      rmspe_folds[v] <- sqrt(mean((ppf_true-ppf)^2))
    }
    
    ## Get mean and standard error over all results and store in a single vector:
    ## AUC
    auc_results <- c(mean(auc_folds), (sd(auc_folds)/(sqrt(V))))
    ## Calibration
    intercept <- c(mean(intercept_folds), (sd(intercept_folds)/(sqrt(V))))
    slope <-  c(mean(slope_folds), (sd(slope_folds)/(sqrt(V))))
    ## Tjur
    tjur_results <- c(mean(tjur_folds), (sd(tjur_folds)/sqrt(V)))
    ## R2 cox snell
    R2 <- c(mean(R2_folds), (sd(R2_folds)/(sqrt(V))))
    ## ECI
    eci <- c(mean(eci_folds), (sd(eci_folds)/(sqrt(V))))
    ## MAPE
    MAPE_results <- c(mean(mape_folds), (sd(mape_folds)/(sqrt(V))))
    ## rMSPE
    rMSPE_results <-  c(mean(rmspe_folds), (sd(rmspe_folds)/(sqrt(V))))
    
 
  # If the it is anything other than 10x10 cv:
  if (x10 == FALSE){
    
    results <- c(paste0(original_V, " fold cross-validation"), auc_results, intercept, slope, tjur_results, R2, 
                 eci,
                 MAPE_results, rMSPE_results, "error_info" = error_info)
    
  } else {
    # 
    # If it is 10x10 cv, only keep the raw results of each fold:  
    results <- list("auc" = auc_folds, 
                    "intercept" = intercept_folds, 
                    "slope" = slope_folds, 
                    "Tjur" = tjur_folds,
                    "R2" = R2_folds,
                    "eci" = eci_folds,
                    "MAPE" = mape_folds,
                    "rMSPE" = rmspe_folds,
                    error_info)
  } # Close else statement (when 10k or 5k CV is used)
  
  }) # Close ErrorsWarnings
  
  # If there are warnings, paste those in the error_info
  if (!is.null(errors_warnings$warning)) {
    
    results['error_info'] <- paste(toString(error_info), toString(errors_warnings$warning), sep = " + ")
    
  }  
  # If there are error messages print this in the error_info
  if ("message" %in% errors_warnings$value) {
    
    results['error_info'] <- paste(toString(error_info), toString(errors_warnings$value$message), sep = " + ")
    
  }
  
  return(results)
  
} # Close function

## Function to apply get_cv_estimands to all datasets in a study:
get_cv_results <- function(study, df, V, studyname) {
  
  # results matrix:
  results_cv <- as.data.frame(matrix(NA, nrow = nrow(study), ncol = length(results_estimands_names), dimnames = list(c(), results_estimands_names)))
  
  # Fill in details of the study
  results_cv$study <- study$study
  # And of each scenario:
  results_cv[, which(colnames(results_cv) %in% study_info)] <- study[, which(colnames(study) %in% study_info)]
  
  
  for (i in 1:nrow(study)) {
    # Paste the correct scenario number:
    results_cv[i, 'scenario'] <- study$scenario[i]
    
    # Check for each scenario whether there are events sampled
    if (any(str_detect(names(df[[i]]),"Error: No events sampled"))){
      
      # If no events were sampled, then the following will be the result:
      results_cv[i, 'error_info'] <- c("Error: No events sampled")
      results_cv[i, 'approach'] <- paste0(V, " fold cross-validation")
      # and all the other columns of the estimands remain NA
      
    } else {
      
      results_cv[i, 'observed events'] <- sum(df[[i]]$y)
      # Else, go along with obtaining the results
      model <- study[i, ]$model
      pred_selection <- study[i, ]$pred_selection
      
      # Obtain dgm, depending on the dimensionality,
      # noise contribution, and therefore dgm settings
      if (study[i, ]$noise == "default"){
        dgm_par <-
          c(study[i, ]$par1, 
            rep(study[i, ]$par2 * 3, round(0.3 * study[i, ]$dim)), 
            rep(study[i, ]$par2, round(0.5 *  study[i, ]$dim)), 
            rep(study[i, ]$par2 * 0, round(0.2 * study[i, ]$dim)))
      }
      
      if (study[i, ]$noise == "half"){
        dgm_par <-
          c(study[i, ]$par1, 
            rep(study[i, ]$par2 * 3, round(1/6 * study[i, ]$dim)), # strong
            rep(study[i, ]$par2, round(2/6 *  study[i, ]$dim)),    # weak
            rep(study[i, ]$par2 * 0, round(3/6 * study[i, ]$dim))) # noise
      } 
      
      if (study[i, ]$noise == "none"){
        dgm_par <-
          c(study[i, ]$par1, 
            rep(study[i, ]$par2 * 3, round(1/6 * study[i, ]$dim)), # strong
            rep(study[i, ]$par2, round(5/6 *  study[i, ]$dim))     # weak
          )
      }
      
      results_cv[i, which(colnames(results_cv) %in% iv_colnames)]  <-
        get_cv_estimands(
          df = df[[i]],
          model = model,
          dgm_par = dgm_par,
          pred_selection = pred_selection,
          V = V,
          x10 = FALSE
        )
      
    } # close else statement
  } # close for loop
  return(results_cv)
}


##############
## 10x10 cv ##
##############

get_10x10_results <- function(study, df, V, studyname){
  
  # results matrix:
  results_cv <- as.data.frame(matrix(NA, nrow = nrow(study), ncol = length(results_estimands_names), dimnames = list(c(), results_estimands_names)))
  
  # Fill in details of the study
  results_cv$study <- study$study
  # And of each scenario:
  results_cv[, which(colnames(results_cv) %in% study_info)] <- study[, which(colnames(study) %in% study_info)]
  
  for (i in 1:length(df)) {
    # Add which scenario we are working on:
    results_cv[i, 'scenario'] <- study$scenario[i]
    
    # Settings for get cv estimands function:
    model <- study[i, ]$model
    pred_selection <- study[i, ]$pred_selection
    data <- df[[i]]
    
    # Check whether there are events sampled
    if (any(str_detect(names(data),"Error: No events sampled") == TRUE)){
      
      # If no events were sampled, then the following will be the result:
      results <- c("Error: No events sampled")
      return(results) 
      
    } else {
      
      results_cv[i, 'observed events'] <- sum(data$y)
      
      # Obtain dgm, depending on the dimensionality,
      # noise contribution, and therefore dgm settings
      if (study[i, ]$noise == "default"){
        dgm_par <-
          c(study[i, ]$par1, 
            rep(study[i, ]$par2 * 3, round(0.3 * study[i, ]$dim)), 
            rep(study[i, ]$par2, round(0.5 *  study[i, ]$dim)), 
            rep(study[i, ]$par2 * 0, round(0.2 * study[i, ]$dim)))
      }
      
      if (study[i, ]$noise == "half"){
        dgm_par <-
          c(study[i, ]$par1, 
            rep(study[i, ]$par2 * 3, round(1/6 * study[i, ]$dim)), # strong
            rep(study[i, ]$par2, round(2/6 *  study[i, ]$dim)),    # weak
            rep(study[i, ]$par2 * 0, round(3/6 * study[i, ]$dim))) # noise
      } 
      
      if (study[i, ]$noise == "none"){
        dgm_par <-
          c(study[i, ]$par1, 
            rep(study[i, ]$par2 * 3, round(1/6 * study[i, ]$dim)), # strong
            rep(study[i, ]$par2, round(5/6 *  study[i, ]$dim))     # weak
          )
      }
      
      # obtain the results of each fold separately 10 times and stored for each replication
      results <- replicate(n = 10, # represents the 10x10 part
                           expr = get_cv_estimands(df = data,
                                                   model = model,
                                                   dgm_par = dgm_par,
                                                   pred_selection = pred_selection,
                                                   V = V,
                                                   x10 = TRUE), 
                           simplify = F)
      
      # Get the results as vectors for each estimand:
      auc_results <- c(sapply(results, "[[", 1))
      intercept_results <- c(sapply(results, "[[", 2))
      slope_results <- c(sapply(results, "[[", 3))
      tjur_results <- c(sapply(results, "[[", 4))
      R2_results <- c(sapply(results, "[[", 5))
      eci_results <- c(sapply(results, "[[", 6))
      MAPE_results <- c(sapply(results, "[[", 7))
      rMSPE_results <- c(sapply(results, "[[", 8))
      
      ## Get mean and standard error over all other results and store in a single vector:
      auc_results <- c(mean(auc_results), (sd(auc_results)/(sqrt(V))))
      
      ## Calibration
      intercept <- c(mean(intercept_results), (sd(intercept_results)/(sqrt(V))))
      slope <-  c(mean(slope_results), (sd(slope_results)/(sqrt(V))))
      
      ## Tjur R2
      tjur_results <- c(mean(tjur_results), (sd(tjur_results)/(sqrt(V))))
      
      ## R2 cox snell
      R2 <- c(mean(R2_results), (sd(R2_results)/(sqrt(V))))
      
      ## ECI
      eci <- c(mean(eci_results), (sd(eci_results)/(sqrt(V))))
      
      ## MAPE
      MAPE_results <- c(mean(MAPE_results), (sd(MAPE_results)/(sqrt(V))))
      
      ## rMSPE
      rMSPE_results <- c(mean(rMSPE_results), (sd(rMSPE_results)/(sqrt(V))))
      
      # Taking care of error messages returned for each repetition. 
      error_info <- paste(c(sapply(results, "[[", 9)), collapse = " + ")
      error_info <- str_remove_all(error_info, "NA \\+ |NA | \\+ NA")
      
      ## Fill results matrix:
      results_cv[i, which(colnames(results_cv) %in% iv_colnames)] <-
        c(
          "10x10 fold cross-validation",
          auc_results,
          intercept,
          slope,
          tjur_results,
          R2,
          eci,
          MAPE_results,
          rMSPE_results,
          error_info
        )
      
    } # close for loop
  } # close else loop
  
  return(results_cv)
}


#########################
## BOOTSTRAP FUNCTIONS ##
#########################
## Obtain lambda
# Bootstrap estimate of optimism: Harrell's
get_lambda <- function(data){
  
  train_results <- data %>% dplyr::select(ends_with("_train"))
  df_results <- data %>% dplyr::select(ends_with("_df"))
  
  lambda <- colMeans(train_results - df_results, na.rm = TRUE)
  names(lambda) <- sub("*_train", "_optimism", names(lambda))
  return(lambda)
}

# ## Obtain gamma
# # Where pred_app is paste0("p_app_", studyname)[[i]]
# get_gamma <- function(data, p_df){
#   p1 <- sum(data$y)/nrow(data)
#   q1_classes <- ifelse(p_df > 0.5, 1, 0)
#   q1 <- sum(q1_classes)/nrow(data)
#   gamma_hat <- p1*(1 - q1) + (1 - p1) * q1
#   return(gamma_hat)
# }


get_bootstrap_estimands <- function(df, model, dgm_par, pred_selection, pred_app, theta_app, nboot){  
  
  error_info <- NA
  errors_warnings <- ErrorsWarnings({  
    
    # create empty matrix for results
    results <- as.data.frame(matrix(NA, nrow = nboot, ncol = 24, dimnames = list(c(), c(
      paste0(estimands_names, c(rep("_df",8), rep("_train",8), rep("_test",8)))
    ))))
    
    # For the eventual results
    results_matrix <- as.data.frame(matrix(NA, ncol = length(iv_colnames), nrow = 3, dimnames = list(c(),
                                                                                                     c(iv_colnames))))
    
    df_matrix <- as.matrix(df[,-which(colnames(df)=="y")])
    df_matrix <- as.matrix(cbind("(Intercept)" = 1, df_matrix))
    
    for (b in 1:nboot) {
      # Obtain indices for the bootstrap sample: "in the bag" = itb
      itb <- sample(nrow(df), size = nrow(df), replace = T) # indices
      
      # Create bootstrap sample and the out sample:
      df_train <- df[itb,] # The bootstrap sample
      df_test <- df[-itb, ] # The test/out sample
      
      # Check whether any events have been sampled
      # within the training sample
      if (var(df_train$y) == 0) {
        error_info <- paste(toString(error_info), paste0(b, " No events sampled in training sample"), sep = " + ")
        next
      }
      
      assign("df_train", df_train, envir = .GlobalEnv)
      assign("df_test", df_test, envir = .GlobalEnv)
      
      ## Create model matrices:
      train_matrix <- as.matrix(df_train[,-which(colnames(df_train)=="y")])
      train_matrix <- as.matrix(cbind("(Intercept)" = 1, train_matrix))
      
      test_matrix <- as.matrix(df_test[,-which(colnames(df_test)=="y")])
      test_matrix <- as.matrix(cbind("(Intercept)" = 1, test_matrix))
      
      # Train/test glm for each fold
      # Fit model depending on scenario
      # And get predicted probabilities and modelmatrix
      if (model == "Firth"| model == "ML" ) {
        
        # Specifiy model formula (necessary for logistf)
        model_form <- as.formula(paste0("y ~", paste(colnames(df)[!colnames(df)%in%"y"], collapse = "+" )))
        
        if (model == "Firth"){
          # Fit Firth model without Backward elimination
          fit <- logistf(formula = model_form, data = df_train, flic = T, firth = T)
        } 
          
        # } else if (model == "Firth" & pred_selection == "<0.157"){
        #   # Fit Firth model with Backward elimination
        #   fit <- logistf(formula = model_form, data = df_train, flic = T, firth = T)
        #   fit <- backward(fit, trace = FALSE)
        #   
        #   # Check whether any predictors have been selected at all. 
        #   if (var(fit$linear.predictors) == 0) {
        #     error_info <- paste(toString(error_info), paste0("No predictors selected -> no calibration slope"), sep = " + ")
        #   }
          
        if (model == "ML"){
          
          if (pred_selection == "none") {
            
            fit <- glm(y ~ ., family = "binomial", data = df_train) 
            
            } else if (pred_selection == "<0.157"){
            
            fit <- glm(y ~ ., family = "binomial", data = df_train) 
            step(fit, trace = FALSE)
            
            if (var(fit$linear.predictors) == 0) {
              error_info <-paste(toString(error_info), paste0("No predictors selected -> no calibration slope"), sep = " + ")
            } # close error handling predictor selection
            
            } # Close for specific models
          
          # Check for separation of ML models
          if(any(sqrt(diag(summary(fit)$cov.unscaled)*summary(fit)$dispersion) > 70)){
            error_info <- paste(toString(error_info), paste0("Data separation might have occured"), sep = " + ")
          } # Close separation detection
        } # close if for ML
        

        # original sample
        df_fit_matrix <- model.matrix(object = fit$formula, data = df)
        p_df <- 1 / (1 + exp(-df_fit_matrix %*% fit$coefficients))
        
        # bootstrapped sample
        train_fit_matrix <- model.matrix(object = fit$formula, data = df_train)
        p_train <- 1 / (1 + exp(-train_fit_matrix %*% fit$coefficients))
        
        # out/test sample
        test_fit_matrix <- model.matrix(object = fit$formula, data = df_test)
        p_test <- 1 / (1 + exp(-test_fit_matrix %*% fit$coefficients))
        
      } else if (model == "Lasso" | model == "Ridge") {
        
        alpha <- ifelse(model == "Lasso", 1, 0) # alpha = 1 for Lasso
        # Make sure that there are at least 8 events or no-events:
        if (sum(df_train$y) < 8 | sum(1-df_train$y) <8 ){
          fit <-  Pen_reg_VC(df = df_train, alpha = alpha, nfolds = nrow(df_train))
          error_info <- paste(toString(error_info), paste0("Too few (non-)events for tuning -> LOOCV"), sep = " + ")
        } else {
          fit <-  Pen_reg_VC(df = df_train, alpha = alpha, nfolds = 10)
        } # close error handling 
        
        # Check whether predictors have been selected
        # Create linear predictor
        lp <- predict(fit, newdata = df_train, s = "lambda.min") 
        if (var(lp) == 0) {
          error_info <- paste(toString(error_info), paste0("No predictors selected -> no calibration slope"), sep = " + ")
        }
        
        # original sample
        p_df <- predict(fit, newdata = df, s = "lambda.min", type = "response")
        # bootstrapped sample
        p_train <- predict(fit, newdata = df_train, s = "lambda.min", type = "response")
        # out/test sample
        p_test <- predict(fit, newdata = df_test, s = "lambda.min", type = "response")
        
      } else if (model == "RF"| model == "CART") {
        
        if (model == "RF"){
          # pre-defining a tune grid
          rf_grid <-
            expand.grid(.mtry = seq(1, ncol(df_train) - 1, by = 2),
                        .splitrule = "gini",
                        .min.node.size = 1)
          
          df_factor <- df
          df_factor$y <-  ifelse(df$y == 1, "yes", "no")
          df_train_factor <- df_train
          df_train_factor$y <- ifelse(df_train$y == 1, "yes", "no")
          df_test_factor <- df_test
          df_test_factor$y <- ifelse(df_test$y == 1, "yes", "no")
          
          # Check number of events:
          if (sum(df_train$y) < 8 | sum(1-df_train$y) <8 ){
            
            # Fitting the appropriate model:
            fit <- caret::train(
              as.factor(y) ~ .,
              data = df_train_factor,
              method = 'ranger',
              num.trees = 250,
              tuneGrid = rf_grid,
              trControl = trainControl(method = "LOOCV", classProbs=TRUE)
            )
            
            # Save the error:
            error_info <- paste(toString(error_info), 
                                paste0("Too few (non-)events for tuning -> LOOCV"),
                                sep = " + ")
            
          } else {
          # Fitting the appropriate model:
          fit <- caret::train(
            as.factor(y) ~ .,
            data = df_train_factor,
            method = 'ranger',
            num.trees = 250,
            tuneGrid = rf_grid,
            trControl = trainControl(method = "cv", classProbs=TRUE)
          )
        
          } # Close events check
          # Get probabilities out:
          p_df <- predict(fit$finalModel, df_factor, type = "response")$predictions[,2]
          p_train <- predict(fit$finalModel, df_train_factor, type = "response")$predictions[,2] 
          p_test <- predict(fit$finalModel, df_test_factor, type = "response")$predictions[,2]
          
        } else {
          # Pre-defining a tuning grid
          rpart_grid <- expand.grid(cp = seq(0, 0.01, by = 0.00034))
          
          # Check number of events:
          if (sum(df_train$y) < 8 | sum(1-df_train$y) <8 ){
            
            # Fitting the appropriate model using loocv:
            fit <- caret::train(
              as.factor(y) ~ .,
              data = df_train,
              method = 'rpart',
              minsplit = 1,
              minbucket = 1,
              tuneGrid = rpart_grid,
              trControl = trainControl(method = "LOOCV")
            )
            
            # Save the error:
            error_info <- paste(toString(error_info), 
                                paste0("Too few (non-)events for tuning -> LOOCV"),
                                sep = " + ")
          } else {
            
          fit <- caret::train(
            as.factor(y) ~ .,
            data = df_train,
            method = 'rpart',
            minsplit = 1,
            minbucket = 1,
            tuneGrid = rpart_grid,
            trControl = trainControl(method = "cv")
          )
          } # Close if else for events check
          
          # Check whether there were more than 1 splits
          if (nrow(fit$finalModel$cptable) == 1) {
            error_info <- paste(toString(error_info), paste0("No predictors selected -> no calibration slope"), sep = " + ")
          } # close check for only 1 splits
          
          # Get probabilities out:
          p_df <- predict(fit$finalModel, newdata = df, type = "prob")[,2]
          p_train <- predict(fit$finalModel, newdata = df_train, type = "prob")[,2] 
          p_test <- predict(fit$finalModel, newdata = df_test, type = "prob")[,2]
        } # Close if else for which specific tree model
        
       
        
      # } else if (model == "SVM"){
      #   # pre-defining a tune grid
      #   svm_grid <- expand.grid(.C = seq(0.001, 10, length.out = 4),
      #                           .sigma = seq(0.001, 1, length.out = 4))
      #   
      #   # Fit svm with Radial Basis Function kernel
      #   invisible(capture.output(fit <- caret::train(as.factor(y) ~.,
      #                                                data = df_train,
      #                                                method = 'svmRadial',
      #                                                tuneGrid = svm_grid,
      #                                                trControl = trainControl(method = "cv"),
      #                                                prob.model = TRUE
      #   )))
      #   p_df <- predict(fit$finalModel, newdata = df[!colnames(df)%in%"y"], type = "prob")[,2]
      #   p_train <- predict(fit$finalModel, newdata = df_train[!colnames(df_train)%in%"y"], type = "prob")[,2]
      #   p_test <- predict(fit$finalModel, newdata = df_test[!colnames(df_test)%in%"y"], type = "prob")[,2]
      #   
      # } else { # Otherwise it is a neural network
      #   nnet_grid <- expand.grid(size = seq(from = 1, to = 10, by = 2),
      #                            decay = seq(from = 0.1, to = 2, by = 0.4))
      #   
      #   # Fit Neural network
      #   fit <- caret::train(as.factor(y) ~.,
      #                       data = df_train, 
      #                       method = 'nnet', 
      #                       tuneGrid = nnet_grid, 
      #                       trace = F, 
      #                       trControl = trainControl(method = "cv"),
      #                       preProcess = c("center", "scale"))
      #   
      #   # Check convergence of model
      #   if (fit$finalModel$convergence == 1){
      #     error_info <- paste(toString(error_info), paste0("nnet: Maximum number of iterations was reached"), sep = " + ")
      #   } # close check for convergence
      #   
      #   # Obtain predictions
      #   p_df <- predict(fit$finalModel, newdata = df)
      #   p_train <- predict(fit$finalModel, newdata = df_train)
      #   p_test <- predict(fit$finalModel, newdata = df_test)
        
      } # close  machine learning models if elses
      
      ## Check whether any predictions are 0:
      if (any(p_df == 0) | any(p_train == 0) | any(p_test == 0)){
        error_info <- paste(toString(error_info), paste0("probabilities of 0 occured"), sep = " + ")
      } # close check for predictions of 0.
      
      if (any(p_df == 1) | any(p_train == 1) | any(p_test == 1)){
        error_info <- paste(toString(error_info), paste0("probabilities of 1 occured"), sep = " + ")
      } # close check for predictions of 0.
      
      p_df <-  ifelse(p_df == 0, 0.000001, p_df)
      p_df <- ifelse(p_df == 1,  0.999999, p_df)
      
      p_train <-  ifelse(p_train == 0, 0.000001, p_train)
      p_train <- ifelse(p_train == 1,  0.999999, p_train)
      
      p_test <-  ifelse(p_test == 0, 0.000001, p_test)
      p_test <- ifelse(p_test == 1,  0.999999, p_test)
      
      # Obtain 'true' predictions:
      p_df_true <- 1 / (1 + exp(-df_matrix %*% dgm_par))
      p_train_true <- 1 / (1 + exp(-train_matrix %*% dgm_par))
      p_test_true <- 1 / (1 + exp(-test_matrix %*% dgm_par))
      
      ## Filling the matrices above withe the respective estimands
      # AUC
      results$auc_df[b] <- fastAUC(p = p_df, y = df$y)
      results$auc_train[b] <- fastAUC(p = p_train, y = df_train$y)
      results$auc_test[b] <- fastAUC(p = p_test, y = df_test$y)
      
      # Calibration slope
      results$calib_int_df[b] <- coef(glm(df$y ~ offset(log(p_df/(1-p_df))), family="binomial"))
      results$calib_int_train[b] <- coef(glm(df_train$y ~ offset(log(p_train/(1-p_train))), family="binomial"))
      results$calib_int_test[b] <- coef(glm(df_test$y ~ offset(log(p_test/(1-p_test))), family="binomial"))
      
      # Calibration intercept
      results$calib_slope_df[b] <- c(coef(glm(df$y ~ log(p_df/(1-p_df)), family="binomial"))[2])
      results$calib_slope_train[b] <- c(coef(glm(df_train$y ~ log(p_train/(1-p_train)), family="binomial"))[2])
      results$calib_slope_test[b] <-  c(coef(glm(df_test$y ~ log(p_test/(1-p_test)), family="binomial"))[2])
      
      # Tjur's R2
      results$Tjur_df[b] <- tjur(p = p_df, y = df$y)
      results$Tjur_train[b] <- tjur(p = p_train, y = df_train$y)
      results$Tjur_test[b] <- tjur(p = p_test, y = df_test$y)
      
      # R2
      results$R2_CS_df[b] <- pseudo_Rsqrs(p = p_df, y = df$y)
      results$R2_CS_train[b] <- pseudo_Rsqrs(p = p_train, y = df_train$y) 
      results$R2_CS_test[b] <- pseudo_Rsqrs(p = p_test, y = df_test$y)
      
      ## ECI
      # DF
      calout <- loess(y ~ log(p_df/(1-p_df)), data = df, span = 10)
      ecir <- mean((p_df - mean(df$y))*(p_df - mean(df$y)))
      ecim <- mean((p_df-fitted(calout))*(p_df-fitted(calout)))
      results$eci_df[b] <- ecim/ecir
      
      # df_train
      calout <- loess(y ~ log(p_train/(1-p_train)), data = df_train, span = 10)
      ecir <- mean((p_train - mean(df_train$y))*(p_train - mean(df_train$y)))
      ecim <- (mean((p_train-fitted(calout))*(p_train-fitted(calout))))
      results$eci_train[b] <- ecim/ecir
      
      # df-test
      calout <- loess(y ~ log(p_test/(1-p_test)), data = df_test, span = 10)
      ecir <- mean((p_test - mean(df_test$y))*(p_test - mean(df_test$y)))
      ecim <- (mean((p_test-fitted(calout))*(p_test-fitted(calout))))
      results$eci_test[b] <- ecim/ecir
      
      # MAPE
      results$mape_df[b] <- mean(abs(p_df_true-p_df))
      results$mape_train[b] <- mean(abs(p_train_true-p_train))
      results$mape_test[b] <- mean(abs(p_test_true-p_test))
      
      #rMSPE
      results$rmspe_df[b]  <- sqrt(mean((p_df_true-p_df)^2))
      results$rmspe_train[b]  <- sqrt(mean((p_train_true-p_train)^2))
      results$rmspe_test[b]  <- sqrt(mean((p_test_true-p_test)^2))
      
    } # Close bootstrap for loop
    
    results_matrix$approach <- c("Harrell's bootstrap", ".632 bootstrap", ".632+ bootstrap")
    
    theta_test <- results %>% dplyr::select(ends_with("_test"))
    test_se <- apply(theta_test, 2, sd)/(sqrt(nboot))
    results_matrix[1 ,which(colnames(results_matrix) %in% estimands_se_names)] <- test_se
    results_matrix[2 ,which(colnames(results_matrix) %in% estimands_se_names)] <- test_se
    results_matrix[3 ,which(colnames(results_matrix) %in% estimands_se_names)] <- test_se
    
    #####################
    # Harrell's results #
    #####################
    lambda <- get_lambda(results)
    Harrell_results <- as.matrix(t(theta_app) - lambda)
    colnames(Harrell_results) <- estimands_names
    
    results_matrix[1,which(colnames(results_matrix) %in% estimands_names)] <- Harrell_results
    
    ##################
    # .632 bootstrap #
    ##################
    ## Average test estimate ##
    theta_test <- colMeans(theta_test, na.rm = TRUE)
    efron_.632_results <- as.matrix(t((0.368 * theta_app) + (0.632 * theta_test)))
    colnames(efron_.632_results) <- estimands_names
    
    results_matrix[2,which(colnames(results_matrix) %in% estimands_names)] <- efron_.632_results
    
    ###################
    # .632+ bootstrap #
    ###################
    gamma <- c(0.5,0.5,0,0,0,1,0.5,0.25) # No-information predictive performance

    # Relative overfitting 
    R <- (theta_test - theta_app) / (gamma - theta_app)
    
    # Exceptional cases:
    # If AUC is smaller than 0.5:
    # https://github.com/nomahi/predboot/blob/master/R/pred.ML 
    if (theta_test['auc_test'] < 0.5) {R[1] <- 1}
    # If the calibration slope is negative
    if (theta_test['calib_slope_test'] < 0) {R[3] <- 1}
    # If one of the pseudo-R^2 is negative:
    if (theta_test['Tjur_test'] < 0) {R[4] <- 1}
    if (theta_test['R2_CS_test'] < 0) {R[5] <- 1}
   
    # In those cases where the apparent performance is smaller than the test performance
    # For AUC, and both Pseudo-R2. 
    app_smaller_test <- which(theta_app[c(1,4,5)] < theta_test[c(1,4,5)])
    R[c(app_smaller_test)] <- 0
    
    # In those cases where the apparent performance is higher than test performance
    # For the ECI, rMSPE and MAPE:
    app_bigger_test <- which(theta_app[c(6,7,8)] > theta_test[c(6,7,8)])
    R[c(app_smaller_test)] <- 0
    
    # In those cases where the absolute distance of the cal_int to 0 is bigger for the apparent performance than the test performance:
    if(abs(theta_test[2]) < abs(theta_app[2])) {R[2] <- 0}
    # In those cases where the absolute distance of the cal_slope to 1 is bigger for the apparent performance than the test performance:
    if(abs(theta_test[3] - 1) <  abs(theta_app[3] - 1)) {R[3] <- 0}
    
    # If R is negative, return absolute value
    R <-  ifelse( R < 0, abs(R), R)
    # Then return all that are bigger than 1 to 1. 
    R <-  ifelse( R > 1, 1, R)
    # weight
    W <-  .632 / (1 - .368 * R)
    
    # results:
    efron_.632_plus_results <- as.matrix(t((1-W) * theta_app) + (W * theta_test)) 
    colnames(efron_.632_plus_results) <- estimands_names
    
    results_matrix[3,which(colnames(results_matrix) %in% estimands_names)] <- efron_.632_plus_results
    
  }) #Close warnings and error check
  
  
  # If there are warnings, paste those in the error_info
  if (!is.null(errors_warnings$warning)) {
    error_info <- paste(toString(error_info), toString(errors_warnings$warning), sep = " + ")
  }  
  # If there are error messages print this in the error_info
  if ("message" %in% errors_warnings$value) {
    error_info <- paste(toString(error_info), toString(errors_warnings$value$message), sep = " + ")
  }
  
  results_matrix$error_info <- error_info 
  return(results_matrix)
} # Close function  


### And to apply this for each scenario in a study and obtaining the results:  
get_bootstrap_results <- function(study, df, nboot, studyname) {
  
  # Matrix for results (all 3 methods so should be three
  # times as large as a single matrix)
  results_boot <- as.data.frame(matrix(NA, nrow =3*nrow(study), 
                                       ncol = length(results_estimands_names), 
                                       dimnames = list(c(), results_estimands_names)))
  
  # Fill in details of the study
  results_boot$study <- study$study
  
  # The results should be pasted in rows of 3 as well
  # Therefore for each scenario, another index is necessary
  # so it corresponds with each 3rd row in the results
  results_i_start <- seq(1, nrow(results_boot), by = 3)
  results_i_end <- seq(3, nrow(results_boot), by = 3)
  
  for (i in 1:nrow(study)){
    # Fill in some details:
    results_boot[c(results_i_start[i]:results_i_end[i]),
                 'scenario'] <- study$scenario[i]
    
    results_boot[c(results_i_start[i]:results_i_end[i]),
                 'observed events'] <- sum(df[[i]]$y)
    
    results_boot[c(results_i_start[i]:results_i_end[i]),
                 which(colnames(results_boot) %in% study_info)] <- study[i, which(colnames(study) %in% study_info)]
    
    # Check whether the data is actually useful
    if (any(str_detect(names(df[[i]]),"Error: No events sampled") == TRUE)) {
      # If no events were sampled, then the following will be
      results_boot[c(results_i_start[i]:results_i_end[i]), 'error_info'] <- c("Error: No events sampled")
      results_boot[c(results_i_start[i]:results_i_end[i]), 'approach'] <- c("Harrell's bootstrap",
                                                                            ".632 bootstrap",
                                                                            ".632+ bootstrap")
    } else {
      
      # Obtain necessary apparent information results
      scen_num <- as.numeric(paste0(str_extract_all(study$scenario[i], "[:digit:]", simplify = T), collapse = ""))
      pred_app <- get(paste0("p_app_", studyname))
      results_app <- as.data.frame(results_app_ext) %>% filter(approach == "Apparent")
      theta_app <- results_app %>% .[,(colnames(results_app) %in% estimands_names)]
      theta_app <- as.numeric(theta_app)
      
      # Obtain necessary information to obtain results:
      results_boot[c(results_i_start[i]:results_i_end[i]),'model'] <- as.character(study[i, ]$model)
      model <- study[i, ]$model
      
      results_boot[c(results_i_start[i]:results_i_end[i]),'pred_selection'] <- as.character(study[i, ]$pred_selection)
      pred_selection <- study[i, ]$pred_selection
      
      # Obtain dgm, depending on the dimensionality,
      # noise contribution, and therefore dgm settings
      if (study[i, ]$noise == "default"){
        dgm_par <-
          c(study[i, ]$par1, 
            rep(study[i, ]$par2 * 3, round(0.3 * study[i, ]$dim)), 
            rep(study[i, ]$par2, round(0.5 *  study[i, ]$dim)), 
            rep(study[i, ]$par2 * 0, round(0.2 * study[i, ]$dim)))
      }
      
      if (study[i, ]$noise == "half"){
        dgm_par <-
          c(study[i, ]$par1, 
            rep(study[i, ]$par2 * 3, round(1/6 * study[i, ]$dim)), # strong
            rep(study[i, ]$par2, round(2/6 *  study[i, ]$dim)),    # weak
            rep(study[i, ]$par2 * 0, round(3/6 * study[i, ]$dim))) # noise
      } 
      
      if (study[i, ]$noise == "none"){
        dgm_par <-
          c(study[i, ]$par1, 
            rep(study[i, ]$par2 * 3, round(1/6 * study[i, ]$dim)), # strong
            rep(study[i, ]$par2, round(5/6 *  study[i, ]$dim))     # weak
          )
      } # Close if statements for dgm_pars
      
      # Paste the results in the results matrix:
      # The rows are per 3, 1:3, 4:6 etc.
      # The columnnames are the iv_colnames
      results_boot[c(results_i_start[i]:results_i_end[i]), 
                   which(colnames(results_boot) %in% iv_colnames)] <- get_bootstrap_estimands(df = df[[i]],
                                                                                              model = model,
                                                                                              dgm_par = dgm_par,
                                                                                              pred_selection = pred_selection,
                                                                                              pred_app = pred_app,
                                                                                              theta_app = theta_app,
                                                                                              nboot = nboot)
      
    }  # Close if else when checking for events in the data
    
  } # close for loop
  return(results_boot)
} # close get_bootstrap_results_function

