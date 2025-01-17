####### Testing study 2 ###########
source("./src/study scenarios.R")
source("./src/data generation functions.R")
source("./src/estimand functions.R")

s3_data <- generate_data(s3, validation = FALSE)
s3_val_data <- readRDS(study_3_val_data)

df <- s3_data[[19]]
df_val <- s3_val_data[[19]]
study <- s3
i <- 19
model <- "SVM"

dgm_par <-
  c(study[i, ]$par1, 
    rep(study[i, ]$par2 * 3, round(0.3 * study[i, ]$dim)), 
    rep(study[i, ]$par2, round(0.5 *  study[i, ]$dim)), 
    rep(study[i, ]$par2 * 0, round(0.2 * study[i, ]$dim)))

pred_selection <- "none"
V <- 5


get_app_ext_estimands <- function(df, df_val, model, dgm_par, pred_selection){
  
  error_info <- NA
  
  errors_warnings <- ErrorsWarnings({
    
    # Fit model depending on scenario
    # And get predicted probabilities and modelmatrix
    
    ## If model is Firth without backwards elimination:
    if (model == "Firth" & pred_selection == "none") {
      
      fit_app <- logistf(y ~ ., data = df, flic = T)
      
      # Apparent
      app_matrix <- model.matrix(object = fit_app$formula, data = df)
      p_app <- 1 / (1 + exp(-app_matrix %*% fit_app$coefficients))
      
      # External
      ext_matrix <- model.matrix(object = fit_app$formula, data = df_val)
      p_ext <- 1 / (1 + exp(-ext_matrix %*% fit_app$coefficients))
      
      # Save coefficients
      coefs <- fit_app$coefficients
      
      # Save dgm 
      dgm_par_app <- dgm_par
      
      ## If model is Firth with backwards elimination:  
    } else if (model == "Firth" & pred_selection == "<0.05") {
      
      df_app <- df
      assign("df_app", as.data.frame(df_app), envir = .GlobalEnv)
      
      model_form <- as.formula(paste0("y ~", paste(colnames(df_app)[!colnames(df_app)%in%"y"], collapse = "+" )))
      fit_app <- logistf(formula = model_form, data = df_app, flic = T)
      fit_app <- backward(fit_app, trace = FALSE)
      
      # Check whether any predictors have been selected.
      if ((var(fit_app$linear.predictors) == 0) == TRUE) {
        error_info <- paste(toString(error_info), paste0("No predictors selected -> no calibration slope"), sep = " + ")
      }
      
      app_matrix <- model.matrix(object = fit_app$formula, data = df_app)
      p_app <- 1 / (1 + exp(-app_matrix %*% fit_app$coefficients))
      
      # External
      ext_matrix <- model.matrix(object = fit_app$formula, data = df_val)
      p_ext <- 1 / (1 + exp(-ext_matrix %*% fit_app$coefficients))
      
      # Save coefficients
      coefs <- fit_app$coefficients
      
      # Get the elements of the design generating mechanism that are
      # belonging to the model after backwards elimination
      # Always get the first element as this is the intercept
      ind <- na.omit(c(1, (as.numeric(str_extract_all(colnames(app_matrix), "(?<=V).*"))
                           # Add 1, because the indices of columns exclude the intercept
                           + 1)))
      
      dgm_par_app <- dgm_par[ind]
      
      ## If model is using ML :  
    } else if (model == "ML") {
      
      fit_app <- glm(y ~ ., family = "binomial", data = df) 
      p_app <- predict(fit_app, type = "response")
      app_matrix <- model.matrix(object = fit_app$formula, data = df)
      
      # External
      ext_matrix <- model.matrix(object = fit_app$formula, data = df_val)
      p_ext <- predict(fit_app, newdata = df_val, type = "response")
      
      # Save coefficients
      coefs <- fit_app$coefficients
      
      # Save dgm
      dgm_par_app <- dgm_par
      
    } else if (model == "Lasso") {
      
      fit_app <-  Pen_reg_VC(df = df, alpha = 1, nfolds = 10)
      
      ## Check whether predictors have been selected
      # Create linear predictor
      lp <- predict(fit_app, newdata = df, s = "lambda.min") 
      
      if (var(lp) == 0) {
        error_info <- paste(toString(error_info), paste0("No predictors selected -> no calibration slope"), sep = " + ")
      }
      
      # Create model matrix:
      # First retrieve the coefficients
      coefs <- as.matrix(coef(fit_app, s = "lambda.min"))
      # Remove those that do are 0
      coefs[coefs == 0] <- NA
      coefs <- na.exclude(coefs)
      
      # Apparent
      p_app <- predict(fit_app, newdata = df, s = "lambda.min", type = "response")
      # Then use the names of the coefficients to get the model matrix:
      df_app_X <- df[colnames(df) %in% rownames(coefs)]
      # Make model matrix by binding a column of 1 to the above
      app_matrix <- as.matrix(cbind("(Intercept)" = 1, df_app_X))
      
      
      # External
      p_ext <- predict(fit_app, newdata = df_val, s = "lambda.min", type = "response")
      # Then use the names of the coefficients to get the model matrix:
      df_ext_X <- df_val[colnames(df_val) %in% rownames(coefs)]
      # Make model matrix by binding a column of 1 to the above
      ext_matrix <- as.matrix(cbind("(Intercept)" = 1, df_ext_X))
      
      # Get the elements of the design generating mechanism that are
      # belonging to the model after backwards elimination
      # Always get the first element as this is the intercept
      ind <- na.omit(c(1, (as.numeric(str_extract_all(colnames(app_matrix), "(?<=V).*"))
                           # Add 1, because the indices of columns exclude the intercept
                           + 1)))
      
      # Save the correct part of dgm_par
      dgm_par_app <- dgm_par[ind]
      
      # If the model is Ridge    
    } else if (model == "Ridge") {
      
      fit_app <-  Pen_reg_VC(df = df, alpha = 0, nfolds = 10)
      
      # First retrieve the coefficients
      coefs <- as.matrix(coef(fit_app, s = "lambda.min"))
      
      # Apparent
      # Get the predictions 
      p_app <- predict(fit_app, newdata = df, s = "lambda.min", type = "response")
      # Then use the names of the coefficients to get the model matrix:
      df_app_X <- df[colnames(df) %in% rownames(coefs)]
      # Make model matrix by binding a column of 1 to the above
      app_matrix <- as.matrix(cbind("(Intercept)" = 1, df_app_X))
      
      
      # External
      p_ext <- predict(fit_app, newdata = df_val, s = "lambda.min", type = "response")
      # Then use the names of the coefficients to get the model matrix:
      df_ext_X <- df_val[colnames(df_val) %in% rownames(coefs)]
      # Make model matrix by binding a column of 1 to the above
      ext_matrix <- as.matrix(cbind("(Intercept)" = 1, df_ext_X))
      
      
      # Save dgm_par to dgm_par_app
      dgm_par_app <- dgm_par
      
    } else if (model == "RF") {
      
      rf_grid <- expand.grid(.mtry = seq(1, ncol(df)-1, by = 1))
      
      fit_app <- caret::train(as.factor(y) ~.,
                          data = df, 
                          method = 'rf',
                          tuneGrid = rf_grid,
                          trControl = trainControl(method = "cv")
      )
      
      
      coefs <- as.matrix(fit_app$coefnames)
      rownames(coefs) <- (fit_app$coefnames)
      
      # Apparent
      p_app <- predict(fit_app$finalModel, newdata = df, type = "prob")[,2]
     
      # Then use the names of the coefficients to get the model matrix:
      df_app_X <- df[colnames(df) %in% rownames(coefs)]
      # Make model matrix by binding a column of 1 to the above
      app_matrix <- as.matrix(cbind("(Intercept)" = 1, df_app_X))
      
      # External 
      p_ext <- predict(fit_app$finalModel, newdata = df_val, type = "prob")[,2]
      # Then use the names of the coefficients to get the model matrix:
      df_ext_X <- df_val[colnames(df_val) %in% rownames(coefs)]
      # Make model matrix by binding a column of 1 to the above
      ext_matrix <- as.matrix(cbind("(Intercept)" = 1, df_ext_X))
      
      
      # Save dgm
      dgm_par_app <- dgm_par
    
    } else if (model == "CART"){

      rpart_grid <- expand.grid(.cp = seq(0, 0.3, by = 0.01))
      
      fit_app <- caret::train(as.factor(y) ~.,
                              data = df,
                              method = 'rpart',
                              tuneGrid = rpart_grid,
                              trControl = trainControl(method = "cv"))
      
      
      coefs <- as.matrix(fit_app$coefnames)
      rownames(coefs) <- (fit_app$coefnames)
      
      # Apparent
      p_app <- predict(fit_app$finalModel, newdata = df, type = "prob")[,2]
      
      # Then use the names of the coefficients to get the model matrix:
      df_app_X <- df[colnames(df) %in% rownames(coefs)]
      # Make model matrix by binding a column of 1 to the above
      app_matrix <- as.matrix(cbind("(Intercept)" = 1, df_app_X))
      
      # External 
      p_ext <- predict(fit_app$finalModel, newdata = df_val, type = "prob")[,2]
      # Then use the names of the coefficients to get the model matrix:
      df_ext_X <- df_val[colnames(df_val) %in% rownames(coefs)]
      # Make model matrix by binding a column of 1 to the above
      ext_matrix <- as.matrix(cbind("(Intercept)" = 1, df_ext_X))
      
      # Save dgm
      dgm_par_app <- dgm_par
      
    } else if (model == "SVM" ){
      
      
      # Pre-defining a tuning grid
      svm_grid <- expand.grid(.C = seq(0.001, 10, length.out = 5),
                              .sigma = seq(0.001, 1, length.out = 5))
      # Fitting the model, making sure that no output
      # is printed!
      invisible(capture.output(fit_app <- caret::train(as.factor(y) ~.,
                                                       data = df,
                                                       method = 'svmRadial',
                                                       tuneGrid = svm_grid,
                                                       trControl = trainControl(method = "cv"),
                                                       prob.model = TRUE
      )))
      
      # Apparent
      p_app <- predict(fit_app$finalModel, newdata = df[,-ncol(df)], type = "prob")[,2]
      # Then use the names of the coefficients to get the model matrix:
      df_app_X <- df[colnames(df)[!colnames(df) %in% "y"]]
      # Make model matrix by binding a column of 1 to the above
      app_matrix <- as.matrix(cbind("(Intercept)" = 1, df_app_X))
      
      # External 
      p_ext <- predict(fit_app$finalModel, newdata = df_val, type = "prob")[,2]
      # Then use the names of the coefficients to get the model matrix:
      df_ext_X <- df_val[colnames(df_val)[!colnames(df_val) %in% "y"]]
      # Make model matrix by binding a column of 1 to the above
      ext_matrix <- as.matrix(cbind("(Intercept)" = 1, df_ext_X))
      
      # Save dgm
      dgm_par_app <- dgm_par
        
    } else if (model == "ANN"){
      
      nnet_grid <- expand.grid(size = seq(from = 1, to = 10, by = 2),
                               decay = seq(from = 0.1, to = 2, by = 0.4))
      
      fit <- caret::train(as.factor(y) ~.,
                          data = df_train, 
                          method = 'nnet', 
                          tuneGrid = nnet_grid, 
                          trace = F, 
                          trControl = trainControl(method = "cv"),
                          preProcess = c("center", "scale"))
      
      if (fit_app$finalModel$convergence == 1){
        error_info <- paste(toString(error_info), paste0("nnet: Maximum number of iterations was reached"), sep = " + ")
      }
      
      coefs <- as.matrix(fit_app$coefnames)
      rownames(coefs) <- (fit_app$coefnames)
      
      # Apparent
      p_app <- predict(fit_app$finalModel, newdata = df, type='raw')
      # Then use the names of the coefficients to get the model matrix:
      df_app_X <- df[colnames(df) %in% rownames(coefs)]
      # Make model matrix by binding a column of 1 to the above
      app_matrix <- as.matrix(cbind("(Intercept)" = 1, df_app_X))
      
      # External 
      p_ext <- predict(fit_app$finalModel, newdata = df_val)
      # Then use the names of the coefficients to get the model matrix:
      df_ext_X <- df_val[colnames(df_val) %in% rownames(coefs)]
      # Make model matrix by binding a column of 1 to the above
      ext_matrix <- as.matrix(cbind("(Intercept)" = 1, df_ext_X))
      
      # Save dgm
      dgm_par_app <- dgm_par
      
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
    MAPE_rMSPE_app <- MAPE_rMSPE(p = p_app, dgm_par = dgm_par_app, iv_matrix = app_matrix) 
    tjur_app <- tjur(p = p_app, y = df$y)
    slope_app <- c(coef(glm(df$y ~ log(p_app/(1-p_app)), family="binomial"))[2])
    intercept_app <- coef(glm(df$y ~ offset(log(p_app/(1-p_app))), family="binomial"))
    calout <- loess(y ~ log(p_app/(1-p_app)), data = df, span = 7)
    eci_app <- (mean((p_app-fitted(calout))*(p_app-fitted(calout))))*(100)

    # obtain external results
    auc_ext <- fastAUC(p = p_ext, y = df_val$y)
    R2_ext <- pseudo_Rsqrs(p = p_ext, y = df_val$y)
    MAPE_rMSPE_ext <- MAPE_rMSPE(p = p_ext, dgm_par = dgm_par_app, iv_matrix = ext_matrix) 
    tjur_ext <- tjur(p = p_ext, y = df_val$y)
    slope_ext <- c(coef(glm(df_val$y ~ log(p_ext/(1-p_ext)), family="binomial"))[2])
    intercept_ext <- coef(glm(df_val$y ~ offset(log(p_ext/(1-p_ext))), family="binomial"))
    calout <- loess(y ~ log(p_ext/(1-p_ext)), data = df_val, span = 7)
    eci_ext <- (mean((p_ext-fitted(calout))*(p_ext-fitted(calout))))*(100)
    
  }) #Close warnings and error check
  
  
  
  # If there are warnings, paste those in the error_info
  if (!is.null(errors_warnings$warning)) {
    
    error_info <- paste(toString(error_info), toString(errors_warnings$warning), sep = " + ")
    
  }  
  # If there are error messages print this in the error_info
  if ("message" %in% errors_warnings$value) {
    
    error_info <- paste(toString(error_info), toString(errors_warnings$value$message), sep = " + ")
    
  }
  
  # Save results
  results <- list(
    c("Apparent",
      auc_app,
      intercept_app,
      slope_app,
      tjur_app,
      R2_app,
      eci_app,
      MAPE_rMSPE_app['MAPE'],
      MAPE_rMSPE_app['rMSPE'],
      error_info
    ),
    c("External",
      auc_ext,
      intercept_ext,
      slope_ext,
      tjur_ext,
      R2_ext,
      eci_ext,
      MAPE_rMSPE_ext['MAPE'],
      MAPE_rMSPE_ext['rMSPE'],
      error_info
    )
  )
  
  
} # End function


################################
## Function to obtain results ##
################################

# Given the study, obtain apparent estimands
# for each scenario

get_app_ext_results <- function(study, df, df_val, studyname) {
  
  # object to store results
  results_app <- as.data.frame(matrix(NA, nrow = nrow(study), 
                                      ncol = length(results_estimands_names), 
                                      dimnames = list(c(), results_estimands_names)))
  results_ext <- as.data.frame(matrix(NA, nrow = nrow(study), 
                                      ncol = length(results_estimands_names), 
                                      dimnames = list(c(), results_estimands_names)))
  
  # Add the information of the study to the respective columns
  results_app[, which(colnames(results_app) %in% study_info)] <- study[, which(colnames(study) %in% study_info)]
  results_ext[, which(colnames(results_ext) %in% study_info)] <- study[, which(colnames(study) %in% study_info)]
  
  # The n for external validation will be overwritten below
  # For each scenario within the study
  for (i in 1:nrow(study)) {
    print(i)
    results_app[i, 'scenario'] <- paste0("Scenario ", i)
    results_ext[i, 'scenario'] <- paste0("Scenario ", i)
    results_ext[i, 'n'] <- nrow(df_val[[i]])
    # Check for each scenario whether there are events sampled
    if (any(str_detect(names(df[[i]]),"Error: No events sampled") == TRUE)){
      
      # If no events were sampled, then the following will be the result:
      results_app[i, 'error_info'] <- c("Error: No events sampled")
      results_app[i, 'approach'] <- c("Apparent")
      results_ext[i, "error_info"] <- c("Error: No events sampled in development set")
      results_ext[i, "approach"] <- c("External")
      results_ext[i, 'observed events'] <- sum(df_val[[i]]$y)
      # and all the other columns of the estimands remain NA
      
    } else {
      
      results_app[i, 'observed events'] <- sum(df[[i]]$y)
      results_ext[i, 'observed events'] <- sum(df_val[[i]]$y)
      # determine which model & pre_selection is specified
      model <- study[i, ]$model
      pred_selection <- study[i, ]$pred_selection
      
      # Obtain dgm, depending on the dimensionality,
      # noise contribution, and therefore dgm settings
      if (study[i, ]$noise == "default"){
        
        dgm_par <-
          c(study[i, ]$par1, 
            rep(study[i, ]$par2 * 3, round(0.3 * study[i, ]$dim)), 
            rep(study[i, ]$par2, round(0.5 *  study[i, ]$dim)), 
            rep(study[i, ]$par2 * 0, round(0.2 * study[i, ]$dim))
            )
        
      } else if (study[i, ]$noise == "half"){
        
        dgm_par <-
          c(study[i, ]$par1, 
            rep(study[i, ]$par2 * 3, round(1/6 * study[i, ]$dim)), # strong
            rep(study[i, ]$par2, round(2/6 *  study[i, ]$dim)),    # weak
            rep(study[i, ]$par2 * 0, round(3/6 * study[i, ]$dim))
            ) # noise
        
      } else if (study[i, ]$noise == "none"){
        
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
      
      results_app[i, which(colnames(results_app) %in% apparent_col_names)] <- results[[1]]
      results_ext[i, which(colnames(results_ext) %in% apparent_col_names)] <- results[[2]]
      
      results_app_ext <- rbind(results_app, results_ext)
      
      # Fill in details of the study
      results_app_ext$study <- studyname
      
    } # end else statement
  } # end for loop
  
  return(results_app_ext)
  
} # end function


get_cv_estimands <- function(df, model, dgm_par, pred_selection, V, x10 = c(FALSE, TRUE)){
  
  
  # If nothing 'bad' happens, the error info remains as NA
  error_info <- NA
  original_V <- V
  
  errors_warnings <- ErrorsWarnings({
    
    # Fit model depending on scenario
    # And get predicted probabilities and modelmatrix
    
    ## If model is Firth without backwards elimination:
    if (model == "Firth" & pred_selection == "none") {
      
      fit_app <- logistf(y ~ ., data = df, flic = T)
      
      # Apparent
      app_matrix <- model.matrix(object = fit_app$formula, data = df)
      p_app <- 1 / (1 + exp(-app_matrix %*% fit_app$coefficients))
      
      # External
      ext_matrix <- model.matrix(object = fit_app$formula, data = df_val)
      p_ext <- 1 / (1 + exp(-ext_matrix %*% fit_app$coefficients))
      
      # Save coefficients
      coefs <- fit_app$coefficients
      
      # Save dgm 
      dgm_par_app <- dgm_par
      
      ## If model is Firth with backwards elimination:  
    } else if (model == "Firth" & pred_selection == "<0.05") {
      
      df_app <- df
      assign("df_app", as.data.frame(df_app), envir = .GlobalEnv)
      
      model_form <- as.formula(paste0("y ~", paste(colnames(df_app)[!colnames(df_app)%in%"y"], collapse = "+" )))
      fit_app <- logistf(formula = model_form, data = df_app, flic = T)
      fit_app <- backward(fit_app, trace = FALSE)
      
      # Check whether any predictors have been selected.
      if ((var(fit_app$linear.predictors) == 0) == TRUE) {
        error_info <- paste(toString(error_info), paste0("No predictors selected -> no calibration slope"), sep = " + ")
      }
      
      app_matrix <- model.matrix(object = fit_app$formula, data = df_app)
      p_app <- 1 / (1 + exp(-app_matrix %*% fit_app$coefficients))
      
      # External
      ext_matrix <- model.matrix(object = fit_app$formula, data = df_val)
      p_ext <- 1 / (1 + exp(-ext_matrix %*% fit_app$coefficients))
      
      # Save coefficients
      coefs <- fit_app$coefficients
      
      # Get the elements of the design generating mechanism that are
      # belonging to the model after backwards elimination
      # Always get the first element as this is the intercept
      ind <- na.omit(c(1, (as.numeric(str_extract_all(colnames(app_matrix), "(?<=V).*"))
                           # Add 1, because the indices of columns exclude the intercept
                           + 1)))
      
      dgm_par_app <- dgm_par[ind]
      
      ## If model is using ML :  
    } else if (model == "ML") {
      
      fit_app <- glm(y ~ ., family = "binomial", data = df) 
      p_app <- predict(fit_app, type = "response")
      app_matrix <- model.matrix(object = fit_app$formula, data = df)
      
      # External
      ext_matrix <- model.matrix(object = fit_app$formula, data = df_val)
      p_ext <- predict(fit_app, newdata = df_val, type = "response")
      
      # Save coefficients
      coefs <- fit_app$coefficients
      
      # Save dgm
      dgm_par_app <- dgm_par
      
    } else if (model == "Lasso") {
      
      fit_app <-  Pen_reg_VC(df = df, alpha = 1, nfolds = 10)
      
      ## Check whether predictors have been selected
      # Create linear predictor
      lp <- predict(fit_app, newdata = df, s = "lambda.min") 
      
      if (var(lp) == 0) {
        error_info <- paste(toString(error_info), paste0("No predictors selected -> no calibration slope"), sep = " + ")
      }
      
      # Create model matrix:
      # First retrieve the coefficients
      coefs <- as.matrix(coef(fit_app, s = "lambda.min"))
      # Remove those that do are 0
      coefs[coefs == 0] <- NA
      coefs <- na.exclude(coefs)
      
      # Apparent
      p_app <- predict(fit_app, newdata = df, s = "lambda.min", type = "response")
      # Then use the names of the coefficients to get the model matrix:
      df_app_X <- df[colnames(df) %in% rownames(coefs)]
      # Make model matrix by binding a column of 1 to the above
      app_matrix <- as.matrix(cbind("(Intercept)" = 1, df_app_X))
      
      
      # External
      p_ext <- predict(fit_app, newdata = df_val, s = "lambda.min", type = "response")
      # Then use the names of the coefficients to get the model matrix:
      df_ext_X <- df_val[colnames(df_val) %in% rownames(coefs)]
      # Make model matrix by binding a column of 1 to the above
      ext_matrix <- as.matrix(cbind("(Intercept)" = 1, df_ext_X))
      
      # Get the elements of the design generating mechanism that are
      # belonging to the model after backwards elimination
      # Always get the first element as this is the intercept
      ind <- na.omit(c(1, (as.numeric(str_extract_all(colnames(app_matrix), "(?<=V).*"))
                           # Add 1, because the indices of columns exclude the intercept
                           + 1)))
      
      # Save the correct part of dgm_par
      dgm_par_app <- dgm_par[ind]
      
      # If the model is Ridge    
    } else if (model == "Ridge") {
      
      fit_app <-  Pen_reg_VC(df = df, alpha = 0, nfolds = 10)
      
      # First retrieve the coefficients
      coefs <- as.matrix(coef(fit_app, s = "lambda.min"))
      
      # Apparent
      # Get the predictions 
      p_app <- predict(fit_app, newdata = df, s = "lambda.min", type = "response")
      # Then use the names of the coefficients to get the model matrix:
      df_app_X <- df[colnames(df) %in% rownames(coefs)]
      # Make model matrix by binding a column of 1 to the above
      app_matrix <- as.matrix(cbind("(Intercept)" = 1, df_app_X))
      
      
      # External
      p_ext <- predict(fit_app, newdata = df_val, s = "lambda.min", type = "response")
      # Then use the names of the coefficients to get the model matrix:
      df_ext_X <- df_val[colnames(df_val) %in% rownames(coefs)]
      # Make model matrix by binding a column of 1 to the above
      ext_matrix <- as.matrix(cbind("(Intercept)" = 1, df_ext_X))
      
      
      # Save dgm_par to dgm_par_app
      dgm_par_app <- dgm_par
      
    } else if (model == "RF") {
      
      rf_grid <- expand.grid(.mtry = seq(1, ncol(df)-1, by = 1))
      
      fit_app <- caret::train(as.factor(y) ~.,
                              data = df, 
                              method = 'rf',
                              tuneGrid = rf_grid,
                              trControl = trainControl(method = "cv")
      )
      
      
      coefs <- as.matrix(fit_app$coefnames)
      rownames(coefs) <- (fit_app$coefnames)
      
      # Apparent
      p_app <- predict(fit_app$finalModel, newdata = df, type = "prob")[,2]
      
      # Then use the names of the coefficients to get the model matrix:
      df_app_X <- df[colnames(df) %in% rownames(coefs)]
      # Make model matrix by binding a column of 1 to the above
      app_matrix <- as.matrix(cbind("(Intercept)" = 1, df_app_X))
      
      # External 
      p_ext <- predict(fit_app$finalModel, newdata = df_val, type = "prob")[,2]
      # Then use the names of the coefficients to get the model matrix:
      df_ext_X <- df_val[colnames(df_val) %in% rownames(coefs)]
      # Make model matrix by binding a column of 1 to the above
      ext_matrix <- as.matrix(cbind("(Intercept)" = 1, df_ext_X))
      
      
      # Save dgm
      dgm_par_app <- dgm_par
      
    } else if (model == "CART"){
      
      rpart_grid <- expand.grid(.cp = seq(0, 0.3, by = 0.01))
      
      fit_app <- caret::train(as.factor(y) ~.,
                              data = df,
                              method = 'rpart',
                              tuneGrid = rpart_grid,
                              trControl = trainControl(method = "cv"))
      
      # get names:
      coefs <- as.matrix(fit_app$coefnames)
      rownames(coefs) <- (fit_app$coefnames)
      
      # Apparent
      p_app <- predict(fit_app$finalModel, newdata = df, type = "prob")[,2]
      
      # Then use the names of the coefficients to get the model matrix:
      df_app_X <- df[colnames(df) %in% rownames(coefs)]
      # Make model matrix by binding a column of 1 to the above
      app_matrix <- as.matrix(cbind("(Intercept)" = 1, df_app_X))
      
      # External 
      p_ext <- predict(fit_app$finalModel, newdata = df_val, type = "prob")[,2]
      # Then use the names of the coefficients to get the model matrix:
      df_ext_X <- df_val[colnames(df_val) %in% rownames(coefs)]
      # Make model matrix by binding a column of 1 to the above
      ext_matrix <- as.matrix(cbind("(Intercept)" = 1, df_ext_X))
      
      # Save dgm
      dgm_par_app <- dgm_par
      
    } else if (model == "SVM" ){
      
      # Use caret::train, predefining a tune grid
      svm_grid <- expand.grid(.C = seq(0.001, 10, length.out = 10),
                              .sigma = seq(0.001, 1, length.out = 10))
      
      
      invisible(capture.output(fit_app <- caret::train(as.factor(y) ~.,
                                                       data = df_train,
                                                       method = 'svmRadial',
                                                       tuneGrid = svm_grid,
                                                       trControl = trainControl(method = "cv"),
                                                       prob.model = TRUE
      )))
      
      coefs <- as.matrix(fit_app$coefnames)
      rownames(coefs) <- (fit_app$coefnames)
      
      # Apparent
      p_app <- predict(fit_app$finalModel, newdata = df, type = "prob")[,2]
      
      # Then use the names of the coefficients to get the model matrix:
      df_app_X <- df[colnames(df) %in% rownames(coefs)]
      # Make model matrix by binding a column of 1 to the above
      app_matrix <- as.matrix(cbind("(Intercept)" = 1, df_app_X))
      
      # External 
      p_ext <- predict(fit_app$finalModel, newdata = df_val, type = "prob")[,2]
      # Then use the names of the coefficients to get the model matrix:
      df_ext_X <- df_val[colnames(df_val) %in% rownames(coefs)]
      # Make model matrix by binding a column of 1 to the above
      ext_matrix <- as.matrix(cbind("(Intercept)" = 1, df_ext_X))
      
      # Save dgm
      dgm_par_app <- dgm_par
      
      
    }# Close model if else statements
      
      results <- list(p, iv_matrix, dgm_par_folds)
      
    } # Close function .doFit
    
    results <- (lapply(seq(V), .doFit, folds = folds, model = model, pred_selection = pred_selection)) # obtain model results for all folds
    p <- c(unlist(sapply(results, "[[", 1))) # Get out pred values as a single vector
    p[which(rownames(df) %in% unlist(folds))] <- p #Re-order pred values, so they are in line with y
    
    
    p_per_fold <- lapply(results, "[[", 1) # getting the predictions per fold
    iv_matrix <- lapply(results, "[[", 2) # iv_matrix per fold
    dgm_par_folds <- lapply(results, "[[", 3) # Dgm_par per fold
    
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
      # create object to shorten future code:
      data <- df[unlist(which(rownames(df) %in% folds[[v]])),]
      ppf <- p_per_fold[[v]]
      dgm_par_ppf <- dgm_par_folds[[v]]
      ppf_true <- 1 / (1 + exp(-iv_matrix[[v]] %*% dgm_par_ppf))
      
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
      calout <- loess(y ~ log(ppf/(1-ppf)), data = data)
      eci_folds[v] <- (mean((ppf-fitted(calout))*(ppf-fitted(calout))))*(100)
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
    
  }) # Close ErrorsWarnings
  
  # If there are warnings, paste those in the error_info
  if (!is.null(errors_warnings$warning)) {
    
    error_info <- paste(toString(error_info), toString(errors_warnings$warning), sep = " + ")
    
  }  
  # If there are error messages print this in the error_info
  if ("message" %in% errors_warnings$value) {
    
    error_info <- paste(toString(error_info), toString(test$value$message), sep = " + ")
    
  }
  
  
  # If the it is anything other than 10x10 cv:
  if (x10 == FALSE){
    
    results <- c(paste0(original_V, " fold cross-validation"), auc_results, intercept, slope, tjur_results, R2, eci, MAPE_results, rMSPE_results, error_info)
    
  } else {
    # 
    # # If it is 10x10 cv, only keep the raw results of each fold and name them:  
    # names(auc_folds) <- c(rep("auc", length(auc_folds)))
    # names(intercept_folds) <- c(rep("intercept", length(intercept_folds)))
    # names(slope_folds) <- c(rep("slope", length(slope_folds)))
    # names(tjur_folds) <- c(rep("Tjur", length(tjur_folds)))
    # names(R2_folds) <- c(rep("R2", length(R2_folds)))
    # names(eci_folds) <- c(rep("eci", length(eci_folds)))
    # names(mape_folds) <- c(rep("MAPE", length(mape_folds)))
    # names(rmspe_folds) <- c(rep("rMSPE", length(rmspe_folds)))
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
  
  return(results)
  
} # Close function

## Function to apply get_cv_estimands to all datasets in a study:
get_cv_results <- function(study, df, V, studyname) {
  
  # results matrix:
  results_cv <- as.data.frame(matrix(NA, nrow = nrow(study), ncol = length(results_estimands_names), dimnames = list(c(), results_estimands_names)))
  
  # Fill in details of the study
  results_cv$study <- studyname
  # And of each scenario:
  results_cv[, which(colnames(results_cv) %in% study_info)] <- study[, which(colnames(study) %in% study_info)]
  
  
  for (i in 1:nrow(study)) {
    print(i) ## REMOVE WHEN SIMULATION COMMENCES
    # Paste the correct scenario number:
    results_cv[i, 'scenario'] <- paste0("Scenario ", i)
    
    # Check for each scenario whether there are events sampled
    if (any(str_detect(names(df[[i]]),"Error: No events sampled") == TRUE)){
      
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


