  assign("df_train", df_train, envir = .GlobalEnv)
  assign("df_test", df_test, envir = .GlobalEnv)
  
  # Train/test glm for each fold
  # Fit model depending on scenario
  # And get predicted probabilities and modelmatrix
  if (model == "Firth"| model == "ML" ) {

    # Specifiy model formula (necessary for logistf)
    model_form <- as.formula(paste0("y ~", paste(colnames(df)[!colnames(df)%in%"y"], collapse = "+" )))
    
    if (model == "Firth" & pred_selection == "none"){
      # Fit Firth model without Backward elimination
      fit <- logistf(formula = model_form, data = df_train, flic = T)
      
    } else if (model == "Firth" & pred_selection == "<0.05"){
      # Fit Firth model with Backward elimination
      fit <- logistf(formula = model_form, data = df_train, flic = T)
      fit <- backward(fit, trace = FALSE)
      
      # Check whether any predictors have been selected at all. 
      if (var(fit$linear.predictors) == 0) {
        error_info <- paste(toString(error_info), paste0("No predictors selected -> no calibration slope"), sep = " + ")
      }
      
    } else if (model == "ML"){
      #Fit logistic regression with Maximum Likelihood
      fit <- glm(y ~ ., family = "binomial", data = df_train)
      
      # Check for separation
      if(any(sqrt(diag(summary(fit)$cov.unscaled)*summary(fit)$dispersion) > 70)){
        error_info <- paste(toString(error_info), paste0("Data separation might have occured"), sep = " + ")
      }
      }
    # original sample
    df_matrix <- model.matrix(object = fit$formula, data = df)
    p_df <- 1 / (1 + exp(-df_matrix %*% fit$coefficients))
    
    # bootstrapped sample
    train_matrix <- model.matrix(object = fit$formula, data = df_train)
    p_train <- 1 / (1 + exp(-train_matrix %*% fit$coefficients))
    
    # out/test sample
    test_matrix <- model.matrix(object = fit$formula, data = df_test)
    p_test <- 1 / (1 + exp(-test_matrix %*% fit$coefficients))
    
    # Get the elements of the design generating mechanism that are
    # belonging to the model after backwards elimination
    # Always get the first element as this is the intercept
    ind <- na.omit(c(1, (as.numeric(str_extract_all(colnames(train_matrix), "(?<=V).*"))
                         # Add 1, because the indices of columns exclude the intercept
                         + 1)))
    
    dgm_par_folds <- dgm_par[ind]
    
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
    
    # Create model matrix:
    # First retrieve the coefficients
    coefs <- as.matrix(coef(fit, s = "lambda.min"))
    # Remove those that are 0
    coefs[coefs == 0] <- NA
    coefs <- na.exclude(coefs)
    
    # original sample
    # Then use the names of the coefficients to get the model matrix:
    df_X <- df[colnames(df) %in% rownames(coefs)]
    # Make model matrix by binding a column of 1 to the above
    df_matrix <- as.matrix(cbind("(Intercept)" = 1, df_X))
    # Obtain predictions
    p_df <- predict(fit, newdata = df, s = "lambda.min", type = "response")
    
    # bootstrapped sample
    train_X <- df_train[colnames(df_train) %in% rownames(coefs)]
    train_matrix <- as.matrix(cbind("(Intercept)" = 1, train_X))
    p_train <- predict(fit, newdata = df_train, s = "lambda.min", type = "response")
    
    # out/test sample
    test_X <- df_test[colnames(df_test) %in% rownames(coefs)]
    test_matrix <- as.matrix(cbind("(Intercept)" = 1, test_X))
    p_test <- predict(fit, newdata = df_test, s = "lambda.min", type = "response")
    
    # Get the elements of the design generating mechanism that are
    # belonging to the model after backwards elimination
    # Always get the first element as this is the intercept
    ind <- na.omit(c(1, (as.numeric(str_extract_all(colnames(train_matrix), "(?<=V).*"))
                         # Add 1, because the indices of columns exclude the intercept
                         + 1)))
    
    dgm_par_boot <- dgm_par[ind]
    
  } else if (model == "RF" | model == "SVM" | model == "ANN" | model == "CART") {
    
    if(model == "RF"| model == "CART"){
      
      if (model == "RF"){
        # pre-defining a tune grid
        rf_grid <- expand.grid(.mtry = seq(1, ncol(df_train)-1, by = 1))
        # Fitting the appropriate model:
        fit <- caret::train(as.factor(y) ~.,
                            data = df_train, 
                            method = 'rf',
                            tuneGrid = rf_grid,
                            trControl = trainControl(method = "cv")
        )
        
      
      } else {
        rpart_grid <- expand.grid(.cp = seq(0, 0.3, by = 0.01))
        fit <- caret::train(as.factor(y) ~.,
                            data = df_train,
                            method = 'rpart',
                            tuneGrid = rpart_grid,
                            parms = list(split = "information"),
                            trControl = trainControl(method = "cv"))
        
        # Check whether there were more than 1 splits
        if (nrow(fit$finalModel$cptable) == 1) {
          # if not, find the complexity parameter with the next highest accuracy:
          ind <-which(fit$results$Accuracy != max(fit$results$Accuracy))
          best_tune <- fit$results[ind,] %>% slice(which.max(Accuracy)) %>% .$cp
          
          # Refit the model
          fit <- caret::train(as.factor(y) ~.,
                              data = df_train,
                              method = 'rpart',
                              parms = list(split = "information"),
                              control = rpart.control(cp = best_tune),
                              trControl = trainControl(method = "none"))
        } # close check for only 1 splits
      } # Close if else for which specific tree model
      
      # Get probabilities out:
      p_df <- predict(fit$finalModel, newdata = df, type = "prob")[,2]
      p_train <- predict(fit$finalModel, newdata = df_train, type = "prob")[,2] 
      p_test <- predict(fit$finalModel, newdata = df_test, type = "prob")[,2]
      
      } else if (model == "SVM"){
      # pre-defining a tune grid
      svm_grid <- expand.grid(.C = seq(0.001, 10, length.out = 5),
                              .sigma = seq(0.001, 1, length.out = 5))
      
      # Fit svm with Radial Basis Function kernel
      invisible(capture.output(fit <- caret::train(as.factor(y) ~.,
                                                   data = df_train,
                                                   method = 'svmRadial',
                                                   tuneGrid = svm_grid,
                                                   trControl = trainControl(method = "cv"),
                                                   prob.model = TRUE
      )))
      p_df <- predict(fit$finalModel, newdata = df[!colnames(df)%in%"y"], type = "prob")[,2]
      p_train <- predict(fit$finalModel, newdata = df_train[!colnames(df_train)%in%"y"], type = "prob")[,2]
      p_test <- predict(fit$finalModel, newdata = df_test[!colnames(df_test)%in%"y"], type = "prob")[,2]
    
      } else { # Otherwise it is a neural network
      nnet_grid <- expand.grid(size = seq(from = 1, to = 10, by = 2),
                               decay = seq(from = 0.1, to = 2, by = 0.4))
      
      # Fit Neural network
      fit <- caret::train(as.factor(y) ~.,
                          data = df_train, 
                          method = 'nnet', 
                          tuneGrid = nnet_grid, 
                          trace = F, 
                          trControl = trainControl(method = "cv"),
                          preProcess = c("center", "scale"))
      
      # Check convergence of model
      if (fit$finalModel$convergence == 1){
        error_info <- paste(toString(error_info), paste0("nnet: Maximum number of iterations was reached"), sep = " + ")
      }
      
      # Obtain predictions
      p_df <- predict(fit$finalModel, newdata = df)
      p_train <- predict(fit$finalModel, newdata = df_train)
      p_test <- predict(fit$finalModel, newdata = df_test)
      
      }
    # Get coef names out
    coefs <- as.matrix(fit$coefnames)
    # Original data
    df_X <- df[colnames(df) %in% coefs]
    df_matrix <- as.matrix(cbind("(Intercept)" = 1, df_X))
    # Boostrapped data
    train_X <- df_train[colnames(df_train) %in% coefs]
    train_matrix <- as.matrix(cbind("(Intercept)" = 1, train_X))
    # Test/out sample
    test_X <- df_test[colnames(df_test) %in% coefs]
    test_matrix <- as.matrix(cbind("(Intercept)" = 1, test_X))
    
    # Get the elements of the design generating mechanism that are
    # belonging to the model after backwards elimination
    # Always get the first element as this is the intercept
    ind <- na.omit(c(1, (as.numeric(str_extract_all(colnames(train_matrix), "(?<=V).*"))
                         # Add 1, because the indices of columns exclude the intercept
                         + 1)))
    
    dgm_par_boot <- dgm_par[ind]
 
}