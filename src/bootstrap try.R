#########################
## BOOTSTRAP FUNCTIONS ##
#########################

########
## [ ] WHAT TO DO IF VALUES BECOME NEGATIVE? LIKE R^2? BECAUSE THE OPTIMISM IS TOO GREAT FOR EXAMPLE?



source("./src/setup.R")
source("./src/estimand functions.R")
source("./src/data generation functions.R")

study <- readRDS(study_1_settings)
s1_data <- generate_data(study, validation = FALSE)

# Getting one dataframe
i <- 1 # or whatever we're interested in1
df <- data[[i]]
model <- study[i, ]$model
pred_selection <- study[i, ]$pred_selection

dgm_par <- c(study[i, ]$par1, 
               rep(study[i, ]$par2 * 3, round(0.3 * study[i, ]$dim)),  # strong
               rep(study[i, ]$par2,     round(0.5 * study[i, ]$dim)),  # medium
               rep(study[i, ]$par2 * 0, round(0.2 * study[i, ]$dim)))  # noise

###########################
###########################
## obtain apparent results:

p_app_study_1 <- readRDS(paste0(estimands_path, "app_preds_for_bootstrap_study1.RDS"))
results_app_ext_s1 <- readRDS(paste0(estimands_path, "app_estimands_for_bootstrap_study1.RDS"))
estimands_app <- results_app %>% filter(approach == "Apparent") %>% .[,(colnames(results_app) %in% estimands_names)]

############
## Lambda ##
############

### Bootstrap estimate of optimism: Harrell's ###
get_lambda <- function(data){
  
  train_results <- data %>% dplyr::select(ends_with("_train"))
  df_results <- data %>% dplyr::select(ends_with("_df"))
  
  lambda <- colMeans(train_results - df_results, na.rm = TRUE)
  names(lambda) <- sub("*_train", "_optimism", names(lambda))
  return(lambda)
}



###############
## for .632+ ##
###############
# Where pred_app is paste0("p_app_", studyname)[[i]]

get_gamma <- function(data, pred_app){
  p1 <- sum(data$y)/nrow(data)
  q1_classes <- ifelse(pred_app > 0.5, 1, 0)
  q1 <- sum(q1_classes)/nrow(data)
  gamma_hat <- p1*(1 - q1) + (1 - p1) * q1
  return(gamma_hat)
}


get_bootstrap_estimands <- function(df, model, dgm_par, pred_selection, pred_app, theta_app, nboot){  
    
  error_info <- NA
    errors_warnings <- ErrorsWarnings({  
    
      # create empty matrix for results
      results <- as.data.frame(matrix(NA, nrow = nboot, ncol = 24, dimnames = list(c(), c(
        paste0(estimands_names, c(rep("_df",8), rep("_train",8), rep("_test",8)))
      ))))
      
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
        
        dgm_par_boot <- dgm_par[ind]
        
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
        
        if (model == "RF"| model == "CART") {
          
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
          } # close check for convergence
          
          # Obtain predictions
          p_df <- predict(fit$finalModel, newdata = df)
          p_train <- predict(fit$finalModel, newdata = df_train)
          p_test <- predict(fit$finalModel, newdata = df_test)
        
          } # close  machine learning models if elses
        
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
      p_df_true <- 1 / (1 + exp(-df_matrix %*% dgm_par_boot))
      p_train_true <- 1 / (1 + exp(-train_matrix %*% dgm_par_boot))
      p_test_true <- 1 / (1 + exp(-test_matrix %*% dgm_par_boot))
      
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
      
      # ECI
      calout <- loess(y ~ log(p_df/(1-p_df)), data = df, span = 10)
      results$eci_df[b] <- (mean((p_df-fitted(calout))*(p_df-fitted(calout))))*(100)
      calout <- loess(y ~ log(p_train/(1-p_train)), data = df_train, span = 10)
      results$eci_train[b] <- (mean((p_train-fitted(calout))*(p_train-fitted(calout))))*(100)
      calout <- loess(y ~ log(p_test/(1-p_test)), data = df_test, span = 10)
      results$eci_test[b] <- (mean((p_test-fitted(calout))*(p_test-fitted(calout))))*(100)
      
      # MAPE
      results$mape_df[b] <- mean(abs(p_df_true-p_df))
      results$mape_train[b] <- mean(abs(p_train_true-p_train))
      results$mape_test[b] <- mean(abs(p_test_true-p_test))
      
      #rMSPE
      results$rmspe_df[b]  <- sqrt(mean((p_df_true-p_df)^2))
      results$rmspe_train[b]  <- sqrt(mean((p_train_true-p_train)^2))
      results$rmspe_test[b]  <- sqrt(mean((p_test_true-p_test)^2))
      
    } # Close bootstrap for loop
      # Matrix to store results for a single scenario:
      results_matrix <- as.data.frame(matrix(NA, ncol = length(iv_colnames), nrow = 3, dimnames = list(c(),
                                                                                                       c(iv_colnames))))
      results_matrix$approach <- c("Harrell's bootstrap", ".632 bootstrap", ".632+ bootstrap")
      
      theta_train <- results %>% dplyr::select(ends_with("_train"))
      train_se <- apply(theta_train, 2, sd)/(sqrt(nboot))
      #####################
      # Harrell's results #
      #####################
      lambda <- get_lambda(results)
      Harrell_results <- as.matrix(t(theta_app) - lambda)
      colnames(Harrell_results) <- estimands_names
      
      results_matrix[1,which(colnames(results_matrix) %in% estimands_names)] <- Harrell_results
      
      ## SE 
      delta_Harrell <- abs(Harrell_results-theta_app)
      SE_Harrell <-  train_se + delta_Harrell
      colnames(SE_Harrell) <- estimands_se_names
      
      results_matrix[1,which(colnames(results_matrix) %in% estimands_se_names)] <- SE_Harrell
      
      ##################
      # .632 bootstrap #
      ##################
      ## Average test estimate ##
      theta_test <- results %>% dplyr::select(ends_with("_test")) %>% colMeans(., na.rm = TRUE)
      efron_.632_results <- as.matrix(t((0.368 * theta_app) + (0.632 * theta_test)))
      colnames(efron_.632_results) <- estimands_names
      
      results_matrix[2,which(colnames(results_matrix) %in% estimands_names)] <- efron_.632_results
      
      ## SE 
      delta_.632 <- abs(efron_.632_results-theta_app)
      SE_.632 <-  train_se + delta_.632
      colnames(SE_.632) <- estimands_se_names
      
      results_matrix[2,which(colnames(results_matrix) %in% estimands_se_names)] <- SE_.632
      
      ###################
      # .632+ bootstrap #
      ###################
      # Obtain gamma:
      gamma <- get_gamma(data = df, pred_app = pred_app)
      
      # Relative overfitting 
      # For almost all the following should give correct estimates
      R <- abs(theta_test - theta_app) / abs(gamma - theta_app)
      
      # However, for the eci, mape and rmspe, the results might be quite different from eachother
      #(R[6:8] <- abs(log(theta_app[6:8]) - log(theta_test[6:8])) / abs(gamma - log(theta_app[6:8])))
      # If R is negative, return absolute value
      R <-  ifelse( R< 0, abs(R), R)
      # Then return all that are bigger than 1 to 1. 
      R <-  ifelse( R> 1, 1, R)
      # weight
      W <-  abs(.632 / (1 - .368 * R))
      
      # results:
      efron_.632_plus_results <- as.matrix(t((1-W) * theta_app) + (W * theta_test)) 
      colnames(efron_.632_plus_results) <- estimands_names
      
      results_matrix[3,which(colnames(results_matrix) %in% estimands_names)] <- efron_.632_plus_results
      
      ## SE 
      delta_.632_plus <- abs(efron_.632_plus_results-theta_app)
      SE_.632_plus <-  train_se + delta_.632_plus
      colnames(SE_.632_plus) <- estimands_se_names
      
      results_matrix[3,which(colnames(results_matrix) %in% estimands_se_names)] <- SE_.632_plus
      
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
}  


### And to apply this for each scenario in a study and obtaining the results:  
get_bootstrap_results <- function(study, df, V, studyname) {

  # Matrix for results (all 3 methods so should be three
  # times as large as a single matrix)
  results_boot <- as.data.frame(matrix(NA, nrow =3*nrow(study), 
                                       ncol = length(results_estimands_names), 
                                       dimnames = list(c(), results_estimands_names)))
  
  # Fill in details of the study
  results_boot$study <- studyname
  
  # The results should be pasted in rows of 3 as well
  # Therefore for each scenario, another index is necessary
  # so it corresponds with each 3rd row in the results
  results_i_start <- seq(1, nrow(results_boot), by = 3)
  results_i_end <- seq(3, nrow(results_boot), by = 3)
  
  for (i in 1:nrow(study)){
    print(i)
    
    # Fill in some details:
    results_boot[c(results_i_start[i]:results_i_end[i]),
                 'scenario'] <- paste0("Scenario ", i)
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
    pred_app <- get(paste0("p_app_", studyname))[[i]]
    results_app <- results_app_ext_s1 %>% filter(approach == "Apparent" & scenario == paste("Scenario", i)) 
    theta_app <- results_app %>% .[,(colnames(results_app) %in% estimands_names)]
    theta_app <- as.numeric(theta_app)
    
    # Go along with obtaining the results
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
                                                                                            nboot = 500)
    
  } # close for loop
  
  return(results_boot)
} # Close if else when checking for events in the data
  
} # close get_bootstrap_results_function