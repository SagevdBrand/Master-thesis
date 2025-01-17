
## Libraries, file paths and functions
source("./src/setup.R")
source("./src/estimand functions.R")
source("./src/data generation functions.R")

library(microbenchmark)

## Load scenario settings
s1 <- readRDS(study_1_settings)
s2 <- readRDS(study_2_settings)
s3 <- readRDS(study_3_settings)

## Load validation data
s1_val_data <- readRDS(study_1_val_data)
s2_val_data <- readRDS(study_2_val_data)
s3_val_data <- readRDS(study_3_val_data)


## Create and load simulation data
set.seed(123)
s1_data <- generate_data(s1, validation = FALSE)
s2_data <- generate_data(s2, validation = FALSE)
s3_data <- generate_data(s3, validation = FALSE)

## create one test set for all studies, shortes scenarios:
df <- c(s1_data[c(7,16)], s2_data[c(10,19)], s3_data[c(1,4,7,10,13,16)])
df_val <- c(s1_val_data[c(7,16)], s2_val_data[c(10,19)], s3_val_data[c(1,4,7,10,13,16)])
study <- rbind(s1[c(7,16),], s2[c(10,19),], s3[c(1,4,7,10,13,16),])


##################################
###### Apparent performance ######
##################################
################################
## Function to obtain results ##
################################

# Given the study, obtain apparent estimands
# for each scenario

get_app_ext_results_test <- function(study, df, df_val, studyname) {
  
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
  
  # Fill in some information for each scenario:
  results_app[, 'scenario'] <- study$scenario
  results_ext[, 'scenario'] <- study$scenario
  # Overwrite n of the development data
  results_ext[, 'n'] <- sapply(df_val, nrow)
  # How many events observed in validation set?
  results_ext[, 'observed events'] <- sapply(df_val, function(x) sum(x$y))
  
  # Check for each scenario whether there are events sampled
  if (any(str_detect(names(df),"Error: No events sampled"))){
    
    ind_no_events  <- which(str_detect(names(df),"Error: No events sampled"))
    # If no events were sampled, then the following will be the result:
    results_app[c(ind_no_events), 'error_info'] <- c("Error: No events sampled")
    results_app[c(ind_no_events), 'approach'] <- c("Apparent")
    # For the external results, there will be the following messages
    results_ext[c(ind_no_events), "error_info"] <- c("Error: No events sampled in development set")
    results_ext[c(ind_no_events), "approach"] <- c("External")
    # and all the other columns of the estimands remain NA
    
  } else {
    
    s <- nrow(study)
    # Fill the columns with apparent results, no SE!
    results <- sapply(seq(s), get_app_ext_estimands_test, df = df, df_val = df_val, study)
    
    # Save the results of each list into the correct matrix
    results_app[i, which(colnames(results_app) %in% apparent_col_names)] <- results[[1]]
    results_ext[i, which(colnames(results_ext) %in% apparent_col_names)] <- results[[2]]
    # Bind the matrices together
    results_app_ext <- rbind(results_app, results_ext)
    # Fill in details of the study
    results_app_ext$study <- studyname
    
    # Save the p_apps to the environment
    p_app_preds[[i]] <- c(results[[3]])
    
  } # end else statement
  
  assign(paste0("p_app_", studyname), p_app_preds, envir = .GlobalEnv)
  return(results_app_ext)
} # end function





get_app_ext_estimands_test <- function(s, df, df_val, study){
  
  # Get correct information of scenario to be used
  df_app <- df[[s]]
  df_val_app <- df_val[[s]]
  
  model <- study$model[s]
  pred_selection <- study$pred_selection[s]
  
  if (study[s, ]$noise == "default"){
    dgm_par <-
      c(study[s, ]$par1, 
        rep(study[s, ]$par2 * 3, round(0.3 * study[s, ]$dim)), 
        rep(study[s, ]$par2, round(0.5 *  study[s, ]$dim)), 
        rep(study[s, ]$par2 * 0, round(0.2 * study[s, ]$dim)))
  }
  
  if (study[s, ]$noise == "half"){
    dgm_par <-
      c(study[s, ]$par1, 
        rep(study[s, ]$par2 * 3, round(1/6 * study[s, ]$dim)), # strong
        rep(study[s, ]$par2, round(2/6 *  study[s, ]$dim)),    # weak
        rep(study[s, ]$par2 * 0, round(3/6 * study[s, ]$dim))) # noise
  } 
  
  if (study[s, ]$noise == "none"){
    dgm_par <-
      c(study[s, ]$par1, 
        rep(study[s, ]$par2 * 3, round(1/6 * study[s, ]$dim)), # strong
        rep(study[s, ]$par2, round(5/6 *  study[s, ]$dim))     # weak
      )
  }
  
  # Create an empty element to store all errors that occur
  # In between and add error_warnings in the end
  error_info <- NA
  
  # Start the actual function
  errors_warnings <- ErrorsWarnings({
    
    # obtaining model matrices
    app_matrix <- as.matrix(df_app[,-which(colnames(df_app)=="y")])
    app_matrix <- as.matrix(cbind("(Intercept)" = 1, app_matrix))
    
    ext_matrix <- as.matrix(df_val_app[, -which(colnames(df_val_app)=="y")])
    ext_matrix <- as.matrix(cbind("(Intercept)" = 1, ext_matrix))
    # Fit model depending on scenario
    # And get predicted probabilities
    
    if (model == "Firth" | model == "ML"){ 
      # Specifiy model formula (necessary for logistf)
      model_form <- as.formula(paste0("y ~", paste(colnames(df_app)[!colnames(df_app)%in%"y"], collapse = "+" )))
      
      ## If model is Firth without backwards elimination:
      if (model == "Firth" & pred_selection == "none") {
        
        fit_app <- logistf(formula = model_form, data = df_app, flic = T, firth = F)
        
        ## If model is Firth with backwards elimination:  
      } else if (model == "Firth" & pred_selection == "<0.05") {
        df_app <- df_app
        assign("df_app", as.data.frame(df_app), envir = .GlobalEnv)
        
        fit_app <- logistf(formula = model_form, data = df_app, flic = T, firth = F)
        fit_app <- backward(fit_app, trace = FALSE)
        
        # Check whether any predictors have been selected.
        if (var(fit_app$linear.predictors) == 0) {
          error_info <- paste(toString(error_info), paste0("No predictors selected -> no calibration slope"), sep = " + ")
        }
        
        ## If model is using ML :  
      } else if (model == "ML") {
        
        fit_app <- glm(y ~ ., family = "binomial", data = df_app) 
        
      }
      # Predictions:
      # Apparent
      fit_app_matrix <- model.matrix(object = fit_app$formula, data = df_app)
      p_app <- 1 / (1 + exp(-fit_app_matrix %*% fit_app$coefficients))
      
      # External
      fit_ext_matrix <- model.matrix(object = fit_app$formula, data = df_val_app)
      p_ext <- 1 / (1 + exp(-fit_ext_matrix %*% fit_app$coefficients))
      
      # If the model uses penalized regression  
    } else if (model == "Lasso" | model == "Ridge") {
      
      alpha <- ifelse(model == "Lasso", 1, 0) # alpha = 1 for Lasso
      
      # Make sure that there are at least 8 events or no-events:
      if (sum(df_app$y) < 8 | sum(1-df_app$y) <8 ){
        # Fit the model depending on alpha with LOOCV
        fit_app <-  Pen_reg_VC(df_app = df_app, alpha = alpha, nfolds = nrow(df_app))
        error_info <- paste(toString(error_info), 
                            paste0("Too few (non-)events for tuning -> LOOCV"),
                            sep = " + ")
      } else {
        # Otherwise just use 10-f-cv for tuning of lambda
        fit_app <-  Pen_reg_VC(df_app = df_app, alpha = alpha, nfolds = 10)
      }  # close error handling too few events
      
      ## Check whether predictors have been selected
      # Create linear predictor
      lp <- predict(fit_app, newdata = df_app, s = "lambda.min") 
      
      if (var(lp) == 0) {
        error_info <- paste(toString(error_info), paste0("No predictors selected -> no calibration slope"), sep = " + ")
      }
      # Apparent
      p_app <- predict(fit_app, newdata = df_app, s = "lambda.min", type = "response")
      # External
      p_ext <- predict(fit_app, newdata = df_val_app, s = "lambda.min", type = "response")
      
    } else  if (model == "RF"| model == "CART") {
      
      if (model == "RF"){
        
        # pre-defining a tune grid
        rf_grid <-
          expand.grid(.mtry = seq(1, ncol(df_app) - 1, by = 2))
        
        fit_app <- caret::train(as.factor(y) ~.,
                                data = df_app, 
                                method = 'rf',
                                tuneGrid = rf_grid,
                                trControl = trainControl(method = "cv")
        )
        # Otherwise fit CART
      } else {
        # Pre-defining a tuning grid
        rpart_grid <- expand.grid(.cp = seq(0, 0.3, by = 0.01))
        # Fitting the model
        fit_app <- caret::train(as.factor(y) ~.,
                                data = df_app,
                                method = 'rpart',
                                tuneGrid = rpart_grid,
                                trControl = trainControl(method = "cv"))
        
        # Check whether there were more than 1 splits
        if (nrow(fit_app$finalModel$cptable) == 1) {
          error_info <- paste(toString(error_info), paste0("No predictors selected -> no calibration slope"), sep = " + ")
        } # close handling if there is only one split in the tree
      } # close which tree
      # Apparent
      p_app <- predict(fit_app$finalModel, newdata = df_app, type = "prob")[,2]
      # External 
      p_ext <- predict(fit_app$finalModel, newdata = df_val_app, type = "prob")[,2]
      
      # } else if (model == "SVM"){
      #   
      #   # Pre-defining a tuning grid
      #   svm_grid <- expand.grid(.C = seq(0.001, 10, length.out = 4),
      #                           .sigma = seq(0.001, 1, length.out = 4))
      #   # Fitting the model, making sure that no output
      #   # is printed!
      #   invisible(capture.output(fit_app <- caret::train(as.factor(y) ~.,
      #                                                    data = df_app,
      #                                                    method = 'svmRadial',
      #                                                    tuneGrid = svm_grid,
      #                                                    trControl = trainControl(method = "cv"),
      #                                                    prob.model = TRUE
      #   )))
      #   
      #   # Apparent
      #   p_app <- predict(fit_app$finalModel, newdata =  df_app[, which(!colnames(df_app) %in% "y")], type = "prob")[,2]
      #   # External 
      #   p_ext <- predict(fit_app$finalModel, newdata = df_val_app[, which(!colnames(df_val_app) %in% "y")], type = "prob")[,2]
      #   
      # } else if (model == "ANN"){
      #   # Pre-defining a tuning grid
      #   nnet_grid <- expand.grid(size = seq(from = 1, to = 10, by = 2),
      #                            decay = seq(from = 0.1, to = 2, by = 0.4))
      #   # Fitting the model
      #   fit_app <- caret::train(as.factor(y) ~.,
      #                           data = df_app, 
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
      #   p_app <- predict(fit_app$finalModel, newdata = df_app, type='raw')
      #   # External 
      #   p_ext <- predict(fit_app$finalModel, newdata = df_val_app, type = 'raw')
      
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
    auc_app <- fastAUC(p = p_app, y = df_app$y)
    R2_app <- pseudo_Rsqrs(p = p_app, y = df_app$y)
    MAPE_rMSPE_app <- MAPE_rMSPE(p = p_app, dgm_par = dgm_par, iv_matrix = app_matrix) 
    tjur_app <- tjur(p = p_app, y = df_app$y)
    slope_app <- c(coef(glm(df_app$y ~ log(p_app/(1-p_app)), family="binomial"))[2])
    intercept_app <- coef(glm(df_app$y ~ offset(log(p_app/(1-p_app))), family="binomial"))
    calout <- loess(y ~ log(p_app/(1-p_app)), data = df_app, span = 10)
    eci_app <- (mean((p_app-fitted(calout))*(p_app-fitted(calout))))*(100)
    
    # obtain external results
    auc_ext <- fastAUC(p = p_ext, y = df_val_app$y)
    R2_ext <- pseudo_Rsqrs(p = p_ext, y = df_val_app$y)
    MAPE_rMSPE_ext <- MAPE_rMSPE(p = p_ext, dgm_par = dgm_par, iv_matrix = ext_matrix) 
    tjur_ext <- tjur(p = p_ext, y = df_val_app$y)
    slope_ext <- c(coef(glm(df_val_app$y ~ log(p_ext/(1-p_ext)), family="binomial"))[2])
    intercept_ext <- coef(glm(df_val_app$y ~ offset(log(p_ext/(1-p_ext))), family="binomial"))
    calout <- loess(y ~ log(p_ext/(1-p_ext)), data = df_val_app, span = 10)
    eci_ext <- (mean((p_ext-fitted(calout))*(p_ext-fitted(calout))))*(100)
    
    # Save all the results in a matrix  
    results <- matrix(nrow = 2, ncol = 11, dimnames = c("Apparent", "External"), c())
    
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
        error_info
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
        error_info
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





