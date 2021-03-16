#######Error handling #######################################################################
## Get libraries, paths and a dataset:
## Libraries, file paths and functions
source("./src/setup.R")
source("./src/estimand functions.R")
source("./src/data generation functions.R")

## Load scenario settings
s1 <- read_rds(study_1_settings)
study <- s1

### Load data as df ###
# Have it run at least once, so that there are files in the folder.
source("./src/validation data generation study 1.R") 

# get the data names
data_files <- list.files(path = study_1_data, recursive = T, full.names = F)

# Generate data, and save temporarily
source("./src/data generation study 1.R")

# Change the names of each element in the list
# to be sure it corresponds to the right scenario
names(s1_data) <- data_files
df <- s1_data[[12]]

#### TESTING SITE FOR ERROR HANDLING
# Set i and V
i <- 12
V <- 5

#Preset which model, predictor selection methods and dgm
model <- study[i, ]$model
pred_selection <- study[i, ]$pred_selection
dgm_par <- c(study[i, ]$par1, 
             rep(study[i, ]$par2 * 3, round(0.3 * study[i, ]$dim)),  # strong
             rep(study[i, ]$par2,     round(0.5 * study[i, ]$dim)),  # weak
             rep(study[i, ]$par2 * 0, round(0.2 * study[i, ]$dim)))  # noise

#################################################################################################
######## Handling errors/warnings #########
## from https://stackoverflow.com/questions/3903157/how-can-i-check-whether-a-function-call-results-in-a-warning

ErrorsWarnings <- function(expr) {
  
  W <- NULL
  
  w_handler <- function(w){
    W <<- w
    invokeRestart("muffleWarning")
  }
  
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e), warning = w_handler), warning = W)
} 

get_cv_estimands_test <- function(df, model, dgm_par, pred_selection, V, x10 = c(FALSE, TRUE)){
  
  # If nothing 'bad' happens, the error info remains as NA
  error_info <- NA
  original_V <- V
  
  errors_warnings <- ErrorsWarnings({
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
    
    #  Check for var(y) == 0 within each fold:
    if (any(lapply(seq(V), function(V) var(df[folds[[V]],]$y)) == 0) == TRUE) {
     
      # Find out which fold has no 1:
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
      error_info <- paste(toString(error_info), paste0("No events sampled in folds ", paste(bad_folds_ind, collapse = " & ")), sep = " + ")
    
    } # Close error handling folds without events 
      # (should not happen because of stratified folds, but just to be sure)
    
    # Getting predictions depending on model used
    .doFit <- function(V, folds, model, pred_selection){ 
      
      # Create the training and test folds
      df_train <- as.data.frame(df[-(which(rownames(df) %in% folds[[V]])),])
      df_test <-  as.data.frame(df[which(rownames(df) %in% folds[[V]]),])
      assign("df_train", df_train, envir = .GlobalEnv)
      assign("df_test", df_test, envir = .GlobalEnv)
      
      # Train/test glm for each fold
      # Fit model depending on scenario
      # And get predicted probabilities and modelmatrix
      if (model == "Firth" & pred_selection == "none") { # If model = Firth , no predictor selection
        
        # Specifiy model formula (necessary for logistf)
        model_form <- as.formula(paste0("y ~", paste(colnames(df)[!colnames(df)%in%"y"], collapse = "+" )))
        
        fit <- logistf(formula = model_form, data = df_train, flic = T)
        iv_matrix <- model.matrix(object = fit$formula, data = df_test)
        p <- 1 / (1 + exp(-iv_matrix %*% fit$coefficients))
        dgm_par_folds <- dgm_par
        
      } else if (model == "Firth" & pred_selection == "<0.05") { # If model = Firth , backwards predictor selection
        
        model_form <- as.formula(paste0("y ~", paste(colnames(df)[!colnames(df)%in%"y"], collapse = "+" )))
        
        fit <- logistf(formula = model_form, data = df_train, flic = T)
        fit <- backward(fit, trace = FALSE)
        
        if (var(fit$linear.predictors) == 0) {
          
          error_info <- paste(toString(error_info), paste0("No predictors selected -> no calibration slope"), sep = " + ")
        }
        
        iv_matrix <- model.matrix(object = fit$formula, data = df_test)
        
        p <- 1 / (1 + exp(-iv_matrix %*% fit$coefficients))
        
        # Get the elements of the design generating mechanism that are
        # belonging to the model after backwards elimination
        # Always get the first element as this is the intercept
        ind <- na.omit(c(1, (as.numeric(str_extract_all(colnames(iv_matrix), "(?<=V).*"))
                             # Add 1, because the indices of columns exclude the intercept
                             + 1)))
        
        dgm_par_folds <- dgm_par[ind]
        
      # If using ML   
      } else if (model == "ML") {
        
        fit <- glm(y ~ ., family = "binomial", data = df_train)
        
        # Check for separation
        if(any(sqrt(diag(summary(fit)$cov.unscaled)*summary(fit)$dispersion) > 70)){
          error_info <- paste(toString(error_info), paste0("Data separation has occured"), sep = " + ")
          }
        
        p <- as.matrix(predict(fit, newdata = df_test , type = "response"))
        iv_matrix <- model.matrix(object = fit$formula, data = df_test)
        dgm_par_folds <- dgm_par
      
      # If the model is Lasso  
      } else if (model == "Lasso") {
       
        # Make sure that there are at least 8 events or no-events:
        if ((sum(df_train$y) > 7 | sum(1-df_train$y) > 7) == FALSE){
          
          fit <-  Pen_reg_VC(df = df_train, alpha = 0, nfolds = nrow(df_train))
          error_info <- paste(toString(error_info), paste0("too few (non-)events for tuning -> LOOCV"), sep = " + ")
        
          } else {
          
          fit <-  Pen_reg_VC(df = df_train, alpha = 0, nfolds = 10)
        
          } # close error handling 
        
        ## Check whether predictors have been selected
        # Create linear predictor
        lp <- predict(fit, newdata = df_train, s = "lambda.min") 
        
        if (var(lp) == 0) {
          error_info <- paste(toString(error_info), paste0("No predictors selected -> no calibration slope"), sep = " + ")
        }
        
        p <- predict(fit, newdata = df_test, s = "lambda.min", type = "response")
        
        # Create model matrix:
        # First retrieve the coefficients
        coefs <- as.matrix(coef(fit, s = "lambda.min"))
        # Remove those that do are 0
        coefs[coefs == 0] <- NA
        coefs <- na.exclude(coefs)
        
        # Then use the names of the coefficients to get the model matrix:
        df_test_X <- df_test[colnames(df_test) %in% rownames(coefs)]
        # Make model matrix by binding a column of 1 to the above
        iv_matrix <- as.matrix(cbind("(Intercept)" = 1, df_test_X))
        
        # Get the elements of the design generating mechanism that are
        # belonging to the model after backwards elimination
        # Always get the first element as this is the intercept
        ind <- na.omit(c(1, (as.numeric(str_extract_all(colnames(iv_matrix), "(?<=V).*"))
                             # Add 1, because the indices of columns exclude the intercept
                             + 1)))
        
        dgm_par_folds <- dgm_par[ind]
    
    # If the model is Ridge    
    } else if (model == "Ridge") {
      
      # Make sure that there are at least 8 events or no-events:
      if ((sum(df_train$y) <= 8 | sum(1-df_train$y) <= 8) == TRUE){
        fit <-  Pen_reg_VC(df = df_train, alpha = 0, nfolds = nrow(df_train))
        error_info <- paste(toString(error_info), paste0("too few (non-)events for tuning -> LOOCV"), sep = " + ")
      
        } else {
        
          fit <-  Pen_reg_VC(df = df_train, alpha = 0, nfolds = 10)
      } # close error handling 
      
      # Get the predictions 
      p <- predict(fit, newdata = df_test, s = "lambda.min", type = "response")
      
      # Create model matrix:
      # First retrieve the coefficients
      coefs <- as.matrix(coef(fit, s = "lambda.min"))
      # Then use the names of the coefficients to get the model matrix:
      df_test_X <- df_test[colnames(df_test) %in% rownames(coefs)]
      # Make model matrix by binding a column of 1 to the above
      iv_matrix <- as.matrix(cbind("(Intercept)" = 1, df_test_X))
      
      # Save dgm_par for all folds
      dgm_par_folds <- dgm_par
      
    } # Close model if else statements
    
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
    MAPE_folds <- c()
    
    # For each fold, calculate calibration intercept & slope, R2, ECI and MAPE
    for (v in 1:V){
      # create object to shorten future code:
      data <- df[unlist(which(rownames(df) %in% folds[[v]])),]
      ppf <- p_per_fold[[v]]
      dgm_par_ppf <- dgm_par_folds[[v]]
      
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
      MAPE_folds[v] <- MAPE(p = ppf, iv_matrix = iv_matrix[[v]], dgm_par = dgm_par_ppf)
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
    MAPE_results <- c(mean(MAPE_folds), (sd(MAPE_folds)/(sqrt(V))))
    
    
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
      
      results <- c(paste0(original_V, " fold cross-validation"), auc_results, intercept, slope, tjur_results, R2, eci, MAPE_results, error_info)
      
    } else {
      
      # If it is 10x10 cv, only keep the raw results of each fold and name them:  
      names(auc_folds) <- c(rep("auc", length(auc_folds)))
      names(intercept_folds) <- c(rep("intercept", length(intercept_folds)))
      names(slope_folds) <- c(rep("slope", length(slope_folds)))
      names(tjur_folds) <- c(rep("Tjur", length(tjur_folds)))
      names(R2_folds) <- c(rep("R2", length(R2_folds)))
      names(eci_folds) <- c(rep("eci", length(eci_folds)))
      names(MAPE_folds) <- c(rep("MAPE", length(MAPE_folds)))
      results <- list("auc" = auc_folds, 
                      "intercept" = intercept_folds, 
                      "slope" = slope_folds, 
                      "Tjur" = tjur_folds,
                      "R2" = R2_folds,
                      "eci" = eci_folds,
                      "MAPE" = MAPE_folds,
                      error_info)
    } # Close else statement (when 10k or 5k CV is used)
  
    return(results)
    
} # Close function


## Function to apply get_cv_estimands to all datasets in a study:
get_cv_results_test <- function(study, df, V, studyname) {
  
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
    dgm_par <- c(study[i, ]$par1, 
                 rep(study[i, ]$par2 * 3, round(0.3 * study[i, ]$dim)),  # strong
                 rep(study[i, ]$par2,     round(0.5 * study[i, ]$dim)),  # weak
                 rep(study[i, ]$par2 * 0, round(0.2 * study[i, ]$dim)))  # noise
    
    
    results_cv[i, which(colnames(results_cv) %in% iv_colnames)]  <-
      get_cv_estimands_test(
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


## Obtain internal validation estimands ##
# 10 fold cross-validation
results_10_cv_test <- get_cv_results_test(study = s1, df = s1_data, V = 10, studyname = "Study 1")

# 5 fold cross-validation
results_5_cv_test <- get_cv_results_test(study = s1, df = s1_data, V = 5, studyname = "Study 1")


get_10x10_results_test <- function(study, df, V, studyname){
  
  # results matrix:
  results_cv <- as.data.frame(matrix(NA, nrow = nrow(study), ncol = length(results_estimands_names), dimnames = list(c(), results_estimands_names)))
  
  # Fill in details of the study
  results_cv$study <- studyname
  # And of each scenario:
  results_cv[, which(colnames(results_cv) %in% study_info)] <- study[, which(colnames(study) %in% study_info)]
  
  for (i in 1:length(df)) {
    print(i)
    
    # Add which scenario we are working on:
    results_cv[i, 'scenario'] <- paste0("Scenario ", i)
    
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
      dgm_par <- c(study[i, ]$par1, 
                   rep(study[i, ]$par2 * 3, round(0.3 * study[i, ]$dim)),  # strong
                   rep(study[i, ]$par2,     round(0.5 * study[i, ]$dim)),  # medium
                   rep(study[i, ]$par2 * 0, round(0.2 * study[i, ]$dim)))  # noise
      
      # obtain the results of each fold separately 10 times and stored for each replication
      results <- replicate(n = 10, # represents the 10x10 part
                           expr = get_cv_estimands_test(df = data, 
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
      
      # Taking care of error messages returned for each repetition. 
      error_info <- paste(c(sapply(results, "[[", 8)), collapse = " + ")
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
          error_info
        )
      
    } # close for loop
  } # close else loop
  
  return(results_cv)
}

# 10X10 fold cross-validation 
results_10x10_cv_test <- get_10x10_results_test(study = s1[c(1:5),], df = s1_data[c(1:5)], V = 10, studyname = "Study 1")

