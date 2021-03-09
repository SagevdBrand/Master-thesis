##########################################################
################## CONSTRUCTION SITE #####################
##########################################################

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
df <- s1_data[[1]]

#### TESTING SITE FOR ERROR HANDLING
# Set i and V
i <- 1
V <- 5

#Preset which model, predictor selection methods and dgm
model <- study[i, ]$model
pred_selection <- study[i, ]$pred_selection
dgm_par <- c(study[i, ]$par1, 
             rep(study[i, ]$par2 * 3, round(0.3 * study[i, ]$dim)),  # strong
             rep(study[i, ]$par2,     round(0.5 * study[i, ]$dim)),  # weak
             rep(study[i, ]$par2 * 0, round(0.2 * study[i, ]$dim)))  # noise



## FIRST CHECK ##
## NO EVENTS SAMPLED IN GENERATED SIMULATION DATA
## ESPECIALLY PROBLEM FOR LASSO AND RIDGE?
# https://stackoverflow.com/questions/12193779/how-to-write-trycatch-in-r
# https://stackoverflow.com/questions/1608130/equivalent-of-throw-in-r 

######## Setting up #########
## set.seed(123)


generate_data_test <- function(scenario, validation = c(TRUE, FALSE)){
  
  s_list <- split(scenario, seq(nrow(scenario))) 
  
  .per_element <- function(s_list){
    
   out <- tryCatch(
      {
        # This is the code of the function to be evaluated: 
        
        # Create covariance matrix to be used as input
        sigma <- matrix(0.2, ncol = s_list$dim, nrow = s_list$dim) 
        
        # set the diagonal to 1
        diag(sigma) <- 1
        
        # Define mu
        mu <- c(rep(0, s_list$dim))
        
        # create candidate predictors
        # If it is the validation set, n = 20 * 100/prev. Otherwise, take n as specified in each scenario
        if (validation == TRUE) {
          X <- mvrnorm(n = 20*(100/s_list$prev), mu = mu, Sigma = sigma)
        }
        else {
          X <-
            mvrnorm(n = s_list$n, mu = mu, Sigma = sigma)
        }
        
        # Putting the above in a data matrix, including intercept
        dm <- cbind(1, X)
        
        # Adding the data generating parameters as specified in each scenario
        dgm_par <- c(s_list$par1, 
                     rep(s_list$par2 * 3, round(0.3 * s_list$dim)),
                     rep(s_list$par2, round(0.5 * s_list$dim)),
                     rep(0, round(0.2 * s_list$dim))
        ) 
        
        # Obtain values for y based on Bernoulli distribution, with input p
        p <- 1/(1+exp(-dm %*% dgm_par))
        y <- as.numeric(p > runif(length(p)))
        
        # UNCOMMENT TO CHECK WHAT HAPPENS:
        y <- rep(0, length(p))
        
        
        # Check whether the variance of y is 0 
        # If it is 0, it means that no events were sampled
        if(var(y) == 0 ) stop("Error: No events sampled")
        
        as.data.frame(cbind(X, y))
      }, # closing expression to be tested
      # What happens
      error = function(e) {
        
        message(e, "\n")
        
        # Create an NA object with the error as its name
        erms <- NA
        names(erms) <- e
        # This should return an NA list instead of a dataframe
        return(erms)
        
            } # Closing error expression
      
      ) # Closing tryCatch function
    
   return(out)
   }
  
  # Apply function above
  data_in_lists <- lapply(s_list, .per_element) 
  return(data_in_lists)
  
}

s1_data <- generate_data_test(s1, validation = FALSE)
df <- s1_data[[1]]
###########################################################################
## Return proper result given s1_data[[x]] == "Error: No events sampled" ##
###########################################################################

get_app_estimands_test <- function(df, model, dgm_par, pred_selection) {
  
  # Check whether the data is actually useful
  if (any(str_detect(names(df),"Error: No events sampled") == TRUE)){
    
  # If no events were sampled, then the following will be the result:
  results <- c("Error: No events sampled")
  return(results) 
    
  } else {
  
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
    c("Apparent",
      auc_app,
      calib_app['intercept'],
      calib_app['slope'],
      R2_app,
      eci_app,
      MAPE_app, 
      NA # Because no error
    )
  
  } # End else statement
  } # End function


get_app_results_test <- function(study, df, studyname) {

  # object to store results
  results_app <- as.data.frame(matrix(NA, nrow = nrow(study), ncol = length(results_estimands_names), dimnames = list(c(), results_estimands_names)))
  
  # Fill in details of the study
  results_app$study <- studyname
  # And of each scenario:
  results_app[, which(colnames(results_app) %in% study_info)] <- study[, which(colnames(study) %in% study_info)]
  
  # For each scenario within the study
  for (i in 1:nrow(study)) {
    results_app[i, 'scenario'] <- paste0("Scenario ", i)
    # determine which model & pre_selection is specified
    model <- study[i, ]$model
    pred_selection <- study[i, ]$pred_selection
    
    # What do the model parameters look like in the
    # Data generating mechanism?
    dgm_par <- c(study[i, ]$par1, 
                 rep(study[i, ]$par2 * 3, round(0.3 * s1[i, ]$dim)),  # strong
                 rep(study[i, ]$par2,     round(0.5 * s1[i, ]$dim)),  # weaker
                 rep(study[i, ]$par2 * 0, round(0.2 * s1[i, ]$dim)))  # noise
    
    # Fill the columns with apparent results, no SE!
    results_app[i, which(colnames(results_app) %in% apparent_col_names)] <-
      get_app_estimands_test(
        df = df[[i]],
        model = model,
        dgm_par = dgm_par,
        pred_selection = pred_selection
      )
    
  } # end for loop

  return(results_app)
  
} # end function

  
  ## Obtain apparent estimands ##
  results_app_test <- get_app_results_test(study = s1, df = s1_data, studyname = "Study 1")

  
get_cv_estimands_test <- function(df, model, dgm_par, pred_selection, V, x10 = c(FALSE, TRUE)){
    
    # Check whether the data is actually useful
  if (any(str_detect(names(df),"Error: No events sampled") == TRUE)){
    
    # If no events were sampled, then the following will be the result:
    results <- c("Error: No events sampled")
    return(results) 
      
    } else {
      
      
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
      
      ####### BUILD CHECK FOR EVENTS SAMPLED IN FOLDS ######
      
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
      
      
      error_info <- NA # If everything went alright this should be NA
      
      # If the it is anything other than 10x10 cv:
      if (x10 == FALSE){
        
        results <- c(paste0(V, " fold cross-validation"), auc_results, intercept, slope, R2, eci, MAPE_results, error_info)
        
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
      } # Close else statement (when 10k or 5k CV is used)
      
      return(results)
      
    } # close else statement (when data is not NA)
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
    results_cv[i, 'scenario'] <- paste0("Scenario ", i)
    print(i)
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
  }
  
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
      
      
      dgm_par <- c(study[i, ]$par1, 
                   rep(study[i, ]$par2 * 3, round(0.3 * study[i, ]$dim)),  # strong
                   rep(study[i, ]$par2,     round(0.5 * study[i, ]$dim)),  # medium
                   rep(study[i, ]$par2 * 0, round(0.2 * study[i, ]$dim)))  # noise
      
      # obtain the results of each fold separately 10 times and stored for each replication
      results <- replicate(n = 10, # represents the 10x10 part
                           expr = get_cv_estimands_test(df = data, model = model, 
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
      
      ###############################################
      ## DONT FORGET TO CHECK AFTER ERROR HANDLING ##
      ###############################################
      error_info <- NA
      
      ## Fill results matrix:
      results_cv[i, which(colnames(results_cv) %in% iv_colnames)] <-
        c(
          "10x10 fold cross-validation",
          auc_results,
          intercept,
          slope,
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
  results_10x10_cv_test <- get_10x10_results_test(study = s1, df = s1_data, V = 10, studyname = "Study 1")
  
  # Bootstrap 3 varieties in one go
  
    
  
  ## Obtain external validation estimands
  
  
  # External validation
  


  