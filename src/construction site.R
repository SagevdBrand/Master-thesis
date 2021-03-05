##########################################################
################## CONSTRUCTION SITE #####################
##########################################################

## Get libraries, paths and a dataset:

## Libraries, file paths and functions
source("./src/setup.R")
source("./src/estimand functions.R")

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

get_cv_estimands(df = df, #df[[i]]
                 model = model,
                 dgm_par = dgm_par,
                 pred_selection = pred_selection,
                 V = V,
                 x10 = FALSE)


## FIRST CHECK ##
## NO EVENTS SAMPLED IN GENERATED SIMULATION DATA
## ESPECIALLY PROBLEM FOR LASSO AND RIDGE?
# https://stackoverflow.com/questions/12193779/how-to-write-trycatch-in-r
# https://stackoverflow.com/questions/1608130/equivalent-of-throw-in-r 

######## Setting up #########
## set.seed(123)
source("./src/data generation functions.R")

s1 <- read_rds(study_1_settings)

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
  if (str_detect(names(df),"Error: No events sampled") == TRUE){
    # If no events were sampled, then the following will be
    results <- list("Error: No events sampled" = NA)
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
    c(
      "AUC_app" = auc_app,
      "calib_int_app" = calib_app['intercept'],
      "calib_slope_app" = calib_app['slope'],
      "R2_cox_snell_app" = R2_app,
      "ECI_app" = eci_app,
      "MAPE_app" = MAPE_app
    )
  
  } # End else statement
  } # End function


get_app_results_test <- function(study, df) {

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
    
    results_app[[i]] <- get_app_estimands_test(df = df[[i]], model = model, dgm_par = dgm_par, pred_selection = pred_selection)
    
  } # end for loop
  
  names(results_app) <- c(1:nrow(study))
  return(results_app)
  
} # end function

  
  ## Obtain apparent estimands ##
  results_app_test <- get_app_results_test(study = s1, df = s1_data)
  
  ## Obtain internal validation estimands ##
  # 10 fold cross-validation
  results_10_cv_test <- get_cv_results(study = s1, df = s1_data, V = 10)
  
  # 5 fold cross-validation
  results_5_cv_test <- get_cv_results(study = s1, df = s1_data, V = 5)
  
  
  
  
  get_10x10_results_test <- function(study, df, V){
    
    results_cv <- list()
    
    for (i in 1:length(df)) {
      print(i)
      # Settings for get cv estimands function:
      model <- study[i, ]$model
      pred_selection <- study[i, ]$pred_selection
      data <- df[[i]]
      
      # Check whether there are events sampled
      if (str_detect(names(data),"Error: No events sampled") == TRUE){
        
        results_cv[[i]] <- list("Error: No events sampled" = NA)
        
      } else {
      
      
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
      
      } # close for loop
      } # close else loop
    
    names(results_cv) <- c(1:length(df))
    return(results_cv)
  }
  
  # 10X10 fold cross-validation 
  results_10x10_cv_test <- get_10x10_results_test(study = s1, df = s1_data, V = 10)
  
  # Bootstrap 3 varieties in one go
  
    
  
  ## Obtain external validation estimands
  
  
  # External validation
  


  