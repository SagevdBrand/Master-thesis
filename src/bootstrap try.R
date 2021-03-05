#########################
## BOOTSTRAP FUNCTIONS ##
#########################

## Create and load simulation data

source("./src/setup.R")
source("./src/estimand functions.R")
source("./src/data generation study 1.R") # Generate data, and save temporarily
data_files <- list.files(path = study_1_data, recursive = T, full.names = F) # get the data names
names(s1_data) <- data_files # Change the names of each element in the list, to be sure it corresponds to the right scenario

# Getting one dataframe
i <- 1 # or whatever we're interested in1
df1 <- s1_data[[i]]
model <- s1[i, ]$model
pred_selection <- s1[i, ]$pred_selection

dgm_par <- c(s1[i, ]$par1, 
               rep(s1[i, ]$par2 * 3, round(0.3 * s1[i, ]$dim)),  # strong
               rep(s1[i, ]$par2,     round(0.5 * s1[i, ]$dim)),  # medium
               rep(s1[i, ]$par2 * 0, round(0.2 * s1[i, ]$dim)))  # noise
  




# Check whether the data is actually useful
if (any(str_detect(names(df),"Error: No events sampled") == TRUE)) {
  # If no events were sampled, then the following will be
  results <- list("Error: No events sampled" = NA)
  return(results) 
  
} else {
  # START ACTUAL BOOTSTRAP FUNCTION:
  #  for (b in 1:nboot) {
  
  # Obtain indices for the bootstrap sample
  ind <- sample(row(df), size = nrow(df), replace = T) # indices
  
  # Create bootstrap sample and the out sample:
  bsamp <- df[ind,] # The bootstrap sample
  out <- df[-ind, ] # The out sample
  
  # Fit model depending on scenario
  # And get predicted probabilities and modelmatrix
  if (model == "Firth"){
    
    # It is either no prediction selection method
    if (pred_selection == "none") {
    
    # Fit model using sampled data
    fit_boot <- logistf(y ~ ., data = bsamp, flic = T)
    
    # Or backward elimination:
    } else {
    
      # Necessary for the backward elimination to work!
      assign("bsamp", as.data.frame(bsamp), envir = .GlobalEnv)
      model_form <- as.formula(paste0("y ~", paste(colnames(bsamp)[!colnames(bsamp)%in%"y"], collapse = "+" )))
      
      # Fit the model and apply backwards elimination:
      fit_boot <- logistf(formula = model_form, data = bsamp, flic = T)
      fit_boot <- backward( fit_boot, trace = FALSE)
    
      # Original data model matrix and predictions
      og_matrix <- model.matrix(object = fit_boot$formula, data = df)
    
      # Get the elements of the design generating mechanism that are
      # belonging to the model after backwards elimination
      # Always get the first element as this is the intercept
      ind <- na.omit(c(1, (as.numeric(str_extract_all(colnames(og_matrix), "(?<=V).*"))
                         # Add 1, because the indices of columns exclude the intercept
                         + 1)))
    
      dgm_par <- dgm_par[ind] 
      
      } # close inner if else statement predictor selection
    
    # After fitting either with or without
    # backwards elimination:
    
    # Original data model matrix and predictions
    og_matrix <- model.matrix(object = fit_boot$formula, data = df)
    og_p <- 1 / (1 + exp(-og_matrix %*% fit_boot$coefficients))
    
    # In bag bootstrap model matrix and predictions
    b_matrix <- model.matrix(object = fit_boot$formula, data = bsamp)
    b_p <- 1 / (1 + exp(-b_matrix %*% fit_boot$coefficients))
    
    # Out of bootstrap sample model matrix and predictions
    out_matrix <- model.matrix(object = fit_boot$formula, data = out)
    out_p <- 1 / (1 + exp(-out_matrix %*% fit_boot$coefficients))
  
  } else if (model == "OLS") {
   
    fit_boot <- glm(y ~ ., family = "binomial", data = bsamp) 
    
    # Original data model matrix and predictions:
    og_p <- predict(fit_boot, newdata = df, type = "response")
    og_matrix <- model.matrix(object = fit_boot$formula, data = df)
    
    # In bag bootstrap model matrix and predictions
    b_p <- predict(fit_boot, newdata = bsamp, type = "response")
    b_matrix <- model.matrix(object = fit_boot$formula, data = bsamp)
    
    # Out of bootstrap sample model matrix and predictions
    out_p <- predict(fit_boot, newdata = out, type = "response")
    out_matrix <- model.matrix(object = fit_boot$formula, data = out)

  } else if (model == "Ridge"){
  
  # Approach by Van Calster et al., 2020     
  fit_boot <- Pen_reg_VC(df = bsamp, alpha = 0) 
  
  og_p <- predict(fit_boot, as.matrix(df[,-ncol(df)]), s = "lambda.min", type = "response")
  ## ADD OTHER MATRICES AND PREDICTIONS    
  } else {
    # This is LASSO
    fit_boot <- Pen_reg_VC(df = bsamp, alpha = 1) 
    
    ## COPY MATRICES AND PREDICTIONS FROM RIDGE
  }
  
  
  
  }
  
  # Performance measures
  auc_ <- fastAUC(p = p_app, y = df$y)
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
} # close else statement
#} # Close function

# Calibration slope for folds
slope_folds[v] <- c(coef(glm(data$y ~ log(p/(1-p)), family="binomial"))[2])
# Calibration intercept for folds
intercept_folds[v] <- coef(glm(data$y ~ offset(log(p/(1-p))), family="binomial"))
# R2 for folds
R2_folds[v] <- pseudo_Rsqrs(p = p, y = data$y)
# ECI for folds
calout <- loess(y ~ log(p/(1-p)), data = data)
eci_folds[v] <- (mean((p-fitted(calout))*(p-fitted(calout))))*(100)
# MAPE for folds
MAPE_folds[v] <- MAPE(p = p, iv_matrix = iv_matrix[[v]], dgm_par = dgm_par)

