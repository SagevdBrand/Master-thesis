
get_app_ext_estimands <- function(df, df_val, model, dgm_par, pred_selection) {
  error_info <- NA
  
  # Check whether the data is actually useful
  if (any(str_detect(names(df),"Error: No events sampled") == TRUE)){
    
    # If no events were sampled, then the following will be the result:
    error_info <- c("Error: No events sampled")
    return(error_info) 
    
  } else {
    
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
        error_info <- paste(error_info, paste0("No predictors selected -> no calibration slope"), sep = " + ")
      }
      
      app_matrix <- model.matrix(object = fit_app$formula, data = df_app)
      p_app <- 1 / (1 + exp(-app_matrix %*% fit_app$coefficients))
      
      # External
      ext_matrix <- model.matrix(object = fit_app$formula, data = df_val)
      p_ext <- 1 / (1 + exp(-ext_matrix %*% fit_app$coefficients))
      
      # Get the elements of the design generating mechanism that are
      # belonging to the model after backwards elimination
      # Always get the first element as this is the intercept
      ind <- na.omit(c(1, (as.numeric(str_extract_all(colnames(app_matrix), "(?<=V).*"))
                           # Add 1, because the indices of columns exclude the intercept
                           + 1)))
      
      dgm_par_app <- dgm_par[ind]
      
      ## If model is using ML :  
    } else {
      
      fit_app <- glm(y ~ ., family = "binomial", data = df) 
      p_app <- predict(fit_app, type = "response")
      app_matrix <- model.matrix(object = fit_app$formula, data = df)
      
      # External
      ext_matrix <- model.matrix(object = fit_app$formula, data = df_val)
      p_ext <- predict(fit_app, newdata = df_val, type = "response")
      
    }

      # Obtain estimands
      auc_app <- fastAUC(p = p_app, y = df$y)
      R2_app <- pseudo_Rsqrs(p = p_app, y = df$y)
      MAPE_app <- MAPE(p = p_app, dgm_par = dgm_par_app, iv_matrix = app_matrix) 
      tjur_app <- tjur(p = p_app, y = df$y)
      
      calib_app <- calib(modelmatrix = app_matrix,
                         data = df,
                         coefs = fit_app$coefficients)
      
      eci_app <- eci_bvc(data = df,
                         modelmatrix = app_matrix,
                         coefs = fit_app$coefficients,
                         preds = p_app)
      
      # obtain external results
     
      auc_ext <- fastAUC(p = p_ext, y = df_val$y)
      R2_ext <- pseudo_Rsqrs(p = p_ext, y = df_val$y)
      MAPE_ext <- MAPE(p = p_ext, dgm_par = dgm_par_app, iv_matrix = ext_matrix) 
      tjur_ext <- tjur(p = p_ext, y = df_val$y)
      
      calib_ext <- calib(modelmatrix = ext_matrix,
                         data = df_val,
                         coefs = fit_app$coefficients)
      
      eci_ext <- eci_bvc(data = df_val,
                         modelmatrix = ext_matrix,
                         coefs = fit_app$coefficients,
                         preds = p_ext)
      
      # Save results
      results <- list(
        c("Apparent",
             auc_app,
             calib_app['intercept'],
             calib_app['slope'],
             tjur_app,
             R2_app,
             eci_app,
             MAPE_app,
             error_info
        ),
        c("External",
             auc_ext,
             calib_ext['intercept'],
             calib_ext['slope'],
             tjur_ext,
             R2_ext,
             eci_ext,
             MAPE_ext,
             error_info
        )
      )
    
  } # End else statement checking for events
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
      # and all the other columns of the estimands remain NA
      
    } else {
      
      results_app[i, 'observed events'] <- sum(df[[i]]$y)
      results_ext[i, 'observed events'] <- sum(df_val[[i]]$y)
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

test <- get_app_ext_results(study = s1, df = s1_data, df_val = s1_val_data, studyname = "s1")
