##################### TEST ENVIRONMENT #######################
## Libraries, file paths and functions
source("scripts/setup.R")
source("scripts/estimand functions.R")
library(e1071)
library(caret)
library(nnet)

s1 <- read_rds(paste0(scenario_1_settings,"s1.Rds"))
data_files <- list.files(path = scenario_1_data, recursive = T, full.names = F)
test <- lapply(paste0(scenario_1_data,data_files),readRDS,.GlobalEnv)
df <- test[[1]]
V = 5


  ####################
  ## Splitting data ##
  ####################
  cvFoldsB <- function(Y, V) {  #Create Balanced CV folds (stratify by outcome)
    Y0 <- split(sample(which(Y=="0")), rep(1:V, length=length(which(Y==0))))
    Y1 <- split(sample(which(Y=="1")), rep(1:V, length=length(which(Y==1))))
    folds <- vector("list", length=V)
    for (v in seq(V)) {folds[[v]] <- c(Y0[[v]], Y1[[v]])}		
    return(folds)
  }
  
  foldsb <- cvFoldsB(Y = df$y, V = V)
  
  ###################################
  ## Testing machine learning code ##
  ###################################
  
  # Use tune alone and then get best model
  fit <- e1071::tune(nnet, y ~., data = df[-foldsb[[V]],], 
                     ranges = list(size = 2^(0:3), decay = 2^(0:0.5)), trace = F)
  
  p <- predict(fit$best.model, newdata = df[foldsb[[V]],])
  fastAUC(p = p, y = df[foldsb[[V]],])
  
  # Use best.nnet, not sure how though
  fit <- e1071::best.nnet(x = df[-foldsb[[V]],-ncol(df)], y = df[-foldsb[[V]], ncol(df)], size = 1, decay = 0.5, tunecontrol = tune.control(sampling = "cross", random = T))
  p <- predict(fit, newdata = df[foldsb[[V]],])
  fastAUC(p = p, y = df[foldsb[[V]],])
  
  # Use caret::train, use tuneLength for random search
  system.time(fit_length <- train(as.factor(y) ~., data = df[-foldsb[[V]],], method = 'nnet', trace = F, tuneLength = 15))
  p <- predict(fit$best.model, newdata = df[foldsb[[V]],])
  fastAUC(p = p, y = df[foldsb[[V]],])
  
  
  # Use caret::train, predefining a tune grid
  tunegrid <- expand.grid(size = seq(from = 1, to = 10, by = 2),
                          decay = seq(from = 0.1, to = 0.5, by = 0.1))
  
  system.time(fit_grid <- train(as.factor(y) ~., data = df[-foldsb[[V]],], method = 'nnet', trace = F, tuneGrid = tunegrid))
  p <- predict(fit_grid, newdata = df[foldsb[[V]],])
  fastAUC(p = p, y = df[foldsb[[V]],])
  
  
  
  
  
  #################################################
  ## Getting predictions depending on model used ##
  #################################################
  doFit <- function(V, folds, model){  #Train/test glm for each fold
    
    if (model == "OLS") {
      fit <- glm(y~., data=df[-foldsb[[V]],], family=binomial) #%>%
      #step(direction = "backward", trace = F) # This should be changeable depending on which thing you're using
      iv_matrix <- model.matrix(object = fit$formula, data = df[foldsb[[V]],]) # model matrix of test fold
      p <- predict(fit, newdata=df[foldsb[[V]],], type = "response")
      
    } else { # If model = Firth
      fit <- logistf(y ~ ., data = df[-foldsb[[V]],], flic = T)
      iv_matrix <- model.matrix(object = fit$formula, data = df[foldsb[[V]],]) # model matrix of test fold
      p <- 1 / (1 + exp(-iv_matrix %*% fit$coefficients))
    }
    
    coefs <- coef(fit) #obtain coefficients of fold
    results <- list(p, coefs, iv_matrix) # only save predictions, coefficients and model matrix
    
  }
  
  results <- (lapply(seq(V), .doFit, folds = foldsb, model = model)) # obtain model results for all folds
  p <- c(unlist(sapply(results, "[[", 1))) # Get out pred values as a single vector
  p[unlist(foldsb)] <- p #Re-order pred values, so they are in line with y
  
  p_per_fold <- lapply(results, "[[", 1) # getting the predictions per fold
  coefs <- (sapply(results, "[[", 2)) # Get out model coefficients
  iv_matrix <- lapply(results, "[[", 3) # get out the test matrices
  
  ######################  
  ## Obtain estimands ##
  ######################
  
  ## Empty objects to store results
  slope_folds <- c()
  intercept_folds <- c()
  R2_folds <- c()
  eci_folds <- c()
  
  # For each fold, calculate calibration intercept & slope, R2 and ECI
  for (v in 1:V){
    # for easy reference
    data <- df[unlist(foldsb[[v]]),]
    ppf <- p_per_fold[[v]]
    
    slope_folds[v] <- c(coef(glm(data$y ~ log(ppf/(1-ppf)),family="binomial"))[2])
    
    intercept_folds[v] <- coef(glm(data$y ~ offset(log(ppf/(1-ppf))),family="binomial"))
    
    R2_folds[v] <- pseudo_Rsqrs(p = ppf, y = data$y)
    
    calout <- loess(y ~ log(ppf/(1-ppf)), data = data)
    eci_folds[v] <- (mean((ppf-fitted(calout))*(ppf-fitted(calout))))*(100)
  }
  
  ## AUC, function by ledell already gives bias corrected se and ci as well
  auc_results <- as.vector(unlist(ci.cvAUC(predictions=p, labels=df$y, folds=foldsb, confidence=0.95)))[-5]
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
  
  results <- c(auc_results, intercept, slope, R2, eci)

