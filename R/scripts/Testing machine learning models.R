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
  
  #####################
  ## Neural networks ##
  #####################
  
  scaled_df <- as.data.frame(scale(df))
  
  # Use tune alone and then get best model
  fit <- e1071::tune(nnet, y ~., data = scaled_df[-foldsb[[V]],], 
                     ranges = list(size = 2^(0:3), decay = 2^(0:0.5)), trace = F)
  
  p <- predict(fit$best.model, newdata = scaled_df[foldsb[[V]],])
  fastAUC(p = p, y = df[foldsb[[V]],])
  
  # Use best.nnet, not sure how though
  fit <- e1071::best.nnet(x = scaled_df[-foldsb[[V]],-ncol(scaled_df)], y = scaled_df[-foldsb[[V]], ncol(scaled_df)], size = 1, decay = 0.5, tunecontrol = tune.control(sampling = "cross", random = T))
  p <- predict(fit, newdata = scaled_df[foldsb[[V]],])
  fastAUC(p = p, y = df[foldsb[[V]],])
  
  # Use caret::train, use tuneLength for random search
  system.time(fit_length <- train(as.factor(y) ~., data = scaled_df[-foldsb[[V]],], method = 'nnet', trace = F, tuneLength = 15))
  p <- predict(fit$best.model, newdata = scaled_df[foldsb[[V]],])
  fastAUC(p = p, y = df[foldsb[[V]],])
  
  
  # Use caret::train, predefining a tune grid
  tunegrid <- expand.grid(size = seq(from = 1, to = 10, by = 2),
                          decay = seq(from = 0.1, to = 0.5, by = 0.1))
  
  system.time(fit_grid <- train(as.factor(y) ~., data = scaled_df[-foldsb[[V]],], method = 'nnet', trace = F, tuneGrid = tunegrid))
  p <- predict(fit_grid, newdata = scaled_df[foldsb[[V]],])
  fastAUC(p = p, y = scaled_df[foldsb[[V]],])
  
  #########
  ## SVM ##
  #########
  
  fit <- e1071::tune(svm, y ~., data = df[-foldsb[[V]],], 
                     ranges = list(size = 2^(0:3), decay = 2^(0:0.5)), trace = F)
  
  p <- predict(fit$best.model, newdata = df[foldsb[[V]],])
  fastAUC(p = p, y = df[foldsb[[V]],])
  
  # Use best.svm, not sure how though
  fit <- e1071::best.svm(x = df[-foldsb[[V]],-ncol(df)], y = df[-foldsb[[V]], ncol(df)], size = 1, decay = 0.5, tunecontrol = tune.control(sampling = "cross", random = T))
  p <- predict(fit, newdata = df[foldsb[[V]],])
  fastAUC(p = p, y = df[foldsb[[V]],])
  
  # Use caret::train, use tuneLength for random search
  system.time(fit_length <- train(as.factor(y) ~., data = df[-foldsb[[V]],], method = 'svm', trace = F, tuneLength = 15))
  p <- predict(fit$best.model, newdata = df[foldsb[[V]],])
  fastAUC(p = p, y = df[foldsb[[V]],])
  
  
  # Use caret::train, predefining a tune grid
  tunegrid <- expand.grid(size = seq(from = 1, to = 10, by = 2),
                          decay = seq(from = 0.1, to = 0.5, by = 0.1))
  
  system.time(fit_grid <- train(as.factor(y) ~., data = df[-foldsb[[V]],], method = 'svm', trace = F, tuneGrid = tunegrid))
  p <- predict(fit_grid, newdata = df[foldsb[[V]],])
  fastAUC(p = p, y = df[foldsb[[V]],])

