##################### TEST ENVIRONMENT #######################
## Libraries, file paths and functions
source("scripts/setup.R")
source("scripts/estimand functions.R")
library(e1071)
library(caret)
library(nnet)
library(rpart)
library(randomForest)
library(kernlab)

s1 <- read_rds(paste0(scenario_1_settings,"s1.Rds"))
data_files <- list.files(path = scenario_1_data, recursive = T, full.names = F)
test <- lapply(paste0(scenario_1_data,data_files),readRDS,.GlobalEnv)
df <- test[[3]]
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
  
  ###################
  ## Random Forest ##
  ###################
  rf_grid <- expand.grid(.mtry = seq(1, 10, by = 1))
  
  system.time(train_rf <- caret::train(as.factor(y) ~., data = df[-foldsb[[V]],], method = 'rf', tuneGrid = rf_grid, trControl = trainControl(method = "cv")))
  p_rf <- predict(train_rf$finalModel, newdata = df[foldsb[[V]],], type = "response")
  fastAUC(p = p_rf, y = df[foldsb[[V]],]$y)
  
  
  ##########
  ## CART ##
  ##########
  rpart_grid <- expand.grid(.cp = seq(0, 0.3, by = 0.01))
  
  system.time(train_rpart <- caret::train(as.factor(y) ~., data = df[-foldsb[[V]],], method = 'rpart', tuneGrid = rpart_grid, trControl = trainControl(method = "cv")))
  p_rpart <- predict(train_rpart$finalModel, newdata = df[foldsb[[V]],], type = "prob")
  fastAUC(p = p_rpart, y = df[foldsb[[V]],]$y)
  
 
  #########
  ## SVM ##
  #########

  # Use caret::train, predefining a tune grid
  svm_grid <- expand.grid(.C = seq(0.001, 10, length.out = 10),
                          .sigma = seq(0.001, 1, length.out = 10))
  
  system.time(fit_svm <- caret::train(as.factor(y) ~., data = df[-foldsb[[V]],], method = 'svmRadial', tuneGrid = svm_grid, trControl = trainControl(method = "cv")))
  
  ##  THIS DOES NOT MAKE SENSE OT ME!! ##
  p_svm <- kernlab::predict(fit_svm$finalModel, newdata = df[foldsb[[V]],-ncol(df)], type = "response")
  fastAUC(p = p_svm, y = df[foldsb[[V]],]$y)
  
  #####################
  ## Neural networks ##
  #####################
  # PRE PROCESSING!
  # preProcValues <- preProcess(training, method = c("center", "scale"))
  # 
  # trainTransformed <- predict(preProcValues, training)
  # testTransformed <- predict(preProcValues, test)
  
  nnet_grid <- expand.grid(size = seq(from = 1, to = 10, by = 2),
                          decay = seq(from = 0.1, to = 2, by = 0.4))
  
  system.time(fit_nnet <- caret::train(as.factor(y) ~., data = df[-foldsb[[V]],], method = 'nnet', tuneGrid = nnet_grid, trace = F, trControl = trainControl(method = "cv")))
  p_nnet <- predict(fit_nnet$finalModel, newdata = df[foldsb[[V]],])
  fastAUC(p = p_nnet, y = df[foldsb[[V]],]$y)
  
 
 
 