##################### TEST ENVIRONMENT #######################
## Libraries, file paths and functions
source("./src/setup.R")
source("./src/estimand functions.R")

## Specific libraries for this study
library(e1071)
library(caret)
library(nnet)
library(rpart)
library(randomForest)
library(kernlab)

s1 <- read_rds(study_1_settings)
data_files <- list.files(path = study_1_data, recursive = T, full.names = F)
test <- lapply(paste0(study_1_data,data_files),readRDS,.GlobalEnv)
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
  
  folds <- cvFoldsB(Y = df$y, V = V)
  
  ###################################
  ## Testing machine learning code ##
  ###################################
machine_learning <- function(folds, V){
  
  df_train <- as.data.frame(df[-(which(rownames(df) %in% folds[[V]])),])
  df_test <-  as.data.frame(df[which(rownames(df) %in% folds[[V]]),])
  
  ###################
  ## Random Forest ##
  ###################
  rf_grid <- expand.grid(.mtry = seq(1, ncol(df_train)-1, by = 1))
  
  system.time(train_rf <- caret::train(as.factor(y) ~.,
                                       data = df_train, 
                                       method = 'rf',
                                       tuneGrid = rf_grid,
                                       trControl = trainControl(method = "cv")
                                       ))
  p_rf <- predict(train_rf$finalModel, newdata = df_test, type = "response")
  print(fastAUC(p = p_rf, y = df_test$y))
  
  
  ##########
  ## CART ##
  ##########
  rpart_grid <- expand.grid(.cp = seq(0, 0.3, by = 0.01))
  
  system.time(train_rpart <- caret::train(as.factor(y) ~.,
                                          data = df_train,
                                          method = 'rpart',
                                          tuneGrid = rpart_grid,
                                          trControl = trainControl(method = "cv")))
  
  p_rpart <- predict(train_rpart$finalModel, newdata = df_test, type = "prob")
  print(fastAUC(p = p_rpart, y = df_test$y))
  
 
  #########
  ## SVM ##
  #########

  # Use caret::train, predefining a tune grid
  svm_grid <- expand.grid(.C = seq(0.001, 10, length.out = 10),
                          .sigma = seq(0.001, 1, length.out = 10))
  
  # Data are scaled internally!
  system.time(fit_svm <- caret::train(as.factor(y) ~., 
                                      data = df_train,
                                      method = 'svmRadial',
                                      tuneGrid = svm_grid,
                                      trControl = trainControl(method = "cv"),
                                      preProcess = c("center", "scale")))
  
  ##  THIS DOES NOT MAKE SENSE OT ME!! ##
  p_svm <- kernlab::predict(fit_svm$finalModel, newdata = df_test[,-ncol(df)], type = "response")
  print(fastAUC(p = p_svm, y = df_test$y))
  
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
  
  system.time(fit_nnet <- caret::train(as.factor(y) ~., 
                                       data = df_train, 
                                       method = 'nnet', 
                                       tuneGrid = nnet_grid, 
                                       trace = F, 
                                       trControl = trainControl(method = "cv"),
                                       preProcess = c("center", "scale")))
  
  p_nnet <- predict(fit_nnet$finalModel, newdata = df_test)
  print(fastAUC(p = p_nnet, y = df_test$y))
  

}

results <- system.time(lapply(seq(V), machine_learning, folds = folds))

 
 