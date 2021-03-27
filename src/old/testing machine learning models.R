##################### TEST ENVIRONMENT #######################
## Libraries, file paths and functions
source("./src/setup.R")
source("./src/estimand functions.R")
source("./src/data generation functions.R")
## Specific libraries for this study
library(e1071)
library(caret)
library(nnet)
library(rpart)
library(randomForest)
library(kernlab)

s1 <- read_rds(study_1_settings)
s3 <- read_rds(study_3_settings)

s1_data <- generate_data(s1, validation = FALSE)
df <- s1_data[[4]]
V = 5
study <- s1
i <- 4
dgm_par <-
  c(study[i, ]$par1, 
    rep(study[i, ]$par2 * 3, round(1/6 * study[i, ]$dim)), # strong
    rep(study[i, ]$par2, round(2/6 *  study[i, ]$dim)),    # weak
    rep(study[i, ]$par2 * 0, round(3/6 * study[i, ]$dim))) # noise


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
  df_train <- as.data.frame(df[-(which(rownames(df) %in% folds[[V]])),])
  df_test <-  as.data.frame(df[which(rownames(df) %in% folds[[V]]),])
  
  #########
  ## SVM ##
  #########

  # Use caret::train, predefining a tune grid
  svm_grid <- expand.grid(.C = seq(0.001, 10, length.out = 10),
                          .sigma = seq(0.001, 1, length.out = 10))
  

  invisible(capture.output(fit <- caret::train(as.factor(y) ~.,
                      data = df_train,
                      method = 'svmRadial',
                      tuneGrid = svm_grid,
                      trControl = trainControl(method = "cv"),
                      prob.model = TRUE
                  )))

  p <- predict(fit$finalModel, newdata = df_test[,-ncol(df)], type = "prob")[,2]

  coefs <- as.matrix(fit$coefnames)
  rownames(coefs) <- fit$coefnames
  # Then use the names of the coefficients to get the model matrix:
  df_test_X <- df_test[colnames(df_test) %in% rownames(coefs)]
  # Make model matrix by binding a column of 1 to the above
  iv_matrix <- as.matrix(cbind("(Intercept)" = 1, df_test_X))

  p_true <- 1 / (1 + exp(-iv_matrix %*% dgm_par))
  
  svm_auc <- fastAUC(p = p, y = df_test$y)
  svm_slope <- c(coef(glm(df_test$y ~ log(p/(1 - p)), family="binomial"))[2])
  svm_intercept <- coef(glm(df_test$y ~ offset(log(p/(1-p))), family="binomial"))
  svm_tjur<- tjur(p = p, y = df_test$y)
  svm_R2 <- pseudo_Rsqrs(p = p, y = df_test$y)
  calout <- loess(y ~ log(p/(1-p)), data = df_test, span = 5)
  svm_eci <- (mean((p-fitted(calout))*(p-fitted(calout))))*(100)
  # MAPE for folds
  svm_mape <- mean(abs(p_true-p))
  # rMSPE for folds
  svm_rmspe <- sqrt(mean((p_true-p)^2))
  
  
  
  
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
  
  fit <- caret::train(as.factor(y) ~.,
                      data = df_train, 
                      method = 'nnet', 
                      tuneGrid = nnet_grid, 
                      trace = F, 
                      trControl = trainControl(method = "cv"),
                      preProcess = c("center", "scale"))
  

  p <- predict(fit$finalModel, newdata = df_test, method = "prob")
  coefs <- as.matrix(fit$coefnames)
  rownames(coefs) <- fit$coefnames
  # Then use the names of the coefficients to get the model matrix:
  df_test_X <- df_test[colnames(df_test) %in% rownames(coefs)]
  # Make model matrix by binding a column of 1 to the above
  iv_matrix <- as.matrix(cbind("(Intercept)" = 1, df_test_X))
  
  p_true <- 1 / (1 + exp(-iv_matrix %*% dgm_par))
  
  nnet_auc <- fastAUC(p = p, y = df_test$y)
  nnet_slope <- c(coef(glm(df_test$y ~ log(p/(1 - p)), family="binomial"))[2])
  nnet_intercept <- coef(glm(df_test$y ~ offset(log(p/(1-p))), family="binomial"))
  nnet_tjur<- tjur(p = p, y = df_test$y)
  nnet_R2 <- pseudo_Rsqrs(p = p, y = df_test$y)
  calout <- loess(y ~ log(p/(1-p)), data = df_test, span = 5)
  nnet_eci <- (mean((p-fitted(calout))*(p-fitted(calout))))*(100)
  # MAPE for folds
  nnet_mape <- mean(abs(p_true-p))
  # rMSPE for folds
  nnet_rmspe <- sqrt(mean((p_true-p)^2))
  


###################
## Random Forest ##
###################
rf_grid <- expand.grid(.mtry = seq(1, ncol(df_train)-1, by = 1))

fit <- caret::train(as.factor(y) ~.,
                    data = df_train, 
                    method = 'rf',
                    tuneGrid = rf_grid,
                    trControl = trainControl(method = "cv")
)


coefs <- as.matrix(train_rf$coefnames)
rownames(coefs) <- (train_rf$coefnames)
# Then use the names of the coefficients to get the model matrix:
df_test_X <- df_test[colnames(df_test) %in% rownames(coefs)]
# Make model matrix by binding a column of 1 to the above
iv_matrix <- as.matrix(cbind("(Intercept)" = 1, df_test_X))



p_rf <- predict(train_rf$finalModel, newdata = df_test, type = "prob")[,2]
p_true <- 1 / (1 + exp(-iv_matrix %*% dgm_par))

rf_auc <- fastAUC(p = p_rf, y = df_test$y)
rf_slope <- c(coef(glm(df_test$y ~ log(p_rf/(1 - p_rf)), family="binomial"))[2])
rf_intercept <- coef(glm(df_test$y ~ offset(log(p_rf/(1-p_rf))), family="binomial"))
rf_tjur<- tjur(p = p_rf, y = df_test$y)
rf_R2 <- pseudo_Rsqrs(p = p_rf, y = df_test$y)
calout <- loess(y ~ log(p_rf/(1-p_rf)), data = df_test)
rf_eci <- (mean((p_rf-fitted(calout))*(p_rf-fitted(calout))))*(100)
# MAPE for folds
rf_mape <- mean(abs(p_true-p_rf))
# rMSPE for folds
rf_rmspe <- sqrt(mean((p_true-p_rf)^2))




##########
## CART ##
##########
rpart_grid <- expand.grid(.cp = seq(0, 0.3, by = 0.01))

fit_app <- caret::train(as.factor(y) ~.,
                        data = df,
                        method = 'rpart',
                        tuneGrid = rpart_grid,
                        maximize = FALSE,
                        parms = list(split = "information"),
                        trControl = trainControl(method = "cv"))

coefs <- as.matrix(fit_app$coefnames)
rownames(coefs) <- fit_app$coefnames

# Check whether there were more than 1 splits
if (nrow(fit_app$finalModel$cptable) == 1) {
  ind <-which(fit_app$results$Accuracy != max(fit_app$results$Accuracy))
  best_tune <- fit_app$results[ind,] %>% slice(which.max(Accuracy)) %>% .$cp
  model_form <- as.formula(paste0("y ~", paste(colnames(df)[!colnames(df)%in%"y"], collapse = "+" )))
  
  fit_app <- caret::train(as.factor(y) ~.,
                          data = df,
                          method = 'rpart',
                          parms = list(split = "information"),
                          control = rpart.control(cp = best_tune),
                          trControl = trainControl(method = "none")
                         )
  
  
  coefs <- colnames(df)[!colnames(df)%in%"y"]
  rownames(coefs) <- fit_app$coefnames
}

# Then use the names of the coefficients to get the model matrix:
df_app_X <- df[colnames(df) %in% rownames(coefs)]
# Make model matrix by binding a column of 1 to the above
app_matrix <- as.matrix(cbind("(Intercept)" = 1, df_app_X))


p_app <- predict(fit_app$finalModel, newdata = df, type = "prob")[,2]
p_true <- 1 / (1 + exp(-app_matrix %*% dgm_par))

rpart_auc <- fastAUC(p = p_app, y = df$y)
rpart_slope <- c(coef(glm(df_test$y ~ log(p/(1 - p)), family="binomial"))[2])
rpart_intercept <- coef(glm(df_test$y ~ offset(log(p/(1-p))), family="binomial"))
rpart_tjur<- tjur(p = p, y = df_test$y)
rpart_R2 <- pseudo_Rsqrs(p = p, y = df_test$y)
calout <- loess(y ~ log(p/(1-p)), data = df_test, span = 5)
rpart_eci <- (mean((p-fitted(calout))*(p-fitted(calout))))*(100)
# MAPE for folds
rpart_mape <- mean(abs(p_true-p))
# rMSPE for folds
rpart_rmspe <- sqrt(mean((p_true-p)^2))




 
 