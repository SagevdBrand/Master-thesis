##################################
##################################
##################################
## get estimands 1 scenario:
set.seed(123)

############ Load necessary stuff ############
source("scripts/setup.R")
source("scripts/libraries.R")
#source("scripts/Scenarios.R") 
source("scripts/estimand functions.R")

#### Load data #####
data_files <- list.files(path = scenario_1_data, recursive = T, full.names = F)
df <- lapply(paste0(scenario_1_data,data_files),readRDS,.GlobalEnv)
s1 <- read_rds(paste0(scenario_1_settings,"s1.Rds"))
system.time(results_app <- get_app_results(scenario = s1, df = df))

###### try-out svm tuning ######
# test <- df[[1]]
# test_svm <- tune.svm(y~., data = test, scale = F, gamma = 2^(-1:1), cost = 2^(2:4), kernel = "linear")
# test_svm$best.model
# p_app <- predict(test_svm$best.model, test, probability = T)
# 
# 
# auc_app <- fastAUC(p = p_app, y = test$y)
# R2_app <- pseudo_Rsqrs(p = p_app, y = test$y)
# 
# 
# model.matrix(object = fit_app$formula, data = df)
# 
# calib_app <-
#   calib(
#     modelmatrix = app_matrix,
#     data = test,
#     coefs = coef(test_svm$best.model)
#   )
# 
# eci_app <-
#   eci_bvc(
#     data = df,
#     modelmatrix = app_matrix,
#     coefs = fit_app$coefficients,
#     preds = p_app
#   )
