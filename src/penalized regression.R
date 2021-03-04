###############################################
###############################################
################ Penalized logistic regression

## Libraries, file paths and functions
source("./src/setup.R")
source("./src/estimand functions.R")

library(glmnet)
library(microbenchmark)
# https://github.com/hongooi73/glmnetUtils maybe for ease of using!

## Getting neceassary files loaded
s1 <- read_rds(paste0(study_1_settings))
data_files <- list.files(path = study_1_data, recursive = T, full.names = F)
test <- lapply(paste0(study_1_data,data_files),readRDS,.GlobalEnv)
df <- test[[1]]
test_matrix <- cbind("(Intercept)"= 1, df[,-which(colnames(df)=="y")])
set.seed(123)

## Van Calster 2020:
# We used 10-fold cross-validation to find the value for k that minimized the deviance,
# using a grid of 251 possible values between 0 (no shrinkage) and 64 (very large shrinkage). 
# The 250 non-null values were equidistant on logarithmic scale.

# To check whether the distance is equal:
test_distance <- function(x){
  x[1]-x[2] == x[2]-x[3]
}

lambda_old <- seq(64, 0, length.out = 251)
test_distance(log(lambda_old)) # It is not equidistant on the log scale

lambda <- c(exp(seq(log(64), log(0.00001), length.out = 250)), 0)
test_distance(log(lambda)) # This is equidistant on the log scale :)

# Ben's approach:
Ben_function <- function(df){
lambda <- c(exp(seq(log(64), log(0.00001), length.out = 250)), 0) # DOUBLE CHECK THIS
fit_Ben <- cv.glmnet(x = as.matrix(df[, -ncol(df)]),
                     y = factor(df[, ncol(df)]),
                     family = "binomial",
                     lambda = lambda,
                     alpha = 0
                     )

return(fit_Ben)
}


# Maarten's approach:
Maarten_function <- function(df){

  nlambda.ridge <- 200
  ridge.defaultlambda <- glmnet(x = as.matrix(df[, -which(colnames(df) == "y")]),
                                y = factor(df[, "y"]),
                                family = "binomial",
                                alpha = 0
                                )$lambda
  steps <- log(ridge.defaultlambda[2]) - log(ridge.defaultlambda[1])
  
  default.lambda <-
    c(
      "min" = min(ridge.defaultlambda),
      "max" = max(ridge.defaultlambda),
      "length" = length(ridge.defaultlambda)
    )
  
  lambda.sequence <-
    exp(seq(
      from = log(max(ridge.defaultlambda)),
      by = steps,
      length.out = nlambda.ridge
    ))
  
  fit_M <- cv.glmnet(x = as.matrix(df[, -which(colnames(df) == "y")]),
                     y = factor(df[, "y"]),
                     family = "binomial",
                     alpha = 0,
                     lambda = lambda.sequence
                     )
return(fit_M)
}

## Results:
fit_results <- function(fit){
p <- predict(fit, as.matrix(df[,-ncol(df)]), s = "lambda.min", type = "response")
cstat <- fastAUC(p = p, y = df$y)
R2 <- pseudo_Rsqrs(p = p, y = df$y)
c_slope <-  c(coef(glm(df$y ~ log(p/(1-p)),family="binomial"))[2])
c_intercept <- coef(glm(df$y ~ offset(log(p/(1-p))),family="binomial"))
calout <- loess(y ~ log(p/(1-p)), data = df)
eci <- (mean((p-fitted(calout))*(p-fitted(calout))))*(100)

c(cstat, R2, c_slope, c_intercept, eci)
}


microbenchmark(Ben_function, Maarten_function, times=30)

fit_M <- replicate(n = 10, Maarten_function(df = df), simplify = F)
results_M <- lapply(fit_M, fit_results)

fit_B <-  replicate(n = 10, Ben_function(df = df), simplify = F)
results_B <- lapply(fit_B, fit_results)

summary(sapply(results_M, "[[",1))
summary(sapply(results_B, "[[",1))

