###############################################
###############################################
################ Penalized logistic regression

## Libraries, file paths and functions
source("./src/setup.R")
source("./src/estimand functions.R")

library(glmnet)
# https://github.com/hongooi73/glmnetUtils maybe for ease of using!

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
system.time({lambda <- seq(64, 0, length.out = 251) # Is this Ben's approach?
fit_Ben <- cv.glmnet(x = as.matrix(df[,-ncol(df)]), y = factor(df[,ncol(df)]), family = "binomial", lambda = lambda, alpha = 0)})
## time: 1.36

## Results:
p_B <- predict(fit_Ben, as.matrix(df[,-ncol(df)]), s = "lambda.min", type = "response")
(c_B <- fastAUC(p = p_B, y = df$y))
(R2_B <- pseudo_Rsqrs(p = p_B, y = df$y))
(slope_B <-  c(coef(glm(df$y ~ log(p_B/(1-p_B)),family="binomial"))[2]))
(intercept_B <- coef(glm(df$y ~ offset(log(p_B/(1-p_B))),family="binomial")))
calout <- loess(y ~ log(p_B/(1-p_B)), data = df)
(eci_b <- (mean((p_B-fitted(calout))*(p_B-fitted(calout))))*(100))


# Maarten's approach:
system.time({
nlambda.ridge <- 200
ridge.defaultlambda <- glmnet(x=as.matrix(df[,-which(colnames(df)=="y")]),y=factor(df[,"y"]),family="binomial",alpha=0)$lambda
steps <- log(ridge.defaultlambda[2])-log(ridge.defaultlambda[1])
default.lambda <- c("min"=min(ridge.defaultlambda),"max"=max( ridge.defaultlambda),"length"=length(ridge.defaultlambda))
lambda.sequence <- exp(seq(from=log(max(ridge.defaultlambda)),by=steps, length.out=nlambda.ridge))
fit_M <- cv.glmnet(x=as.matrix(df[,-which(colnames(df)=="y")]),
                                       y=factor(df[,"y"]),family="binomial",
                     alpha = 0,lambda=lambda.sequence)})

## Results:
p_M <- predict(fit_M, as.matrix(df[,-ncol(df)]), s = "lambda.min", type = "response")
(c_M <- fastAUC(p = p_M, y = df$y))
(R2_M <- pseudo_Rsqrs(p = p_M, y = df$y))
(slope_M <-  c(coef(glm(df$y ~ log(p_M/(1-p_M)),family="binomial"))[2]))
(intercept_M <- coef(glm(df$y ~ offset(log(p_M/(1-p_M))),family="binomial")))
calout <- loess(y ~ log(p_M/(1-p_M)), data = df)
(eci_M <- (mean((p_M-fitted(calout))*(p_M-fitted(calout))))*(100))

## Time:
# 1.47