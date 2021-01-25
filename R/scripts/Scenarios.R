############### Scenario 1 #################
## Varying data gen mechanism
AUC1 <- 0.75
dim1 <- 10
n1 <- as.factor(c("at", "below"))
prev1 <- c(0.05, 0.2)

## Models used:
models1 <-  as.factor(c("OLS",  "Firth", "SVM"))

## all scenarios:
s1 <- expand.grid(AUC = AUC1, dim = dim1, n = n1, prev = prev1, model = models1)

############### Scenario 2 #################
## Varying data gen mechanism
AUC2 <- 0.75
dim2 <- c(6, 33, 60)
n2 <- as.factor(c("at", "below"))
prev2 <- 0.2

## Models used:
models2 <-  as.factor(c("OLS",  "Firth", "SVM"))

## all scenarios:
s2 <- expand.grid(AUC = AUC2, prev = prev2, dim = dim2, n = n2, model = models2)


############### Scenario 3 ##################
## Varying data gen mechanism
AUC3 <- 0.75
dim3 <- 10
n3 <- as.factor(c("at", "below"))
prev3 <- 0.2

## Models used:
models3 <-  as.factor(c("OLS", "Ridge", "Firth", "LASSO", "Elastic Net", "SVM", "ANN", "RF"))


## all scenarios:
s3 <- expand.grid(AUC = AUC3, dim = dim3, n = n3, prev = prev3, model = models3)

