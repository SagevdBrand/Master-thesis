###############################################
###############################################
###############################################
### Define scenario's

#### To do:
## Add smart solution to define sample sizes!

### Setting up ###
library(tidyverse)
library(pmsampsize)
set.seed(123)

############### Scenario 1 #################
## Varying data gen mechanism
AUC1 <- 0.75
dim1 <- 10
n1 <- as.factor(c("at", "below"))
prev1 <- c(0.05, 0.2)

## Models used:
models1 <-  as.factor(c("OLS",  "Firth", "SVM"))

## all combinations:
s1 <- expand.grid(AUC = AUC1, dim = dim1, n_state = n1, prev = prev1, model = models1)

############### Scenario 2 #################
## Varying data gen mechanism
AUC2 <- 0.75
dim2 <- c(6, 32, 60)
n2 <- as.factor(c("at", "below"))
prev2 <- 0.2

## Models used:
models2 <-  as.factor(c("OLS",  "Firth", "SVM"))

## all combinations:
s2 <- expand.grid(AUC = AUC2, prev = prev2, dim = dim2, n_state = n2, model = models2)


############### Scenario 3 ##################
## Varying data gen mechanism
AUC3 <- 0.75
dim3 <- 10
n3 <- as.factor(c("at", "below"))
prev3 <- 0.2

## Models combinations:
models3 <-  as.factor(c("OLS", "Ridge", "Firth", "LASSO", "Elastic Net", "SVM", "ANN", "RF"))

## all combinations:
s3 <- expand.grid(AUC = AUC3, prev = prev3, dim = dim3, n_state = n3, model = models3)

#### Expected R^2 ####

source("scripts/Data generation functions.R")
# Within the script above, a desired level of R^2 Cs is approximated for a given
# prevalence and C-statistic. The same script also allows for the calculation of
# the minimum sample size in each simulation.

# # If more changes:
# AUC <- c(0.6, 0.75, 0.9)
# prevalence <- c(0.05, 0.2)
# expand.grid(AUC, prevalence)

## Determine which R2 belongs to each situation:
R2 <- c("prev_0.05" = approximate_R2(auc = AUC1, prev = 0.05)$R2.coxsnell,
        "prev_0.2" = approximate_R2(auc = AUC1, prev = 0.2)$R2.coxsnell)

s1 <- s1 %>% mutate(R2 = case_when(prev == prev1[1] ~ R2[1], prev == prev1[2] ~ R2[2]))
s2$R2 <- R2[2]
s3$R2 <- R2[2]


#### Actual sample size ####
at_n_10x_e.2 <- pmsampsize(type = "b", parameters = 10, prevalence = 0.2, rsquared = R2[2])$sample_size
at_n_10x_e.05 <- pmsampsize(type = "b", parameters = 10, prevalence = 0.05, rsquared = R2[1])$sample_size
# at_n_32x_e.2 <- pmsampsize(type = "b", parameters = 33, prevalence = 0.2, rsquared = R2[2])$sample_size
# at_n_60x_e.2 <- pmsampsize(type = "b", parameters = 60, prevalence = 0.2, rsquared = R2[2])$sample_size
# 
below_n_10x_e.2 <- ceiling(0.8 * at_n_10x_e.2)
below_n_10x_e.05<- ceiling(0.8 * at_n_10x_e.05)
# below_n_32x_e.2<- ceiling(0.8 * at_n_32x_e.2)
# below_n_60x_e.2<- ceiling(0.8 * at_n_60x_e.2)

#

s1 <- s1 %>%
  mutate(
    n = case_when(
      prev == 0.2 & n_state == "at" ~ at_n_10x_e.2,
      prev == 0.2 & n_state == "below" ~ below_n_10x_e.2,
      prev == 0.05 & n_state == "at" ~ at_n_10x_e.05,
      prev == 0.05 & n_state == "below" ~ below_n_10x_e.05,
      TRUE ~ NA_real_
      )
    )


# Remove everything except the three scenario matrices
rm(list=ls()[! ls() %in% c("s1","s2", "s3", "R2")])

####################################
####################################
####################################
## END SCRIPT
