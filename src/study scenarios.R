###############################################
###############################################
###############################################
### Define studies


#######################################################################################################################
#######################################################################################################################
############################################### WORK IN PROGRESS ######################################################
#######################################################################################################################
#######################################################################################################################

## Look at simulation protocol to define the final studies and scenarios!



### Setting up ###
set.seed(123)

source("./src/setup.R")
source("./src/data generation functions.R")

#####################
###### Study 1 ######
#####################

# Settings that do not depend on calculations
AUC1 <- 0.75
dim1 <- 10
es1 <- c(0.9,0.6,0.3)
prev1 <- c(0.05, 0.2, 0.5)
model1 <- "Firth"
pred_sel1 <- c("none", "<0.05")

## all combinations:
s1 <-
  expand.grid(
    AUC = AUC1,
    dim = dim1,
    shrinkage = es1,
    prev = prev1,
    model = model1,
    pred_selection = pred_sel1,
    KEEP.OUT.ATTRS = F
  )


### Expected R^2 ###

# Determine which R2 belongs to each situation:
# Here, a desired level of R^2 Cs is approximated for a given
# prevalence and C-statistic. 
# This can then be used to determine the minimal sample size needed 
# and to optimize beta coefficients.

R2 <- c(approximate_R2(auc = AUC1, prev = prev1[1])$R2.coxsnell,
        approximate_R2(auc = AUC1, prev = prev1[2])$R2.coxsnell, 
        approximate_R2(auc = AUC1, prev = prev1[3])$R2.coxsnell)

# Add the respective value of R2 depending on prev
s1 <- s1 %>% mutate(R2 = case_when(prev == prev1[1] ~ R2[1], 
                                   prev == prev1[2] ~ R2[2],
                                   prev == prev1[3] ~ R2[3]))

### Study scenario sample size ###

# Calculate sample size 
# Using pmsampsize package
n_1 <- c("e0.05_s0.9" = pmsampsize(type = "b", parameters = 10, prevalence = prev1[1], rsquared = R2[1])$sample_size,
                "e0.2_s0.9" = pmsampsize(type = "b", parameters = 10, prevalence = prev1[2], rsquared = R2[2])$sample_size,
                "e0.5_s0.9" = pmsampsize(type = "b", parameters = 10, prevalence = prev1[3], rsquared = R2[3])$sample_size,
                "e0.05_s0.6" = pmsampsize(type = "b", parameters = 10, prevalence = prev1[1], rsquared = R2[1], shrinkage = es1[2])$sample_size,
                "e0.2_s0.6" = pmsampsize(type = "b", parameters = 10, prevalence = prev1[2], rsquared = R2[2], shrinkage = es1[2])$sample_size,
                "e0.5_s0.6" = pmsampsize(type = "b", parameters = 10, prevalence = prev1[3], rsquared = R2[3], shrinkage = es1[2])$sample_size,
                "e0.05_s0.3" = pmsampsize(type = "b", parameters = 10, prevalence = prev1[1], rsquared = R2[1], shrinkage = es1[3])$sample_size,
                "e0.2_s0.3" = pmsampsize(type = "b", parameters = 10, prevalence = prev1[2], rsquared = R2[2], shrinkage = es1[3])$sample_size,
                "e0.5_s0.3" = pmsampsize(type = "b", parameters = 10, prevalence = prev1[3], rsquared = R2[3], shrinkage = es1[3])$sample_size)

# Bind the n to the study scenarios
s1 <- s1 %>%
  mutate(
    n = case_when(
      prev == 0.05 & shrinkage == 0.9 ~ n_1['e0.05_s0.9'],
      prev == 0.05 & shrinkage == 0.6 ~ n_1['e0.05_s0.6'],
      prev == 0.05 & shrinkage == 0.3 ~ n_1['e0.05_s0.3'],
      prev == 0.2 & shrinkage == 0.9 ~ n_1['e0.2_s0.9'],
      prev == 0.2 & shrinkage == 0.6 ~ n_1['e0.2_s0.6'],
      prev == 0.2 & shrinkage == 0.3 ~ n_1['e0.2_s0.3'],
      prev == 0.5 & shrinkage == 0.9 ~ n_1['e0.5_s0.9'],
      prev == 0.5 & shrinkage == 0.6 ~ n_1['e0.5_s0.6'],
      prev == 0.5 & shrinkage == 0.3 ~ n_1['e0.5_s0.3'],
      TRUE ~ NA_real_
      )
    )

# Remove attributes of 'n'
attr(s1$n, "ATT") <- NULL

### Coefficients ###

## Development data for coefficients
# setting n to develop betas on. 
n <- 30000 #                                                        

# For study 1, the number of candidate
# predictors is the same across all scenarios.
n_pred <- s1[1,]$dim

# create covariance matrix to be used as input
# Correlation between predictors is set to 0.2.
sigma <- matrix(0.2, ncol = n_pred, nrow = n_pred)
diag(sigma) <- 1 # set the diagonal to 1                                                          
mu <- c(rep(0, n_pred)) # provide a vector of values for mu -> standard normal                  
X <- mvrnorm(n = n, mu = mu, Sigma = sigma) # create predictor columns                        

# Model matrix (or data matrix, dm)
dm <- cbind(1, X)                                                                               
                                                                                                
## Validation data for coefficients                                                                              
# setting n  
n_val <- 100000                                                                    
X_val <- mvrnorm(n = n_val, mu = mu, Sigma = sigma)

# Putting the above in a data matrix, including intercept 
dm_val <- cbind(1, X_val)            

## Obtain regression coefficients
# Prevalence at 0.05 and AUC at 0.75
beta_0.05 <- optim_beta(prev_scenario = 0.05, R2_scenario= R2[1])
checking_val(par = beta_0.05)

# Prevalence at 0.2 and AUC at 0.75
beta_0.2 <- optim_beta(prev_scenario = 0.2, R2_scenario = R2[2])
checking_val(par = beta_0.2)

# Prevalence at 0.5 and AUC at 0.75
beta_0.5 <- optim_beta(prev_scenario = 0.5, R2_scenario = R2[3])
checking_val(par = beta_0.5)

## Bind parameters to study 1 matrix
# par1 represents intercept
# par2 represents regression coefficients
s1 <- s1 %>% 
  mutate(par1 = case_when(prev == 0.05 ~ beta_0.05[1],
                          prev == 0.2 ~ beta_0.2[1],
                          prev == 0.5 ~ beta_0.5[1],
                          TRUE ~ NA_real_
  )) %>%
  mutate(par2 = case_when(prev == 0.05 ~ beta_0.05[2],
                          prev == 0.2 ~ beta_0.2[2],
                          prev == 0.5 ~ beta_0.5[2],
                          TRUE ~ NA_real_))

### Save the study scenarios in the settings folder
write_rds(s1, file = study_1_settings)

#####################
###### Study 2 ######
#####################

## Data gen mechanism
AUC2 <- 0.75
dim2 <- c(6, 30, 60)
es2 <- c(0.3,0.6,0.9)
prev2 <- 0.2

## Models used:
model2 <-  "Firth"

## all combinations:
s2 <- expand.grid(AUC = AUC2, prev = prev2, dim = dim2, shrinkage = es2, model = models2, KEEP.OUT.ATTRS = F)

#################
## Expected R2 ##
#################

s2$R2 <- R2[2]

########################
## Actual sample size ##
########################

# at_n_10x_e.2 <- pmsampsize(type = "b", parameters = 10, prevalence = 0.2, rsquared = R2[2])$sample_size
# at_n_32x_e.2 <- pmsampsize(type = "b", parameters = 33, prevalence = 0.2, rsquared = R2[2])$sample_size
# at_n_60x_e.2 <- pmsampsize(type = "b", parameters = 60, prevalence = 0.2, rsquared = R2[2])$sample_size
# below_n_10x_e.2 <- ceiling(0.8 * at_n_10x_e.2)
# below_n_32x_e.2 <- ceiling(0.8 * at_n_32x_e.2)
# below_n_60x_e.2 <- ceiling(0.8 * at_n_60x_e.2)


#####################
###### Study 3 ######
#####################

## Data gen mechanism
AUC3 <- 0.75
dim3 <- 10
n3 <- as.factor(c("at", "below"))
prev3 <- 0.2

## Models combinations:
models3 <-  as.factor(c("OLS", "Ridge", "Firth", "LASSO", "Elastic Net", "SVM", "ANN", "RF"))

## all combinations:
s3 <- expand.grid(AUC = AUC3, prev = prev3, dim = dim3, n_state = n3, model = models3, KEEP.OUT.ATTRS = F)

#################
## Expected R2 ##
#################

s3$R2 <- R2[2]


########################
## Actual sample size ##
########################

# at_n_10x_e.2 <- pmsampsize(type = "b", parameters = 10, prevalence = 0.2, rsquared = R2[2])$sample_size
# below_n_10x_e.2 <- ceiling(0.8 * at_n_10x_e.2)


#####################
###### Study 4 ######
#####################

## Data gen mechanism
AUC4 <- 0.75
dim4 <- 30
n4 <- as.factor(c("60%", "80%", "at (100%)"))
prev4 <- 0.2
predictor_effects4 <- c("50% strong 50% noise")


## Models combinations:
models4 <-  as.factor(c("OLS", "Firth"))

## all combinations:
s4 <- expand.grid(AUC = AUC4, dim = dim4, n = n4,  prevalence = prev4, model = models4, pred_effect = predictor_effects4, KEEP.OUT.ATTRS = F)

#################
## Expected R2 ##
#################

s4$R2 <- R2[2]


########################
## Actual sample size ##
########################

# at_n_10x_e.2 <- pmsampsize(type = "b", parameters = 10, prevalence = 0.2, rsquared = R2[2])$sample_size
# below_n_10x_e.2 <- ceiling(0.8 * at_n_10x_e.2)


###################################
## Rounding the script to an end ##
###################################

# Remove everything except the three study matrices
rm(list=ls()[! ls() %in% c("s1","s2", "s3", "R2", "actual_n_1")])

####################################
####################################
####################################
## END SCRIPT


# # If more changes, for example:
# AUC <- c(0.6, 0.75, 0.9)
# prevalence <- c(0.05, 0.2)
# test <-  expand.grid(AUC = AUC, prevalence = prevalence)
# 
# test$R2 <- 0
# for (i in 1:nrow(test)){
#   test[i,3] <- approximate_R2(auc = test[i,1], prev = test[i,2])$R2.coxsnell
# }


