###############################################
###############################################
###############################################
### Define scenario's

### Setting up ###

set.seed(123)

source("scripts/Setup.R")
source("scripts/Data generation functions.R")

# Within the script above, the function "approximate_R2" is called. Here, a desired level of R^2 Cs is approximated for a given
# prevalence and C-statistic. This can then be used to optimize beta coefficients and to determine the minimal sample size needed. 

########################
###### Scenario 1 ######
########################

## Data gen mechanism
AUC1 <- 0.75
dim1 <- 10
es1 <- c(0.3,0.6,0.9)
prev1 <- c(0.05, 0.2, 0.5)
model1 <- "Firth"
## all combinations:
s1 <- expand.grid(AUC = AUC1, dim = dim1, Shrinkage = es1, prev = prev1, model = model1, KEEP.OUT.ATTRS = F)


##################
## Expected R^2 ##
##################

## Determine which R2 belongs to each situation:
R2 <- c(approximate_R2(auc = AUC1, prev = prev1[1])$R2.coxsnell,
        approximate_R2(auc = AUC1, prev = prev1[2])$R2.coxsnell, 
        approximate_R2(auc = AUC1, prev = prev1[3])$R2.coxsnell)

s1 <- s1 %>% mutate(R2 = case_when(prev == prev1[1] ~ R2[1], 
                                   prev == prev1[2] ~ R2[2],
                                   prev == prev1[3] ~ R2[3]))

########################
## Actual sample size ##
########################
###################################################
############# EXTEND FOR SHRINKAGE!!! #############
###################################################

actual_n_1 <- c("at_e_.05" = pmsampsize(type = "b", parameters = 10, prevalence = 0.05, rsquared = R2[1])$sample_size,
                "at_e_0.2" = pmsampsize(type = "b", parameters = 10, prevalence = 0.2, rsquared = R2[2])$sample_size,
                "at_e_0.5" = pmsampsize(type = "b", parameters = 10, prevalence = 0.5, rsquared = R2[3])$sample_size)

s1 <- s1 %>%
  mutate(
    n = case_when(
      prev == 0.05 & n_state == "at" ~ actual_n_1[1],
      prev == 0.05 & n_state == "below" ~ ceiling(0.8 * actual_n_1[1]),
      prev == 0.2 & n_state == "at" ~ actual_n_1[2],
      prev == 0.2 & n_state == "below" ~ ceiling(0.8 * actual_n_1[2]),
      prev == 0.5 & n_state == "at" ~ actual_n_1[3],
      prev == 0.5 & n_state == "below" ~ ceiling(0.8 * actual_n_1[3]),
      TRUE ~ NA_real_
      )
    )

##################
## coefficients ##
##################
beta_0.05 <- optim_beta(prev_scenario = 0.05, R2_scenario = R2[1])
beta_0.2 <- optim_beta(prev_scenario = 0.2, R2_scenario = R2[2])
beta_0.5 <- optim_beta(prev_scenario = 0.5, R2_scenario = R2[3])

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

write_rds(s1, file = paste0(scenario_1_settings, "s1.Rds"))


#######################################################################################################################
#######################################################################################################################
############################################### WORK IN PROGRESS ######################################################
#######################################################################################################################
#######################################################################################################################

########################
###### Scenario 2 ######
########################

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


########################
###### Scenario 3 ######
########################

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


########################
###### Scenario 4 ######
########################

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

# Remove everything except the three scenario matrices
rm(list=ls()[! ls() %in% c("s1","s2", "s3", "R2")])

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


