###############################################
###############################################
###############################################
### Define studies
source("./src/setup.R")

## Look at simulation protocol to define the final studies and scenarios!

#####################
###### Study 1 ######
#####################

# Settings that do not depend on calculations
AUC1 <- 0.75
dim1 <- 10
n_setting1 <- c("n/2", "n", "n*2")
prev1 <- c(0.05, 0.2, 0.5)
model1 <- "ML"
pred_sel1 <- c("none", "<0.157")
noise_1 <- c("default")

## all combinations:
s1 <-
  expand.grid(
    AUC = AUC1,
    dim = dim1,
    n_setting = n_setting1,
    noise = noise_1,
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

# Within the testing site optimization script,
# The following values were obtained:
R2 <- c(0.04131983, 0.12658143, 0.18407039) 

# Add the respective value of R2 depending on prev
s1 <- s1 %>% mutate(R2 = case_when(prev == prev1[1] ~ R2[1], 
                                   prev == prev1[2] ~ R2[2],
                                   prev == prev1[3] ~ R2[3]))

# Calculate sample size 
# Using pmsampsize package
n_1 <- c("e0.05" = pmsampsize(type = "b", parameters = 10, prevalence = prev1[1], rsquared = R2[1])$sample_size,
                "e0.2" = pmsampsize(type = "b", parameters = 10, prevalence = prev1[2], rsquared = R2[2])$sample_size,
                "e0.5" = pmsampsize(type = "b", parameters = 10, prevalence = prev1[3], rsquared = R2[3])$sample_size
         )
                
# Bind the n to the study scenarios
s1 <- s1 %>%
  mutate(
    n = ceiling(case_when(
      prev == 0.05 & n_setting1 == "n/2" ~ n_1["e0.05"]/2,
      prev == 0.05 & n_setting1 == "n" ~ n_1["e0.05"],
      prev == 0.05 & n_setting1 == "n*2" ~ n_1["e0.05"]*2,
      prev == 0.2 & n_setting1 == "n/2" ~ n_1["e0.2"]/2,
      prev == 0.2 & n_setting1 == "n" ~ n_1["e0.2"],
      prev == 0.2 & n_setting1 == "n*2" ~ n_1["e0.2"]*2,
      prev == 0.5 & n_setting1 == "n/2" ~ n_1["e0.5"]/2,
      prev == 0.5 & n_setting1 == "n" ~ n_1["e0.5"],
      prev == 0.5 & n_setting1 == "n*2" ~ n_1["e0.5"]*2,
      TRUE ~ NA_real_
      )
    )
    )

# Remove attributes of "n"
attr(s1$n, "ATT") <- NULL

####################
### Coefficients ###
####################
### Using the testing site optimization script
## The following values for the parameters were 
## found during optimization:

# Prevalence at 0.05 and AUC at 0.75
beta_0.05 <- c(-3.3705418, 0.1222722)

# Prevalence at 0.2 and AUC at 0.75
beta_0.2 <- c(-1.6740826, 0.1280781)

# Prevalence at 0.5 and AUC at 0.75
beta_0.5 <- c(0.01403259, 0.13144766)

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


s1$scenario <- c(paste0("Scenario_", 1:nrow(s1)))
### Save the study scenarios in the settings folder
write_rds(s1, file = study_1_settings)

#####################
###### Study 2 ######
#####################

## Data gen mechanism
AUC2 <- 0.75
dim2 <- c(6, 30)
n_setting2 <- c("n/2", "n", "n*2")
prev2 <- 0.2
model2 <- "ML"
pred_sel2 <- c("none", "<0.157")
noise_2 <- c("none", "half")

## all combinations:
s2 <- expand.grid(AUC = AUC2, 
                  dim = dim2,
                  n_setting = n_setting2, 
                  noise = noise_2, 
                  prev = prev2, 
                  model = model2, 
                  pred_selection = pred_sel2, 
                  KEEP.OUT.ATTRS = F)


#################
## Expected R2 ##
#################

s2$R2 <- R2[2]

########################
## Actual sample size ##
########################

n_2 <- c("dim_6" = pmsampsize(type = "b", parameters = 6, prevalence = prev2, rsquared = R2[2])$sample_size,
         "dim_30" = pmsampsize(type = "b", parameters = 30, prevalence = prev2, rsquared = R2[2])$sample_size
         )
         

s2 <- s2 %>%
  mutate(
    n = ceiling(case_when(
      n_setting == "n/2" & dim == 6 ~ n_2["dim_6"]/2,
      n_setting == "n"   & dim == 6 ~ n_2["dim_6"],
      n_setting == "n*2" & dim == 6 ~ n_2["dim_6"]*2,
      n_setting == "n/2" & dim == 30 ~ n_2["dim_30"]/2,
      n_setting == "n"   & dim == 30 ~ n_2["dim_30"],
      n_setting == "n*2" & dim == 30 ~ n_2["dim_30"]*2,
      TRUE ~ NA_real_
    )
    )
  )


# Remove attributes of "n"
attr(s2$n, "ATT") <- NULL

####################
### Coefficients ###
####################
### Using the testing site optimization script
## The following values for the parameters were 
## found during optimization:

# 6 candidate predictors, with half noise:
beta_6_half <- c(-1.671603, 0.27602491)
# 30 candidate predictors, with half noise:
beta_30_half <- c(-1.623487, 0.07781071)


# 6 candidate predictors, with no noise:
beta_6_none <- c(-1.670588, 0.20894297)
# 6 candidate predictors, with no noise:
beta_30_none <- c(-1.667019,  -0.05344094)

## Bind parameters to study 1 matrix
# par1 represents intercept
# par2 represents regression coefficients
s2 <- s2 %>% 
  mutate(par1 = case_when(dim == 6 & noise == "half" ~ beta_6_half[1],
                          dim == 30 & noise == "half" ~ beta_30_half[1],
                          dim == 6 & noise == "none" ~ beta_6_none[1],
                          dim == 30 & noise == "none" ~ beta_30_none[1],
                          TRUE ~ NA_real_
  )) %>%
  mutate(par2 = case_when(dim == 6 & noise == "half" ~ beta_6_half[2],
                          dim == 30 & noise == "half" ~ beta_30_half[2],
                          dim == 6 & noise == "none" ~ beta_6_none[2],
                          dim == 30 & noise == "none" ~ beta_30_none[2],
                          TRUE ~ NA_real_
  ))


s2$scenario <- c(paste0("Scenario_", 1:nrow(s2)))

### Save the study scenarios in the settings folder
write_rds(s2, file = study_2_settings)


#####################
###### Study 3 ######
#####################

## Data gen mechanism
AUC3 <- 0.75
dim3 <- 20
n_setting3 <- c("n/2", "n", "n*2")
noise_3 <- c("default")
prev3 <- 0.2
pred_sel3 <- c("none")
models3 <-
  c("ML", 
    "Firth",
    "Ridge",
    "Lasso",
    "CART",
    #"ANN",
    #"SVM",
    "RF")

## all combinations:
s3 <- expand.grid(AUC = AUC3,
                  dim = dim3,
                  n_setting = n_setting3,
                  noise = noise_3,
                  prev = prev3, 
                  model = models3,
                  pred_selection = pred_sel3,
                  KEEP.OUT.ATTRS = F)

#################
## Expected R2 ##
#################

s3$R2 <- R2[2]

########################
## Actual sample size ##
########################
n_3 <- pmsampsize(type = "b", parameters = 20, prevalence = 0.2, rsquared = R2[2])$sample_size


s3 <- s3 %>%
  mutate(
    n = ceiling(case_when(
      n_setting == "n/2" ~ n_3/2,
      n_setting == "n" ~ n_3,
      n_setting == "n*2" ~ n_3*2,
      TRUE ~ NA_real_
    )
    )
  )

####################
### Coefficients ###
####################
### Using the testing site optimization script
## The following values for the parameters were 
## found during optimization:
beta_20 <- c(-1.64975835, -0.07082455)

s3$par1 <- beta_20[1]
s3$par2 <- beta_20[2]

s3$scenario <- c(paste0("Scenario_", 1:nrow(s3)))

### Save the study scenarios in the settings folder
write_rds(s3, file = study_3_settings)

# Remove everything except the three study matrices
#rm(list=ls()[! ls() %in% c("s1","s2", "s3")])


########################
## All studies in one ##
########################

s1$study <- "Study_1"
s2$study <- "Study_2"
s3$study <- "Study_3"

studies <- rbind(s1,s2,s3)
saveRDS(studies, paste0(setting_path, "studies.RDS"))



####################################
####################################
####################################
## END SCRIPT

