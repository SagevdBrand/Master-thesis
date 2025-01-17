#################################
######## Visualization ##########
#################################
################
## Setting up ##
################
source("./src/setup.R")
scenarios <- readRDS(paste0(setting_path, "studies.RDS"))

## Loading the separate datafiles and binding them togetgher saving it for easy reference 
#df_all <- list.files(path = estimands_path, pattern = "*.Rds", full.names = T) %>%
#     map_dfr(readRDS)

# First save it, just to be sure!
#saveRDS(df_all, file = paste0(estimands_general_path, "all_estimands_part_2.RDS"))
df_all <- readRDS(paste0(estimands_general_path, "all_estimands_part_2.RDS"))

### Pre-processing ###
# Convert all "NA" to actual NA 
df_all$error_info[df_all$error_info == "NA"] <- NA

# Make sure that all estimands and estimand_se are numeric
df_all[, estimands_names] <- lapply(estimands_names, function(x) as.numeric(df_all[[x]]))
df_all[, estimands_se_names] <- lapply(estimands_se_names, function(x) as.numeric(df_all[[x]]))


##########################################################
############## Check which results are present ###########
##########################################################

# How many iterations have been analyzed?
counts <- df_all %>% group_by(study, scenario) %>% summarise(count = n())
# # For 1000 iterations there should be 8000 results per scenario!
ind <- seq(1000)
 
# # If it is not 8000, check what is missing by:
df_s3_s18 <- df_all %>% filter(study == "Study_3", scenario == "Scenario_18")
which(!ind %in% df_s3_s18$iteration)

##############################
## ESTIMAND SPECIFIC ERRORS ##
##############################

#########
## AUC ##
#########

summary(df_all$auc)
auc_na <- df_all %>% filter(is.na(auc))

####################
## Slope problems ##
####################
summary(df_all$calib_slope)
slope_problems <- df_all %>% filter(is.na(calib_slope)|calib_slope<0)
negative_slope <- slope_problems %>% filter

big_slope <-   slope_problems %>% filter(calib_slope > 10)
na_slope <- slope_problems %>% filter(is.na(calib_slope))

########################
## intercept problems ##
########################
summary(df_all$calib_int)
intercept_problems <- df_all %>% filter(is.na(calib_int)|calib_int > 5|calib_int< -5)
very_negative_intercept <- intercept_problems %>% filter(calib_int < -5) # All CART S3_s13!
very_positive_intercept <-   intercept_problems %>% filter(calib_int > 5)
na_intercept <- intercept_problems %>% filter(is.na(calib_int))

##################
## ECI problems ##
##################
summary(df_all$eci)
eci_problems <- df_all %>% filter(is.na(eci)|is.infinite(eci)| eci>1| eci<0)
infinite_eci <- eci_problems %>% filter(eci == Inf)
negative_eci <- eci_problems %>% filter(eci <0)
bigger_1_eci <- eci_problems %>% filter(eci > 1 & eci != Inf)
eci_na <- eci_problems %>% filter(is.na(eci))


###################
## R2CS problems ##
###################
summary(df_all$R2_CS)
R2_CS_problems <- df_all %>% filter(is.na(R2_CS)|R2_CS>1| R2_CS<0)
negative_R2_CS <- R2_CS_problems %>% filter(R2_CS <0)
na_R2_CS <-   R2_CS_problems %>% filter(is.na(R2_CS))

#####################
## R2Tjur problems ##
#####################

summary(df_all$Tjur)
R2_Tjur_problems <- df_all %>% filter(is.na(Tjur)|Tjur>1| Tjur<0)
negative_R2_Tjur <- R2_Tjur_problems %>% filter(Tjur <0)
na_R2_Tjur <-   R2_Tjur_problems %>% filter(is.na(Tjur))

###########
## rMSPE ##
###########

summary(df_all$rmspe)
na_rmspe <-  df_all %>% filter(is.na(rmspe))

##########
## MAPE ##
##########

summary(df_all$mape)
na_mape <-  df_all %>% filter(is.na(mape))

##################################
## What are the error messages? ##
##################################

### Count all the errors ###
error_counting <- df_all %>% group_by(study,scenario) %>% count(error_info)

# ## Specific errors ##
no_pred <- "No predictors selected -> no calibration slope"
separation <- "Data separation might have occured"
separation_warning <- "simpleWarning: glm.fit: fitted probabilities numerically"
few_events <- "Too few (non-)events for tuning -> LOOCV"
max_probs <- "probabilities of 0 or 1 occured"
loess_warning <- "simpleWarning in simpleLoess"
no_events <- "Error: No events sampled"
no_events_folds <- "No events sampled in folds "
no_events_training_samp <- " No events sampled in training sample"



## Count how many were present per error and scenario ##
per_scenario <- df_all %>% group_by(study,scenario) %>% summarise(no_pred_count = sum(str_count(error_info, no_pred), na.rm = T),
                                                                  no_pred_app = n_distinct(which(is.na(auc))), #This were the only times when this error occurred
                                                                  negative_slopes = n_distinct(which(calib_slope <0)),
                                                                  large_slope = n_distinct(which(calib_slope >10)),
                                                                  NA_slopes = n_distinct(which(is.na(calib_slope))),
                                                                  very_negative_intercepts = n_distinct(which(calib_int < -5)),
                                                                  very_positive_intercepts = n_distinct(which(calib_int > 5)),
                                                                  negative_R2_CS = n_distinct(which(is.na(R2_CS))),
                                                                  negative_R2_Tjur = n_distinct(which(is.na(Tjur))),
                                                                  infinite_eci = n_distinct(c(which(is.infinite(eci)), which(eci < 0))),
                                                                  eci_bigger_1 = n_distinct(which(eci > 1)),
                                                                  NA_ECI = n_distinct(which(is.na(eci))),
                                                                  max_probs_count = sum(str_count(error_info, max_probs), na.rm = T),
                                                                  sep_count = sum(str_count(error_info, paste0(separation, " | ",separation_warning)), na.rm = T),
                                                                  few_events_count = sum(str_count(error_info, few_events), na.rm = T),
                                                                  loess_warning_count = sum(str_count(error_info, loess_warning), na.rm = T),
                                                                  no_events_count = sum(str_count(error_info, no_events), na.rm = T),
                                                                  no_events_folds_count = sum(str_count(error_info, no_events_folds), na.rm = T),
                                                                  no_events_training_samp_count = sum(str_count(error_info, no_events_training_samp), na.rm = T),
                                                                  .groups = "keep")

per_scenario$scenario <- as.numeric(gsub("Scenario_", "", per_scenario$scenario))
per_scenario <- per_scenario %>% arrange(study, scenario)
colnames(per_scenario) <- c("Study",
"Scenario",
"No predictors selected",
"NA .632+ bootstrap results",
"Negative Slopes",
"Slopes > 10",
"NA slopes",
"Intercepts < -5",
"Intercepts > 5",
"Negative R2 CS",
"Negative R2 Tjur",
"(-)Infinite ECI",
"ECI > 1",
"NA ECI",
"Prob. of 0 or 1",
"Separation",
"< 8 events",
"eci: LOESS warning",
"No events",
"No events in fold",
"No events in bootstrap sample"
)

#per_scenario[57,14] <- 55-52 # 29 out of 32 errors were already covered by the NA .632+ bootstrap results!
 
#saveRDS(per_scenario, file = paste0(errors_path,"all_errors_per_scenario.Rds"))
 
#all_together <- as.matrix(colSums(per_scenario[,c(3:21)]))
 
# saveRDS(all_together, file = paste0(errors_path,"all_errors_together.Rds"))

per_scenario <- readRDS(paste0(errors_path,"all_errors_per_scenario.Rds"))
all_together <- readRDS(paste0(errors_path,"all_errors_together.Rds"))

kable(per_scenario, format = "markdown")

# How many models in total?
num_models_total <- (1 + 5 + 10 + 100 + 500) * 60 * 1000 # Not accounting for hyperparameter tuning
num_models_pred_sel <- (1 + 5 + 10 + 100 + 500) * 30 * 1000
num_models_tuning <- (1 + 5 + 10 + 100 + 500) * 12 * 1000 # Not accounting for hyperparameter tuning
num_models_cv <- (5 + 10 + 100) * 60 * 1000
num_models_boot <- 500*60*1000
num_models_tree_based <- (1 + 5 + 10 + 100 + 500) * 6 * 1000
num_datasets <- 60*1000

all_together <- t(all_together)
percentage_errors <- c(all_together[,1]/num_models_pred_sel *100,
                       all_together[,13]/num_models_tree_based *100,
                       all_together[,14]/num_models_total *100,
                       all_together[,15]/num_models_tuning *100,
                       all_together[,16]/num_models_total *100,
                       all_together[,17]/num_datasets *100,
                       all_together[,18]/num_models_cv *100,
                       all_together[,19]/num_models_boot *100
                         )
percentage_errors <- lapply(percentage_errors, round, 4)  

#########################################################
# For those situations were no predictors were selected #
# #########################################################
## do this first to have a check later:
# Find the max values for calibration slope for each scenario
max_slopes <- df_all %>% group_by(study, scenario) %>% summarise(max_slope = max(calib_slope, na.rm = T))
# Which rows have the relevant warning? "No predictors selected -> no calibration slope"
no_pred_ind <- which(str_detect(df_all$error_info, "No predictors selected -> no calibration slope"))


## For all those results which have the error of having no calibration slope, replace them with the
## Highest calibration slope within the scenario.
df_all <- df_all %>% group_by(study, scenario) %>%
  mutate(calib_slope =
           case_when(
             str_detect(error_info,
                        "No predictors selected -> no calibration slope") == TRUE ~ max(calib_slope, na.rm = T),
             TRUE ~ calib_slope
           ),
         eci = case_when(eci > 1 ~ 1,
                         is.infinite(eci) ~ 1,
                         eci < 0 ~ 1, 
                         TRUE ~ eci)
         )

# Check
# df_all$calib_slope[no_pred_ind]
saveRDS(df_all, file = paste0(estimands_general_path, "all_estimands_part_2.RDS"))

