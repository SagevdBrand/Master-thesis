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
#  map_dfr(readRDS)

df_all <- readRDS(paste0(estimands_general_path, "all_estimands_batch_6.RDS"))

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
#For 500 iterations there should be 4000 results per scenario!
ind <- seq(500)
# 
# If it is not 4000, check what is missing by:
df_s3_s18 <- df_all %>% filter(study == "Study_3", scenario == "Scenario_18")
which(!ind %in% df_s3_s18$iteration)



##################################
## What are the error messages? ##
################################## 

### Count all the errors ###
error_counting <- df_all %>% group_by(study,scenario) %>% count(error_info)

## Specific errors ##
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
# per_scenario <- df_all %>% group_by(study,scenario) %>% summarise(no_pred_count = sum(str_count(error_info, no_pred), na.rm = T),
#                                                           max_probs_count = sum(str_count(error_info, max_probs), na.rm = T),
#                                                           sep_count = sum(str_count(error_info, paste0(separation, " | ",separation_warning)), na.rm = T),
#                                                           few_events_count = sum(str_count(error_info, few_events), na.rm = T),
#                                                           loess_warning_count = sum(str_count(error_info, loess_warning), na.rm = T),
#                                                           no_events_count = sum(str_count(error_info, no_events), na.rm = T),
#                                                           no_events_folds_count = sum(str_count(error_info, no_events_folds), na.rm = T),
#                                                           no_events_training_samp_count = sum(str_count(error_info, no_events_training_samp), na.rm = T),
#                                                           .groups = "keep")

# per_scenario$scenario <- as.numeric(gsub("Scenario_", "", per_scenario$scenario))
# per_scenario <- per_scenario %>% arrange(study, scenario)
# saveRDS(per_scenario, file = paste0(performance_general_path,"all_errors_per_scenario.Rds"))
# 
# all_together <- df_all %>% summarise(no_pred_count = sum(str_count(.$error_info, no_pred), na.rm = T),
#                                                           max_probs_count = sum(str_count(.$error_info, max_probs), na.rm = T),
#                                                           sep_count = sum(str_count(.$error_info, paste0(separation, " | ",separation_warning)), na.rm = T),
#                                                           few_events_count = sum(str_count(.$error_info, few_events), na.rm = T),
#                                                           loess_warning_count = sum(str_count(.$error_info, loess_warning), na.rm = T),
#                                                           no_events_count = sum(str_count(.$error_info, no_events), na.rm = T),
#                                                           no_events_folds_count = sum(str_count(.$error_info, no_events_folds), na.rm = T),
#                                                           no_events_training_samp_count = sum(str_count(.$error_info, no_events_training_samp), na.rm = T),
#                                                           )
# 
# all_together <- all_together[1, c(3:10)]
# saveRDS(all_together, file = paste0(performance_general_path,"all_errors_together.Rds"))
per_scenario <- readRDS(paste0(performance_general_path,"all_errors_per_scenario.Rds"))
all_together <- readRDS(paste0(performance_general_path,"all_errors_together.Rds"))

t(all_together)

# How many models in total?
num_models_total <- (1 + 5 + 10 + 100 + 500) * 60 * 500
num_models_pred_sel <- (1 + 5 + 10 + 100 + 500) * 30 * 500
num_models_tuning <- (1 + 5 + 10 + 100 + 500) * 12 * 500
num_models_cv <- (5 + 10 + 100) * 60 * 500
num_models_boot <- 500*60*500
num_models_tree_based <- (1 + 5 + 10 + 100 + 500) * 6 * 500
num_datasets <- 60*500


percentage_errors <- c(all_together[,1]/num_models_pred_sel,
                       all_together[,2]/num_models_tree_based,
                       all_together[,3]/num_models_total,
                       all_together[,4]/num_models_tuning,
                       all_together[,5]/num_models_total,
                       all_together[,6]/num_datasets,
                       all_together[,7]/num_models_cv,
                       all_together[,8]/num_models_boot
                         )
percentage_errors <- lapply(percentage_errors, round, 4)  

#########################################################
# For those situations were no predictors were selected #
# #########################################################
# ## do this first to have a check later:
# # Find the max values for calibration slope for each scenario
# max_slopes <- df_all %>% group_by(study, scenario) %>% summarise(max_slope = max(calib_slope, na.rm = T))
# # Which rows have the relevant warning? "No predictors selected -> no calibration slope"
# no_pred_ind <- which(str_detect(df_all$error_info, "No predictors selected -> no calibration slope"))
# 
# 
# ## For all those results which have the error of having no calibration slope, replace them with the
# ## Highest calibration slope within the scenario.
# df_all <- df_all %>% group_by(study, scenario) %>%
#   mutate(calib_slope = 
#            case_when(str_detect(error_info, 
#                                 "No predictors selected -> no calibration slope") == TRUE ~ max(calib_slope, na.rm = T),
#                      TRUE ~ calib_slope))
# 
# # s3_scen_14 <- df_all %>% filter(study == "Study_3", scenario == "Scenario_14")
# # test <- which(is.na(s3_scen_14$calib_slope))
# # test_df <- s3_scen_14[test, ]
# 
# # Check

# df_all$calib_slope[no_pred_ind]
saveRDS(df_all, file = paste0(estimands_general_path, "all_estimands_batch_6.RDS"))

