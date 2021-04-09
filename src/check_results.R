################
## Setting up ##
################
source("./src/setup.R")
scenarios <- readRDS(paste0(setting_path, "studies.RDS"))

## Loading the separate datafiles and binding them togetgher saving it for easy reference 
#df_all <- list.files(path = estimands_path, pattern = "*.Rds", full.names = T) %>%
#   map_dfr(readRDS)
# 
# Saving it for easy reference 
#saveRDS(df_all, file = paste0(estimands_general_path, "all_estimands_batch_1_2_3_4.RDS"))

## Loading the above created file
df_all <- readRDS(paste0(estimands_general_path, "all_estimands_batch_1_2_3_4.RDS"))
df_all$error_info[df_all$error_info == "NA"] <- NA

## Make sure that all estimands and estimand_se are numeric
df_all[, estimands_names] <- lapply(estimands_names, function(x) as.numeric(df_all[[x]]))
df_all[, estimands_se_names] <- lapply(estimands_se_names, function(x) as.numeric(df_all[[x]]))

##########################################################
############## Check which results are present ###########
##########################################################

# How many iterations have been analyzed?
counts <- df_all %>% group_by(study, scenario) %>% summarise(count = n())
# For 500 iterations there should be 4000 results per scenario!
# ind <- seq(500)
# 
# # # If it is not 4000, check what is missing by:
# df_s2_s24 <- df_all %>% filter(study == "Study_2", scenario == "Scenario_24")
# which(!ind %in% df_s2_s24$iteration)

error_counting <- df_all %>% group_by(study,scenario) %>% count(error_info)

##########################################################
# For those situations were no predictors were selected #
##########################################################
## do this first to have a check later:
# Find the max values for calibration slope for each scenario
max_slopes <- df_all %>% group_by(study, scenario) %>% summarise(max_slope = max(calib_slope))
# Which rows have the relevant warning? "No predictors selected -> no calibration slope"
no_pred_ind <- str_which(df_all$error_info, "No predictors selected -> no calibration slope")

df_all <- df_all %>% group_by(study, scenario) %>%
  mutate(calib_slope = 
           case_when(str_detect(error_info, 
                                "No predictors selected -> no calibration slope") == TRUE ~ max(calib_slope),
                     TRUE ~ calib_slope))

# Check
df_all$calib_slope[no_pred_ind]

# test <- which(is.na(df_all$calib_slope))
# test_df <- df_all[test, ]
# union(which(is.na(df_all$calib_slope)), no_pred_ind)

#####################################
####### Performance measures ########
#####################################

# The simulation results will be summarized by the average difference 
# between external validation and internal validation and the root mean 
# squared difference (RMSD). Empirical standard errors of the mean of the difference
# and RMSD will also be calculated. For the calibration slope, median and IQR are 
# calculated instead of mean and SE, due to expected skewness in smaller datasets 
# (van Smeden et al., SMMR, 2019). 

# Split the external from the rest so it can be merged later based on: 
# study
# scenario
# iteration

external_results <- df_all %>% filter(approach == "External")
# Only keep study, scenario, iteration and value (change value to "external value")
external_results <- external_results[, (colnames(external_results) %in% c("study", "scenario", "iteration", "approach", estimands_names))]
colnames(external_results) <- c("iteration", "study", "scenario", "approach", paste0(estimands_names, "_ext"))

# Use everything but the external values as a base
# for the column of "ext_value" to be merged on:
iv_app_results <- df_all%>% filter(approach != "External")

perform_results <- iv_app_results %>% group_by(approach) %>% merge(., external_results, by = c("study", "scenario", "iteration"))

### Function to obtain the performance measures

performance_measures <- function(ext, int, name){
  
  estimand_results <- perform_results[, (colnames(perform_results) %in% c("study", "scenario", "approach.x", ext, int))]
  
  if (name == "calib_slope"){
    
    results <- estimand_results %>% group_by(study, scenario, approach.x) %>% summarise(md = median(get(ext) - get(int) , na.rm = T),
                                                                            md_var = IQR(get(ext) - get(int), na.rm = T),
                                                                            rmsd = sqrt(median((get(ext) - get(int))^2, na.rm = T)),
                                                                            rmsd_var = sqrt(IQR((get(ext) - get(int))^2, na.rm = T)))
    
    } else {
      
      results <- estimand_results %>% group_by(study, scenario, approach.x) %>% summarise(md = mean(get(ext) - get(int) , na.rm = T),
                                                                     md_var = sd(get(ext) - get(int))/sqrt(nrow(.)),
                                                                     rmsd = sqrt(mean((get(ext) - get(int))^2, na.rm = T)),
                                                                     rmsd_var = sqrt(sd((get(ext) - get(int))^2/sqrt(nrow(.)))))
    }
  results$scenario <- as.numeric(gsub("Scenario_", "", results$scenario))
  results$estimand <- name
  colnames(results) <-c("study", "scenario", "approach", "md", "md_dist", "rmsd", "rmsd_dist", "estimand")  
  results <- results %>% arrange(study, scenario)
  return(results)

}

auc_p <- performance_measures(ext = "auc_ext", int = "auc", name = "auc")
calib_slope_p <- performance_measures(ext = "calib_slope_ext", int = "calib_slope", name = "calib_slope")
calib_int_p <- performance_measures(ext = "calib_int_ext", int = "calib_int", name = "calib_int")
tjur_p <- performance_measures(ext = "Tjur_ext", int = "Tjur", name = "Tjur")
R2_CS_p <- performance_measures(ext = "R2_CS_ext", int = "R2_CS", name = "R2_CS")
eci_p <- performance_measures(ext = "eci_ext", int = "eci", name = "eci")
mape_p <- performance_measures(ext = "mape_ext", int = "mape", name = "mape")
rmspe_p <- performance_measures(ext = "rmspe_ext", int = "rmspe", name = "rmspe")

m_perform_results <-  rbind(auc_p, calib_slope_p, calib_int_p, tjur_p, R2_CS_p,eci_p, mape_p, rmspe_p)



