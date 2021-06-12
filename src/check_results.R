################
## Setting up ##
################
source("./src/setup.R")
scenarios <- readRDS(paste0(setting_path, "studies.RDS"))

## Loading the file, created after error handling. 
df_all <- readRDS(paste0(estimands_general_path, "all_estimands_part_2.RDS"))

scenarios_n_setting <- scenarios[, (colnames(scenarios) %in% c("study", "scenario", "n_setting", "noise"))]
df_all <- merge(df_all, scenarios_n_setting, by = c("study", "scenario"))

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
colnames(external_results) <- c("study", "scenario", "iteration", "approach", paste0(estimands_names, "_ext"))


# Get the performance results in a format fit for 
# the nested loop plot. Here, we can then use the external validation as a reference.
perform_results <- df_all %>% group_by(approach) %>% merge(., external_results, by = c("study", "scenario", "iteration"))
saveRDS(perform_results, paste0(performance_general_path, "internal_external_by_columns.RDS"))

# For the other plots:
# Use everything but the external values as a base
# for the column of "ext_value" to be merged on:
iv_app_results <- df_all%>% filter(approach != "External")
perform_results <- iv_app_results %>% group_by(approach) %>% merge(., external_results, by = c("study", "scenario", "iteration"))


###############################
## Obtain median differences ##
###############################


performance_measures <- function(ext, int, name){
  
  estimand_results <- perform_results[, (colnames(perform_results) %in% c("study", "scenario", "approach.x", ext, int))]
  
  if (name == "calib_slope"){
    
    results <- estimand_results %>% group_by(study, scenario, approach.x) %>% summarise(md = median(log(get(int)) - log(get(ext)), na.rm = T),
                                                                                        md_q1 = summary(ifelse(get(int) == 0,
                                                                                                               log(0.01), 
                                                                                                               log(get(int))) - 
                                                                                                          ifelse(get(ext) == 0,
                                                                                                                 log(0.01),
                                                                                                                 log(get(ext))))[2],
                                                                                        md_q3 = summary(ifelse(get(int) == 0,
                                                                                                               log(0.01), 
                                                                                                               log(get(int))) - 
                                                                                                          ifelse(get(ext) == 0,
                                                                                                                 log(0.01),
                                                                                                                 log(get(ext))))[5]) #,
    #rmsd = sqrt(median((get(ext) - get(int))^2, na.rm = T)),
    #rmsd_var = sqrt(IQR((get(ext) - get(int))^2, na.rm = T)))
    results$scenario <- as.numeric(gsub("Scenario_", "", results$scenario))
    results$estimand <- name
    colnames(results) <-c("study", "scenario", "approach", "md", "md_q1", "md_q3", "estimand") 
    
  } else {
    
    results <- estimand_results %>% group_by(study, scenario, approach.x) %>% summarise(md = median(get(int) - get(ext), na.rm = T),
                                                                                        md_var = sd(get(int) - get(ext)))#,
    #rmsd = sqrt(mean((get(ext) - get(int))^2, na.rm = T)),
    #rmsd_var = sqrt(sd((get(ext) - get(int))^2)))
    
    
    
    results$scenario <- as.numeric(gsub("Scenario_", "", results$scenario))
    results$estimand <- name
    colnames(results) <-c("study", "scenario", "approach", "md", "md_dist", "estimand") 
  }
  
  
  
  results <- results %>% arrange(study, scenario)
  return(results)
  
}


#################################
### OBTAIN THE ACTUAL RESULTS ###
#################################

auc_p <- performance_measures(ext = "auc_ext", int = "auc", name = "auc")
calib_slope_p <- performance_measures(ext = "calib_slope_ext", int = "calib_slope", name = "calib_slope")
calib_int_p <- performance_measures(ext = "calib_int_ext", int = "calib_int", name = "calib_int")
tjur_p <- performance_measures(ext = "Tjur_ext", int = "Tjur", name = "Tjur")
R2_CS_p <- performance_measures(ext = "R2_CS_ext", int = "R2_CS", name = "R2_CS")
eci_p <- performance_measures(ext = "eci_ext", int = "eci", name = "eci")
mape_p <- performance_measures(ext = "mape_ext", int = "mape", name = "mape")
rmspe_p <- performance_measures(ext = "rmspe_ext", int = "rmspe", name = "rmspe")

m_perform_results <-  rbind(auc_p, calib_slope_p, calib_int_p, tjur_p, R2_CS_p,eci_p, mape_p, rmspe_p)

saveRDS(m_perform_results, file = paste0(performance_general_path, "all_pm_part_2.RDS"))


##########################
## CURIOUS CASE OF CART ##
##########################

# Comparing the results of ML with CART
ML_CART_results <- 
  df_all %>% 
  filter(study == "Study_3" & scenario %in% c("Scenario_3", "Scenario_15")) %>% 
  group_by(scenario, approach, study) %>% 
  summarise(auc_mean = mean(auc, na.rm = T),
            auc_median = median(auc, na.rm = T),
            auc_sd = sd(auc, na.rm = T),
            slope_median = median(calib_slope, na.rm = T),
            slope_lower = summary(calib_slope)[2],
            slope_upper = summary(calib_slope)[5],
            intercept_mean = mean(calib_int, na.rm = T),
            intercept_median = median(calib_int, na.rm = T),
            intercept_sd = sd(calib_int, na.rm = T),
            Tjur_mean = mean(Tjur, na.rm = T),
            Tjur_median = median(Tjur, na.rm = T),
            Tjur_sd = sd(Tjur, na.rm = T),
            R2_CS_mean = mean(R2_CS, na.rm = T),
            R2_CS_median = median(R2_CS, na.rm = T),
            R2_CS_sd = sd(R2_CS, na.rm = T),
            eci_mean = mean(eci, na.rm = T),
            eci_median = median(eci, na.rm = T),
            eci_sd = sd(eci, na.rm = T),
            mape_mean = mean(mape, na.rm = T),
            mape_median = median(mape, na.rm = T),
            mape_sd = sd(mape, na.rm = T),
            rmspe_mean = mean(rmspe, na.rm = T),
            rmspe_median = median(rmspe, na.rm = T),
            rmspe_sd = sd(rmspe, na.rm = T)
            )


ML_CART_table <- ML_CART_results %>% group_by(scenario) %>%
  dplyr::select(approach, scenario, study, auc_mean, auc_sd, slope_median, slope_lower, slope_upper, R2_CS_mean, R2_CS_sd, eci_mean, eci_sd, rmspe_mean, rmspe_sd) %>% 
  mutate(model = case_when(scenario == "Scenario_15" ~ "CART",
                           scenario == "Scenario_3" ~ "ML"),
         approach = fct_relevel(approach,
                                "5 fold cross-validation",
                                "10 fold cross-validation",
                                "10x10 fold cross-validation",
                                "Harrell's bootstrap",
                                ".632 bootstrap",
                                ".632+ bootstrap",
                                "Apparent")) %>%
  dplyr::select(model, everything()) %>%
  arrange(scenario, approach, model)

ML_CART_table <- ML_CART_table[,-c(3,4)] 
ML_CART_table[,c(3:13)] <- round(ML_CART_table[, c(3:13)], 3)

ML_CART_table$auc <- paste0(ML_CART_table$auc_mean," (", ML_CART_table$auc_sd,")")
ML_CART_table$slope <- paste0(ML_CART_table$slope_median," (", ML_CART_table$slope_lower," - ", ML_CART_table$slope_upper,")")
ML_CART_table$R2_CS <- paste0(ML_CART_table$R2_CS_mean," (", ML_CART_table$R2_CS_sd,")")
ML_CART_table$eci <- paste0(ML_CART_table$eci_mean," (", ML_CART_table$eci_sd,")")
ML_CART_table$rmspe <- paste0(ML_CART_table$rmspe_mean," (", ML_CART_table$rmspe_sd,")")

kable(ML_CART_table[, c(2, 14:18)], 
      digits = 3, 
      format = "latex",
      row.names = NA, 
      booktabs = TRUE) %>% 
  pack_rows(index = c("CART" = 8, "ML" = 8))


