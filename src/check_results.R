################
## Setting up ##
################
source("./src/setup.R")
scenarios <- readRDS(paste0(setting_path, "studies.RDS"))

## Loading the file, created after error handling. 
df_all <- readRDS(paste0(estimands_general_path, "all_estimands_batch_6.RDS"))

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

saveRDS(m_perform_results, file = paste0(performance_general_path, "all_pm_batch_5.RDS"))


