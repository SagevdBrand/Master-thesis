################
## Setting up ##
################
source("./src/setup.R")
scenarios <- readRDS(paste0(setting_path, "studies.RDS"))

## Loading the separate datafiles and binding them togetgher saving it for easy reference 
# df_all <- list.files(path = estimands_path, pattern = "*.Rds", full.names = T) %>%
#   map_dfr(readRDS)
# 
# Saving it for easy reference 
# saveRDS(df_all, file = paste(estimands_path, "all_estimands_batch_1.RDS"))

## Loading the above created file
df_all <- readRDS(paste(estimands_path, "all_estimands_batch_1.RDS"))
df_all$error_info[df_all$error_info == "NA"] <- NA

## Make sure that all estimands and estimand_se are numeric
df_all[, estimands_names] <- lapply(estimands_names, function(x) as.numeric(df_all[[x]]))
df_all[, estimands_se_names] <- lapply(estimands_se_names, function(x) as.numeric(df_all[[x]]))


##########################################################
############## Check which results are present ###########
##########################################################

# How many iterations have been analyzed?
counts <- df_all %>% group_by(study, scenario) %>% summarise(count = n())
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

external_results <- df_all %>% filter(approach == "External")
iv_results <- df_all %>% filter(approach != "External", approach != "Apparent")
app_results <- df_all %>% filter(approach == "Apparent")







#################################
######## Visualization ##########
#################################
# Turn the results into long format:
df_all_long <- df_all %>% pivot_longer(., cols = estimands_names, names_to = "estimand", values_to = "value")
df_all_long$value <- as.numeric(df_all_long$value)

# Specify color palette:
colors <- c("#02c39a", #Mountain meadow
            "#f28482", #Light coral
            "#84a59d", # Morning blue
            "#f6bd60", # Maximum yellow red
            "white",   # White
            "#457b9d", # Celadon blue
            "#f4a261", # Sandy brown
            "#e63946" # imperial red
)



# Results of study 1:
df_s1_long<- df_all_long %>% filter(study == "Study_1")


s1_s1_4_7 <- df_s1_long %>% filter(scenario %in% c("Scenario_1", "Scenario_4", "Scenario_7")) %>% dplyr::select(c(estimand, study, scenario,model, approach, value, prev)) 
s1_s1_4_7 %>% ggplot(data = ., mapping = aes(x = as.factor(prev), y = value, fill = approach)) +
  geom_boxplot(position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
  theme_set(theme_bw(base_size = 11)) +
  scale_fill_manual(values = colors) +
  labs(y = " ",
       x = "Prevalence",
       fill = "Validation approach") +
  facet_wrap(~estimand, nrow = 3, scales = "free") +
  theme(legend.position = c(0.8, 0.12))





