#################################
######## Visualization ##########
#################################
################
## Setting up ##
################
source("./src/setup.R")
scenarios <- readRDS(paste0(setting_path, "studies.RDS"))

## Loading the above estimand files:
df_all <- readRDS(paste0(estimands_general_path, "all_estimands_batch_5.RDS"))

################################
## Pre-processing of the data ##
################################

# As I missed this when programming the study, need to add it here:
# Add n_setting to the columns:
scenarios_n_setting <- scenarios[, (colnames(scenarios) %in% c("study", "scenario", "n_setting", "noise"))]
df_all_2 <- merge(df_all, scenarios_n_setting, by = c("study", "scenario"))

# Then for all of the next steps, a long format is more useful:
# Turn the results into long format:
df_all_long <- df_all_2 %>% pivot_longer(., cols = all_of(estimands_names), names_to = "estimand", values_to = "value")

############
## Colors ##
############

# Specify color palette:
colors <- c("#02c39a", #Mountain meadow
            "#f28482", #Light coral
            "#84a59d", # Morning blue
            "#f6bd60", # Maximum yellow red
            "#457b9d", # Celadon blue
            "#f4a261", # Sandy brown
            "white",   # White
            "#e63946" # imperial red
)

############################
############################
########## STUDY 1 #########
############################
############################

# Results of study 1:
df_s1_long<- df_all_long %>% filter(study == "Study_1")
df_s1_long$scenario <- as.numeric(gsub("Scenario_", "", df_s1_long$scenario))
df_s1_long <- df_s1_long %>% mutate(value = case_when(value > 5 ~ 5, TRUE ~ value))

# For each 3 scenarios a plot can be created
start_scen <- seq(1,18, by = 3)
end_scen <- seq(3, 18, by = 3)

for (i in 1:6){
preval <- df_s1_long %>% filter(scenario == start_scen[i]) %>% dplyr::select(prev) %>% .[1,1]
model <- df_s1_long %>% filter(scenario == start_scen[i]) %>% dplyr::select(model) %>% .[1,1]
predictor_selection <- df_s1_long %>% filter(scenario == start_scen[i]) %>% dplyr::select(pred_selection) %>% .[1,1]
  
p1 <- df_s1_long %>% 
    filter(scenario %in% c(start_scen[i]:end_scen[i])) %>% 
    ggplot(data = ., mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
    geom_boxplot(position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
    theme_set(theme_bw(base_size = 11)) +
    scale_fill_manual(values = colors) +
    labs(y = " ",
         x = "Prevalence",
         fill = "Validation approach",
         title = paste("Estimands for study 1: Scenario", start_scen[i], "until scenario", end_scen[i] ),
         subtitle = paste0("Prevalence = ", preval,", ",
                          "Model = ", model$model, ", ",
                          "Predictor Selection = ", predictor_selection$pred_selection)) +
    facet_wrap(~estimand, nrow = 3, scales = "free") +
    theme(legend.position = c(0.8, 0.12))
  
assign(paste0("Study_1_set_", i), p1)
}


############################
############################
########## STUDY 2 #########
############################
############################

# Results of study 2:
df_s2_long<- df_all_long %>% filter(study == "Study_2")
df_s2_long$scenario <- as.numeric(gsub("Scenario_", "", df_s2_long$scenario))


df_s2_long <- df_s2_long %>% mutate(value = case_when(value > 5 ~ 5, TRUE ~ value))

# For each 2 scenarios a plot can be created
start_scen <- seq(1,24, by = 2)
end_scen <- seq(2, 24, by = 2)

for (i in 1:12){
  noise <- df_s2_long %>% filter(scenario == start_scen[i]) %>% dplyr::select(noise) %>% .[1,1]
  model <- df_s2_long %>% filter(scenario == start_scen[i]) %>% dplyr::select(model) %>% .[1,1]
  predictor_selection <- df_s2_long %>% filter(scenario == start_scen[i]) %>% dplyr::select(pred_selection) %>% .[1,1]
  n_setting <- df_s2_long %>% filter(scenario == start_scen[i]) %>% dplyr::select(n_setting) %>% .[1,1]
  
  p1 <- df_s2_long %>% 
    filter(scenario %in% c(start_scen[i]:end_scen[i])) %>% 
    ggplot(data = ., mapping = aes(x = as.factor(dim), y = value, fill = approach)) +
    geom_boxplot(position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
    theme_set(theme_bw(base_size = 11)) +
    scale_fill_manual(values = colors) +
    labs(y = " ",
         x = "Number of Candidate Predictors",
         fill = "Validation approach",
         title = paste("Estimands for study 2: Scenario", start_scen[i], "and scenario", end_scen[i] ),
         subtitle = paste0("Noise = ", noise$noise,", ",
                           "model = ", model$model, ", ",
                           "predictor selection = ", predictor_selection$pred_selection,", ",
                           "n-setting = ", n_setting$n_setting)) +
    facet_wrap(~estimand, nrow = 3, scales = "free") +
    theme(legend.position = c(0.8, 0.12))
  
  assign(paste0("Study_2_set_", i), p1)
}


############################
############################
########## STUDY 3 #########
############################
############################

# Results of study 2:
df_s3_long<- df_all_long %>% filter(study == "Study_3")
df_s3_long$scenario <- as.numeric(gsub("Scenario_", "", df_s3_long$scenario))


df_s3_long <- df_s3_long %>% mutate(value = case_when(value > 5 ~ 5, TRUE ~ value))

# For each 3 scenarios a plot can be created
start_scen <- seq(1,18, by = 3)
end_scen <- seq(3, 18, by = 3)

for (i in 1:6){
  model <- df_s3_long %>% filter(scenario == start_scen[i]) %>% dplyr::select(model) %>% .[1,1]
  preval <- df_s3_long %>% filter(scenario == start_scen[i]) %>% dplyr::select(prev) %>% .[1,1]
  predictor_selection <- df_s3_long %>% filter(scenario == start_scen[i]) %>% dplyr::select(pred_selection) %>% .[1,1]
  n_setting <- df_s3_long %>% filter(scenario == start_scen[i]) %>% dplyr::select(n_setting) %>% .[1,1]
  
  p1 <- df_s3_long %>% 
    filter(scenario %in% c(start_scen[i]:end_scen[i])) %>% 
    ggplot(data = ., mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
    geom_boxplot(position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
    theme_set(theme_bw(base_size = 11)) +
    scale_fill_manual(values = colors) +
    labs(y = " ",
         x = "Number of Candidate Predictors",
         fill = "Validation approach",
         title = paste("Estimands for study 3: Scenario", start_scen[i], "until scenario", end_scen[i] ),
         subtitle = paste0("Model = ", model$model
                           )) +
    facet_wrap(~estimand, nrow = 3, scales = "free") +
    theme(legend.position = c(0.8, 0.12))
  
  assign(paste0("Study_3_set_", i), p1)
}

#############################
## VIEW THE VISUALIZATIONS ##
#############################

Study_1_set_1
Study_1_set_2
Study_1_set_3
Study_1_set_4
Study_1_set_5
Study_1_set_6

Study_2_set_1
Study_2_set_2
Study_2_set_3
Study_2_set_4
Study_2_set_5
Study_2_set_6
Study_2_set_7
Study_2_set_8
Study_2_set_9
Study_2_set_10
Study_2_set_11
Study_2_set_12

Study_3_set_1
Study_3_set_2
Study_3_set_3
Study_3_set_4
Study_3_set_5
Study_3_set_6

