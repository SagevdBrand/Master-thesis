########################################
########################################
###### Performance measure visualization

################
## Setting up ##
################
source("./src/setup.R")
scenarios <- readRDS(paste0(setting_path, "studies.RDS"))

df_perf <- readRDS(paste0(performance_general_path, "all_pm_part_2.RDS"))

################################
## Pre-processing of the data ##
################################

## Add some scenario information to the columns:
scenarios_settings <- scenarios[, (colnames(scenarios) %in% c("study", "scenario", "n_setting", "noise", "dim", "model", "pred_selection", "prev"))]
df_perf$scenario <- as.character(paste0("Scenario_", df_perf$scenario))
df_perf <- full_join(df_perf, scenarios_settings, by = c("study", "scenario"))


df_perf <- df_perf %>% mutate(estimand = case_when(estimand == "auc" ~ "AUC",
                                                   estimand == "calib_slope" ~  "log(Slope)",
                                                   estimand == "calib_int" ~ "Calibration Intercept",
                                                   estimand == "Tjur" ~ "R2 Tjur",
                                                   estimand == "R2_CS" ~ "R2 Cox Snell",
                                                   estimand == "eci" ~ "ECI",
                                                   estimand == "mape" ~ "MAPE",
                                                   estimand == "rmspe" ~ "rMSPE",
                                                   TRUE ~ estimand),
                              estimand = fct_relevel(estimand,
                                                     "AUC", 
                                                     "Calibration Intercept",
                                                     "log(Slope)",
                                                     "R2 Tjur",
                                                     "R2 Cox Snell",
                                                     "ECI",
                                                     "MAPE",
                                                     "rMSPE"),
                              approach = fct_relevel(approach,
                                                     "5 fold cross-validation",
                                                     "10 fold cross-validation",
                                                     "10x10 fold cross-validation",
                                                     "Harrell's bootstrap",
                                                     ".632 bootstrap",
                                                     ".632+ bootstrap",
                                                     "Apparent"),
                              prev = case_when(prev == 0.05 ~ "event-fraction = 0.05",
                                               prev == 0.2 ~ "event-fraction = 0.2",
                                               prev == 0.5 ~ "event-fraction = 0.5",
                                               TRUE ~ as.character(prev)),
                              dim = case_when(dim == 6 ~ "6 Candidate predictors",
                                              dim == 30 ~"30 Candidate predictors",
                                              TRUE ~ as.character(dim)),
                              noise = case_when(noise == "none" ~ "No noise predictors",
                                                noise == "half" ~"50% noise predictors",
                                                noise == "default" ~ "20% noise predictors",
                                                TRUE ~ as.character(dim)),
                              n_setting = case_when(n_setting == "n/2" ~ "N/2",
                                                    n_setting == "n" ~ "N",
                                                    n_setting == "n*2" ~ "N*2",
                                                    TRUE ~ as.character(n_setting))
)

# step1 <- df_perf %>% pivot_longer(cols = c(md), names_to = "performance measures", values_to = "value")
# step2 <- step1 %>% dplyr::select(!md_dist & !rmsd_dist)
# step2b <- df_perf 
# 
# step3 <- step2b %>% 
#   pivot_longer(cols = c(md_dist, rmsd_dist), names_to = "performance measures", values_to = "value_se") %>% 
#   mutate(`performance measures` = case_when(`performance measures` == "md_dist" ~ "md",
#                                             `performance measures` == "rmsd_dist" ~ "rmsd",
#                                             TRUE ~ `performance measures`))
# 
# df_perf <- full_join(step3, step2) %>%
#   mutate(`performance measures` = case_when(`performance measures` == "md" ~ "Mean Difference",
#                                             `performance measures` == "rmsd" ~ "Root Mean Squared Difference",
#                                             TRUE ~ `performance measures`))
# 


# Specify color palette:
colors_perf <- c(
  "#90be6d", # green
  "#f28482", #Light coral
  "#84a59d", # Morning blue
  "#2ec4b6", # Blue
  "#a21c1c", # Dark red
  "#ff9f1c", # Orange
  "#e2a0ff", # purple
  "black" # black
)


# Performance measures for the new estimand names:
estimand_names_thesis <- c(
  "AUC",
  "Calibration Intercept",
  "log(Slope)",
  "R2 Tjur",
  "R2 Cox Snell",
  "ECI",
  "MAPE",
  "rMSPE"
)

####################
## visualizations ##
####################
p1_all <- 
    ggplot(data = df_perf %>% filter(study == "Study_1", 
                                     estimand != "log(Slope)"),
           mapping = aes(x = as.factor(n_setting),
                         y = md,
                         ymin = md-md_dist,
                         ymax = md+md_dist,
                         color = approach,
                         shape = pred_selection
                         )) +
    geom_point(position = position_dodge(width = 0.99), size = 1.5) +
    geom_linerange(position = position_dodge(width = 0.99), size = 0.1) +
    geom_point(data = df_perf %>% filter(study == "Study_1", 
                                         estimand == "log(Slope)"), 
               mapping = aes(x = as.factor(n_setting),
                             y = md,
                             color = approach,
                             shape = pred_selection),
               position = position_dodge(width = 0.99), size = 1.5) +    
  geom_linerange(data = df_perf %>% filter(study == "Study_1",
                                           estimand == "log(Slope)"), 
                 mapping = aes(x = as.factor(n_setting),
                               y = md,
                               ymin = md_q1,
                               ymax = md_q3,
                               color = approach
                               ),
                 position = position_dodge(width = 0.99), size = 0.1) +
    scale_color_manual(values = colors_perf) +
    labs(y = "Median difference",
         x = "Sample size setting",
         color = "Validation approach",
         shape = "Predictor selection"
    ) +
    facet_grid(rows = vars(estimand), cols = vars(prev), scales = "free")+
    theme_set(theme_bw(base_size = 11)) +
    #theme(legend.position = c(0.75, 0.23)) +
    guides(color = guide_legend(nrow=4, ncol=2),
           shape = guide_legend(nrow = 2)) +
    theme(legend.position="bottom")+
  theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))

  

p1_thesis <- 
  ggplot(data = df_perf %>% filter(study == "Study_1",
                                   estimand %in% c("AUC", "rMSPE", "R2 Tjur")), 
         mapping = aes(x = as.factor(n_setting),
                       y = md,
                       ymin = md-md_dist,
                       ymax = md+md_dist,
                       color = approach,
                       shape = pred_selection
                       )) +
  geom_point(position = position_dodge(width = 0.99), size = 1.5) +
  geom_linerange(position = position_dodge(width = 0.99), size = 0.1) +
  geom_point(data = df_perf %>% filter(study == "Study_1",
                                       estimand == "log(Slope)"),
             mapping = aes(x = as.factor(n_setting),
                           y = md,
                           color = approach,
                           shape = pred_selection),
             position = position_dodge(width = 0.99), size = 1.5) +
  geom_linerange(data = df_perf %>% filter(study == "Study_1",
                                           estimand == "log(Slope)"),
                 mapping = aes(x = as.factor(n_setting),
                               y = md,
                               ymin = md_q1,
                               ymax = md_q3,
                               color = approach
                               ),
                 position = position_dodge(width = 0.99),
                 size = 0.1) +
  scale_color_manual(values = colors_perf) +
  labs(y = "Median difference",
       x = "Sample size setting",
       color = "Validation approach",
       shape = "Predictor selection"
  ) +
  facet_grid(rows = vars(estimand), cols = vars(prev), scales = "free")+
  theme_set(theme_bw(base_size = 11)) +
  #theme(legend.position = c(0.75, 0.23)) +
  guides(color = guide_legend(nrow=4, ncol=2),
         shape = guide_legend(nrow = 2)) +
  theme(legend.position="bottom")+
  theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))



p2_all <- 
    ggplot(data = df_perf %>% filter(study == "Study_2",
                                     estimand != "log(Slope)"),
                            mapping = aes(x = as.factor(n_setting),
                                y = md,
                                ymin = md-md_dist,
                                ymax = md+md_dist,
                                color = approach,
                                shape = pred_selection
                  )) +
             geom_point(position = position_dodge(width = 0.99), size = 1.5) +
             geom_linerange(position = position_dodge(width = 0.99), size = 0.1) +
             geom_point(data = df_perf %>% filter(study == "Study_2",
                                                  estimand == "log(Slope)"), 
                        mapping = aes(x = as.factor(n_setting),
                                      y = md,
                                      color = approach,
                                      shape = pred_selection),
                        position = position_dodge(width = 0.99), size = 1.5) +
             geom_linerange(data = df_perf %>% filter(study == "Study_2", 
                                                      estimand == "log(Slope)"),
                            mapping = aes(x = as.factor(n_setting),
                                          y = md,
                                          ymin = md_q1,
                                          ymax = md_q3,
                                          color = approach
                            ),
                            position = position_dodge(width = 0.99), size = 0.1) +
             scale_color_manual(values = colors_perf) +
             labs(y = "Median difference",
                  x = "Sample size setting",
                  color = "Validation approach",
                  shape = "Predictor selection"
             ) +
             facet_grid(rows = vars(estimand), cols = vars(dim, noise), scales = "free")+
             theme_set(theme_bw(base_size = 11)) +
             #theme(legend.position = c(0.75, 0.23)) +
             guides(color = guide_legend(nrow=4, ncol=2),
                    shape = guide_legend(nrow = 2)) +
             theme(legend.position="bottom")+
  theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))

    
p2_thesis <- 
  ggplot(data = df_perf %>% filter(study == "Study_2",
                                   estimand %in% c("AUC", "ECI"),
                                   noise == "50% noise predictors"),
         mapping = aes(x = as.factor(n_setting),
                       y = md,
                       ymin = md-md_dist,
                       ymax = md+md_dist,
                       color = approach,
                       shape = pred_selection
                       )) +
  geom_point(position = position_dodge(width = 0.99), size = 1.5) +
  geom_linerange(position = position_dodge(width = 0.99), size = 0.1) +
  geom_point(data = df_perf %>% filter(study == "Study_2",
                                       estimand == "log(Slope)",
                                       noise == "50% noise predictors"), 
             mapping = aes(x = as.factor(n_setting),
                           y = md,
                           color = approach,
                           shape = pred_selection),
             position = position_dodge(width = 0.99), size = 1.5) +
  geom_linerange(data = df_perf %>% filter(study == "Study_2", 
                                           estimand == "log(Slope)",
                                           noise == "50% noise predictors"),
                 mapping = aes(x = as.factor(n_setting),
                               y = md,
                               ymin = md_q1,
                               ymax = md_q3,
                               color = approach
                               ),
  position = position_dodge(width = 0.99), size = 0.1) +
  scale_color_manual(values = colors_perf) +
  labs(y = "Median difference",
       x = "Sample size setting",
       color = "Validation approach",
       shape = "Predictor selection"
  ) +
  facet_grid(rows = vars(estimand), cols = vars(dim), scales = "free")+
  theme_set(theme_bw(base_size = 11)) +
  #theme(legend.position = c(0.75, 0.23)) +
  guides(color = guide_legend(nrow=4, ncol=2),
         shape = guide_legend(nrow = 2)) +
  theme(legend.position="bottom")+
  theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))



p3_all_models <- 
  ggplot(data = df_perf %>% 
           filter(study == "Study_3", 
                  estimand != "log(Slope)", 
                  model %in% c("Firth", "Ridge", "Lasso", "RF")),
         mapping = aes(x = as.factor(n_setting),    
                       y = md,
                       ymin = md-md_dist,
                       ymax = md+md_dist,
                       color = approach,
                       shape = pred_selection
         )) +
  geom_point(position = position_dodge(width = 0.99), size = 1.5) +
  geom_linerange(position = position_dodge(width = 0.99), size = 0.1) +
  geom_point(data = df_perf %>% 
               filter(study == "Study_3", 
                      estimand == "log(Slope)", 
                      model %in% c("Firth", "Ridge", "Lasso", "RF")), 
             mapping = aes(x = as.factor(n_setting),
                           y = md,
                           color = approach,
                           shape = pred_selection),
             position = position_dodge(width = 0.99), size = 1.5) +
  geom_linerange(data = df_perf %>%
                   filter(study == "Study_3",
                          estimand == "log(Slope)", 
                          model %in% c("Firth", "Ridge", "Lasso", "RF")), 
                 mapping = aes(x = as.factor(n_setting),
                               y = md,
                               ymin = md_q1,
                               ymax = md_q3,
                               color = approach
                               ),
                 position = position_dodge(width = 0.99), size = 0.1) +
  scale_color_manual(values = colors_perf) +
  labs(y = "Median difference",
       x = "Sample size setting",
       color = "Validation approach",
       shape = "Predictor selection"
  ) +
  facet_grid(rows = vars(estimand), cols = vars(model), scales = "free")+
  theme_set(theme_bw(base_size = 11)) +
  #theme(legend.position = c(0.75, 0.23)) +
  guides(color = guide_legend(nrow=4, ncol=2),
         shape = guide_legend(nrow = 2)) +
  theme(legend.position="bottom")+
  theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))




p3_thesis_models <- 
  ggplot(data = df_perf %>% filter(study == "Study_3",  
                                   estimand %in% c("AUC", "R2 Tjur"),
                                   model %in% c("Firth", "Ridge", "Lasso", "RF")
                                   ),
         mapping = aes(x = as.factor(n_setting),
                       y = md,
                       ymin = md-md_dist,
                       ymax = md+md_dist,
                       color = approach
                       )) +
  geom_point(position = position_dodge(width = 0.99), size = 1.5) +
  geom_linerange(position = position_dodge(width = 0.99), size = 0.1) +
  geom_point(data = df_perf %>% filter(study == "Study_3",
                                       estimand == "log(Slope)", 
                                       model %in% c("Firth", "Ridge", "Lasso", "RF")
                                       ),
             mapping = aes(x = as.factor(n_setting),
                           y = md,
                           color = approach
                           ),
             position = position_dodge(width = 0.99), size = 1.5) +
  geom_linerange(data = df_perf %>% filter(study == "Study_3", 
                                           estimand == "log(Slope)", 
                                           model %in% c("Firth", "Ridge", "Lasso", "RF")
                                           ),
                 mapping = aes(x = as.factor(n_setting),
                               y = md,
                               ymin = md_q1,
                               ymax = md_q3,
                               color = approach
                               ),
  position = position_dodge(width = 0.99), size = 0.1) +
  scale_color_manual(values = colors_perf) +
  labs(y = "Median difference",
       x = "Sample size setting",
       color = "Validation approach"
  ) +
  facet_grid(rows = vars(estimand), 
             cols = vars(model), 
             scales = "free")+
  theme_set(theme_bw(base_size = 11)) +
  #theme(legend.position = c(0.75, 0.23)) +
  guides(color = guide_legend(nrow=4, ncol=2),
         shape = guide_legend(nrow = 2)) +
  theme(legend.position="bottom")+
  theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))






p3_all_ML_CART <- 
  ggplot(data = df_perf %>% 
           filter(study == "Study_3", 
                  estimand != "log(Slope)", 
                  model %in% c("ML", "CART")),
         mapping = aes(x = as.factor(n_setting),    
                       y = md,
                       ymin = md-md_dist,
                       ymax = md+md_dist,
                       color = approach,
                       shape = pred_selection
         )) +
  geom_point(position = position_dodge(width = 0.99), size = 1.5) +
  geom_linerange(position = position_dodge(width = 0.99), size = 0.1) +
  geom_point(data = df_perf %>% 
               filter(study == "Study_3", 
                      estimand == "log(Slope)", 
                      model %in%c("ML", "CART")), 
             mapping = aes(x = as.factor(n_setting),
                           y = md,
                           color = approach,
                           shape = pred_selection),
             position = position_dodge(width = 0.99), size = 1.5) +
  geom_linerange(data = df_perf %>%
                   filter(study == "Study_3",
                          estimand == "log(Slope)", 
                          model %in%c("ML", "CART")), 
                 mapping = aes(x = as.factor(n_setting),
                               y = md,
                               ymin = md_q1,
                               ymax = md_q3,
                               color = approach
                 ),
                 position = position_dodge(width = 0.99), size = 0.1) +
  scale_color_manual(values = colors_perf) +
  labs(y = "Median difference",
       x = "Sample size setting",
       color = "Validation approach",
       shape = "Predictor selection"
  ) +
  facet_grid(rows = vars(estimand), cols = vars(model), scales = "free")+
  theme_set(theme_bw(base_size = 11)) +
  #theme(legend.position = c(0.75, 0.23)) +
  guides(color = guide_legend(nrow=4, ncol=2),
         shape = guide_legend(nrow = 2)) +
  theme(legend.position="bottom")+
  theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))




p3_thesis_ML_CART <- 
  ggplot(data = df_perf %>% filter(study == "Study_3",  
                                   estimand %in% c("AUC", "R2 Tjur"),
                                   model %in%c("ML", "CART")),
         mapping = aes(x = as.factor(n_setting),
                       y = md,
                       ymin = md-md_dist,
                       ymax = md+md_dist,
                       color = approach
         )) +
  geom_point(position = position_dodge(width = 0.99), size = 1.5) +
  geom_linerange(position = position_dodge(width = 0.99), size = 0.1) +
  geom_point(data = df_perf %>% filter(study == "Study_3",
                                       estimand == "log(Slope)", 
                                       model %in%c("ML", "CART")),
             mapping = aes(x = as.factor(n_setting),
                           y = md,
                           color = approach
             ),
             position = position_dodge(width = 0.99), size = 1.5) +
  geom_linerange(data = df_perf %>% filter(study == "Study_3", 
                                           estimand == "log(Slope)", 
                                           model %in% c("ML", "CART")),
                 mapping = aes(x = as.factor(n_setting),
                               y = md,
                               ymin = md_q1,
                               ymax = md_q3,
                               color = approach
                 ),
                 position = position_dodge(width = 0.99), size = 0.1) +
  scale_color_manual(values = colors_perf) +
  labs(y = "Median difference",
       x = "Sample size setting",
       color = "Validation approach"
  ) +
  facet_grid(rows = vars(estimand), cols = vars(model), scales = "free")+
  theme_set(theme_bw(base_size = 11)) +
  #theme(legend.position = c(0.75, 0.23)) +
  guides(color = guide_legend(nrow=4, ncol=2),
         shape = guide_legend(nrow = 2)) +
  theme(legend.position="bottom")+
  theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))



p3_all_ML_Firth <- 
  ggplot(data = df_perf %>% filter(study == "Study_3",  
                                   estimand != "log(Slope)", 
                                   model %in%c("ML", "Firth")),
         mapping = aes(x = as.factor(n_setting),
                       y = md,
                       ymin = md-md_dist,
                       ymax = md+md_dist,
                       color = approach
         )) +
  geom_point(position = position_dodge(width = 0.99), size = 1.5) +
  geom_linerange(position = position_dodge(width = 0.99), size = 0.1) +
  geom_point(data = df_perf %>% filter(study == "Study_3",
                                       estimand == "log(Slope)", 
                                       model %in%c("ML", "Firth")),
             mapping = aes(x = as.factor(n_setting),
                           y = md,
                           color = approach
             ),
             position = position_dodge(width = 0.99), size = 1.5) +
  geom_linerange(data = df_perf %>% filter(study == "Study_3", 
                                           estimand == "log(Slope)", 
                                           model %in% c("ML", "Firth")),
                 mapping = aes(x = as.factor(n_setting),
                               y = md,
                               ymin = md_q1,
                               ymax = md_q3,
                               color = approach
                 ),
                 position = position_dodge(width = 0.99), size = 0.1) +
  scale_color_manual(values = colors_perf) +
  labs(y = "Median difference",
       x = "Sample size setting",
       color = "Validation approach"
  ) +
  facet_grid(rows = vars(estimand), cols = vars(model), scales = "free")+
  theme_set(theme_bw(base_size = 11)) +
  #theme(legend.position = c(0.75, 0.23)) +
  guides(color = guide_legend(nrow=4, ncol=2),
         shape = guide_legend(nrow = 2)) +
  theme(legend.position="bottom")+
  theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))



  ggsave(paste0(full_performance_plots,"S7_Performance_measures_study_1.pdf"), plot = p1_all,  width = 20, height = 25, units = "cm")
  #ggsave(paste0(performance_measures_plots,"performance_study_1_thesis.pdf"), plot = p1_thesis, width = 19, height = 15, units = "cm")
  
  
  ggsave(paste0(full_performance_plots,"S8_Performance_measures_study_2.pdf"), plot = p2_all, width = 25, height = 30, units = "cm")
  #ggsave(paste0(performance_measures_plots,"performance_study_2_thesis.pdf"), plot = p2_thesis, width = 21, height = 15, units = "cm")
  
  ggsave(paste0(full_performance_plots,"S9_Performance_measures_study_3_F_Ridge_Lasso_RF.pdf"), plot = p3_all_models, width = 25, height = 30, units = "cm")
  #ggsave(paste0(performance_measures_plots,"performance_study_3_thesis_models.pdf"), plot = p3_thesis_models, width = 19, height = 15, units = "cm")
  ggsave(paste0(full_performance_plots,"S10_Performance_measures_study_3_ML_CART.pdf"), plot = p3_all_ML_CART, width = 25, height = 30, units = "cm")
  #ggsave(paste0(performance_measures_plots,"performance_study_3_thesis_ML_CART.pdf"), plot = p3_thesis_ML_CART, width = 19, height = 15, units = "cm")
  #ggsave(paste0(full_performance_plots,"performance_study_3_full_Firth_ML.pdf"), plot = p3_all_ML_Firth, width = 25, height = 30, units = "cm")
  


  #####################
  # Nested loop plots #
  #####################
  source("./src/setup.R")
  scenarios <- readRDS(paste0(setting_path, "studies.RDS"))
  
  df_perf <- readRDS(paste0(performance_general_path, "all_pm_part_2.RDS"))
  
  ################################
  ## Pre-processing of the data ##
  ################################
  
  ## Add some scenario information to the columns:
  scenarios_settings <- scenarios[, (colnames(scenarios) %in% c("study", "scenario", "n_setting", "noise", "dim", "model", "pred_selection", "prev"))]
  df_perf$scenario <- as.character(paste0("Scenario_", df_perf$scenario))
  df_perf <- full_join(df_perf, scenarios_settings, by = c("study", "scenario"))
  
  
  df_perf <- df_perf %>% mutate(estimand = case_when(estimand == "auc" ~ "AUC",
                                                     estimand == "calib_slope" ~  "log(Slope)",
                                                     estimand == "calib_int" ~ "Calibration Intercept",
                                                     estimand == "Tjur" ~ "R2 Tjur",
                                                     estimand == "R2_CS" ~ "R2 Cox Snell",
                                                     estimand == "eci" ~ "ECI",
                                                     estimand == "mape" ~ "MAPE",
                                                     estimand == "rmspe" ~ "rMSPE",
                                                     TRUE ~ estimand),
                                estimand = fct_relevel(estimand,
                                                       "AUC", 
                                                       "Calibration Intercept",
                                                       "log(Slope)",
                                                       "R2 Tjur",
                                                       "R2 Cox Snell",
                                                       "ECI",
                                                       "MAPE",
                                                       "rMSPE"),
                                approach = fct_relevel(approach,
                                                       "5 fold cross-validation",
                                                       "10 fold cross-validation",
                                                       "10x10 fold cross-validation",
                                                       "Harrell's bootstrap",
                                                       ".632 bootstrap",
                                                       ".632+ bootstrap",
                                                       "Apparent"),
                                dim = case_when(dim == 6 ~ "6 Candidate predictors",
                                                dim == 30 ~"30 Candidate predictors",
                                                TRUE ~ as.character(dim)),
                                noise = case_when(noise == "none" ~ "No noise predictors",
                                                  noise == "half" ~"50% noise predictors",
                                                  noise == "default" ~ "20% noise predictors",
                                                  TRUE ~ as.character(dim)),
                                n_setting = case_when(n_setting == "n/2" ~ "N/2",
                                                      n_setting == "n" ~ "N",
                                                      n_setting == "n*2" ~ "N*2",
                                                      TRUE ~ as.character(n_setting)),
                                apparent_result = ifelse(approach == "Apparent", "1", "0")
  )
  
  df_perf$scenario <- as.numeric(gsub("Scenario_", "", df_perf$scenario))
  
  
  ###########
  # Study 1 #
  ###########
  plot_data_s1 <- df_perf %>% filter(study == "Study_1") %>% 
    mutate(n_setting_plot = case_when(n_setting == "N/2" ~ 0.33,
                                      n_setting == "N" ~ 0.35,
                                      n_setting == "N*2" ~ 0.37,
                                      TRUE ~ 0),
           prev_plot = case_when(prev == 0.05 ~ 0.43,
                                 prev == 0.2 ~ 0.45,
                                 prev == 0.5 ~ 0.47,
                                 TRUE ~ 0),
           pred_sel_plot = case_when(pred_selection == "none" ~ 0.53,
                                     pred_selection == "<0.157" ~ 0.56,
                                     TRUE ~ 0)
    )
  
  
  
  p1_thesis <- ggplot(df_perf %>% filter(study == "Study_1", estimand == "log(Slope)"),
                      mapping = aes(x = scenario,
                                    y = md,
                                    color = approach)) +
    geom_step(direction = "hv", size = 0.7, aes(linetype = apparent_result)) +
    geom_step(data = plot_data_s1, mapping = aes(x = scenario, y = pred_sel_plot), color = "black",
              size = 0.4) +
    geom_label(x = 1, y = 0.55, 
               label = "Predictor Selection: None, < 0.157",
               hjust = "inward",
               vjust = "outward",
               color = "black",
               label.size = 0) +
    geom_step(data = plot_data_s1, mapping = aes(x = scenario, y = prev_plot),color = "black", 
              size = 0.4) +
    geom_label(x = 1, y = 0.48, 
               label = "Event fraction: 0.05, 0.2, 0.5",
               hjust = "inward",
               vjust = "outward",
               color = "black",
               label.size = 0) +
    geom_step(data = plot_data_s1, mapping = aes(x = scenario, y = n_setting_plot), color = "black", 
              size = 0.4) +
    geom_label(x = 1, y = 0.38, 
               label = "Sample size setting: N/2, N, N*2",
               hjust = "inward",
               vjust = "outward",
               color = "black",
               label.size = 0) +
    geom_hline(yintercept =  0,linetype="dotted") +
    theme_bw(base_size = 11) + 
    theme(legend.position="bottom")+
    scale_color_manual(values = colors_perf) +
    scale_linetype_manual(values = c("solid", "dashed"))  +
    labs(y = "Median difference log(Calibration Slope)",
         x = "Scenarios study 1",
         color = "Validation approach",
         linetype = NULL
    ) +
    scale_x_continuous(breaks = 1:18) +
    scale_y_continuous(breaks = seq(-0.2, 0.3, by = 0.05)) +
    guides(color = guide_legend(nrow=2, ncol=4),
           linetype = FALSE) +
    theme(plot.margin = unit(c(0.1, 0.6, 0.1, 0.1), "cm"),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"))
  
  p1_thesis
  
  ggsave(paste0(performance_measures_plots,"nlp_study_1_thesis.pdf"), plot = p1_thesis, width = 19, height = 15, units = "cm")
  
  
  ###########
  # Study 2 #
  ###########
  plot_data_s2 <- df_perf %>% filter(study == "Study_2") %>% 
    mutate(dim_plot = case_when(dim == "6 Candidate predictors" ~ 0.64,
                                dim == "30 Candidate predictors" ~ 0.68,
                                TRUE ~ 0),
           
           n_setting_plot = case_when(n_setting == "N/2" ~ 0.77,
                                      n_setting == "N" ~ 0.8,
                                      n_setting == "N*2" ~ 0.83,
                                      TRUE ~ 0),
           
           noise_plot = case_when(noise == "No noise predictors" ~ 0.94,
                                  noise == "50% noise predictors" ~ 0.97,
                                  noise == "20% noise predictors" ~ 1.00,
                                  TRUE ~ 0),
           
           pred_sel_plot = case_when(pred_selection == "none" ~ 1.09,
                                     pred_selection == "<0.157" ~ 1.12,
                                     TRUE ~ 0)
    )
  
  
  
  p2_thesis <- ggplot(df_perf %>% filter(study == "Study_2", estimand == "ECI"),
                      mapping = aes(x = scenario,
                                    y = md,
                                    color = approach)) +
    geom_step(direction = "hv", size = 0.7, aes(linetype = apparent_result)) +
    geom_step(data = plot_data_s2, mapping = aes(x = scenario, y = pred_sel_plot), color = "black",
              size = 0.4) +
    geom_label(x = 1, y = 1.13,
               label = "Predictor Selection: None, < 0.157",
               hjust = "inward",
               vjust = "outward",
               color = "black",
               label.size = 0) +
    geom_step(data = plot_data_s2, mapping = aes(x = scenario, y = dim_plot),color = "black",
              size = 0.4) +
    geom_label(x = 1, y = 0.69,
               label = "Candidate predictors: 6, 30",
               hjust = "inward",
               vjust = "outward",
               color = "black",
               label.size = 0) +
    geom_step(data = plot_data_s2, mapping = aes(x = scenario, y = noise_plot), color = "black",
              size = 0.4) +
    geom_label(x = 1, y = 0.84,
               label = "Sample size setting: N/2, N, N*2",
               hjust = "inward",
               vjust = "outward",
               color = "black",
               label.size = 0) +
    geom_step(data = plot_data_s2, mapping = aes(x = scenario, y = n_setting_plot), color = "black",
              size = 0.4) +
    geom_label(x = 1, y = 1,
               label = "Noise contribution: No noise, 50% noise",
               hjust = "inward",
               vjust = "outward",
               color = "black",
               label.size = 0) +
    geom_hline(yintercept =  0,linetype="dotted") +
    theme_bw(base_size = 11) + 
    theme(legend.position="bottom")+
    scale_color_manual(values = colors_perf) +
    labs(y = "Median difference ECI",
         x = "Scenarios study 2",
         color = "Validation approach",
         linetype = NULL
    ) +
    scale_x_continuous(breaks = 1:24) +
    scale_y_continuous(breaks = seq(-0.2, 0.6, by = 0.1),
                       limits = c(-0.1, 1.15)) +
    guides(color = guide_legend(nrow=2),
          linetype = FALSE) +
    theme(plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm"),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"))
  
  ggsave(paste0(performance_measures_plots,"nlp_study_2_thesis.pdf"), plot = p2_thesis, width = 19, height = 13, units = "cm")
  
