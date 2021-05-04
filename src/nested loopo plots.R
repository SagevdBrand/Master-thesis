#####################
# Nested loop plots #
#####################
source("./src/setup.R")
scenarios <- readRDS(paste0(setting_path, "studies.RDS"))

df_perf <- readRDS(paste0(performance_general_path, "all_pm_batch_7.RDS"))

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
                                                TRUE ~ as.character(dim))
)

df_perf$scenario <- as.numeric(gsub("Scenario_", "", df_perf$scenario))


###########
# Study 1 #
###########
plot_data_s1 <- df_perf %>% filter(study == "Study_1") %>% 
  mutate(n_setting_plot = case_when(n_setting == "n/2" ~ 0.033,
                                    n_setting == "n" ~ 0.035,
                                    n_setting == "n*2" ~ 0.037,
                                    TRUE ~ 0),
         prev_plot = case_when(prev == 0.05 ~ 0.04,
                               prev == 0.2 ~ 0.042,
                               prev == 0.5 ~ 0.044,
                               TRUE ~ 0),
         pred_sel_plot = case_when(pred_selection == "none" ~ 0.047,
                                   pred_selection == "<0.157" ~ 0.05,
                                   TRUE ~ 0)
  )



p1_thesis <- ggplot(df_perf %>% filter(study == "Study_1", estimand == "rMSPE"),
       mapping = aes(x = scenario,
                     y = md,
                     color = approach)) +
  geom_step(direction = "hv", size = 1) +
  geom_step(data = plot_data_s1, mapping = aes(x = scenario, y = pred_sel_plot), color = "black",
            size = 0.5) +
  geom_label(x = 1, y = 0.05, 
            label = "Predictor Selection: None, < 0.157",
            hjust = "inward",
            vjust = "outward",
            color = "black",
            label.size = 0) +
  geom_step(data = plot_data_s1, mapping = aes(x = scenario, y = prev_plot),color = "black", 
            size = 0.5) +
  geom_label(x = 1, y = 0.044, 
             label = "Event fraction: 0.05, 0.2, 0.5",
             hjust = "inward",
             vjust = "outward",
             color = "black",
             label.size = 0) +
  geom_step(data = plot_data_s1, mapping = aes(x = scenario, y = n_setting_plot), color = "black", 
            size = 0.5) +
  geom_label(x = 1, y = 0.037, 
             label = "Sample size setting: n/2, n, n*2",
             hjust = "inward",
             vjust = "outward",
             color = "black",
             label.size = 0) +
  geom_hline(yintercept =  0,linetype="dotted") +
  theme_bw(base_size = 11) + 
  theme(legend.position="bottom")+
  scale_color_manual(values = colors_perf) +
  labs(y = "Median difference rMSPE",
       x = "Scenarios study 1",
       color = "Validation approach"
  ) +
  scale_x_continuous(breaks = 1:18) +
  guides(color = guide_legend(nrow=3, ncol=3)) +
  theme(plot.margin = unit(c(0.1, 0.4, 0.1, 0.1), "cm"))

ggsave(paste0(performance_measures_plots,"nlp_study_1_thesis.pdf"), plot = p1_thesis, width = 18, height = 15, units = "cm")


###########
# Study 2 #
###########
plot_data_s2 <- df_perf %>% filter(study == "Study_2") %>% 
  mutate(n_setting_plot = case_when(n_setting == "n/2" ~ 0.63,
                                    n_setting == "n" ~ 0.66,
                                    n_setting == "n*2" ~ 0.69,
                                    TRUE ~ 0),
         dim_plot = case_when(dim == "6 Candidate predictors" ~ 0.75,
                              dim == "30 Candidate predictors" ~ 0.78,
                              TRUE ~ 0),
         pred_sel_plot = case_when(pred_selection == "none" ~ 0.84,
                                   pred_selection == "<0.157" ~ 0.87,
                                   TRUE ~ 0),
         noise_plot = case_when(noise == "No noise predictors" ~ 0.94,
                                noise == "50% noise predictors" ~ 0.97,
                                noise == "20% noise predictors" ~ 1,
                                TRUE ~ 0)
  )



p2_thesis <- ggplot(df_perf %>% filter(study == "Study_2", estimand == "ECI"),
                    mapping = aes(x = scenario,
                                  y = md,
                                  color = approach)) +
  geom_step(direction = "hv", size = 1) +
  geom_step(data = plot_data_s2, mapping = aes(x = scenario, y = pred_sel_plot), color = "black",
            size = 0.5) +
  geom_label(x = 1, y = 0.87,
             label = "Predictor Selection: None, < 0.157",
             hjust = "inward",
             vjust = "outward",
             color = "black",
             label.size = 0) +
  geom_step(data = plot_data_s2, mapping = aes(x = scenario, y = dim_plot),color = "black",
            size = 0.5) +
  geom_label(x = 1, y = 0.78,
             label = "Candidate predictors: 6, 30",
             hjust = "inward",
             vjust = "outward",
             color = "black",
             label.size = 0) +
  geom_step(data = plot_data_s2, mapping = aes(x = scenario, y = noise_plot), color = "black",
            size = 0.5) +
  geom_label(x = 1, y = 0.69,
             label = "Sample size setting: n/2, n, n*2",
             hjust = "inward",
             vjust = "outward",
             color = "black",
             label.size = 0) +
  geom_step(data = plot_data_s2, mapping = aes(x = scenario, y = n_setting_plot), color = "black",
            size = 0.5) +
  geom_label(x = 1, y = 0.97,
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
       color = "Validation approach"
  ) +
  scale_x_continuous(breaks = 1:24) +
  guides(color = guide_legend(nrow=3, ncol=3)) +
  theme(plot.margin = unit(c(0.1, 0.4, 0.1, 0.1), "cm"))

ggsave(paste0(performance_measures_plots,"nlp_study_2_thesis.pdf"), plot = p2_thesis, width = 18, height = 15, units = "cm")

              

