#################################
######## Visualization ##########
#################################
################
## Setting up ##
################
source("./src/setup.R")
scenarios <- readRDS(paste0(setting_path, "studies.RDS"))

## Loading the above estimand files:
df_all <- readRDS(paste0(estimands_general_path, "all_estimands_part_2.RDS"))

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

################################################################
## Rename estimands for plots and change the order of factors ##
################################################################
df_all_long <- df_all_long %>% mutate(estimand = case_when(estimand == "auc" ~ "AUC",
                                                           estimand == "calib_slope" ~  "Calib. Slope",
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
                                                             "Calib. Slope",
                                                             "R2 Tjur",
                                                             "R2 Cox Snell",
                                                             "ECI",
                                                             "MAPE",
                                                             "rMSPE"),
                                      approach = fct_relevel(approach,
                                                             "External",
                                                             "Apparent",
                                                             "5 fold cross-validation",
                                                             "10 fold cross-validation",
                                                             "10x10 fold cross-validation",
                                                             "Harrell's bootstrap",
                                                             ".632 bootstrap",
                                                             ".632+ bootstrap"),
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

############
## Colors ##
############

# Specify color palette:
colors <- c("white",   # white
            "#e2a0ff", # purple
            "#90be6d", # green
            "#f28482", #Light coral
            "#84a59d", # Morning blue
            "#2ec4b6", # Blue
            "#a21c1c", # Dark red
            "#ff9f1c" # Orange
)


############################
############################
########## STUDY 1 #########
############################
############################

# Results of study 1:
df_s1_long<- df_all_long %>% filter(study == "Study_1")
df_s1_long$scenario <- as.numeric(gsub("Scenario_", "", df_s1_long$scenario))
df_s1_long <- df_s1_long %>% mutate(value = case_when(estimand == "Calib. Slope" & value > 10 ~ 10,
                                                      TRUE ~ value))

## Estimands without predictor selection, but only a selection
p1 <- 
  ggplot(data = df_s1_long %>% 
           filter(pred_selection == "none",
                  estimand %in% c("AUC", "rMSPE", "Calib. Slope", "R2 Tjur")),
         mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
  geom_boxplot(position = position_dodge(width = 0.95), 
               size = 0.005, 
               outlier.shape = NA) +
  theme_set(theme_bw(base_size = 11)) +
  scale_fill_manual(values = colors) +
  labs(y = "Estimand value",
       x = "Sample size setting",
       fill = "Validation approach"
       ) +
  facet_grid(rows = vars(estimand), cols = vars(prev), scales = "free")+ 
  theme(legend.position="bottom")+
  guides(color = guide_legend(nrow=4, ncol=3)) +
  theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))

### Figure to be included in thesis
## Manually setting the scales
scales_y_1 <- list(
  `AUC` = scale_y_continuous(limits = c(0.5, 0.9)),
  `Calib. Slope` = scale_y_continuous(trans = "log10", limits = c(0.45, 1.1)),
  `R2 Tjur` = scale_y_continuous(limits = c(0,0.3)),
  `rMSPE` = scale_y_continuous(limits = c(0,0.25))
)


p2 <- 
  ggplot(data = df_s1_long %>% 
           filter(pred_selection == "<0.157",
                  estimand %in% c("AUC", "rMSPE", "Calib. Slope",  "R2 Tjur")),
         mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
  geom_boxplot(position = position_dodge(width = 0.95), 
               outlier.shape = NA,
               size = 0.005) +
  theme_set(theme_bw(base_size = 11)) +
  scale_fill_manual(values = colors) +
  labs(y = "Estimand value",
       x = "Sample size setting",
       fill = "Validation approach") +
  facet_grid_sc(rows = vars(estimand), cols = vars(prev), scales = list(y = scales_y_1)) +
  theme(legend.position="bottom")+
  guides(color = guide_legend(nrow=4, ncol=3))+
  theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))




# Figures of all estimands (For supplementary material)
p3 <- 
  ggplot(data = df_s1_long %>% 
           filter(pred_selection == "none"), 
         mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
  geom_boxplot(position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
  theme_set(theme_bw(base_size = 11)) +
  scale_fill_manual(values = colors) +
  labs(y = "Estimand value",
       x = "Sample size setting",
       fill = "Validation approach") +
  facet_grid(rows = vars(estimand), cols = vars(prev), scales = "free")+ 
  theme(legend.position="bottom")+
  guides(color = guide_legend(nrow=4, ncol=3))+
  theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))


p4 <- 
  ggplot(data = df_s1_long %>% 
           filter(pred_selection == "<0.157"),
         mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
  geom_boxplot(position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
  theme_set(theme_bw(base_size = 11)) +
  scale_fill_manual(values = colors) +
  labs(y = "Estimand value",
       x = "Sample size setting",
       fill = "Validation approach") +
  facet_grid(rows = vars(estimand), cols = vars(prev), scales = "free")+ 
  theme(legend.position="bottom")+
  guides(color = guide_legend(nrow=4, ncol=3))+
  theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))


## Saving the plots ##

### No predictor selection 
# original scale
# ggsave(paste0(estimand_plots,"estimands_study_1_no_pred_sel.pdf"), plot = p1, width = 19, height = 15, units = "cm")
# All estimands
ggsave(paste0(full_estimand_plots,"S1_All_metrics_study_1_no_pred_sel.pdf"), plot = p3, width = 20, height = 25, units = "cm")

### Predictor selection
# original scale
ggsave(paste0(estimand_plots,"estimands_study_1_pred_sel.pdf"), plot = p2, width = 19, height = 18, units = "cm")
# All estimands
ggsave(paste0(full_estimand_plots,"S2_All_metrics_study_1_pred_sel.pdf"), plot = p4, width = 20, height = 25, units = "cm")


############################
############################
########## STUDY 2 #########
############################
############################

# Results of study 2:
df_s2_long<- df_all_long %>% filter(study == "Study_2")
df_s2_long$scenario <- as.numeric(gsub("Scenario_", "", df_s2_long$scenario))
df_s2_long <- df_s2_long %>% mutate(value = case_when(estimand == "Calib. Slope" & value > 10 ~ 10,
                                                      TRUE ~ value))
# 50% noise, part of the estimands, no predictor selection
 p1 <-  
    ggplot(data = df_s2_long %>%
             filter(pred_selection == "none",
                    noise == "50% noise predictors",
                    estimand %in% c("AUC", "Calib. Slope", "ECI")), 
           mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
    geom_boxplot(position = position_dodge(width = 0.95), size = 0.005, outlier.shape = NA) +
    theme_set(theme_bw(base_size = 11)) +
    scale_fill_manual(values = colors) +
    labs(y = "Estimand value",
         x = "Sample size setting",
         fill = "Validation approach") +
    facet_grid(rows = vars(estimand), cols = vars(dim), scales = "free")+ 
    theme(legend.position="bottom")+
    guides(color = guide_legend(nrow=4, ncol=3))+
    theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))
  
 ### Figure to be included in thesis
 # 50% noise, part of the estimands, predictor selection
 
 ## Manually setting the scales
 scales_y_2 <- list(
   `AUC` = scale_y_continuous(limits = c(0.5, 0.9)),
   `Calib. Slope` = scale_y_continuous(trans = "log10", limits = c(0.2, 10)),
   `ECI` = scale_y_continuous(limits = c(0,1))
 )
 
 p2 <- ggplot(data = df_s2_long %>%
                 filter(pred_selection == "<0.157",
                        noise == "50% noise predictors",
                        estimand %in% c("AUC", "ECI", "Calib. Slope")),
                        mapping = aes(x = as.factor(n_setting),  y =value, fill = approach)) +
    geom_boxplot(position = position_dodge(width = 0.95),size = 0.005,  outlier.shape = NA) +
    theme_set(theme_bw(base_size = 11)) +
    scale_fill_manual(values = colors) +
    labs(y = "Estimand value",
         x = "Sample size setting",
         fill = "Validation approach") +
   facet_grid_sc(rows = vars(estimand), cols = vars(dim), scales = list(y = scales_y_2)) +
    theme(legend.position="bottom")+
    guides(color = guide_legend(nrow=4, ncol=3))+
    theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))
  
  
  ##### ALL ESTIMANDS ######
  p3 <- ggplot(data = df_s2_long %>%
                 filter(pred_selection == "none"), 
               mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
    geom_boxplot(position = position_dodge(width = 0.95),size = 0.005,  outlier.shape = NA) +
    theme_set(theme_bw(base_size = 11)) +
    scale_fill_manual(values = colors) +
    labs(y = "Estimand value",
         x = "Sample size setting",
         fill = "Validation approach") +
    facet_grid(rows = vars(estimand), cols = vars(noise, dim), scales = "free")+ 
    theme(legend.position="bottom")+
    guides(color = guide_legend(nrow=4, ncol=3))+
    theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))
  
  
  
  p4 <-  ggplot(data = df_s2_long %>%
                 filter(pred_selection == "<0.157"), 
               mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
    geom_boxplot(position = position_dodge(width = 0.95),size = 0.005,  outlier.shape = NA) +
    theme_set(theme_bw(base_size = 11)) +
    scale_fill_manual(values = colors) +
    labs(y = "Estimand value",
         x = "Sample size setting",
         fill = "Validation approach") +
    facet_grid(rows = vars(estimand), cols = vars(noise, dim), scales = "free")+ 
    theme(legend.position="bottom")+
    guides(color = guide_legend(nrow=4, ncol=3))+
    theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))
  
  
  # Saving the plots:
  ### 50% noise and no predictor selection
  # Original scale
  #ggsave(paste0(estimand_plots,"estimands_study_2_no_pred_sel.pdf"), plot = p1, width = 19, height = 15, units = "cm")
  # All estimands: both with and without noise, no predictor selection
  ggsave(paste0(full_estimand_plots,"S3_All_metrics_study_2_no_pred_sel.pdf"), plot = p3, width = 20, height = 25, units = "cm")
  
  
  ### 50% noise and predictor selection
  # Original scale
  ggsave(paste0(estimand_plots,"estimands_study_2_pred_sel.pdf"), plot = p2, width = 19, height = 18, units = "cm")
  # All estimands: both with and without noise, predictor selection
  ggsave(paste0(full_estimand_plots,"S4_All_metrics_study_2_pred_sel.pdf"), plot = p4, width = 20, height = 25, units = "cm")
  
  
############################
############################
########## STUDY 3 #########
############################
############################

# Results of study 2:
df_s3_long<- df_all_long %>% filter(study == "Study_3")
df_s3_long$scenario <- as.numeric(gsub("Scenario_", "", df_s3_long$scenario))
df_s3_long <- df_s3_long %>% mutate(value = case_when(estimand == "Calib. Slope" & value > 10 ~ 10,
                                                      TRUE ~ value))

# For thesis:
## Manually setting the scales
scales_y_3 <- list(
  `AUC` = scale_y_continuous(limits = c(0.5, 1)),
  `Calib. Slope` = scale_y_continuous(trans = "log10", limits = c(0.03, 10)),
  `R2 Tjur`= scale_y_continuous(limits = c(0,0.7)),
  `ECI` = scale_y_continuous(limits = c(0,0.8))
)



p1 <-  
  ggplot(data = df_s3_long %>% 
                filter(model %in% c("Firth", "Ridge",  "RF"), 
                       estimand %in% c("AUC", "Calib. Slope", "R2 Tjur", "ECI")), 
              mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
  geom_boxplot(position = position_dodge(width = 0.95),  size = 0.005, outlier.shape = NA) +
  theme_set(theme_bw(base_size = 11)) +
  scale_fill_manual(values = colors) +
  labs(y = "Estimand value",
       x = "Sample size setting",
       fill = "Validation approach") +
  facet_grid_sc(rows = vars(estimand), cols = vars(model), scales = list(y =scales_y_3))+ 
  theme(legend.position="bottom")+
  guides(color = guide_legend(nrow=4, ncol=3))+
  theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))
   

 
  # All estimands
  p2 <-ggplot(data = df_s3_long %>% 
                filter(model %in% c("Firth", "Ridge", "RF")), 
              mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
    geom_boxplot(position = position_dodge(width = 0.95), size = 0.005, outlier.shape = NA) +
    theme_set(theme_bw(base_size = 11)) +
    scale_fill_manual(values = colors) +
    labs(y = "Estimand value",
         x = "Sample size setting",
         fill = "Validation approach") +
    facet_grid(rows = vars(estimand), cols = vars(model), scales = "free")+ 
    theme(legend.position="bottom")+
    guides(color = guide_legend(nrow=4, ncol=3))+
    theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))
  
  
  
  # ML AND CART
  p3 <- 
    ggplot(data = df_s3_long %>% 
             filter(model %in% c("ML","CART"), 
                    estimand %in% c("AUC", "Calib. Slope", "R2 Tjur")),
           mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
    geom_boxplot(position = position_dodge(width = 0.95), size = 0.005, outlier.shape = NA)+
    theme_set(theme_bw(base_size = 11)) +
    scale_fill_manual(values = colors) +
    labs(y = "Estimand value",
         x = "Sample size setting",
         fill = "Validation approach") +
    facet_grid(rows = vars(estimand), cols = vars(model), scales = "free")+ 
    theme(legend.position="bottom")+
    guides(color = guide_legend(nrow=4, ncol=3))+
    theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))
  

  # All estimands
  p4 <-ggplot(data = df_s3_long %>% 
                filter(model %in% c("ML", "CART")), 
              mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
    geom_boxplot(position = position_dodge(width = 0.95), size = 0.005, outlier.shape = NA) +
    theme_set(theme_bw(base_size = 11)) +
    scale_fill_manual(values = colors) +
    labs(y = "Estimand value",
         x = "Sample size setting",
         fill = "Validation approach") +
    facet_grid(rows = vars(estimand), cols = vars(model), scales = "free")+ 
    theme(legend.position="bottom")+
    guides(color = guide_legend(nrow=4, ncol=3))+
    theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))
  
  
  # ML and Firth 
  p5 <-ggplot(data = df_s3_long %>% 
                filter(model %in% c("ML", "Firth", "CART")), 
              mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
    geom_boxplot(position = position_dodge(width = 0.95), size = 0.005, outlier.shape = NA) +
    theme_set(theme_bw(base_size = 11)) +
    scale_fill_manual(values = colors) +
    labs(y = "Estimand value",
         x = "Sample size setting",
         fill = "Validation approach") +
    facet_grid(rows = vars(estimand), cols = vars(model), scales = "free")+ 
    theme(legend.position="bottom")+
    guides(color = guide_legend(nrow=4, ncol=3))+
    theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))
  
  
  # Ridge and LASSO
  p6 <-ggplot(data = df_s3_long %>% 
                filter(model %in% c("Ridge", "Lasso", "RF")), 
              mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
    geom_boxplot(position = position_dodge(width = 0.95), size = 0.005, outlier.shape = NA) +
    theme_set(theme_bw(base_size = 11)) +
    scale_fill_manual(values = colors) +
    labs(y = "Estimand value",
         x = "Sample size setting",
         fill = "Validation approach") +
    facet_grid(rows = vars(estimand), cols = vars(model), scales = "free")+ 
    theme(legend.position="bottom")+
    guides(color = guide_legend(nrow=4, ncol=3))+
    theme(plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))
  
 
  
  ### Firth, ridge, lasso and rf
  # original scale 
  ggsave(paste0(estimand_plots,"estimands_study_3_thesis.pdf"), plot = p1, width = 19, height = 17, units = "cm")
  # all estimands
  # ggsave(paste0(full_estimand_plots,"full_estimands_study_3_thesis.pdf"), plot = p2, width = 20, height = 20, units = "cm")
  
  
  ### ML and CART
  # original scale
  #ggsave(paste0(estimand_plots,"estimands_study_3_ML_CART.pdf"), plot = p3, width = 19, height = 15, units = "cm")
  # all estimands
  #ggsave(paste0(full_estimand_plots,"full_estimands_study_3_ML_CART.pdf"), plot = p4, width = 20, height = 20, units = "cm")
  
  ## ML, Firth and RF
  ggsave(paste0(full_estimand_plots,"S5_All_metrics_study_3_ML_Firth_CART.pdf"), plot = p5, width = 20, height = 20, units = "cm")
  
  ## Ridge, LASSO and CART
  ggsave(paste0(full_estimand_plots,"S6_All_metrics_study_3_Ridge_Lasso_RF.pdf"), plot = p6, width = 20, height = 20, units = "cm")
