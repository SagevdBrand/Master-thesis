#################################
######## Visualization ##########
#################################
################
## Setting up ##
################
source("./src/setup.R")
scenarios <- readRDS(paste0(setting_path, "studies.RDS"))

## Loading the above estimand files:
df_all <- readRDS(paste0(estimands_general_path, "all_estimands_batch_6.RDS"))

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
                                                      TRUE ~ as.character(dim))
                                      )

############
## Colors ##
############

# Specify color palette:
colors <- c("white",   # White
            "#e63946", # imperial red
            "#02c39a", #Mountain meadow
            "#f28482", #Light coral
            "#84a59d", # Morning blue
            "#f6bd60", # Maximum yellow red
            "#457b9d", # Celadon blue
            "#f4a261" # Sandy brown
)

############################
############################
########## STUDY 1 #########
############################
############################

# Results of study 1:
df_s1_long<- df_all_long %>% filter(study == "Study_1")
df_s1_long$scenario <- as.numeric(gsub("Scenario_", "", df_s1_long$scenario))
df_s1_long <- df_s1_long%>% mutate(value = case_when(estimand == "log(Slope)" & value > 5 ~ 5,
                                                     estimand == "log(Slope)" & value < 0 ~ 0,
                                                     estimand == "log(Slope)" & is.infinite(value) ~ 5,
                                                     TRUE ~ value))

### Figures to be included in thesis
p1 <- 
  ggplot(data = df_s1_long %>% 
           filter(pred_selection == "none",
                  estimand %in% c("AUC", "rMSPE")),
         mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
  geom_boxplot(data = df_s1_long %>% 
                 filter(pred_selection == "none",
                        estimand == "log(Slope)"),
               mapping = aes(x = as.factor(n_setting), y =ifelse(value == 0, 0, log(value)), fill = approach),
               position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
  #coord_cartesian(ylim = c(-1, 2.5)) +
  geom_boxplot(position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
  theme_set(theme_bw(base_size = 11)) +
  scale_fill_manual(values = colors) +
  labs(y = "Estimand value",
       x = "Sample size setting",
       fill = "Validation approach"
       ) +
  facet_grid(rows = vars(estimand), cols = vars(prev), scales = "free")+ 
  theme(legend.position="bottom")+
  guides(color = guide_legend(nrow=4, ncol=3))

p2 <- 
  ggplot(data = df_s1_long %>% 
           filter(pred_selection == "<0.157",
                  estimand %in% c("AUC", "rMSPE")),
         mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
  geom_boxplot(data = df_s1_long %>% 
                 filter(pred_selection == "<0.157",
                        estimand == "log(Slope)"),
               mapping = aes(x = as.factor(n_setting), y =ifelse(value == 0, 0, log(value)), fill = approach),
               position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
  geom_boxplot(position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
  # stat_summary(fun=mean, geom="point",position = position_dodge(width = 0.95), shape=1, size= 1.5) +
  theme_set(theme_bw(base_size = 11)) +
  scale_fill_manual(values = colors) +
  labs(y = "Estimand value",
       x = "Sample size setting",
       fill = "Validation approach") +
  facet_grid(rows = vars(estimand), cols = vars(prev), scales = "free")+ 
  theme(legend.position="bottom")+
  guides(color = guide_legend(nrow=4, ncol=3))


# Figures of all estimands (For supplementary material)
p3 <- 
  ggplot(data = df_s1_long %>% 
           filter(pred_selection == "none", 
                  estimand != "log(Slope)"), 
         mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
  geom_boxplot(position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
  geom_boxplot(data = df_s1_long %>% 
                 filter(pred_selection == "none",
                        estimand == "log(Slope)"),
               mapping = aes(x = as.factor(n_setting), y =ifelse(value == 0, 0, log(value)), fill = approach),
               position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
  theme_set(theme_bw(base_size = 11)) +
  scale_fill_manual(values = colors) +
  labs(y = "Estimand value",
       x = "Sample size setting",
       fill = "Validation approach") +
  facet_grid(rows = vars(estimand), cols = vars(prev), scales = "free")+ 
  theme(legend.position="bottom")+
  guides(color = guide_legend(nrow=4, ncol=3))

p4 <- 
  ggplot(data = df_s1_long %>% 
           filter(pred_selection == "<0.157",
                  estimand != "log(Slope)"),
         mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
  geom_boxplot(position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
  geom_boxplot(data = df_s1_long %>% 
                 filter(pred_selection == "<0.157",
                        estimand == "log(Slope)"),
               mapping = aes(x = as.factor(n_setting), y =ifelse(value == 0, 0, log(value)), fill = approach),
               position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
  theme_set(theme_bw(base_size = 11)) +
  scale_fill_manual(values = colors) +
  labs(y = "Estimand value",
       x = "Sample size setting",
       fill = "Validation approach") +
  facet_grid(rows = vars(estimand), cols = vars(prev), scales = "free")+ 
  theme(legend.position="bottom")+
  guides(color = guide_legend(nrow=4, ncol=3))

## Saving the plots
assign(paste0("thesis_estimands_study_1_no_pred_sel"), p1)
ggsave(paste0(estimand_plots,"estimands_study_1_no_pred_sel.pdf"), plot = p1, width = 20, height = 25, units = "cm")
assign(paste0("thesis_estimands_study_1_pred_sel"), p2)
ggsave(paste0(estimand_plots,"estimands_study_1_pred_sel.pdf"), plot = p2, width = 20, height = 25, units = "cm")
ggsave(paste0(full_estimand_plots,"full_estimands_study_1_no_pred_sel.pdf"), plot = p3, width = 20, height = 25, units = "cm")
ggsave(paste0(full_estimand_plots,"full_estimands_study_1_pred_sel.pdf"), plot = p4, width = 20, height = 25, units = "cm")


############################
############################
########## STUDY 2 #########
############################
############################

# Results of study 2:
df_s2_long<- df_all_long %>% filter(study == "Study_2")
df_s2_long$scenario <- as.numeric(gsub("Scenario_", "", df_s2_long$scenario))
df_s2_long <- df_s2_long %>% mutate(value = case_when(estimand == "log(Slope)" & value > 5 ~ 5,
                                                      estimand == "log(Slope)" & value < 0 ~ 0,
                                                      estimand == "log(Slope)" & is.infinite(value) ~ 5,
                                                      TRUE ~ value))

  p1 <-  
    ggplot(data = df_s2_long %>%
             filter(pred_selection == "none",
                    estimand %in% c("AUC", "rMSPE", "eci")), 
           mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
    geom_boxplot(position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
    geom_boxplot(data = df_s2_long %>% 
                   filter(pred_selection == "none",
                          estimand == "log(Slope)"),
                 mapping = aes(x = as.factor(n_setting), y =ifelse(value == 0, 0, log(value)), fill = approach),
                 position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
    theme_set(theme_bw(base_size = 11)) +
    scale_fill_manual(values = colors) +
    labs(y = "Estimand value",
         x = "Sample size setting",
         fill = "Validation approach") +
    facet_grid(rows = vars(estimand), cols = vars(noise, dim), scales = "free")+ 
    theme(legend.position="bottom")+
    guides(color = guide_legend(nrow=4, ncol=3))
  
  
  p2 <- ggplot(data = df_s2_long %>%
                 filter(pred_selection == "<0.157",
                        estimand %in% c("AUC", "rMSPE", "eci")), 
               mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
    geom_boxplot(position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
    geom_boxplot(data = df_s2_long %>% 
                   filter(pred_selection == "<0.157",
                          estimand == "log(Slope)"),
                 mapping = aes(x = as.factor(n_setting), y =ifelse(value == 0, 0, log(value)), fill = approach),
                 position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
    theme_set(theme_bw(base_size = 11)) +
    scale_fill_manual(values = colors) +
    labs(y = "Estimand value",
         x = "Sample size setting",
         fill = "Validation approach") +
    facet_grid(rows = vars(estimand), cols = vars(noise, dim), scales = "free")+ 
    theme(legend.position="bottom")+
    guides(color = guide_legend(nrow=4, ncol=3))
  
  p3 <- ggplot(data = df_s2_long %>%
                 filter(pred_selection == "none",
                        estimand != "log(Slope)"), 
               mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
    geom_boxplot(position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
    geom_boxplot(data = df_s2_long %>% 
                   filter(pred_selection == "none",
                          estimand == "log(Slope)"),
                 mapping = aes(x = as.factor(n_setting), y =ifelse(value == 0, 0, log(value)), fill = approach),
                 position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
    theme_set(theme_bw(base_size = 11)) +
    scale_fill_manual(values = colors) +
    labs(y = "Estimand value",
         x = "Sample size setting",
         fill = "Validation approach") +
    facet_grid(rows = vars(estimand), cols = vars(noise, dim), scales = "free")+ 
    theme(legend.position="bottom")+
    guides(color = guide_legend(nrow=4, ncol=3))
  
  
  p4 <- ggplot(data = df_s2_long %>%
                 filter(pred_selection == "<0.157",
                        estimand != "log(Slope)"), 
               mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
    geom_boxplot(position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
    geom_boxplot(data = df_s2_long %>% 
                   filter(pred_selection == "<0.157",
                          estimand == "log(Slope)"),
                 mapping = aes(x = as.factor(n_setting), y =ifelse(value == 0, 0, log(value)), fill = approach),
                 position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
    theme_set(theme_bw(base_size = 11)) +
    scale_fill_manual(values = colors) +
    labs(y = "Estimand value",
         x = "Sample size setting",
         fill = "Validation approach",
         title = "Estimands for study 2") +
    facet_grid(rows = vars(estimand), cols = vars(noise, dim), scales = "free")+ 
    theme(legend.position="bottom")+
    guides(color = guide_legend(nrow=4, ncol=3))
  
  assign("thesis_estimands_study_2_no_pred_sel", p1)
  ggsave(paste0(estimand_plots,"estimands_study_2_no_pred_sel.pdf"), plot = p1, width = 20, height = 25, units = "cm")
  assign("thesis_estimands_study_2_pred_sel", p2)
  ggsave(paste0(estimand_plots,"estimands_study_2_pred_sel.pdf"), plot = p2, width = 20, height = 25, units = "cm")
  ggsave(paste0(full_estimand_plots,"full_estimands_study_2_no_pred_sel.pdf"), plot = p3, width = 20, height = 25, units = "cm")
  ggsave(paste0(full_estimand_plots,"full_estimands_study_2_pred_sel.pdf"), plot = p4, width = 20, height = 25, units = "cm")
  

############################
############################
########## STUDY 3 #########
############################
############################

# Results of study 2:
df_s3_long<- df_all_long %>% filter(study == "Study_3")
df_s3_long$scenario <- as.numeric(gsub("Scenario_", "", df_s3_long$scenario))
df_s3_long <- df_s3_long %>% mutate(value = case_when(estimand == "log(Slope)" & value > 5 ~ 5,
                                                      estimand == "log(Slope)" & value < 0 ~ 0,
                                                      estimand == "log(Slope)" & is.infinite(value) ~ 5,
                                                      TRUE ~ value))

  p1 <- 
    ggplot(data = df_s3_long %>% 
             filter(model %in% c("Firth", "Ridge", "Lasso"), 
                    !(estimand %in% c( "log(Slope)", "R2 Cox Snell", "MAPE")) ), 
           mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
    geom_boxplot(position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
    geom_boxplot(data = df_s2_long %>% 
                   filter(pred_selection == "<0.157",
                          model %in% c("Firth", "Ridge", "Lasso"), 
                          estimand == "log(Slope)"),
                 mapping = aes(x = as.factor(n_setting), y = ifelse(value == 0, 0, log(value)), fill = approach),
                 position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
    theme_set(theme_bw(base_size = 11)) +
    scale_fill_manual(values = colors) +
    labs(y = "Estimand value",
         x = "Sample size setting",
         fill = "Validation approach") +
    facet_grid(rows = vars(estimand), cols = vars(model), scales = "free")+ 
    theme(legend.position="bottom")+
    guides(color = guide_legend(nrow=4, ncol=3))
  
 
  
  p2 <-ggplot(data = df_s3_long %>% 
                filter(model %in% c("Firth", "Ridge", "Lasso"), 
                       estimand != "log(Slope)"), 
              mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
    geom_boxplot(position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
    geom_boxplot(data = df_s2_long %>% 
                   filter(pred_selection == "<0.157",
                          model %in% c("Firth", "Ridge", "Lasso"), 
                          estimand == "log(Slope)"),
                 mapping = aes(x = as.factor(n_setting), y = ifelse(value == 0, 0, log(value)), fill = approach),
                 position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
    theme_set(theme_bw(base_size = 11)) +
    scale_fill_manual(values = colors) +
    labs(y = "Estimand value",
         x = "Sample size setting",
         fill = "Validation approach") +
    facet_grid(rows = vars(estimand), cols = vars(model), scales = "free")+ 
    theme(legend.position="bottom")+
    guides(color = guide_legend(nrow=4, ncol=3))
  
  
  
  p3 <- 
    ggplot(data = df_s3_long %>% 
             filter(model %in% c("ML", "CART", "RF"), 
                    !(estimand %in% c("log(Slope)", "R2 Cox Snell", "MAPE")) ), 
           mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
    geom_boxplot(position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
    geom_boxplot(data = df_s2_long %>% 
                   filter(pred_selection == "<0.157",
                          model %in% c("ML", "CART", "RF"), 
                          estimand == "log(Slope)"),
                 mapping = aes(x = as.factor(n_setting), y = ifelse(value == 0, 0, log(value)), fill = approach),
                 position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
    theme_set(theme_bw(base_size = 11)) +
    scale_fill_manual(values = colors) +
    labs(y = "Estimand value",
         x = "Sample size setting",
         fill = "Validation approach") +
    facet_grid(rows = vars(estimand), cols = vars(model), scales = "free")+ 
    theme(legend.position="bottom")+
    guides(color = guide_legend(nrow=4, ncol=3))
  
  p4 <-ggplot(data = df_s3_long %>% 
                filter(model %in% c("ML", "CART", "RF"), 
                       estimand  != "log(Slope)"), 
              mapping = aes(x = as.factor(n_setting), y = value, fill = approach)) +
    geom_boxplot(position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
    geom_boxplot(data = df_s2_long %>% 
                   filter(pred_selection == "<0.157",
                          model %in% c("ML", "CART", "RF"), 
                          estimand == "log(Slope)"),
                 mapping = aes(x = as.factor(n_setting), y = ifelse(value == 0, 0, log(value)), fill = approach),
                 position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
    theme_set(theme_bw(base_size = 11)) +
    scale_fill_manual(values = colors) +
    labs(y = "Estimand value",
         x = "Sample size setting",
         fill = "Validation approach") +
    facet_grid(rows = vars(estimand), cols = vars(model), scales = "free")+ 
    theme(legend.position="bottom")+
    guides(color = guide_legend(nrow=4, ncol=3))
  
  
  assign("thesis_estimands_study_3_penalized", p1)
  ggsave(paste0(estimand_plots,"estimands_study_3_penalized.pdf"), plot = p1, width = 20, height = 25, units = "cm")
  ggsave(paste0(full_estimand_plots,"full_estimands_study_3_penalized.pdf"), plot = p2, width = 20, height = 20, units = "cm")
  assign("full_estimands_study_3_treebased", p3)
  ggsave(paste0(estimand_plots,"estimands_study_3_treebased.pdf"), plot = p3, width = 20, height = 25, units = "cm")
  ggsave(paste0(full_estimand_plots,"full_estimands_study_3_treebased.pdf"), plot = p4, width = 20, height = 20, units = "cm")
