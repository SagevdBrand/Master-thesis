##########################################################
############## Check which results are present ###########
##########################################################
source("./src/setup.R")

## Loading the separate datafiles and binding them togetgher saving it for easy reference 
# df_all <- list.files(path = estimands_path, pattern = "*.Rds", full.names = T) %>%
#   map_dfr(readRDS)
# 
# Saving it for easy reference 
# saveRDS(df_all, file = paste(estimands_path, "all_estimands_batch_1.RDS"))

### Loading the above created file
df_all <- readRDS(paste(estimands_path, "all_estimands_batch_1.RDS"))

# How many iterations have been analyzed?
counts <- df_all %>% group_by(study, scenario) %>% summarise(count = n())

# Check results of study 1:
df_s1 <- df_all %>% filter(study == "Study_1")

df_s1 <- df_s1 %>% pivot_longer(., cols = estimands_names, names_to = "estimand", values_to = "value")
df_s1$value <- as.numeric(df_s1$value)

colors <- c("#02c39a", #Mountain meadow
            "#f28482", #Light coral
            "#84a59d", # Morning blue
            "#f6bd60", # Maximum yellow red
            "white",   # White
            "#457b9d", # Celadon blue
            "#f4a261", # Sandy brown
            "#e63946" # imperial red
            )

s1_s1_4_7 <- df_s1 %>% filter(scenario %in% c("Scenario_1", "Scenario_4", "Scenario_7")) %>% dplyr::select(c(estimand, study, scenario,model, approach, value, prev)) 
s1_s1_4_7 %>% ggplot(data = ., mapping = aes(x = as.factor(prev), y = value, fill = approach)) +
  geom_boxplot(position = position_dodge(width = 0.95), size = 0.005, outlier.size = 0.00) +
  theme_set(theme_bw(base_size = 11)) +
  scale_fill_manual(values = colors) +
  labs(y = " ",
       x = "Prevalence",
       fill = "Validation approach") +
  facet_wrap(~estimand, nrow = 3, scales = "free") +
  theme(legend.position = c(0.8, 0.12))


