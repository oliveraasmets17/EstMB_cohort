

# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load the packages
library("tidyverse")
library("stringr")
library("xlsx")

# Gather the outputs of the single ML models
gather_outputs <- FALSE
if (gather_outputs == TRUE){
  output_df <- data.frame()
  files_list <- list.files("RData/ML")
  files_list_clean <- files_list[str_detect(files_list, "_output")]
  for (i in files_list_clean){
    file = readRDS(file.path("RData/ML", i))
    output_df <- dplyr::bind_rows(output_df, file)
  }
  saveRDS(output_df, "RData archive/tidymodels_output_summarized.rds")
} else{
  output_df <- readRDS("RData archive/tidymodels_output_summarized.rds")
}







# Summarize and visualize modelling results ----
#------------------------------------------------#
#                                                #
#    SUMMARIZE AND VISUALIZE MODELLING REULTS    # 
#                                                #
#------------------------------------------------#

# Clean the output data
ML_output <- output_df %>% 
  dplyr::left_join(readRDS("RData archive/Factor_groups_names.rds"), by = c("analysis_factor" = "factor")) %>%   # Merge factor names - ICD names
  dplyr::mutate(predictorset_name = case_when(predictors == "SET0" ~ "MB only",
                                              predictors == "SET2" ~ "Null",
                                              predictors == "SET2_MB" ~ "Null + MB",
                                              predictors == "SET3" ~ "Null + AB history",
                                              predictors == "SET3_MB" ~ "Null + AB history + MB"),
                predictorset_name = factor(predictorset_name, levels = c("MB only", "Null", "Null + MB", "Null + AB history", "Null + AB history + MB")))


# Aggregate the results over 10 splits
output_aggregated <- ML_output %>% 
  dplyr::filter(set == "test") %>%   # Consider performances on the test set
  dplyr::group_by(factor_name, predictors) %>% 
  dplyr::summarize(mean_value = mean(value),
                   median_value = median(value), 
                   n = n()) %>% 
  dplyr::ungroup() 

# AUC values of the "baseline" model - model with age, BMI, gender and Bristol stool scale as predictors
baseline <- output_aggregated %>% 
  dplyr::filter(predictors == "SET2") %>% 
  dplyr::rename("baseline" = "mean_value") %>% 
  dplyr::select(factor_name, baseline)

# Apply arbitrary threshold to divide the diseases into groups - whether the microbiome had some additional effect on the preditions etc
n_p = 0.01
output_grouped <- output_aggregated %>% 
  dplyr::left_join(baseline, by = "factor_name") %>%
  dplyr::mutate(change = mean_value-baseline) %>% #(mean_value - baseline)/baseline*100) %>% 
  dplyr::select(factor_name, predictors, change) %>% 
  tidyr::spread(predictors, change) %>% 
  dplyr::mutate(change_group = case_when(SET2_MB < n_p & SET3_MB < n_p & SET3 < n_p ~ "No additional effect", 
                                         SET2_MB >= n_p & SET3_MB >= n_p & SET3 < n_p ~ "MB only effect",
                                         SET2_MB >= n_p & SET3_MB >= n_p & SET3_MB-SET3 >= n_p ~ "MB + AB effect",
                                         SET2_MB >= n_p & SET3_MB >= n_p & SET3 >= n_p ~ "MB effect covered in AB",
                                         SET2_MB < n_p & SET3_MB >= n_p & SET3 >= n_p ~ "AB usage effect", 
                                         SET2_MB < n_p & SET3_MB < n_p & SET3 >= n_p ~ "AB usage effect", 
                                         TRUE ~ "Something else"), 
                change_group = factor(change_group, levels = c("No additional effect","MB only effect", "AB usage effect", 
                                                               "MB effect covered in AB", "Something else")))


# Save the aggregated results for the prediction models
ML_output_final <- output_aggregated %>% 
  dplyr::select(factor_name, predictors, mean_value) %>% 
  tidyr::spread(predictors, mean_value) %>% 
  dplyr::select(factor_name, SET0, SET2, SET2_MB, SET3, SET3_MB) %>%
  dplyr::rename("Disease" = "factor_name", 
                "MB only" = "SET0",
                "Null" = "SET2",
                "Null+MB" = "SET2_MB",
                "Null+AB" = "SET3",
                "Null+MB+AB" = "SET3_MB") %>%
  as.data.frame() %>% 
  dplyr::arrange(desc(`MB only`))

write.xlsx(ML_output_final, file = "Results/Supplementary Table 7. Predictive modelling.xlsx", 
           sheetName = "AUROC", col.names = TRUE, row.names = FALSE, append = FALSE)

# Determine factor order for plotting
factor_order <- ML_output %>% 
  dplyr::left_join(output_grouped[ ,c("factor_name", "change_group")], by = "factor_name") %>% 
  dplyr::filter(set == "test" & predictors == "SET0") %>% 
  dplyr::group_by(factor_name, change_group) %>% 
  dplyr::summarise(mean_AUC = mean(value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(change_group, desc(mean_AUC))

# Visualize the modelling results
p1 <- ggplot(ML_output %>% dplyr::filter(set == "test"), 
             aes(x = factor(factor_name, levels = factor_order$factor_name), 
                 y = value, fill = predictorset_name)) + 
  geom_boxplot() + 
  geom_abline(intercept = 0.5, slope = 0, linetype = 1, size = 1, color = "black") + 
  geom_vline(xintercept = c(13.45, 14.6, 26.5), color = "gray60", linetype = 2, size = 1) + 
  annotate("label", x = 6.5, y = 0.98, label = "Group 1", size = 5) + 
  annotate("label", x = 14, y = 0.98, label = "Group 2", size = 5) + 
  annotate("label", x = 21, y = 0.98, label = "Group 3", size = 5) + 
  annotate("label", x = 27.2, y = 0.98, label = "Group 4", size = 5) + 
  scale_fill_manual(name = "Predictor set", values = c("gray30", "#fb6a4a", "#a50f15",
                                                       "#74a9cf","#045a8d")) + 
  ylab("Area Under the Curve (AUC)") + 
  xlab("") + 
  theme_bw() + 
  theme(legend.position = "right") +
  expand_limits(x = c(0, 28)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14)) 
p1

ggsave(plot = p1, filename = "Figures/LASSO_predictions_boxplot.png", width = 20, height = 10)
ggsave(plot = p1, filename = "Figures/LASSO_predictions_boxplot.svg", width = 20, height = 10)



