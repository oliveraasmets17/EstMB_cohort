


# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load libraries
library("tidyverse")
library("viridis")
library("ggpubr")






# Visualize diversity analysis results ----
#------------------------------------------------#
#                                                #
#               VISUALIZING RESULTS              # 
#                                                #
#------------------------------------------------#

# Define the function
# @input - PERMANOVA_allFactors_file - beta-diversity analysis results file
# @input - alpha_diversity_associations_file - alpha-diversity analysis results file
# @input - output_file - output file name for the figures
# Function saves the diversity figures

visualize_diversity = function(PERMANOVA_allFactors_file, alpha_diversity_associations_file, output_file){
  
  # Read results files
  PERMANOVA_output_all <- readRDS(PERMANOVA_allFactors_file) %>% 
    dplyr::mutate(FDR = p.adjust(`Pr(>F)`, method = "BH")) %>% 
    dplyr::select(analysis_factor, Df, SumsOfSqs, MeanSqs, F.Model, R2, `Pr(>F)`, adjusted, FDR) 
  
  alpha_diversity_summary_output <- readRDS(alpha_diversity_associations_file) %>%
    dplyr::filter(method == "Spearman" & adjusted == "NONE")
  
  # Read help files
  factor_groups <- readRDS(file = "RData archive/Factor_groups_names.rds")

  # Define factors to plot - at least one significant association with either alpha or beta diversity
  Has_alpha_associations <- alpha_diversity_summary_output %>% 
    dplyr::group_by(method, diversity_metric) %>% 
    dplyr::mutate(FDR = p.adjust(p.value, method = "BH")) %>% 
    dplyr::ungroup() %>%
    dplyr::filter(FDR <= 0.05) %>% 
    dplyr::distinct(analysis_factor) %>% 
    dplyr::pull(analysis_factor)
  
  Has_beta_associations <- PERMANOVA_output_all %>%
    dplyr::filter(FDR <= 0.05) %>% 
    dplyr::distinct(analysis_factor) %>%
    dplyr::pull(analysis_factor)
  
  # Consider only continuous measures of antibiotics and antidepressants usage
  Has_associations <- setdiff(union(Has_alpha_associations, Has_beta_associations),
                              c("antidepressants_history_quartile", "antibiotics_history_quartile"))
  
  # Arrange factors according to R2
  factor_ordering <- PERMANOVA_output_all %>%
    dplyr::left_join(factor_groups, by = c("analysis_factor" = "factor")) %>%
    dplyr::filter(analysis_factor %in% Has_associations) %>% 
    dplyr::mutate(factor_group = factor(factor_group, levels = c("Procedures", "Medications", "Diseases", 
                                                                 "Dietary", "Intrinsic", "Other"))) %>%
    dplyr::arrange(factor_group, R2) %>%
    dplyr::pull(factor_name)
  
  # Prepare data for plotting
  alpha_diversity_plotData <- alpha_diversity_summary_output %>%
    dplyr::group_by(method, adjusted, diversity_metric) %>% 
    dplyr::mutate(FDR = p.adjust(p.value, method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::filter(analysis_factor %in% Has_associations) %>% 
    dplyr::left_join(factor_groups, by = c("analysis_factor" = "factor")) %>% 
    dplyr::filter(is.nan(estimate) == FALSE) %>%
    dplyr::mutate(factor_name = factor(factor_name, levels = factor_ordering), 
                  setting = ifelse(diversity_metric == "diversity_observed", "Observed", "Shannon"),
                  setting = factor(setting, levels = c("Shannon", "Observed"))) 
  
  PERMANOVA_plotData <- PERMANOVA_output_all  %>%
    dplyr::left_join(factor_groups, by = c("analysis_factor" = "factor")) %>%
    dplyr::filter(analysis_factor %in% Has_associations) %>% 
    dplyr::mutate(factor_name = factor(factor_name, levels = factor_ordering))
  
  # Visualize Alpha diversity output - Spearman
  p_Alpha_Spearman <- ggplot(alpha_diversity_plotData, 
                             aes(x = factor_name, y = setting, fill = estimate)) + 
    geom_tile(color = "black", size = 1, width = 0.84) + 
    geom_text(aes(label = paste(round(estimate, 2), ifelse(FDR <= 0.05, "*", ""), sep = "")), size = 4.5, vjust = 0.45) + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-0.28,0.28), guide = F) +  
    theme_classic() + 
    theme(axis.ticks = element_blank(),
          axis.line = element_blank(), 
          axis.text.x = element_text(hjust = 0.5, vjust = 0.25),
          text = element_text(size = 18),
          plot.margin=unit(c(-0.1,1,1,1), "cm")) + 
    xlab("") +
    ylab ("") + 
    scale_y_discrete(position = "right")  + 
    coord_flip()
  
  # Visualize PERMANOVA output
  p_PERMANOVA <- ggplot(PERMANOVA_plotData, aes(x = factor_name, y = R2*100, fill = factor_group)) + 
    geom_bar(stat = "identity", color = "black", size = 1, width = 0.84) + 
    geom_text(aes(y = R2*100 + 0.01, label = ifelse(FDR <= 0.05, "*", "")), size = 7, vjust = 0.75) + 
    scale_fill_manual(name = "", values = c("#8ed0d1ff",  "#282554ff", "#de2524ff",
                                            "#ffcb4dff", "#6b9064ff","#ddddd6ff")) + 
    theme_classic() +
    scale_y_continuous(position = "right", expand = c(0.02, 0)) +
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          legend.position = c(0.6, 0.03),
          legend.text = element_text(size = 12),
          legend.direction = "horizontal",
          text = element_text(size = 18),
          axis.line.x = element_line(),
          plot.margin=unit(c(1,1,-0.5,1), "cm")) + 
    xlab("") +
    ylab("") + 
    coord_flip()
  
  # Merge plots
  p_overall <- ggpubr::ggarrange(p_Alpha_Spearman, p_PERMANOVA, ncol = 2, align = "h", widths = c(2, 1))
  p_overall
  
  ggsave(plot = p_overall, filename = paste(output_file, ".png", sep = ""), height = 20, width = 15)    
  ggsave(plot = p_overall, filename = paste(output_file, ".svg", sep = ""), height = 20, width = 15)    
  
}


# Apply the function
visualize_diversity(PERMANOVA_allFactors_file = "RData archive/PERMANOVA_KEGGko_allFactors_NONE.rds",
                    alpha_diversity_associations_file = "RData archive/Diversity_associations_KEGGko.rds",
                    output_file = "Figures/Diversity_overall_KEGGko")

visualize_diversity(PERMANOVA_allFactors_file = "RData archive/PERMANOVA_species_allFactors_NONE.rds",
                    alpha_diversity_associations_file = "RData archive/Diversity_associations_species.rds",
                    output_file = "Figures/Diversity_overall_species")

