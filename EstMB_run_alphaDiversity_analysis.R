

# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load necessary packages
library("dplyr")
library("ppcor")
library("stringr")

# Read skood vkood link data - for linking with EHR
skood_vkood_link <- readRDS("RData/Skood_Vkood_link.rds") 

# Samples with extremely low reads to be filtered out
low_read_scodes <- readRDS("RData/low_read_scodes.rds")

# Read the phenotype data
full_factor_data <- readRDS(file = "RData/Data_master.rds") %>% 
  as.data.frame()






# Function that outputs spearman correlations with alpha diversity for the factor analyzed ----
#------------------------------------------------#
#                                                #
#        CORRELATIONS WITH ALPHA DIVERSITY       # 
#                                                #
#------------------------------------------------#

# Define the function
# @input - factor_chosen - character name of the column in the phenotype dataframe 
# @input - diversity_data - dataframe countaining phenotype data and alpha-diversity metrics
# Function returns data frame with Spearman correlation estimates, p-values and relevant metadata

find_diversity_associations = function(factor_chosen, diversity_data){

  # Create temporary variable for analysis
  metadata <- diversity_data %>%
    dplyr::mutate(analysis_factor = diversity_data[ ,factor_chosen]) 
  
  # Diversity metrics
  diversity_factors <- metadata %>% dplyr::select(contains("diversity")) %>% colnames(.)
  
  # Loop over all diversity metrics
  run_output <- data.frame()
  for (i in diversity_factors){
    
    # Run Spearman correlation
    if (!(factor_chosen %in% c("samplingSeason", "eatingHabit_name", "antibiotics_history_quartile", "antidepressants_history_quartile"))){
      Spearman <- cor(x = metadata %>% dplyr::pull(factor_chosen), 
                      y = metadata %>% dplyr::pull(i), 
                      method = "spearman", use = "complete.obs")
      
      Spearman_p <- cor.test(x = metadata %>% dplyr::pull(factor_chosen), 
                             y = metadata %>% dplyr::pull(i), 
                             method = "spearman", use = "complete.obs")$p.value
    }else{
      Spearman = as.numeric(NA)
      Spearman_p = as.numeric(NA)
    }  
    
    # Run Spearman partial correlation - age, sex and BMI
    if (!(factor_chosen %in% c("gender", "Age_at_MBsample", "BMI", "samplingSeason", "eatingHabit_name"))){
      metadata_clean <- metadata %>%
        dplyr::select(all_of(factor_chosen), gender, Age_at_MBsample, BMI, all_of(i)) %>%
        dplyr::filter(complete.cases(.))
      
      Spearman_partial1 <- spcor.test(x = metadata_clean %>% dplyr::pull(factor_chosen), 
                                      y = metadata_clean %>% dplyr::pull(i), 
                                      z = metadata_clean %>% dplyr::select(gender, Age_at_MBsample, BMI), 
                                      method = "spearman") %>%
        dplyr::mutate(method = "Spearman", 
                      adjusted = "SET1") %>%
        dplyr::select(estimate, p.value, method, adjusted)
    } else{
      Spearman_partial1 = NULL
    }  
    
    # Run Spearman partial correlation
    if (!(factor_chosen %in% c("gender",  "Age_at_MBsample", "BMI", "usualStoolType_category",
                               "samplingSeason", "eatingHabit_name"))){
      metadata_clean <- metadata %>%
        dplyr::select(all_of(factor_chosen), gender, Age_at_MBsample, usualStoolType_category, BMI, all_of(i)) %>%
        dplyr::filter(complete.cases(.))
      
      Spearman_partial2 <- spcor.test(x = metadata_clean %>% dplyr::pull(factor_chosen), 
                                      y = metadata_clean %>% dplyr::pull(i), 
                                      z = metadata_clean %>% dplyr::select(gender, usualStoolType_category, Age_at_MBsample, 
                                                                           BMI), 
                                      method = "spearman") %>%
        dplyr::mutate(method = "Spearman", 
                      adjusted = "SET2") %>%
        dplyr::select(estimate, p.value, method, adjusted)
    } else{
      Spearman_partial2 = NULL
    }
    
    # Output the results
    loop_output <- data.frame(estimate = Spearman, 
                              p.value = Spearman_p, 
                              method = "Spearman", 
                              adjusted = "NONE",
                              stringsAsFactors = FALSE) %>%
      dplyr::bind_rows(Spearman_partial1, Spearman_partial2) %>%
      dplyr::mutate(diversity_metric = i, 
                    analysis_factor = factor_chosen)
    
    run_output = rbind(run_output, loop_output)
  }
  return(run_output)
}





# Function that applies the previous function to all factors analyzed ----
#------------------------------------------------#
#                                                #
#            ANALYZE MULTIPLE FACTORS            # 
#                                                #
#------------------------------------------------#

# Define the function
# @input - diversity_data_file - character name of file that contains alpha-diversity metrics (dataframe)
# @input - factors_analyzed_file - character name of file that contains the names of the factors to be analyzed (vector)
# @input - output_file - character name of the output file name
# Function saves the final output

find_diversity_associated_factors = function(diversity_data_file, factors_analyzed_file, output_file){
  
  # Read diversity indices data
  diversity <- readRDS(diversity_data_file)
  
  # Read the factors analyzed
  factors_analyzed <- readRDS(factors_analyzed_file)
  
  # Merge data containing scodes - linking with registry data
  diversity_data <- diversity %>%
    dplyr::mutate(vkood_new = str_replace(vkood, coll("."), "_"),
                  vkood_new = ifelse(substr(vkood_new, nchar(vkood_new)-1, nchar(vkood_new)) == "_1", 
                                     substr(vkood_new, 1, nchar(vkood_new)-2), vkood_new)) %>% 
    dplyr::left_join(skood_vkood_link, by = c("vkood_new" = "vkood")) %>%
    dplyr::filter(!(skood %in% low_read_scodes)) %>% 
    dplyr::select(skood, contains("diversity")) %>%
    dplyr::filter(complete.cases(.)) %>%
    dplyr::left_join(full_factor_data, by = "skood") %>%
    # Recoding factor levels
    dplyr::mutate(usualStoolType_category = as.numeric(factor(usualStoolType_category, levels = c("Consipation", "Normal", "Diarrhea"))),
                  alcohol_frequency_category = as.numeric(factor(alcohol_frequency_category, levels = c("Low", "Average", "High"))))
  
  # Run the function for all factors
  diversity_summary_output <- data.frame()
  
  for (i in factors_analyzed){
    output_run = find_diversity_associations(factor_chosen = i, diversity_data = diversity_data)
    diversity_summary_output = rbind(diversity_summary_output, output_run)
  }
  
  # Save the output
  saveRDS(diversity_summary_output, file = output_file)
}


# Apply the function
find_diversity_associated_factors(diversity_data_file = "RData/Diversity_KEGGko.rds", 
                                  factors_analyzed_file = "RData/All_factors_analyzed_alpha.rds",
                                  output_file = "RData/Test_diversity_associations_KEGGko.rds")

find_diversity_associated_factors(diversity_data_file = "RData/Diversity_species.rds", 
                                  factors_analyzed_file = "RData/All_factors_analyzed_alpha.rds",
                                  output_file = "RData/Test_diversity_associations_species.rds")






