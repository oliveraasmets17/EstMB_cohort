



# Command line arguments ----
#------------------------------------------------#
#                                                #
#          SET UP COMMAND LINE ARGUMENTS         # 
#                                                #
#------------------------------------------------#

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)             
# Expected input:
# args = c(distance_matrix_file, metadata_file, factors_analyzed_file, output_file, included_covariates)
# args = c("RData/Abundance_species.rds", "RData/Data_master.rds", "RData/Disease_factors_analyzed.rds", "RData/Test_MWAS_species_NONE.rds", "SET2") 

cat(args, sep = "\n")

# Test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file)", call. = FALSE)
} 

# Declare count matrix input
count_matrix_file <- args[1]

# Declare phenotype data used by model - dataframe
metadata_file <- args[2]

# Declare factors to be analyzed - character vector
factors_analyzed_file <- args[3]

# Declare output file name
output_file <- args[4]

# Declare which covariates to use - character
include_covariates <- args[5]






# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load packages
library("dplyr")
library("ALDEx2")
library("tibble")
library("stringr")

# Read count matrix used as input (genes, species etc)
count_data <- readRDS(file = count_matrix_file)  %>%  # sample x feature
  t() %>%
  as.data.frame()

# Read the phenotype data
full_factor_data <- readRDS(file = metadata_file) %>% 
  as.data.frame() 

# Read factors analyzed
factors_analyzed <- readRDS(file = factors_analyzed_file)

# Read skood vkood link data
skood_vkood_link <- readRDS("RData/Skood_Vkood_link.rds") 

# Antibiotics users (last 180 days)
q_antibiotic_users <- readRDS("RData/q_antibiotic_users.rds")

# Samples with extremely low reads to be filtered out
low_read_scodes <- readRDS("RData/low_read_scodes.rds")

# Merge data containing scodes - linking with registry data
corresponding_scodes <- data.frame(vkood = rownames(count_data), stringsAsFactors = FALSE) %>%
  dplyr::mutate(vkood_new = str_replace(vkood, coll("."), "_"),
                vkood_new = ifelse(substr(vkood_new, nchar(vkood_new)-1, nchar(vkood_new)) == "_1", 
                                   substr(vkood_new, 1, nchar(vkood_new)-2), vkood_new)) %>% 
  dplyr::left_join(skood_vkood_link, by = c("vkood_new" = "vkood")) %>%
  dplyr::filter(complete.cases(.)) 

# Subset the count matrix and rename with scodes
count_data_scoded <- count_data %>% 
  dplyr::select(-contains("Other")) %>%
  tibble::rownames_to_column(var = "vkood") %>%
  dplyr::left_join(corresponding_scodes[ ,c("vkood", "skood")], by = "vkood") %>%
  dplyr::select(-one_of("vkood")) %>%
  dplyr::filter(is.na(skood) == FALSE) 






# Function that outputs ALDEx2 output for the ICD10 code ----
#------------------------------------------------#
#                                                #
#       DIFFERENTIAL ABUNDANCE BETWEEN MGC       # 
#                                                #
#------------------------------------------------#

# Define the function
# @input - factor_chosen - character name of the column in the phenotype dataframe 
# @input - use_covariates - character indicationg, which covariates should be used. Options:
#                           NONE - no covariates
#                           SET1 - age, BMI and gender as covariates
#                           SET2 - age, BMI, gender and Bristol stool scale as covariates
#                           SET3 - age, BMI, gender, Bristol stool scale and history of antibiotics usage as covariates
# Function returns data frame of the modelling results

find_differentially_abundant_features = function(factor_chosen, use_covariates){
  
  # Subset complete observations for the given factor
  if (use_covariates == "NONE"){
    
    metadata = full_factor_data %>%
      dplyr::mutate(analysis_factor = full_factor_data[, factor_chosen]) %>%
      dplyr::filter(!(skood %in% q_antibiotic_users)) %>%   # Exclude current antibiotics users (482)
      dplyr::filter(!(skood %in% low_read_scodes)) %>%   # Exclude samples with extremely low reads (3)
      dplyr::select(skood, analysis_factor) %>%
      dplyr::filter(complete.cases(.))
    
  } else if (use_covariates == "SET1"){
    
    metadata = full_factor_data %>%
      dplyr::mutate(analysis_factor = full_factor_data[ ,factor_chosen]) %>%
      dplyr::filter(!(skood %in% q_antibiotic_users)) %>%   # Exclude current antibiotics users (482)
      dplyr::filter(!(skood %in% low_read_scodes)) %>%   # Exclude samples with extremely low reads (3)
      dplyr::select(skood, analysis_factor, gender, Age_at_MBsample, BMI) %>%
      dplyr::filter(complete.cases(.))
    
  } else if (use_covariates == "SET2"){
    
    metadata = full_factor_data %>%
      dplyr::mutate(analysis_factor = full_factor_data[ ,factor_chosen]) %>%
      dplyr::filter(!(skood %in% q_antibiotic_users)) %>%   # Exclude current antibiotics users (482)
      dplyr::filter(!(skood %in% low_read_scodes)) %>%   # Exclude samples with extremely low reads (3)
      dplyr::select(skood, analysis_factor, gender, Age_at_MBsample, BMI, 
                    usualStoolType_category) %>%
      dplyr::filter(complete.cases(.))
    
  } else if (use_covariates == "SET3"){
    
    metadata = full_factor_data %>%
      dplyr::mutate(analysis_factor = full_factor_data[ ,factor_chosen]) %>%
      dplyr::filter(!(skood %in% q_antibiotic_users)) %>%   # Exclude current antibiotics users (264)
      dplyr::filter(!(skood %in% low_read_scodes)) %>%   # Exclude samples with extremely low reads (3)
      dplyr::select(skood, analysis_factor, gender, Age_at_MBsample, BMI, 
                    usualStoolType_category, antibiotics_history_cont) %>%
      dplyr::filter(complete.cases(.))
  } 
  
  # Merge metadata with count data
  full_data <- metadata %>% 
    dplyr::left_join(count_data_scoded, by = "skood")
  
  # Create formulas for the analysis
  if (use_covariates == "NONE"){
    
    if(is.numeric(metadata$analysis_factor) == TRUE){
      formula = as.formula(paste(" ~ analysis_factor", sep = ""))
    }else{
      formula = as.formula(paste(" ~ factor(analysis_factor)", sep = ""))
    }
    start_index = 3
    
  } else if(use_covariates == "SET1"){
    
    if(is.numeric(metadata$analysis_factor) == TRUE){
      formula = as.formula(paste(" ~ factor(gender) + Age_at_MBsample + BMI + analysis_factor", sep = ""))
    }else{
      formula = as.formula(paste(" ~ factor(gender) + Age_at_MBsample + BMI + factor(analysis_factor)", sep = ""))
    }
    start_index = 6
    
  } else if(use_covariates == "SET2"){
    
    if(is.numeric(metadata$analysis_factor) == TRUE){
      formula = as.formula(paste(" ~ factor(gender) + Age_at_MBsample + BMI + 
                               factor(usualStoolType_category) + analysis_factor", sep = ""))
    }else{
      formula =  as.formula(paste(" ~ factor(gender) + Age_at_MBsample + BMI + 
                               factor(usualStoolType_category) + factor(analysis_factor)", sep = ""))
    }
    start_index = 7
  }  else if(use_covariates == "SET3"){
    
    if(is.numeric(metadata$analysis_factor) == TRUE){
      formula = as.formula(paste(" ~ factor(gender) + Age_at_MBsample + BMI + 
                               factor(usualStoolType_category) + antibiotics_history_cont + analysis_factor", sep = ""))
    }else{
      formula =  as.formula(paste(" ~ factor(gender) + Age_at_MBsample + BMI + 
                               factor(usualStoolType_category) + antibiotics_history_cont + factor(analysis_factor)", sep = ""))
    }
    start_index = 8
  } 
  
  # Run ALDEx2 to find differentially abundant microbialfeatures
  n_cases <- metadata %>% 
    dplyr::filter(analysis_factor == 1) %>% 
    nrow()
  
  if (n_cases == 0){
    output <- NULL
  } else{
    model_matrix <- model.matrix(formula, full_data)
    set.seed(1)
    
    suppressMessages( 
      aldex_clr <- aldex.clr(t(full_data[ ,start_index:ncol(full_data)]), 
                             model_matrix, 
                             mc.samples = 128, 
                             denom = "all")
    )
    
    aldex_glm <- aldex.glm(aldex_clr)
    
    # Clean output
    output <- aldex_glm %>% 
      tibble::rownames_to_column(var = "feature") %>%
      dplyr::mutate(analysis_factor = factor_chosen, 
                    adjusted = use_covariates)
  }
  return(output)
}






# Apply the function for all factors analyzed ----
#------------------------------------------------#
#                                                #
#            ANALYZE MULTIPLE FACTORS            # 
#                                                #
#------------------------------------------------#
# Run MWAS on all the codes chosen
MWAS_summary_output <- data.frame()

run_all_MGC <- TRUE
if (run_all_MGC == TRUE){
  for (i in factors_analyzed){
    output_run = find_differentially_abundant_features(factor_chosen = i, use_covariates = include_covariates)
    MWAS_summary_output = rbind(MWAS_summary_output, output_run)
  }
  saveRDS(MWAS_summary_output, file = output_file)
}


