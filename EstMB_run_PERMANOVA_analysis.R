

# Command line arguments ----
#------------------------------------------------#
#                                                #
#          SET UP COMMAND LINE ARGUMENTS         # 
#                                                #
#------------------------------------------------#

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)     
# Expected input:
# args = c(distance_matrix_file, metadata_file, factors_analyzed_file, output_file, use_covariates)
# args = c("RData/CLR_distance_species_clean.rds", "RData/Data_master.rds", "RData/Disease_factors_analyzed.rds", "RData/Test_PERMANOVA_species_NONE.rds", "NONE") 

cat(args, sep = "\n")

# Test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file)", call. = FALSE)
} 

# Declare count matrix input - dist object
distance_matrix_file <- args[1]

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
library("vegan")
library("stringr")

# Read distance object to be used as input (genes, species etc)
distance_matrix_raw <- readRDS(file = distance_matrix_file) # subject x subject

# Read the phenotype data. Data frame. 
full_factor_data <- readRDS(file = metadata_file) %>% 
  as.data.frame() 

# Read factors to be analyzed. Character vector of column names of phenotype data.
factors_analyzed <- readRDS(file = factors_analyzed_file)

# Read skood vkood link data - for linking with EHR
skood_vkood_link <- readRDS("RData/Skood_Vkood_link.rds") 

# Samples with extremely low reads to be filtered out
low_read_scodes <- readRDS("RData/low_read_scodes.rds")





# Work with the distance matrix ----
#------------------------------------------------#
#                                                #
#           PREPROCESS DISTANCE MATRIX           # 
#                                                #
#------------------------------------------------#

# Create matrix object from the dist object
distance_matrix_mat <- distance_matrix_raw %>%
  as.matrix() 

# Merge data containing scodes - linking with registry data
corresponding_scodes <- data.frame(vkood = colnames(distance_matrix_mat), stringsAsFactors = FALSE) %>%
  dplyr::mutate(vkood_new = str_replace(vkood, fixed("."), "_"),
                vkood_new = ifelse(substr(vkood_new, nchar(vkood_new)-1, nchar(vkood_new)) == "_1", 
                                   substr(vkood_new, 1, nchar(vkood_new)-2), vkood_new)) %>% 
  dplyr::left_join(skood_vkood_link, by = c("vkood_new" = "vkood")) %>%
  dplyr::filter(complete.cases(.)) 

# Subset the distance matrix and rename with scodes
distance_matrix_scoded <- distance_matrix_mat[corresponding_scodes$vkood, corresponding_scodes$vkood]
colnames(distance_matrix_scoded) <- rownames(distance_matrix_scoded) <- corresponding_scodes$skood






# Function that runs PERMANOVA for the factor of interest ----
#------------------------------------------------#
#                                                #
#       PERMANOVA FOR COMMUNITY DIFFERENCES      # 
#                                                #
#------------------------------------------------#

# Define the function
# @input - factor_chosen - character name of the column in the phenotype dataframe 
# @input - use_covariates - character indicationg, which covariates should be used. Options:
#                           NONE - no covariates
#                           SET1 - age, BMI and gender as covariates
#                           SET2 - age, BMI, gender and Bristol stool scale as covariates
# Function returns PERMANOVA aov table

find_MB_associated_factors = function(factor_chosen, use_covariates){
  
  # Subset complete observations for the given factor
  if (use_covariates == "NONE"){
    
    metadata = full_factor_data %>%
      dplyr::mutate(analysis_factor = full_factor_data[, factor_chosen]) %>%
      dplyr::filter(!(skood %in% low_read_scodes)) %>%   # Exclude samples with extremely low reads (3)
      dplyr::select(skood, analysis_factor) %>%
      dplyr::filter(complete.cases(.))
    
  } else if (use_covariates == "SET1"){
    
    metadata = full_factor_data %>%
      dplyr::mutate(analysis_factor = full_factor_data[ ,factor_chosen]) %>%
      dplyr::filter(!(skood %in% low_read_scodes)) %>%   # Exclude samples with extremely low reads (3)
      dplyr::select(skood, analysis_factor, gender, Age_at_MBsample, BMI) %>%
      dplyr::filter(complete.cases(.))
    
  } else if (use_covariates == "SET2"){
    
    metadata = full_factor_data %>%
      dplyr::mutate(analysis_factor = full_factor_data[ ,factor_chosen]) %>%
      dplyr::filter(!(skood %in% low_read_scodes)) %>%   # Exclude samples with extremely low reads (3)
      dplyr::select(skood, analysis_factor, gender, Age_at_MBsample, BMI, 
                    alcohol_frequency_category, usualStoolType_category) %>%
      dplyr::filter(complete.cases(.))
    
  } 
  
  # Subset distance matrix - keep rows and cols without missing values
  distance_matrix_subset <- distance_matrix_scoded[metadata$skood, metadata$skood]
  
  # Create formulas for the analysis
  if (use_covariates == "NONE"){
    
    if(is.numeric(metadata$analysis_factor) == TRUE){
      formula = as.formula(paste("distance_matrix_subset ~ analysis_factor", sep = ""))
    }else{
      formula = as.formula(paste("distance_matrix_subset ~ factor(analysis_factor)", sep = ""))
    }
    
  } else if(use_covariates == "SET1"){
    
    if(is.numeric(metadata$analysis_factor) == TRUE){
      formula = as.formula(paste("distance_matrix_subset ~ factor(gender) + Age_at_MBsample + BMI + analysis_factor", sep = ""))
    }else{
      formula = as.formula(paste("distance_matrix_subset ~ factor(gender) + Age_at_MBsample + BMI + factor(analysis_factor)", sep = ""))
    }
    
  } else if(use_covariates == "SET2"){
    if(is.numeric(metadata$analysis_factor) == TRUE){
      formula = as.formula(paste("distance_matrix_subset ~ factor(gender) + Age_at_MBsample + BMI + 
                               factor(usualStoolType_category) + analysis_factor", sep = ""))
    }else{
      formula =  as.formula(paste("distance_matrix_subset ~ factor(gender) + Age_at_MBsample + BMI + 
                               factor(usualStoolType_category) + factor(analysis_factor)", sep = ""))
    }
  } 
  
  # Run PERMANOVA
  set.seed(1)
  permanova <- vegan::adonis(formula, 
                             data = metadata,
                             permutations = 10000, 
                             parallel = 20)
  
  # Clean output
  permanova_output <- permanova$aov.tab %>% 
    as.data.frame() %>%
    dplyr::filter(str_detect(rownames(.), coll("analysis_factor")) == TRUE) %>%
    dplyr::mutate(analysis_factor = factor_chosen,
                  adjusted = use_covariates) %>%
    dplyr::select(analysis_factor, everything())
  
  return(permanova_output) 
}





# Apply the function for all factors analyzed ----
#------------------------------------------------#
#                                                #
#            ANALYZE MULTIPLE FACTORS            # 
#                                                #
#------------------------------------------------#

# Run PERMANOVA on all the codes chosen
PERMANOVA_summary_output <- data.frame()

run_all_PERMANOVA <- TRUE
if (run_all_PERMANOVA == TRUE){
  # Run the function for all factors 
  for (i in factors_analyzed){
    output_run = find_MB_associated_factors(factor_chosen = i, use_covariates = include_covariates)
    PERMANOVA_summary_output = rbind(PERMANOVA_summary_output, output_run)
  }
  
  # Save the output as AOV table
  saveRDS(PERMANOVA_summary_output, file = output_file)
}



