

# Command line arguments ----
#------------------------------------------------#
#                                                #
#          SET UP COMMAND LINE ARGUMENTS         # 
#                                                #
#------------------------------------------------#

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)   
# Expected input:
# args = c(count_matrix_file, metadata_file, model_name, preprocessing, seed, grid_size, predictorSet, output_file_prefix)
# args = c("RData/Abundance_species.rds", "RData/Data_master.rds", "LASSO", "CLR", "1", "5", "SET2", "Rdata/ML") 

cat(args, sep = "\n")

# Test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file)", call. = FALSE)
} 

# Declare count matrix input
count_matrix_file <- args[1]

# Declare phenotype data used, dataframe
metadata_file <- args[2]

# Declare ML algorithm name, character (currently RF and LASSO supported)
model_name <- args[3]

# Declare microbiome data preprocessing step used (Currently CLR and TSS)
preprocessing <- args[4]

# Declare seed used for data splitting
seed <- as.numeric(args[5])

# Declare grid size
grid_size <- as.numeric(args[6])

# Declare predictor set used (SET0, SET2, SET2_MB, SET3, SET3_MB)
predictorSet <- args[7]

# Declare prefix for output files
output_file_prefix <- args[8]






# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load packages
library("tidymodels")
library("dplyr")
library("stringr")

# Read count matrix used as input (genes, species etc)
count_data <- readRDS(file = count_matrix_file)  %>%  # feature x sample
  t() %>%
  as.data.frame()

# Read the source data
full_factor_data <- readRDS(file = metadata_file) 

# Read factors used in analysis - factors with associations with beta-diversity with > 50 cases
factors_analyzed <- c("I10", "M10", "E11", "I11", "J20", "F41", "K58", "K21", "K35", "J35", "E78", "K80", "K29", 
                      "J45", "H25", "F32", "N20", "F33", "N10", "J30", "K35", "D50", "M15", "B37", "N95", "N81") 

# Read skood vkood link data
skood_vkood_link <- readRDS("RData/Skood_Vkood_link.rds") 

# Anitbiotics users (last 180 days)
q_antibiotic_users <- readRDS("RData/q_antibiotic_users.rds")

# Scodes with low read numbers left aside
low_read_scodes <- readRDS("RData/low_read_scodes.rds")







# Data transformations ----
#------------------------------------------------#
#                                                #
#              DATA TRANSFORMATION               # 
#                                                #
#------------------------------------------------#

if (preprocessing == "TSS"){
  
  count_data_transformed <- count_data/rowSums(count_data)
  
} else if (preprocessing == "CLR"){
  
  # Apply pseudocount
  pseudocount = 0.5
  
  count_data_help <- count_data
  count_data_help[count_data_help == 0] <- pseudocount
  
  # CLR_transformation
  count_data_transformed <- apply(count_data_help, 1, function(x) log(x) - mean(log(x))) %>%
    t() %>%
    as.data.frame()
}





# Data preprocessing ----
#------------------------------------------------#
#                                                #
#                DATA PREPROCESSING              # 
#                                                #
#------------------------------------------------#

# Merge data containing scodes - linking with registry data
corresponding_scodes <- data.frame(vkood = rownames(count_data_transformed), stringsAsFactors = FALSE) %>%
  dplyr::mutate(vkood_new = stringr::str_replace(string = vkood, pattern = coll("."), replacement = "_"),
                vkood_new = ifelse(substr(vkood_new, nchar(vkood_new)-1, nchar(vkood_new)) == "_1", 
                                   substr(vkood_new, 1, nchar(vkood_new)-2), vkood_new)) %>% 
  dplyr::left_join(skood_vkood_link, by = c("vkood_new" = "vkood")) %>%
  dplyr::filter(complete.cases(.)) 

# Subset the distance matrix and rename with scodes
count_data_scoded <- count_data_transformed %>% 
  tibble::rownames_to_column(var = "vkood") %>%
  dplyr::left_join(corresponding_scodes[ ,c("vkood", "skood")], by = "vkood") %>%
  dplyr::select(-one_of("vkood")) %>%
  dplyr::filter(is.na(skood) == FALSE)

# Merge data for modelling
data_merged <- full_factor_data %>%
  dplyr::filter(!(skood %in% q_antibiotic_users)) %>%
  dplyr::filter(!(skood %in% low_read_scodes)) %>%
  dplyr::left_join(count_data_scoded, by = "skood")






# ML modelling ----
#------------------------------------------------#
#                                                #
#               TIDYMODELS FUNCTION              # 
#                                                #
#------------------------------------------------#

# Define the function
# @input - ICD_chosen - character name of the disease (column in the phenotype dataframe)
# @input - p_predictorSet - character indicationg, which predictors should be used. Options:
#                           SET0 - MB features only
#                           SET2 - age, BMI, gender and Bristol stool scale
#                           SET2_MB - age, BMI, gender and Bristol stool scale + MB features
#                           SET3 - age, BMI, gender and Bristol stool scale + antibiotics usage 
#                           SET3_MB - age, BMI, gender and Bristol stool scale + antibiotics usage + MB features
# Function returns PERMANOVA aov table

find_predictive_model = function(ICD_chosen, p_predictorSet){
  
  # Initial data manipulations 
  if (p_predictorSet == "SET0"){
    metadata <- data_merged %>% 
      dplyr::mutate(analysis_factor = data_merged %>% dplyr::pull(ICD_chosen), 
                    analysis_factor = factor(analysis_factor)) %>%
      dplyr::select(analysis_factor, all_of(colnames(count_data_scoded))) %>% 
      dplyr::select(-one_of("skood")) %>%       
      dplyr::filter(complete.cases(.))
    
  } else if (p_predictorSet == "SET2"){
    metadata <- data_merged %>% 
      dplyr::mutate(analysis_factor = data_merged %>% dplyr::pull(ICD_chosen), 
                    analysis_factor = factor(analysis_factor)) %>%
      dplyr::select(gender, Age_at_MBsample, BMI, usualStoolType_category, analysis_factor) %>% 
      dplyr::filter(complete.cases(.))
    
  } else if (p_predictorSet == "SET2_MB"){
    metadata <- data_merged %>% 
      dplyr::mutate(analysis_factor = data_merged %>% dplyr::pull(ICD_chosen), 
                    analysis_factor = factor(analysis_factor)) %>%
      dplyr::select(gender, Age_at_MBsample, BMI, usualStoolType_category, 
                    analysis_factor, all_of(colnames(count_data_scoded))) %>% 
      dplyr::select(-one_of("skood")) %>% 
      dplyr::filter(complete.cases(.))
    
  }else if (p_predictorSet == "SET3"){
    metadata <- data_merged %>% 
      dplyr::mutate(analysis_factor = data_merged %>% dplyr::pull(ICD_chosen), 
                    analysis_factor = factor(analysis_factor)) %>%
      dplyr::select(gender, Age_at_MBsample, BMI, usualStoolType_category, antibiotics_history_cont, 
                    analysis_factor) %>% 
      dplyr::filter(complete.cases(.))
    
  } else if (p_predictorSet == "SET3_MB"){
    metadata <- data_merged %>% 
      dplyr::mutate(analysis_factor = data_merged %>% dplyr::pull(ICD_chosen), 
                    analysis_factor = factor(analysis_factor)) %>%
      dplyr::select(gender, Age_at_MBsample, BMI, usualStoolType_category, antibiotics_history_cont, 
                    analysis_factor, all_of(colnames(count_data_scoded))) %>% 
      dplyr::select(-one_of("skood")) %>% 
      dplyr::filter(complete.cases(.))
  }
  
  
  
  # Initial data split ----
  #------------------------------------------------#
  
  # Create random data split - fixed seed
  set.seed(seed)
  data_split <- rsample::initial_split(metadata, prop = .75, strata = analysis_factor) 
  
  # Split the data
  data_train <- rsample::training(data_split)
  data_test <- rsample::testing(data_split)
  
  # Create CV sets
  set.seed(seed)
  data_cv <- rsample::vfold_cv(data_train, v = 5, repeats = 1, strata = analysis_factor) 
  
  
  
  # Preprocessing steps for the data ----
  #------------------------------------------------#
  
  # Data Preprocessing
  data_preprocessing <- 
    # Specify the outcome variable and predictor variables
    recipes::recipe(analysis_factor ~ ., data = data_train) %>%
    recipes::step_dummy(all_nominal(), -all_outcomes()) %>%
    # Apply preprocessing steps to the data (either one by one or for numeric/character etc)
    recipes::step_normalize(all_predictors()) %>%
    recipes::step_nzv(all_predictors())
  
  
  
  # Define the model used ----
  #------------------------------------------------#
  
  # Define the model - no training done yet
  if (model_name == "RF"){
    
    model <- parsnip::rand_forest() %>%
      parsnip::set_args(
        trees = 1000,    # Can set either specific values or leave them to be tuned
        mtry = tune::tune(),  
        min_n = tune::tune()) %>%
      set_engine("ranger") %>% # Engine = underlying package
      set_mode("classification")
    
  }else if (model_name == "LASSO"){
    
    model <- parsnip::logistic_reg() %>%
      parsnip::set_args(
        penalty = tune::tune(),  
        mixture = tune::tune()) %>%
      set_engine("glmnet") %>% 
      set_mode("classification")
  }
  
  
  # Construct a workflow that combines your recipe and your model
  ml_workflow <- workflows::workflow() %>%
    workflows::add_recipe(data_preprocessing) %>%
    workflows::add_model(model)
  
  
  # Start tuning the model ----
  #------------------------------------------------#
  
  
  # # Tune the models
  set.seed(seed)
  tuned_model <- ml_workflow %>%
    tune::tune_grid(resamples = data_cv,
                    grid = grid_size, 
                    metrics = yardstick::metric_set(yardstick::roc_auc))

  # Return best parameters
  best_params <-
    tuned_model %>%
    tune::select_best(metric = "roc_auc")
  
  
  
  # Evaluate model performance ----
  #------------------------------------------------#
  
  # Describe the best models according to the best parameters gathered earlier
  best_model <-
    ml_workflow %>%
    tune::finalize_workflow(best_params) 
  
  # Fit the models
  model_final <- best_model %>%
    tune::last_fit(data_split)
  
  # Collect CV metrics
  CV_performance <- tuned_model %>% tune::show_best(n = 1, metric = "roc_auc") 
  
  # Collect test metrics
  test_performance <- model_final %>% tune::collect_metrics(metric = "roc_auc") %>% dplyr::filter(.metric == "roc_auc")
  
  
  
  # Output ----
  #------------------------------------------------#
  
  # Collect model performance data
  model_performance <- dplyr::bind_rows(CV_performance, 
                                        test_performance) %>%
    dplyr::mutate(value = ifelse(is.na(mean), .estimate, mean), 
                  metric = .metric, 
                  set = c("CV", "test"), 
                  analysis_factor = ICD_chosen, 
                  algorithm = model_name,
                  predictors = p_predictorSet,
                  MB_preprocessing = preprocessing, 
                  seed) %>%
    dplyr::select(analysis_factor, algorithm, predictors, MB_preprocessing,  set, metric, value, seed)

  # Output
  if (model_name == "LASSO"){
    predictors <- best_model %>%
      fit(data = data_train) %>% 
      pull_workflow_fit() %>% 
      tidy() %>% 
      dplyr::filter(estimate != 0) 
    
    return(list(best_params, model_performance, predictors))
  } else{
    return(list(model_final, model_performance))
  }
}








# Apply the function for numerous ICD codes ----
#------------------------------------------------#
#                                                #
#              ANALYSE ALL DISEASES              # 
#                                                #
#------------------------------------------------#

# Run tidymodels on all the codes chosen
run_all_models <- TRUE
if (run_all_models == TRUE){
  for (i in factors_analyzed){
    
    # Create names for the outputs
    model_output_file <- paste("tidymodels_", model_name, "_", i, "_", predictorSet, "_", preprocessing, "_seed", seed, "_model.rds", sep = "")
    params_output_file <- paste("tidymodels_", model_name, "_", i, "_", predictorSet, "_", preprocessing, "_seed", seed, "_params.rds", sep = "")
    results_output_file <- paste("tidymodels_", model_name, "_", i, "_", predictorSet, "_", preprocessing, "_seed", seed, "_output.rds", sep = "")
    estimates_output_file <- paste("tidymodels_", model_name, "_", i, "_", predictorSet, "_", preprocessing, "_seed", seed, "_estimates.rds", sep = "")
    
    output_run = find_predictive_model(ICD_chosen = i, p_predictorSet = predictorSet)
    model = output_run[[1]] # Model object
    modelling_results = output_run[[2]] # Model performance
    
    if (model_name == "LASSO"){
      estimates = output_run[[3]] # LASSO estimates
      
      saveRDS(modelling_results, file = file.path(output_file_prefix, results_output_file))
      saveRDS(estimates, file = file.path(output_file_prefix, estimates_output_file))
      saveRDS(model, file = file.path(output_file_prefix, model_output_file))
      
    } else{
      saveRDS(modelling_results, file = file.path(output_file_prefix, results_output_file))
      saveRDS(model, file = file.path(output_file_prefix, model_output_file))
    }
  }
}
