

# Command line arguments ----
#------------------------------------------------#
#                                                #
#          SET UP COMMAND LINE ARGUMENTS         # 
#                                                #
#------------------------------------------------#

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)            
# Expected input:
# args = c(count_matrix_file, metadata_file, factors_analyzed_file, output_file)
# args = c("RData/Abundance_species.rds", "RData/Data_master.rds", "RData/RDA_intrinsic_factors_analyzed.rds", "RData/Test_RDA_output.rds") 

cat(args, sep = "\n")

# Test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file)", call. = FALSE)
} 

# Declare count matrix input
count_matrix_file <- args[1]

# Declare metadata used by model
metadata_file <- args[2]

# Declare facors analyzed
factors_analyzed_file <- args[3]

# Declare prefix for output files
output_file <- args[4]






# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load packages
library("vegan")
library("dplyr")
library("stringr")

# Read count matrix used as input (genes, species etc)
count_data <- readRDS(file = count_matrix_file)  %>% # feature x subject
  as.data.frame()

# Read the phenotype data
full_factor_data <- readRDS(file = metadata_file) 

# Read factors used in analysis
factors_analyzed <- readRDS(file = factors_analyzed_file)

# Read skood vkood link data
skood_vkood_link <- readRDS("RData/Skood_Vkood_link.rds") 

# Samples with extremely low reads to be filtered out
low_read_scodes <- readRDS("RData/low_read_scodes.rds")





# Data transformations ----
#------------------------------------------------#
#                                                #
#              DATA TRANSFORMATION               # 
#                                                #
#------------------------------------------------#

# Apply pseudocount
pseudocount = 0.5

count_data_help <- count_data
count_data_help[count_data_help == 0] <- pseudocount

# CLR_transformation
count_data_transformed <- apply(count_data_help, 2, function(x) log(x) - mean(log(x))) %>%
  t() %>%
  as.data.frame()







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

# Merge data containing scodes - linking with registry data
count_data_scoded <- count_data_transformed %>% 
  tibble::rownames_to_column(var = "vkood") %>%
  dplyr::left_join(corresponding_scodes[ ,c("vkood", "skood")], by = "vkood") %>%
  dplyr::select(-one_of("vkood")) %>%
  dplyr::filter(is.na(skood) == FALSE) %>% 
  tibble::column_to_rownames(var = "skood")

# Merge data for modelling
metadata <- full_factor_data %>%
  dplyr::filter(!(skood %in% low_read_scodes)) %>%   # Exclude samples with extremely low reads (3)
  dplyr::select(skood, all_of(factors_analyzed)) %>% 
  dplyr::filter(complete.cases(.)) %>% 
  tibble::column_to_rownames(var = "skood")

# Subset count data
count_data_ordered <- count_data_scoded[rownames(metadata), ]






# Run RDA ----
#------------------------------------------------#
#                                                #
#             RUN REDUNCANCY ANALYSIS            # 
#                                                #
#------------------------------------------------#

# Run RDA model
mod0 <- vegan::rda(count_data_ordered ~ 1, metadata)  # Model with intercept only
mod1 <- vegan::rda(count_data_ordered ~ ., metadata)  # Model with all explanatory variables

set.seed(1)
mod <- ordistep(mod0, scope = formula(mod1), direction = "forward", steps = 300, permutations = 200, parallel = 20)

# Save output
saveRDS(list(RsquareAdj(mod)$r.squared*100,
             mod$anova), output_file)



