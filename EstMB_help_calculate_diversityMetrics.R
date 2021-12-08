

# Command line arguments ----
#------------------------------------------------#
#                                                #
#          SET UP COMMAND LINE ARGUMENTS         # 
#                                                #
#------------------------------------------------#

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)  
# Expected input:
# args = c(abundance_matrix_file, calculate_shannon, output_file)
# args = c("RData/Abundance_species.rds", "TRUE", "RData/Test_diversity_species.rds")

cat(args, sep = "\n")

# Test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file)", call. = FALSE)
} 


# Declare count matrix file
input_file <- args[1]

# Calculate shannon index?
calculate_shannon <- args[2]

# Declare output file name
output_file <- args[3]






# Calculate diversity metrics ----
#------------------------------------------------#
#                                                #
#          CALCULATE DIVERSITY METRICS           # 
#                                                #
#------------------------------------------------#

# Load packages
library("dplyr")
library("tibble")
library("data.table")
library("vegan")

# Read data
count_data_raw <- readRDS(input_file) # feature x subject

# Calculate diversity metrics
diversity_observed <- apply(count_data_raw > 0, 2, sum)

if (calculate_shannon == TRUE){
  diversity_shannon <- vegan::diversity(t(count_data_raw), index = "shannon", MARGIN = 1)
  
  diversity_data <- data.frame(vkood = colnames(count_data_raw), 
                               diversity_observed, 
                               diversity_shannon,
                               input_file) 
} else{
  diversity_data <- data.frame(vkood = colnames(count_data_raw), 
                               diversity_observed,
                               input_file)
}

# Save data
saveRDS(diversity_data, output_file)

