

# Command line arguments ----
#------------------------------------------------#
#                                                #
#          SET UP COMMAND LINE ARGUMENTS         # 
#                                                #
#------------------------------------------------#

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)       
# Expected input:
# args = c(abundance_matrix_file, output_file)
# args = c("RData/Abundance_species.rds", "RData/Test_CLR_distance_species.rds")

cat(args, sep = "\n")

# Test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file)", call. = FALSE)
} 


# Declare count matrix input
abundance_matrix_file <- args[1]

# Declare output file name
output_file <- args[2]






# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load packages
library("dplyr")

# Read data
count_data_raw <- readRDS(abundance_matrix_file) 

# Initial data manipulations
count_data <- count_data_raw %>%
  t()

# Define pseudocount
pseudocount <- 0.5

# Impute zero
count_data_zeroReplaced <- count_data
count_data_zeroReplaced[count_data_zeroReplaced == 0] <- pseudocount

# CLR transformation
count_data_CLR <- t(apply(count_data_zeroReplaced, 1, function(x) log(x) - mean(log(x))))

# Calculate Euclidean distances
Euclidean_CLR <- dist(count_data_CLR, method = "euclidean")

# Save the distance matrix 
save_dist_matrix <- TRUE
if (save_dist_matrix == TRUE){
  saveRDS(Euclidean_CLR, output_file)
}

