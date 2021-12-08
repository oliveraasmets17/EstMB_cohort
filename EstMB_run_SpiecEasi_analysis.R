



# Command line arguments ----
#------------------------------------------------#
#                                                #
#          SET UP COMMAND LINE ARGUMENTS         # 
#                                                #
#------------------------------------------------#

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)     
# Expected input:
# args = c(count_matrix_file, output_file)
# args = c("RData/Abundance_species.rds", "RData/SpiecEasi_species.rds") 

cat(args, sep = "\n")

# Test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file)", call. = FALSE)
} 

# Declare count matrix input
count_matrix_file <- args[1]

# Declare prefix for output files
output_file <- args[2]







# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load libraries
library("dplyr")
library("SpiecEasi")
library("stringr")

# Read skood vkood link data
skood_vkood_link <- readRDS("RData/Skood_Vkood_link.rds") 

# Samples with extremely low reads to be filtered out
low_read_scodes <- readRDS("RData/low_read_scodes.rds")






# Create co-abundance network ----
#------------------------------------------------#
#                                                #
#               VISUALIZING RESULTS              # 
#                                                #
#------------------------------------------------#

# Read count matrix used as input (genes, species etc)
count_data <- readRDS(file = count_matrix_file)  %>%  # feature x subject
  t() %>%
  as.data.frame()

# Merge data containing scodes - linking with registry data
corresponding_scodes <- data.frame(vkood = rownames(count_data), stringsAsFactors = FALSE) %>%
  dplyr::mutate(vkood_new = stringr::str_replace(string = vkood, pattern = coll("."), replacement = "_"),
                vkood_new = ifelse(substr(vkood_new, nchar(vkood_new)-1, nchar(vkood_new)) == "_1", 
                                   substr(vkood_new, 1, nchar(vkood_new)-2), vkood_new)) %>% 
  dplyr::left_join(skood_vkood_link, by = c("vkood_new" = "vkood")) %>%
  dplyr::filter(complete.cases(.)) 

# Subset the distance matrix and rename with scodes
count_data_scoded <- count_data %>% 
  tibble::rownames_to_column(var = "vkood") %>%
  dplyr::left_join(corresponding_scodes[ ,c("vkood", "skood")], by = "vkood") %>%
  dplyr::select(-one_of("vkood")) %>%
  dplyr::filter(!(skood %in% low_read_scodes)) %>% 
  dplyr::filter(is.na(skood) == FALSE) %>%
  tibble::column_to_rownames(var = "skood") %>%
  as.matrix()

# Run SpiecEasi
spiecEasi <- SpiecEasi::spiec.easi(count_data_scoded, 
                                   method = "glasso", 
                                   sel.criterion = "bstars", 
                                   pulsar.select = TRUE,
                                   verbose = FALSE, 
                                   lambda.min.ratio = 1e-3, 
                                   nlambda = 20,
                                   list(rep.num = 20, 
                                        thresh = 0.05))

# Save the SpiecEasi object
saveRDS(spiecEasi, output_file)

