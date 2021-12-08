

# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

#Load packages
library("tidyverse")
library("Matrix")

# Read distance matrix used as input (genes, species etc)
distance_matrix_raw <- readRDS(file = "C:/Users/oliver17/Desktop/Doktorantuur/Projekt_EstMiBiom/RData/CLR_distance_species_clean.rds") # subject x subject





# Work with the distance matrix ----
#------------------------------------------------#
#                                                #
#           PREPROCESS DISTANCE MATRIX           # 
#                                                #
#------------------------------------------------#

# Distance object to matrix
distance_matrix_mat <- distance_matrix_raw %>%
  as.matrix() 

# Data for plotting
replicates_plotData <- Matrix::tril(distance_matrix_mat) %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "sample1") %>% 
  tidyr::gather(sample2, value, -sample1) %>%
  dplyr::filter(value != 0) %>% 
  dplyr::mutate(isReplicate = case_when(substr(sample1, 1, 6) == "VR4UUJ" & substr(sample2, 1, 6) == "VR4UUJ"~ "Biological replicates", 
                                     TRUE ~ "Random pairs"))

# Visualize replicates
replicates_plot <- ggplot(replicates_plotData, aes(x = isReplicate, y = value)) + 
  geom_boxplot() + 
  ylab("CLR-distance") +
  xlab("") + 
  theme_bw() + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

ggsave(replicates_plot, file = "Figures/Replicates.png", height = 8, width = 8)
