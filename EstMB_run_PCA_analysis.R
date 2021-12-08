

# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load libraries
library("tidyverse")
library("stringr")
library("ggthemes")
library("gridExtra")

# Read distance matrix used as input (genes, species etc)
distance_matrix_raw <- readRDS(file = "RData/CLR_distance_species_clean.rds") # subject x subject

# Read the phenotype data
full_factor_data <- readRDS(file = "RData/Data_master.rds") %>% 
  as.data.frame() 

# Read skood vkood link data
skood_vkood_link <- readRDS("RData/Skood_Vkood_link.rds") 

# Samples with extremely low reads to be filtered out
low_read_scodes <- readRDS("RData/low_read_scodes.rds")

# Dominant genus for each sample
dominant_genus <- readRDS("RData/Dominant_genera_perSample.rds")







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
  dplyr::mutate(vkood_new = str_replace(vkood, coll("."), "_"),
                vkood_new = ifelse(substr(vkood_new, nchar(vkood_new)-1, nchar(vkood_new)) == "_1", 
                                   substr(vkood_new, 1, nchar(vkood_new)-2), vkood_new)) %>% 
  dplyr::left_join(skood_vkood_link, by = c("vkood_new" = "vkood")) %>%
  dplyr::filter(!(skood %in% low_read_scodes)) %>% 
  dplyr::filter(complete.cases(.)) 

# Subset the distance matrix and rename with scodes
distance_matrix_scoded <- distance_matrix_mat[corresponding_scodes$vkood, corresponding_scodes$vkood]
colnames(distance_matrix_scoded) <- rownames(distance_matrix_scoded) <- corresponding_scodes$skood







# Perform multidimensional scaling ----
#------------------------------------------------#
#                                                #
#                   PERFORM MDS                  # 
#                                                #
#------------------------------------------------#

# Genus data - for visualizing the relative abundance of most dominant genera
genus_data_TSS <- readRDS("Data/Abundance_genus.rds") %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "vkood") %>% 
  dplyr::left_join(corresponding_scodes, by = "vkood") %>% 
  dplyr::select(-one_of("vkood", "vkood_new")) %>%
  tidyr::pivot_longer(-skood, names_to = "key", values_to = "value") %>%
  dplyr::group_by(skood) %>% 
  dplyr::mutate(n_reads = sum(value),
                rel_value = value/n_reads) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(is.na(skood) == FALSE) %>%
  dplyr::filter(str_detect(key, "g__Prevotella") | str_detect(key, "g__Bacteroides") | str_detect(key, "g__Clostridium")) %>% 
  dplyr::select(skood, key, rel_value) %>% 
  tidyr::spread(key, rel_value)

# Perform MDS
MDS <- cmdscale(distance_matrix_scoded, k = 50, eig = TRUE) # Calculate first 50 principal components
plot(MDS$eig[1:50]/sum(MDS$eig)*100)
cumsum(MDS$eig[1:50]/sum(MDS$eig)*100)

# PCA score data
MDS_data <- MDS$points %>%
  as.data.frame() %>%
  stats::setNames(paste("PC", seq(1, ncol(.)), sep = "")) %>%
  tibble::rownames_to_column(var = "skood")

saveRDS(MDS_data, "RData/PCA_components_data.rds")

# Merge factor data with MDS data
factor_data_MDS <- full_factor_data %>% 
  dplyr::left_join(MDS_data, by = "skood") %>%
  dplyr::left_join(corresponding_scodes, by = "skood") %>% 
  dplyr::left_join(dominant_genus, by = "vkood") %>% 
  dplyr::left_join(genus_data_TSS, by = "skood") %>% 
  as.data.frame()


# Visualize the biplot according to the dominant genera relative abundances
p1 <- ggplot(factor_data_MDS %>% dplyr::arrange(`k__Bacteria / p__Bacteroidetes / c__Bacteroidia / o__Bacteroidales / f__Prevotellaceae / g__Prevotella`), 
             aes(x = PC1, y = PC2, color = `k__Bacteria / p__Bacteroidetes / c__Bacteroidia / o__Bacteroidales / f__Prevotellaceae / g__Prevotella`)) + 
  geom_point(size = 5) + 
  ggtitle("Prevotella") + 
  scale_color_gradient2_tableau(name = "Relative\nabundance", trans = "reverse") +
  theme_bw() + 
  theme(axis.title = element_text(size = 12))

p2 <- ggplot(factor_data_MDS %>% dplyr::arrange(`k__Bacteria / p__Bacteroidetes / c__Bacteroidia / o__Bacteroidales / f__Bacteroidaceae / g__Bacteroides`), 
             aes(x = PC1, y = PC2, color = `k__Bacteria / p__Bacteroidetes / c__Bacteroidia / o__Bacteroidales / f__Bacteroidaceae / g__Bacteroides`)) + 
  geom_point(size = 5) + 
  ggtitle("Bacteroides") + 
  scale_color_gradient2_tableau(name = "Relative\nabundance", trans = "reverse") +
  theme_bw() + 
  theme(axis.title = element_text(size = 12))


p3 <- ggplot(factor_data_MDS %>% dplyr::arrange(`k__Bacteria / p__Firmicutes / c__Clostridia / o__Clostridiales / f__Clostridiaceae / g__Clostridium`), 
             aes(x = PC1, y = PC2, color = `k__Bacteria / p__Firmicutes / c__Clostridia / o__Clostridiales / f__Clostridiaceae / g__Clostridium`)) + 
  geom_point(size = 5) + 
  ggtitle("Clostridium") + 
  scale_color_gradient2_tableau(name = "Relative\nabundance", trans = "reverse") +
  theme_bw() + 
  theme(axis.title = element_text(size = 12))

p_PCA <- gridExtra::grid.arrange(p1, p2, p3, ncol = 1)

ggsave(plot = p_PCA, filename = "Figures/PCA_biplot_genera.png", height = 12, width = 6)
ggsave(plot = p_PCA, filename = "Figures/PCA_biplot_genera.svg", height = 12, width = 6)






