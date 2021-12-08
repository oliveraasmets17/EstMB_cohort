

# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load packages
library("tidyverse")
library("viridis")
library("ggpubr")
library("SpiecEasi")
library("igraph")
library("GGally")
library("ggthemes")
library("reshape2")
library("xlsx")

# Read raw results files
MWAS_species_diseases <- readRDS(file = "RData archive/MWAS_species_p10d001_diseaseFactors_SET2.rds")
MWAS_species_diseases_ABadjusted <- readRDS(file = "RData archive/MWAS_species_p10d001_diseaseFactors_SET3.rds")
MWAS_species_diseases_medicationAdjusted <- readRDS(file = "RData archive/MWAS_species_p10d001_diseaseFactors_medicationAdjusted.rds")

MWAS_species_dietary <- readRDS(file = "RData archive/MWAS_species_p10d001_dietaryFactors_SET2.rds")
MWAS_species_medications <- dplyr::bind_rows(readRDS(file = "RData archive/MWAS_species_p10d001_medicationFactorsSet1_SET2.rds"),
                                             readRDS(file = "RData archive/MWAS_species_p10d001_medicationFactorsSet2_SET2.rds"),
                                             readRDS(file = "RData archive/MWAS_species_p10d001_medicationFactorsSet3_SET2.rds"))

MWAS_KEGGko_diseases <- readRDS(file = "RData archive/MWAS_KEGGko_p10d001_diseaseFactors_SET2.rds")
MWAS_KEGGko_diseases_ABadjusted <- readRDS(file = "RData archive/MWAS_KEGGko_p10d001_diseaseFactors_SET3.rds")
MWAS_KEGGko_diseases_medicationAdjusted <- readRDS(file = "RData archive/MWAS_KEGGko_p10d001_diseaseFactors_medicationAdjusted.rds")

MWAS_KEGGko_dietary <- readRDS(file = "RData archive/MWAS_KEGGko_p10d001_dietaryFactors_SET2.rds")
MWAS_KEGGko_medications <- dplyr::bind_rows(readRDS(file = "RData archive/MWAS_KEGGko_p10d001_medicationFactorsSet1_SET2.rds"),
                                            readRDS(file = "RData archive/MWAS_KEGGko_p10d001_medicationFactorsSet2_SET2.rds"),
                                            readRDS(file = "RData archive/MWAS_KEGGko_p10d001_medicationFactorsSet3_SET2.rds"))

# Gather all results
MWAS_species_raw <- dplyr::bind_rows(MWAS_species_diseases, MWAS_species_diseases_ABadjusted, MWAS_species_dietary, MWAS_species_medications)
MWAS_KEGGko_raw <- dplyr::bind_rows(MWAS_KEGGko_diseases, MWAS_KEGGko_diseases_ABadjusted, MWAS_KEGGko_dietary, MWAS_KEGGko_medications)

# Calculate CLR median values for species network node size
count_data <- readRDS(file = "RData/Abundance_filtered_species_p10d001.rds")  %>%  # sample x eafture
  t() %>%
  as.data.frame()

pseudocount = 0.5

count_data_help <- count_data
count_data_help[count_data_help == 0] <- pseudocount

count_data_transformed <- apply(count_data_help, 1, function(x) log(x) - mean(log(x))) %>%
  t() %>%
  as.data.frame()

CLR_mean_value <- data.frame(taxa = colnames(count_data_transformed), 
                             meanCLR = colMeans(count_data_transformed))

# Read help datasets
KEGG_db_clean <- readRDS("RData/KEGG_db_clean.rds")
ICD10_name_data <- readRDS("RData/ICD10_name_data.rds")





# Process data for plotting ----
#------------------------------------------------#
#                                                #
#               DATA PREPROCESSING               # 
#                                                #
#------------------------------------------------#

# Clean data and keep the most important columns
MWAS_species_data <- MWAS_species_raw %>% 
  dplyr::filter(`model.analysis_factor Pr(>|t|).BH` <= 0.05) %>%
  dplyr::rename("Estimate" = "model.analysis_factor Estimate",
                "Std" = "model.analysis_factor Std. Error",
                "BH" = "model.analysis_factor Pr(>|t|).BH") %>%
  dplyr::select(Estimate, Std, BH, analysis_factor, feature, adjusted)

MWAS_KEGGko_data <- MWAS_KEGGko_raw %>% 
  dplyr::filter(`model.analysis_factor Pr(>|t|).BH` <= 0.05) %>%
  dplyr::rename("Estimate" = "model.analysis_factor Estimate",
                "Std" = "model.analysis_factor Std. Error",
                "BH" = "model.analysis_factor Pr(>|t|).BH") %>%
  dplyr::select(Estimate, Std, BH, analysis_factor, feature, adjusted) 

# Filter disease associations
MWAS_species_diseaseAssociations <- MWAS_species_data %>% 
  dplyr::filter(analysis_factor %in% readRDS("RData archive/Disease_factors_analyzed.rds")) %>% 
  dplyr::filter(!(analysis_factor %in% c("definedHealthy", "health_status_ok", "mental_health_status_ok")))

MWAS_KEGGko_diseaseAssociations <- MWAS_KEGGko_data %>% 
  dplyr::filter(analysis_factor %in% readRDS("RData archive/Disease_factors_analyzed.rds")) %>% 
  dplyr::filter(!(analysis_factor %in% c("definedHealthy", "health_status_ok", "mental_health_status_ok")))







# Analyze associations ----
#------------------------------------------------#
#                                                #
#               ANALYZE ASSOCIATIONS             # 
#                                                #
#------------------------------------------------#

# Summarize the number of associations by phylum/genus
associations_by_phylum <- MWAS_species_diseaseAssociations %>% 
  tidyr::separate(feature, into = c("k", "p", "c", "o", "f", "g", "s"), sep = ";", remove = FALSE) %>% 
  dplyr::group_by(k, p, adjusted) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  tidyr::spread(adjusted, n, fill = 0) %>% 
  dplyr::arrange(desc(SET2)) %>% 
  dplyr::mutate(prop = SET3/SET2*100)

associations_by_genus <- MWAS_species_diseaseAssociations %>% 
  tidyr::separate(feature, into = c("k", "p", "c", "o", "f", "g", "s"), sep = ";", remove = TRUE) %>% 
  dplyr::group_by(k, p, c, o, f, g, adjusted) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  tidyr::spread(adjusted, n, fill = 0) %>% 
  dplyr::arrange(desc(SET2)) %>% 
  dplyr::mutate(prop = SET3/SET2*100)

# Summarize the number of associations by KEGG categories
associations_by_pathway <- MWAS_KEGGko_diseaseAssociations %>% 
  dplyr::left_join(KEGG_db_clean, by = c("feature" = "KO")) %>% 
  dplyr::group_by(Pathway, analysis_factor, adjusted) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  tidyr::spread(adjusted, n, fill = 0) %>% 
  dplyr::arrange(analysis_factor, desc(SET2)) %>% 
  dplyr::mutate(prop = SET3/SET2*100)

associations_by_subPathway <- MWAS_KEGGko_diseaseAssociations %>% 
  dplyr::left_join(KEGG_db_clean, by = c("feature" = "KO")) %>% 
  dplyr::group_by(`Sub-pathway`,analysis_factor,  adjusted) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  tidyr::spread(adjusted, n, fill = 0) %>% 
  dplyr::arrange(analysis_factor, desc(SET2)) %>% 
  dplyr::mutate(prop = SET3/SET2*100)

associations_by_metabolicPathway <- MWAS_KEGGko_diseaseAssociations %>% 
  dplyr::left_join(KEGG_db_clean, by = c("feature" = "KO")) %>% 
  dplyr::group_by(`Metabolic pathway`, analysis_factor, adjusted) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  tidyr::spread(adjusted, n, fill = 0) %>% 
  dplyr::arrange(analysis_factor, desc(SET2)) %>% 
  dplyr::mutate(prop = SET3/SET2*100)






# Visualize overlap in associations between diseases ----
#------------------------------------------------#
#                                                #
#       PAIRWISE OVERLAP IN ASSOCIATIONS         # 
#                                                #
#------------------------------------------------#

# Focus on diseases with at least 2 associations
at_least_2associations <- MWAS_species_diseaseAssociations %>% 
  dplyr::filter(adjusted == "SET2") %>%
  dplyr::group_by(analysis_factor) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::arrange(desc(n)) %>% 
  dplyr::filter(n > 0) %>% 
  dplyr::pull(analysis_factor)

# Create blank dataset
help_diseases <- unique(MWAS_species_diseaseAssociations$analysis_factor)
help_df2 <- data.frame(one = rep(help_diseases, each = length(help_diseases)),
                       two = rep(help_diseases, length(help_diseases)))

# Identify overlap of associations
find_overlapping_associations <- function(p_adjusted){
  
  associations_diseases_sub <- MWAS_species_diseaseAssociations %>% 
    dplyr::filter(adjusted == p_adjusted)
  
  help_df1 <- associations_diseases_sub %>% 
    dplyr::mutate(help = 1) %>%
    dplyr::select(feature, analysis_factor, help) %>% 
    tidyr::spread(analysis_factor, help, fill = 0)
  
  for (i in 1:nrow(help_df2)){
    one = help_df2[i, "one"]
    two = help_df2[i, "two"]
    
    if (one %in% colnames(help_df1) & two %in% colnames(help_df1)){
      nrows <- help_df1 %>% 
        dplyr::mutate(one = help_df1[ ,one],
                      two = help_df1[ ,two]) %>% 
        dplyr::filter(one == 1 & two == 1) %>% 
        nrow()
    } else{
      nrows = 0
    }
    
    help_df2[i, "n"] <- nrows
  }
  return(help_df2)
}

overlapping_associations_SET2 <- find_overlapping_associations(p_adjusted = "SET2") %>% 
  dplyr::left_join(ICD10_name_data, by = c("one" = "ICD10_category")) %>% 
  dplyr::mutate(one_name = paste(one, " ", ICD10_name, sep = ""))

overlapping_associations_SET3 <- find_overlapping_associations(p_adjusted = "SET3") %>% 
  dplyr::left_join(ICD10_name_data, by = c("one" = "ICD10_category")) %>% 
  dplyr::mutate(one_name = paste(one, " ", ICD10_name, sep = ""))

# Define order of the variables
factor_order <- overlapping_associations_SET2 %>% 
  dplyr::filter(one == two) %>% 
  dplyr::group_by(one, one_name) %>% 
  dplyr::summarise(s = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(s))

# Creat triangular matrices
create_triangular_matrix <- function(input_data, side = "upper"){
  data_wide <- input_data %>% 
    dplyr::select(one_name, two, n) %>% 
    tidyr::spread(two, n, fill = 0) %>% 
    tibble::column_to_rownames(var = "one_name")
  
  data_wide_ordered <- data_wide[factor_order$one_name, factor_order$one]
  data_wide_ordered_lowerTriangle <- data_wide_ordered 
  
  if (side == "upper"){
    data_wide_ordered_lowerTriangle[upper.tri(data_wide_ordered_lowerTriangle)] <- NA
  } else{
    data_wide_ordered_lowerTriangle[lower.tri(data_wide_ordered_lowerTriangle)] <- NA
  }
  
  output_data <- data_wide_ordered_lowerTriangle %>% 
    tibble::rownames_to_column(var = "one") %>% 
    tidyr::gather(two, n, -one) %>% 
    dplyr::filter(is.na(n) == FALSE) %>% 
    dplyr::filter(substr(one, 1, 3) %in% at_least_2associations & 
                    two %in% at_least_2associations) %>% 
    dplyr::mutate(one = factor(one, levels = factor_order$one_name),
                  two = factor(two, levels = factor_order$one))
  
  return(output_data)
}

association_overlap_data_SET2 <- create_triangular_matrix(input_data = overlapping_associations_SET2, side = "lower")
association_overlap_data_SET3 <- create_triangular_matrix(input_data = overlapping_associations_SET3, side = "lower")

# Visualize the pairwise associations
# Before AB adjustment
p1 <- ggplot(association_overlap_data_SET2, aes(x = two, y = one, fill = log(n))) + 
  geom_tile(color = "white") + 
  geom_text(aes(label = n, size = 6)) + 
  scale_fill_gradientn(colors = rev(c("#9E3D22", "#C25122", "#DF6D27", "#F69035", "#E9BE99",
                                      "#BBC7CD", "#77ACD3", "#5C8FBB", "#4374A3", "#2B5C8A")), limits = c(0, log(248))) + 
  theme_bw() + 
  xlab("") + 
  ylab("") + 
  theme_classic() + 
  scale_y_discrete(position = "right") + 
  theme(axis.text.y = element_text(size = 14),
        legend.position = "NONE", 
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_text(size = 14))
p1

ggsave(plot = p1, filename = "Figures/MWAS_associations_overlap.png", width = 15, height = 10)
ggsave(plot = p1, filename = "Figures/MWAS_associations_overlap.svg", width = 15, height = 10)

# After AB adjustment
p2 <- ggplot(association_overlap_data_SET3, aes(x = two, y = one, fill = log(n))) + 
  geom_tile(color = "white") + 
  geom_text(aes(label = n, size = 6)) + 
  scale_fill_gradientn(colors = rev(c("#9E3D22", "#C25122", "#DF6D27", "#F69035", "#E9BE99",
                                      "#BBC7CD", "#77ACD3", "#5C8FBB", "#4374A3", "#2B5C8A")), limits = c(0, log(248))) + 
  theme_bw() + 
  xlab("") + 
  ylab("") + 
  theme_classic() + 
  scale_y_discrete(position = "right") + 
  theme(axis.text.y = element_text(size = 14),
        legend.position = "NONE", 
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_text(size = 14))

p2

ggsave(plot = p2, filename = "Figures/MWAS_associations_overlap_ABadjusted.png", width = 15, height = 10)
ggsave(plot = p2, filename = "Figures/MWAS_associations_overlap_ABadjusted.svg", width = 15, height = 10)








# Look for common dysbiosis in species level data ----
#------------------------------------------------#
#                                                #
#                COMMON DYSBIOSIS                # 
#                                                #
#------------------------------------------------#

# Summarize the number of associations by phylum/ genus
associations_by_phylum <- MWAS_species_diseaseAssociations %>% 
  tidyr::separate(feature, into = c("k", "p", "c", "o", "f", "g", "s"), sep = ";", remove = FALSE) %>% 
  dplyr::group_by(k, p, adjusted) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  tidyr::spread(adjusted, n, fill = 0) %>% 
  dplyr::arrange(desc(SET2)) %>% 
  dplyr::mutate(prop = SET3/SET2*100)

associations_by_genus <- MWAS_species_diseaseAssociations %>% 
  tidyr::separate(feature, into = c("k", "p", "c", "o", "f", "g", "s"), sep = ";", remove = TRUE) %>% 
  dplyr::group_by(k, p, c, o, f, g, adjusted) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  tidyr::spread(adjusted, n, fill = 0) %>% 
  dplyr::arrange(desc(SET2)) %>% 
  dplyr::mutate(prop = SET3/SET2*100)

# Summarize number of diseases bacteria are associated with 
n_shared_associations <- MWAS_species_diseaseAssociations %>% 
  dplyr::select(analysis_factor, feature, adjusted) %>% 
  dplyr::group_by(feature, adjusted) %>% 
  dplyr::summarise(n = n(),
                   diseases = paste(analysis_factor, collapse = ", ")) %>% 
  dplyr::ungroup()

n_shared_associations_wide <- n_shared_associations %>% 
  dplyr::mutate(adjusted = ifelse(adjusted == "SET2", "Excl_AB", "Incl_AB")) %>% 
  tidyr::pivot_wider(names_from = adjusted, values_from = c("n", "diseases"), values_fill = list(n = 0, diseases = "")) %>%
  dplyr::arrange(desc(n_Excl_AB, n_Incl_AB)) %>% 
  dplyr::mutate(common_dysbiosis_ExclAB = ifelse(n_Excl_AB >= 3, "Yes", "-"),
                common_dysbiosis_InclAB = ifelse(n_Incl_AB >= 3, "Yes", "-"))  %>%
  tidyr::separate(feature, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", remove = FALSE) %>%
  dplyr::mutate(kingdom = substring(kingdom, 4),
                phylum = substring(phylum, 4),
                class = substring(class, 4),
                order = substring(order, 4),
                family = substring(family, 4),
                genus = substring(genus, 4), 
                species = substring(species, 4)) %>% 
  dplyr::select(kingdom, phylum, class, order, family, genus, species, n_Excl_AB, diseases_Excl_AB, common_dysbiosis_ExclAB, n_Incl_AB, diseases_Incl_AB, common_dysbiosis_InclAB) %>% 
  as.data.frame()

write.xlsx(n_shared_associations_wide, file = "Results/Supplementary Table 6. Shared dysbiosis.xlsx", sheetName = "Shared species", col.names = TRUE, row.names = FALSE, append = FALSE)


# Summarize the common dysbiosis by phylum and genera
n_shared_associations_final %>% 
  dplyr::filter(common_dysbiosis_SET2 == "Yes") %>% 
  dplyr::group_by(phylum) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::arrange(desc(n))

n_shared_associations_final %>% 
  dplyr::filter(common_dysbiosis_SET3 == "Yes") %>% 
  dplyr::group_by(phylum) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::arrange(desc(n))







# Create co-abundance network ----
#------------------------------------------------#
#                                                #
#              VISUALIZING NETWORKS              # 
#                                                #
#------------------------------------------------#

# Read igraph object
igraph_raw <- readRDS("RData archive/SpiecEasi_igraph_species_p10d001.rds")

# Remove nodes with no connections 
igraph_clean <- delete.vertices(igraph_raw, which(igraph::degree(igraph_raw, mode = "all") == 0))

# Metadata
ggnet_data <- ggnet2(igraph_clean)$data %>%
  dplyr::left_join(CLR_mean_value, by = c("label" = "taxa")) %>%
  dplyr::mutate(CLR_mean_value_normalized = (meanCLR-min(meanCLR))/(max(meanCLR)-min(meanCLR))) 

# Function to visualize the associations in a co-occurrence network
return_network_plot <- function(factor_of_interest, factor_name = "", p_adjusted = "SET2", 
                                positive_estimate = "Enriched in disease", negative_estimate = "Decreased in disease"){
  
  # Subset MWAS data
  MWAS_subset <- MWAS_data_clean %>%
    dplyr::filter(analysis_factor == factor_of_interest, 
                  adjusted == p_adjusted) %>%
    dplyr::select(feature, Estimate, Node_color)

  # Modifivations to the input data
  ggnet_data_disease <- ggnet_data %>%
    dplyr::left_join(MWAS_subset, by = c("label" = "feature")) %>%
    dplyr::mutate(effect_color = case_when(is.na(Estimate) ~ "No association", 
                                           Estimate >= 0 ~ positive_estimate, 
                                           Estimate < 0 ~ negative_estimate),
                  effect_color = factor(effect_color, levels = c(negative_estimate, "No association", positive_estimate)), 
                  Node_color = ifelse(is.na(Node_color) == TRUE, "gray75", Node_color)) %>%
    tidyr::separate(label, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>%
    dplyr::mutate(CLR_mean_value_normalized = CLR_mean_value_normalized,
                  phylum = substring(phylum, 4),
                  phylum = ifelse(phylum %in% c("Bacteroidetes", "Firmicutes", "Proteobacteria", "Actinobacteria", 
                                                "Verrucomicrobia"), phylum, "Other"),
                  phylum = factor(phylum, levels = c("Bacteroidetes", "Firmicutes", "Proteobacteria", "Actinobacteria", 
                                                     "Verrucomicrobia", "Other")))
  
  factor_title_name <- ifelse(factor_name == "", factor_of_interest, factor_name)
  p_adjusted_name <- ifelse(p_adjusted == "SET2", "Not adjusted to AB usage", "Adjusted to AB usage")
  
  # Visualize the network
  set.seed(8)
  disease_net <- ggnet2(igraph_clean, 
                        mode = "fruchtermanreingold",
                        node.color = ggnet_data_disease$effect_color,
                        edge.size = 0.1,
                        node.size = ggnet_data_disease$CLR_mean_value_normalized,
                        edge.color = "gray50",
                        size.cut = 5,
                        max_size = 3, 
                        #edge.alpha = 0.7,
                        node.alpha = 0.8) + 
    scale_color_manual(name = "", values = c("navy blue", "gray75", "#ff0303")) + 
    ggtitle(paste("Associations with ", factor_title_name, " (", p_adjusted_name, ")", sep = "")) + 
    theme(legend.position = "bottom",
          legend.text = element_text(size = 20),
          aspect.ratio = 1, 
          plot.title = element_text(size = 20, hjust = 0.5)) +
    guides(size = FALSE,
           color = guide_legend(override.aes = list(size = 4)))

  # Change color scheme
  if (sum(ggnet_data_disease$effect_color == negative_estimate) == 0){
    disease_net <- disease_net +
      scale_color_manual(name = "", values = c("gray75", "#ff0303"))
  } else if (sum(ggnet_data_disease$effect_color == positive_estimate) == 0){
    disease_net <- disease_net +
      scale_color_manual(name = "", values = c("navy blue",  "gray75"))
  } else{
    disease_net <- disease_net +
      scale_color_manual(name = "", values = c("navy blue",  "gray75", "#ff0303"))
  }

  disease_net_name_png <- paste("Figures/Networks/Species_Network_", factor_title_name, "_", p_adjusted, ".png", sep = "")
  ggsave(plot = disease_net, filename = disease_net_name_png, width = 10, height = 10)
  
  return(disease_net)
}

# Disease networks without AB adjustment
M10_SET2 <- return_network_plot("M10", p_adjusted = "SET2") 
I10_SET2 <- return_network_plot("I10", p_adjusted = "SET2") 
I11_SET2 <- return_network_plot("I11", p_adjusted = "SET2") 
F41_SET2 <- return_network_plot("F41", p_adjusted = "SET2") 
K50_SET2 <- return_network_plot("K50", p_adjusted = "SET2") 
K21_SET2 <- return_network_plot("K21", p_adjusted = "SET2") 
K51_SET2 <- return_network_plot("K51", p_adjusted = "SET2") 
E11_SET2 <- return_network_plot("E11", p_adjusted = "SET2") 
K58_SET2 <- return_network_plot("K58", p_adjusted = "SET2") 
C18_SET2 <- return_network_plot("C18", p_adjusted = "SET2") 


# Disease networks with AB adjustment
M10_SET3 <- return_network_plot("M10", p_adjusted = "SET3") 
I10_SET3 <- return_network_plot("I10", p_adjusted = "SET3") 
I11_SET3 <- return_network_plot("I11", p_adjusted = "SET3") 
F41_SET3 <- return_network_plot("F41", p_adjusted = "SET3") 
K50_SET3 <- return_network_plot("K50", p_adjusted = "SET3") 
K21_SET3 <- return_network_plot("K21", p_adjusted = "SET3") 
K51_SET3 <- return_network_plot("K51", p_adjusted = "SET3") 
E11_SET3 <- return_network_plot("E11", p_adjusted = "SET3") 
K58_SET3 <- return_network_plot("K58", p_adjusted = "SET3") 
C18_SET3 <- return_network_plot("C18", p_adjusted = "SET3") 







