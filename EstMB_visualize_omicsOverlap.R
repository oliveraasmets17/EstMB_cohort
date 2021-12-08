

library("readr")
library("ComplexUpset")

# Read the phenotype data for scodes
full_factor_data <- readRDS(file = "RData/Data_master.rds") 

# Read omics data and create omic-based code sets
omics_samples <- read_delim("Data/EGCUT_DataLayers_2021-03-02.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
CNV_WGS <- unique(omics_samples$`CNV-WGSbased_Maarja_Available_scode_31.01.2020_N=2269`)
CNV_GSA <- unique(omics_samples$`CNV-GSAbased_Maarja_Available_scode_31.01.2020_N=49426`)
WGS <- unique(omics_samples$`WGS-set1_Mart_PostQC_scode_28.01.2020_N=2240`, 
              omics_samples$`WGS-set2_Mart_PostQC_scode_28.01.2020_N=2420`)
WES <- unique(omics_samples$`WES_Mart_PostQC_scode_28.01.2020_N=2445`)
MB <- unique(full_factor_data$skood)

# Gather the sets 
omics_overlap <- data.frame(skood = MB, 
                            stringsAsFactors = FALSE) %>%
  dplyr::mutate(MB = ifelse(skood %in% MB, 1, 0),
                NMR = MB,
                CNV_WGS = ifelse(skood %in% CNV_WGS, 1, 0),
                CNV_GSA = ifelse(skood %in% CNV_GSA, 1, 0),
                WGS = ifelse(skood %in% WGS, 1, 0),
                WES = ifelse(skood %in% WES, 1, 0))

# Create upset-plot of the interaction
overlap_plot <- upset(omics_overlap[, 2:ncol(omics_overlap)], 
                      intersect = colnames(omics_overlap)[2:ncol(omics_overlap)],
                      queries = list(upset_query(set = 'NMR', fill = "#2b7aa1")),
                      stripes = "white", 
                      themes = upset_default_themes(panel.grid = element_blank()),
                      wrap = TRUE) 

ggsave(overlap_plot, file = "Figures/Omics_data_overlap.svg", width = 12, height = 5)
ggsave(overlap_plot, file = "Figures/Omics_data_overlap.png", width = 12, height = 6)





