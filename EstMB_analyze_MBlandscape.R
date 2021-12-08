

# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load thepackages
library("tidyverse")
library("microbiome")
library("ggthemes")
library("readxl")
library("readr")
library("sf")
library("xlsx")

# Samples with biological replicates filtered out
used_vcodes <- readRDS("RData/vcodes_replicates_removed.rds")

# Samples with extremely low reads to be filtered out
low_read_vcodes <- readRDS("RData/low_read_vcodes.rds")

# Color palette
Geva_colors <- c("#8a2889", "#2484c6", "#e2e7ea") # EstBB colors








# Find the most dominant genera per sample ----
#------------------------------------------------#
#                                                #
#               FIND DOMINANT GENERA             # 
#                                                #
#------------------------------------------------#

# Find dominant genus for each sample
dominant_genus_data <- readRDS("Data/Abundance_genus.rds") %>% 
  tibble::rownames_to_column(var = "genus") %>% 
  tidyr::gather(vkood, value, -genus) %>% 
  dplyr::filter(vkood %in% setdiff(used_vcodes, low_read_vcodes)) %>% 
  dplyr::group_by(vkood) %>% 
  dplyr::mutate(n_total = sum(value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(rel_value = value/n_total*100, 
                help = 1) %>% 
  dplyr::arrange(vkood, desc(rel_value)) %>% 
  dplyr::group_by(vkood) %>% 
  dplyr::mutate(counter = cumsum(help)) %>% 
  dplyr::ungroup() %>%
  dplyr::filter(counter == 1)

# Summarize the proportion of dominant genera
dominant_genus_summarized <- dominant_genus_data %>% 
  dplyr::group_by(genus) %>% 
  dplyr::summarise(p_dominant = n()/length(unique(dominant_genus_data$vkood))*100)

# Top 3 are "g__Prevotella", "g__Clostridium" and "g__Bacteroides"
# Rename all else as "Other"
Dominant_genera_perSample <- dominant_genus_data %>% 
  dplyr::mutate(genus = ifelse(genus %in% c("g__Prevotella", "g__Clostridium", "g__Bacteroides"), genus, "g__Other")) %>% 
  dplyr::rename("dominant_genus" = "genus") %>% 
  dplyr::select(vkood, dominant_genus)

saveRDS(Dominant_genera_perSample, "RData/Dominant_genera_perSample.rds")










# Calcuate summary statistics for various taxonomic levels ----
#------------------------------------------------#
#                                                #
#           CALCULATE SUMMARY STATISTICS         # 
#                                                #
#------------------------------------------------#

# Summary statistics of kingdom-level data
kingdom_data <- readRDS("Data/Abundance_kingdom.rds") %>%
  tibble::rownames_to_column(var = "kingdom") %>%
  tidyr::gather(vkood, value, -kingdom) %>%
  dplyr::filter(vkood %in% setdiff(used_vcodes, low_read_vcodes)) %>% 
  dplyr::group_by(vkood) %>% 
  dplyr::mutate(n_total = sum(value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(rel_value = value/n_total*100)

kingdom_data_summarized <- kingdom_data %>% 
  dplyr::group_by(kingdom) %>% 
  dplyr::summarise(Prevalence = sum(rel_value != 0)/length(unique(phylum_data$vkood))*100, 
                   `Mean Relative Abundance (%)` = mean(rel_value), 
                   `Min Relative Abundance (%)` = min(rel_value), 
                   `Max Relative Abundance (%)` = max(rel_value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(str_detect(kingdom, "Unclassified") == FALSE) %>%
  dplyr::mutate(kingdom_print = substring(kingdom, 4)) %>% 
  dplyr::select(kingdom_print, Prevalence, `Mean Relative Abundance (%)`,  `Min Relative Abundance (%)`,  `Max Relative Abundance (%)`) %>% 
  dplyr::arrange(desc(`Mean Relative Abundance (%)`)) %>% 
  as.data.frame()

kingdom_data_summarized 
# k__Bacteria        98.5   
# k__Viruses          0.547 
# k__Eukaryota        0.354 
# k__Archaea          0.0834

# Save the summary statistics for kingdom
write.xlsx(kingdom_data_summarized, file = "Results/Supplementary Table 2. Taxa abundance summary statistics.xlsx", sheetName = "Kingdom_summary_statistics", 
           col.names = TRUE, row.names = FALSE, append = FALSE)

# Summary statistics of phylum-level data
phylum_data <- readRDS("Data/Abundance_phylum.rds") %>% 
  tidyr::gather(vkood, value, -c("kingdom", "phylum")) %>% 
  dplyr::group_by(vkood) %>% 
  dplyr::mutate(n_total = sum(value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(rel_value = value/n_total*100)  %>% 
  dplyr::filter(vkood %in% setdiff(used_vcodes, low_read_vcodes))

phylum_data_summarized <- phylum_data %>% 
  dplyr::group_by(kingdom, phylum) %>% 
  dplyr::summarise(Prevalence = sum(rel_value != 0)/length(unique(phylum_data$vkood))*100, 
                   `Mean Relative Abundance (%)` = mean(rel_value), 
                   `Min Relative Abundance (%)` = min(rel_value), 
                   `Max Relative Abundance (%)` = max(rel_value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(phylum_print = paste(substring(phylum, 4), " (", substring(kingdom, 4), ")", sep = "")) %>% 
  dplyr::select(phylum_print, Prevalence, `Mean Relative Abundance (%)`,  `Min Relative Abundance (%)`,  `Max Relative Abundance (%)`) %>% 
  dplyr::arrange(desc(`Mean Relative Abundance (%)`)) %>% 
  as.data.frame()

write.xlsx(phylum_data_summarized, file = "Results/Supplementary Table 2. Taxa abundance summary statistics.xlsx", sheetName = "Phylum_summary_statistics", 
           col.names = TRUE, row.names = FALSE, append = TRUE)







# Identify core microbiome on genus level ----
#------------------------------------------------#
#                                                #
#                 CORE MICROBIOME                # 
#                                                #
#------------------------------------------------#

# Calulate summary statistics for each genus
calculate_genus_statistics <- FALSE
if (calculate_genus_statistics == TRUE){
  genus_data <- readRDS("Data/Abundance_genus.rds") %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "vkood") %>%
    tidyr::gather(key, value, -vkood) %>%
    dplyr::group_by(vkood) %>%
    dplyr::mutate(n_reads = sum(value)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(rel_value = value/n_reads*100) %>% 
    dplyr::filter(vkood %in% setdiff(used_vcodes, low_read_vcodes))
  
  genus_statistics <- genus_data %>%
    dplyr::group_by(key) %>%
    dplyr::summarise(prevalence_FIN = sum(rel_value >= 0.1)/length(unique(genus_data$vkood))*100, 
                     prevalence = sum(rel_value > 0)/length(unique(genus_data$vkood))*100, 
                     median_abundance = median(rel_value), 
                     mean_abundance = mean(rel_value))
  
  saveRDS(genus_statistics, "RData/Genus_abundance_statistics.rds")
} else{
  genus_statistics <- readRDS("RData/Genus_abundance_statistics.rds")
}

# Core genera statistics
core_genera_data <- genus_statistics %>%
  dplyr::filter(prevalence_FIN >= 10) %>%
  dplyr::arrange(desc(prevalence_FIN), desc(mean_abundance)) %>% 
  tidyr::separate(key, into = c("k", "p", "c", "o", "f", "g"), sep = " / ", remove = FALSE) %>% 
  dplyr::mutate(genus_print = paste(substring(g, 4), " (", substring(k, 4), ")", sep = "")) %>% 
  dplyr::left_join(dominant_genus_summarized, by = c("key" = "genus"))

core_genera_output <- core_genera_data %>% 
  dplyr::mutate(p_dominant = ifelse(is.na(p_dominant), 0, p_dominant)) %>% 
  dplyr::rename("Mean Relative Abundance (%)" = "mean_abundance", 
                "Prevalence at 0.1 Percent Detection (%)" = "prevalence_FIN",
                "Prevalence" = "prevalence",
                "Proportion of samples the genus is most dominant in (%)" = "p_dominant") %>%
  dplyr::select(genus_print,  Prevalence, `Prevalence at 0.1 Percent Detection (%)`, `Mean Relative Abundance (%)`, `Proportion of samples the genus is most dominant in (%)`) %>% 
  as.data.frame()

# Compare with finnish core
finnish_cohort_core <- read_delim("Data/Comparison/Finnish_genus_core.txt", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::filter(`Present in 1 Percent Core (%)`  == "Yes") %>%
  dplyr::mutate(EstMBcore = ifelse(Genus %in% core_genera$genus_print, "Y", "")) %>%
  dplyr::select(Genus, `Mean Relative Abundance (%)`, EstMBcore)

sum(finnish_cohort_core$Genus %in% core_genera_output$genus_print) # 43 overlapping genera

# Save the summary statistics for core genera
write.xlsx(core_genera_output, file = "Results/Supplementary Table 2. Taxa abundance summary statistics.xlsx", sheetName = "Core_summary_statistics", 
           col.names = TRUE, row.names = FALSE, append = TRUE)








# Identify core microbiome on species level ----
#------------------------------------------------#
#                                                #
#                 CORE MICROBIOME                # 
#                                                #
#------------------------------------------------#


# Calulate summary statistics for each species
calculate_species_statistics <- FALSE
if (calculate_species_statistics == TRUE){
  species_data <- readRDS("RData/Abundance_species.rds") %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "vkood") %>%
    tidyr::gather(key, value, -vkood) %>%
    dplyr::group_by(vkood) %>%
    dplyr::mutate(n_reads = sum(value)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(rel_value = value/n_reads*100)  %>% 
    dplyr::filter(vkood %in% setdiff(used_vcodes, low_read_vcodes))
  
  species_statistics <- species_data %>% 
    dplyr::group_by(key) %>%
    dplyr::summarise(prevalence_FIN = sum(rel_value >= 0.1)/length(unique(species_data$vkood))*100, 
                     prevalence = sum(rel_value > 0)/length(unique(species_data$vkood))*100, 
                     median_abundance = median(rel_value), 
                     mean_abundance = mean(rel_value))
  
  saveRDS(species_statistics, "RData/Species_abundance_statistics.rds")
} else{
  species_statistics <- readRDS("RData archive/Species_abundance_statistics.rds")
}

core_species_data <- species_statistics %>%
  dplyr::filter(prevalence_FIN >= 10) %>%
  dplyr::arrange(desc(prevalence_FIN), desc(mean_abundance)) %>% 
  tidyr::separate(key, into = c("k", "p", "c", "o", "f", "g"), sep = " / ", remove = FALSE) %>% 
  dplyr::mutate(genus_print = paste(substring(g, 4), " (", substring(k, 4), ")", sep = "")) %>% 
  dplyr::left_join(dominant_genus_summarized, by = c("key" = "genus"))









# See distribution of reads mapped to various categories ----
#------------------------------------------------#
#                                                #
#          ANALYZE READS BY CATEGORIES           # 
#                                                #
#------------------------------------------------#

# Core coverage
core_genera_coverage <- genus_data %>%
  dplyr::filter(vkood %in% setdiff(used_vcodes, low_read_vcodes)) %>% 
  dplyr::filter(key %in% core_genera_data$key) %>%
  dplyr::group_by(vkood, n_reads) %>%
  dplyr::summarise(core_reads = sum(value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(core_coverage = core_reads/n_reads*100)

median(core_genera_coverage$core_coverage) # 85.68%
mean(core_genera_coverage$core_coverage) # 84.42%

core_species_coverage <- species_data %>%
  dplyr::filter(vkood %in% setdiff(used_vcodes, low_read_vcodes)) %>% 
  dplyr::filter(key %in% core_species_data$key) %>%
  dplyr::group_by(vkood, n_reads) %>%
  dplyr::summarise(core_reads = sum(value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(core_coverage = core_reads/n_reads*100)

median(core_species_coverage$core_coverage) # 78.78%
mean(core_species_coverage$core_coverage) # 77.83%


# Coverage of chosen species
species_filtered <- rownames(readRDS("RData/Abundance_filtered_species_p10d001.rds")) # 1231

chosen_species_coverage <- readRDS("RData/Abundance_species.rds") %>%
  tibble::rownames_to_column(var = "taxonomy") %>% 
  tidyr::gather(vkood, value, -taxonomy) %>% 
  dplyr::mutate(is_chosen = ifelse(taxonomy %in% species_filtered, "chosen", "exc")) %>% 
  dplyr::filter(vkood %in% used_vcodes) %>% 
  dplyr::group_by(vkood, is_chosen) %>% 
  dplyr::summarize(n_reads = sum(value)) %>% 
  dplyr::ungroup() %>% 
  tidyr::spread(is_chosen, n_reads) %>% 
  dplyr::mutate(chosenSpecies_coverage = chosen/(chosen + exc)*100)

mean(chosen_species_coverage$chosenSpecies_coverage) # 94.69%
median(chosen_species_coverage$chosenSpecies_coverage) # 95.16%









# Visualize MB characteristics in Estonian context ----
#------------------------------------------------#
#                                                #
#    VISUALIZE COUNTIES BY MB CHARACTERISTICS    # 
#                                                #
#------------------------------------------------#

# Date MB sample was taken
MB_sampleTaken_date <- read_excel("Data/query-time_sample_taken.xlsx") %>%
  dplyr::filter(is.na(`CONCATSTR(BiomaterialTaken drawDate)`) == FALSE) %>% # Keep rows that have sample date available
  dplyr::mutate(MB_sample_date = as.Date(`CONCATSTR(BiomaterialTaken drawDate)`)) %>%
  dplyr::select(skood, MB_sample_date)

# County
query_MB_counties <- read_excel("Data/query-BirthPlace.xlsx", sheet = "elukoht 2016-2020") %>%
  dplyr::left_join(MB_sampleTaken_date, by = c("SKOOD" = "skood")) %>% 
  dplyr::mutate(sampleYear = as.numeric(substr(MB_sample_date, 1, 4))) %>% 
  dplyr::filter(aasta == sampleYear) %>% 
  dplyr::rename("skood" = "SKOOD", 
                "county_name" = "maakond") %>% 
  dplyr::select(skood, county_name)

# Aggregate samples collected by county
county_samples_aggregated <- query_MB_counties %>% 
  dplyr::group_by(county_name) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(prop = n/nrow(query_MB_counties)*100)

# Visalize samples by county
download_shFiles <- FALSE
if (download_shFiles == TRUE){
  download.file("https://geoportaal.maaamet.ee/docs/haldus_asustus/maakond_shp.zip", destfile = "maakond_shp.zip")
  unzip("maakond_shp.zip")
}
list.files(pattern = ".shp")

Estonian_county <- st_read("maakond_20210601.shp") %>%
  dplyr::left_join(county_samples_aggregated, by = c("MNIMI" = "county_name"))

# Visualize participants by county
p1 <- ggplot(Estonian_county, aes(fill = prop)) +
  geom_sf() + 
  geom_sf_text(aes(label =  paste(round(prop, 2), " %", sep = ""))) + 
  scale_fill_gradient2(low = Geva_colors[3], mid = Geva_colors[2], high = Geva_colors[1], midpoint = mean(Estonian_county$prop),  guide = FALSE) + 
  theme_void() 

ggsave(plot = p1, filename = "Figures/EstonianMap_sampleCollection.svg", height = 10, width = 10, units = "cm")    

