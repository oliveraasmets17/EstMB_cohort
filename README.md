# EstMB_cohort

The following scripts were used for the results of "Gut metagenome associations with extensive digital health data in a volunteer-based Estonian microbiome cohort".

Each script loads the necessary R packages individually. 
The scripts should be run in the following order. 

If questions, please contact oliver.aasmets@ut.ee. 


Scripts to be run before analyzing microbiome – phenotype associations    
1) EstMB_help_calculate_diversityMetrics.R – script that calculates alpha diversity metrics (observed richness and Shannon diversity index) based on the given count table
2) EstMB_help_calculate_CLR_distance.R - script that calculates the Aitchison distance matrix based on the given count table.

Scripts that analyzing microbiome – phenotype associations    
3) EstMB_run_PCA_analysis.R – Runs PCA on the distance matrix obtained by the script EstMB_help_calculate_CLR_distance.R. Visualizes the PCA biplot according to the most dominant genera.     
4) EstMB_run_PERMANOVA_analysis.R – Runs PERMANOVA on the distance matrix obtained by the script EstMB_help_calculate_CLR_distance.R. Arguments include a character vector of names of factors to be analyzed, data frame of phenotype data and character indicating, which features are to be used as covariates.    
5) EstMB_run_alphaDiversity_analysis.R - Calculates Spearman correlations with the phenotypes of interest using the diversity metrics obtained by the script EstMB_help_calculate_diversityMetrics.R.     
6) EstMB_run_MWAS_analysis.R – Runs ALDEx2 models to detect differentially abundant microbial features. Arguments include a count table of microbial features, a character vector of names of factors to be analyzed, data frame of phenotype data and character indicating, which features are to be used as covariates.    
7) EstMB_run_stepwiseRDA_analysis.R – Runs stepwise RDA analysis. Takes raw counts, phenotype data and character vector of factors of interest as input.     
8) EstMB_run_SpiecEasi_analysis.R – Runs SpiecEasi function. Takes read count matrix as input.    
9) EstMB_run_tidymodels.R –Runs machine learning models using the tidymodels framework. Arguments include a count table of microbial features, data frame of phenotype data, ML algorithm name, seed for train-test split, grid_size for model tuning, preprocessing step for microbial features and character indicating, which features are to be used as predictors.    

Scripts to be run after analyzing microbiome – phenotype associations    
10) EstMB_visualize_predictiveModellingResults.R – Script that visualizes the AUROC results of the predictive models produced by EstMB_run_tidymodels.R.    
11) EstMB_visualize_MWASresults.R – script that visualizes the overlap in disease-taxa associations between various diseases. Additionally visualizes the co-occurrence networks according to the identified associations.     
12) EstMB_visualize_diversityResults.R - Script that visualizes the results of the diversity-analysis results produced by EstMB_run_PERMANOVA_analysis.R and EstMB_run_alphaDiversity_analysis.R.    

Other scripts for visualizations    
13) EstMB_analyze_MBlandscape.R – Script that identifies core genera, summary statistics for taxonomic abundances and visualized the participants by county.     
14) EstMB_visualize_omicsOverlap.R – Script that visualized the overlap of different omics data in EstBB cohort.    
15) EstMB_visualize_replicates.R - Script that visualizes the Aitchison distance between biological replicates and between unrelated subjects.    



