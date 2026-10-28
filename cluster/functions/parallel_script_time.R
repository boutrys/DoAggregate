#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=26) {
  stop("26 arguments must be supplied ()",
   call.=FALSE)
}


####### Parameter of main function

path_to_patient <- args[1]
path_to_control <- args[2]
KindOfRegion <- args[3] #gene, pathway, GO
know_target_region <- args[4] #if list of gene/pathway/GO term to be tested
#Perform different analysis based on different filtering criteria which each criteria having possibly several different levels
filtering_criteria <- args[5] #filtering_criteria <- list(c("gnomad_wes_af", "consensus_prediction"))
filtering_level <- args[6] #filtering_level <- list(c(0.03,3), c(0.01, 10))
filtering_direction <- args[7] #filtering_direction <- list(c("-", "+"))
min_var <- args[8] #unique variant to be present within a genetic region for the region to be analyzed
max_var <- args[9] #max nbr unique variant to be within a genetic region for the region to be analyzed
rare_maf_threshold <- args[10] #thresold to be used to separate rare VS common variants
#covariates for aggregation model
var_covariate <- args[11] #put your covariates here
var_covariate_pca <- args[12] #If set to true, covariates used will be the PCA
importation_method <- args[13]
#quality filters criterias for data preparation
mapqv_threshold <- args[14]
min_coverage <- args[15]
remove_individual <- args[16]
remove_position <- args[17]
path_dir_results <- args[18] #directory where to store results of the analysis
coverage_matrix_path <- args[19]
do_data_preparation <- args[20]
path_data_annotate <- args[21]
DO_over_representation <- args[22]
only_over_represented <- args[23]
version <- args[24]
path_Excalibur_function <- args[25]
RLIBRARIES <- args[26]

if(KindOfRegion == 0){
  KindOfRegion <- c("gene", "pathway", "GO")
}else if(KindOfRegion == 1){
  KindOfRegion <- "gene"
}else if(KindOfRegion == 2){
  KindOfRegion <- "pathway"
}else if(KindOfRegion == 3){
  KindOfRegion <- "GO"
}else{
  stop("KindOfRegion argument badly initialized, should be (0=ALL, 1=gene, 2=pathway, 3=GO)")
}

#Perform different analysis based on different filtering criteria which each criteria having possibly several different levels
if(filtering_criteria == 0){
  filtering_criteria <- list() #filtering_criteria <- list(c("gnomad_wes_af", "consensus_prediction"))
}
if(filtering_level == 0){
  filtering_level <- list() #filtering_level <- list(c(0.03,3), c(0.01, 10))
}
if(filtering_direction == 0){
  filtering_direction <- list() #filtering_direction <- list(c("-", "+"))
}

min_var <- as.numeric(min_var)
max_var <- as.numeric(max_var)
rare_maf_threshold <- as.numeric(rare_maf_threshold)

if(var_covariate == 0){
  var_covariate <- 0
}else{
  stop("still to be implemented")
}
if(var_covariate_pca == 0){
  var_covariate_pca <- FALSE
}else{
  var_covariate_pca <- TRUE
}

if(importation_method == 0){
  importation_method <- "default"
}else{
  stop("No other importation method implemented")
}

mapqv_threshold <- as.numeric(mapqv_threshold)
min_coverage <- as.numeric(min_coverage)

if(remove_individual == 0){
  remove_individual <- FALSE
}else{
  remove_individual <- TRUE
}

if(remove_position == 0){
  remove_position <- FALSE
}else{
  remove_position <- TRUE
}

path_dir_results <- paste(path_dir_results, "/Results/", sep="")
dir.create(path_dir_results)

if(do_data_preparation == 0){
  do_data_preparation <- FALSE
}else{
  do_data_preparation <- TRUE
}

if(DO_over_representation == 0){
  DO_over_representation <- FALSE
}else{
  DO_over_representation <- TRUE
}

if(only_over_represented == 0){
  only_over_represented <- FALSE
}else{
  only_over_represented <- TRUE
}

if(version == 0){
  version <- "optimal"
}else{
  stop("this version is currently not available, please use version == 0 (optimal)")
}


#path_data_annotate <- "/data/gehu/sboutry/Excalibur_cluster/Data/"
#version <- "optimal"
#path_Excalibur_function <- "/data/gehu/sboutry/Excalibur_cluster/functions/"


#TO DO add gene_list arg in main function and option to load this list and filter data

# Let's be sure Excalibur use correct R libraries
setwd("/tmp/")
.libPaths(RLIBRARIES,FALSE)
print(.libPaths())


library(RCurl)

library(writexl)
library(readxl)
library(data.table) 
library(GOSemSim) #Bioconductor
library(clusterProfiler) #Bioconductor
library(org.Hs.eg.db) #Bioconductor
library(DOSE) #Bioconductor
library(biomaRt) #Bioconductor
library(ReactomePA) #Bioconductor
library(stringr)
library(ggdendro)
library(KATSP) #install_github("baolinwu/KATSP")
library(minqa)
library(CompQuadForm)
library(scatterpie)
library(tidyverse)
library(SKAT)
library(AssotesteR)
library(DoEstRare)

#Load excalibur functions
source(paste(path_Excalibur_function, "load_input.R", sep = ""))
source(paste(path_Excalibur_function, "checkFormat.R", sep = ""))
source(paste(path_Excalibur_function, "genotypeMatrix.R", sep = ""))
source(paste(path_Excalibur_function, "coverage_matrix_creation.R", sep = ""))
source(paste(path_Excalibur_function, "genotype_covered_matrix.R", sep = ""))
source(paste(path_Excalibur_function, "fullAnnotator.R", sep = ""))
source(paste(path_Excalibur_function, "data_preparation.R", sep = ""))
source(paste(path_Excalibur_function, "filtering_option.R", sep = ""))
source(paste(path_Excalibur_function, "gene_set_analysis.R", sep = ""))
source(paste(path_Excalibur_function, "hypotheses_generation.R", sep = ""))
source(paste(path_Excalibur_function, "single_analysis.R", sep = ""))
source(paste(path_Excalibur_function, "region_selection.R", sep = ""))
source(paste(path_Excalibur_function, "stat_Framework.R", sep = ""))
source(paste(path_Excalibur_function, "ADA_test.R", sep = ""))
source(paste(path_Excalibur_function, "main_function.R", sep = ""))

listOklm <- main_function(path_to_patient = path_to_patient,
                          path_to_control = path_to_control,
                          KindOfRegion = KindOfRegion,
                          know_target_region = know_target_region,
                          filtering_criteria = filtering_criteria,
                          filtering_level = filtering_level,
                          filtering_direction = filtering_direction,
                          min_var = min_var,
                          max_var = max_var,
                          rare_maf_threshold = rare_maf_threshold,
                          var_covariate = var_covariate,
                          var_covariate_pca = var_covariate_pca,
                          importation_method = importation_method,
                          mapqv_threshold = mapqv_threshold,
                          min_coverage = min_coverage,
                          remove_individual = remove_individual,
                          remove_position = remove_position,
                          path_dir_results = path_dir_results,
                          coverage_matrix_path = coverage_matrix_path,
                          do_data_preparation = do_data_preparation,
                          path_data_annotate = path_data_annotate,
                          DO_over_representation = DO_over_representation,
                          only_over_represented = only_over_represented,
                          version = version)


