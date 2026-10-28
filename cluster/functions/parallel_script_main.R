#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=27) {
  stop("27 arguments must be supplied ()",
   call.=FALSE)
}


####### Parameter of main function

path_to_patient <- args[1]
path_to_control <- args[2]
know_target_region <- args[3] #if list of gene/pathway/GO term to be tested
importation_method <- args[4]
#quality filters criterias for data preparation
mapqv_threshold <- args[5]
min_coverage <- args[6]
remove_individual <- args[7]
remove_position <- args[8]
path_dir_results <- args[9] #directory where to store results of the analysis
coverage_matrix_path <- args[10]
do_data_preparation <- args[11]
path_data_annotate <- args[12]
path_Excalibur_function <- args[13]
tmp_folder_data_preparation <- args[14]
RLIBRARIES <- args[15]
path_GVCF <- args[16]
sample <- args[17]
chr <- args[18]
pos <- args[19]
reference <- args[20] 
alternative <- args[21] 
zygosity <- args[22]
gene_symbol <- args[23]
gnomad_wes_af <- args[24]
gt_GT <- args[25]
gt_DP <- args[26]
rare_maf_threshold <- args[27]



if(importation_method == 0){
  importation_method <- "default"
}else if(importation_method == 1){
  importation_method <- "gvcf"
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


if(do_data_preparation == 0){
  do_data_preparation <- FALSE
}else{
  do_data_preparation <- TRUE
}


if(zygosity == 0){
  zygosity <- ""
}
mandatory_column_mapping <- list(sample = sample, 
                                  chr = chr, 
                                  pos = pos, 
                                  reference = reference, 
                                  alternative = alternative, 
                                  zygosity = zygosity,
                                  gene_symbol = gene_symbol, 
                                  gnomad_wes_af = gnomad_wes_af,
                                  gt_GT = gt_GT, 
                                  gt_DP = gt_DP)

rare_maf_threshold <- as.numeric(rare_maf_threshold)


# Let's be sure Excalibur use correct R libraries
setwd("/tmp/")
.libPaths(RLIBRARIES,FALSE)
print(.libPaths())
print(getRversion())

library(RCurl)
library(writexl)
library(readxl)
library(data.table) 
library(ggdendro)
library(tidyverse)
library(vcfR)

print("WOOOOP WOOOOOOPPPP")
#Load excalibur functions
source(paste(path_Excalibur_function, "load_input.R", sep = ""))
source(paste(path_Excalibur_function, "checkFormat.R", sep = ""))
source(paste(path_Excalibur_function, "genotypeMatrix.R", sep = ""))
source(paste(path_Excalibur_function, "coverage_matrix_creation.R", sep = ""))
source(paste(path_Excalibur_function, "genotype_covered_matrix.R", sep = ""))
source(paste(path_Excalibur_function, "fullAnnotator.R", sep = ""))
source(paste(path_Excalibur_function, "data_preparation.R", sep = ""))
source(paste(path_Excalibur_function, "parallel_main_function.R", sep = ""))
source(paste(path_Excalibur_function, "gvcf_pretreatment.R", sep = ""))

output_R_function <- parallel_main_function(path_to_patient = path_to_patient,
                          path_to_control = path_to_control,
                          know_target_region = know_target_region,
                          importation_method = importation_method,
                          mapqv_threshold = mapqv_threshold,
                          min_coverage = min_coverage,
                          remove_individual = remove_individual,
                          remove_position = remove_position,
                          path_dir_results = path_dir_results,
                          coverage_matrix_path = coverage_matrix_path,
                          do_data_preparation = do_data_preparation,
                          path_data_annotate = path_data_annotate,
                          tmp_folder_data_preparation = tmp_folder_data_preparation,
                          path_gvcf = path_GVCF,
                          mandatory_column_mapping = mandatory_column_mapping,
                          MAF = rare_maf_threshold)


