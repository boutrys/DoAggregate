#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=11) {
  stop("11 arguments must be supplied ()",
   call.=FALSE)
}


####### Parameter of main function
id <- args[1]
path_data <- args[2]
max_var <- args[3]
rare_maf_threshold <- args[4] #thresold to be used to separate rare VS common variants
#covariates for aggregation model
var_covariate <- args[5] #put your covariates here
var_covariate_pca <- args[6] #If set to true, covariates used will be the PCA
version <- args[7]
path_Excalibur_function <- args[8]
path_output <- args[9]
time_path <- args[10]
RLIBRARIES <- args[11]


max_var <- as.integer(max_var)
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


# Let's be sure Excalibur use correct R libraries
setwd("/tmp/")
.libPaths(RLIBRARIES,FALSE)
print(.libPaths())

library(KATSP) #install_github("baolinwu/KATSP")
library(minqa)
library(CompQuadForm)
library(tidyverse)
library(SKAT)
library(AssotesteR)
library(CATT)
library(iGasso)
library(WGScan)
library(REBET)
library(RVtests)

#Load excalibur functions
source(paste(path_Excalibur_function, "stat_Framework.R", sep = ""))
source(paste(path_Excalibur_function, "ADA_test.R", sep = ""))
source(paste(path_Excalibur_function, "parallel_region_analysis.R", sep = ""))


if(version != "optimal" && version != "all"){
  name_test <- stat_Framework(get_test = TRUE,
                              version = version)
  if(length(name_test) != length(version)){
    stop("this version is currently not available, please use version == 0 (optimal)")
  }
}


list_id <- strsplit(split = "/", id)[[1]]
for(i in 1:length(list_id)){
  res <- readRDS(paste(path_data, list_id[i], "_data.rds", sep = ""))
  focus_data_patient <- res[[1]]
  focus_data_control <- res[[2]]
  focus_genotype_matrix <- res[[3]]
  region_of_interest <- res[[4]]
  Nbr_patient <- res[[5]]
  Nbr_control <- res[[6]]




  parallel_region_analysis(focus_data_patient = focus_data_patient,
                            focus_data_control = focus_data_control,
                            Nbr_patient = Nbr_patient,
                            Nbr_control = Nbr_control,
                            Binary = TRUE,
                            focus_genotype_matrix = focus_genotype_matrix,
                            max_var = max_var,
                            var_covariate = var_covariate,
                            var_covariate_pca = var_covariate_pca,
                            rare_maf_threshold = rare_maf_threshold,
                            pvalue_threshold = 0.05,
                            version = version,
                            path_output = path_output,
                            region_id = list_id[i],
                            region_of_interest = region_of_interest,
                            time_path = time_path)
  print(paste("FINISH ANALYSIS OF ", list_id[i], "_data.rds", sep = ""))
}



