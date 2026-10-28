#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)<10) {
  stop("10 or 11 arguments must be supplied ()",
   call.=FALSE)
}


####### Parameter of main function

tmp_folder_input <- args[1]
KindOfRegion <- args[2] 
#Perform different analysis based on different filtering criteria which each criteria having possibly several different levels
filtering_criteria <- args[3] #filtering_criteria <- list(c("gnomad_wes_af", "consensus_prediction"))
filtering_level <- args[4] #filtering_level <- list(c(0.03,3), c(0.01, 10))
filtering_direction <- args[5] #filtering_direction <- list(c("-", "+"))
tmp_folder_output <- args[6]
path_Excalibur_function <- args[7]
max_var <- args[8]
remove_private <- args[9]
RLIBRARIES <- args[10]


if(KindOfRegion == 0){
  KindOfRegion <- c("gene", "pathway")
}else if(KindOfRegion == 1){
  KindOfRegion <- "gene"
}else if(KindOfRegion == 2){
  KindOfRegion <- "pathway"
  tmp_path_pathway <- args[11]
}else if(KindOfRegion == 4){
  KindOfRegion <- c("over_represented_pathway")
  tmp_path_pathway <- args[11]
  name_file <- list.files(tmp_path_pathway)
  if(length(name_file) > 0){
      id_pathway <- grep("pathway", name_file)
    if(length(id_pathway)>0){
      path_pathway <- name_file[id_pathway]
    }
  }else{
    stop("useless use of parallel_script_analisi_prep.R, please contact Simon Boutry")
  }
}else if(KindOfRegion == 5){
  KindOfRegion <- "sig_aggregate_gene"
  tmp_path_pathway <- args[11]
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

remove_private <- as.integer(remove_private)
if(remove_private == 0){
  remove_private <- FALSE
}else if(remove_private == 1){
  remove_private <- TRUE
}

# Let's be sure Excalibur use correct R libraries
setwd("/tmp/")
.libPaths(RLIBRARIES,FALSE)
print(.libPaths())


library(RCurl)
library(data.table) 
library(ggdendro)
library(tidyverse)


#Load excalibur functions
source(paste(path_Excalibur_function, "parallel_filtering_option.R", sep = ""))
source(paste(path_Excalibur_function, "hypotheses_generation.R", sep = ""))
source(paste(path_Excalibur_function, "region_selection.R", sep = ""))
source(paste(path_Excalibur_function, "parallel_analysis_prep.R", sep = ""))


#Preparing analysis
output <- parallel_analysis_prep(tmp_folder_input = tmp_folder_input,
                          KindOfRegion = KindOfRegion,
                          filtering_criteria = filtering_criteria,
                          filtering_level = filtering_level,
                          filtering_direction = filtering_direction,
                          max_var = max_var,
                          remove_private = remove_private,
                          tmp_path_pathway = tmp_path_pathway,
                          tmp_folder_output = tmp_folder_output)


