#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=8) {
  stop("8 arguments must be supplied ()",
   call.=FALSE)
}


####### Parameter of main function
path_to_results <- args[1]
path_store_results <- args[2]
version <- args[3]
level <- args[4]
path_general <- args[5]
path_Excalibur_function <- args[6]
path_Excalibur_Data <- args[7]
RLIBRARIES <- args[8]


# Let's be sure Excalibur use correct R libraries
setwd("/tmp/")
.libPaths(RLIBRARIES,FALSE)
print(.libPaths())

library(writexl)
library(tidyverse)
library(readxl)
if(level == "gene"){
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(GOSemSim)
  library(ReactomePA)
  library(DOSE)
  library(ComplexHeatmap)
  require(circlize)
  library(Cairo)
}

if(level == "gene"){
  tmp_file <- list.files(path_to_results)
  tmp_gene <- tmp_file[which(tmp_file == "gene")]
  file <- c()
  if(length(tmp_gene) > 0){
    file <- tmp_gene
  }
}else if(level == "pathway"){
  tmp_file <- list.files(path_to_results)
  tmp_pathway <- tmp_file[which(tmp_file %in% c("pathway", "over_represented_pathway"))]
  print(tmp_pathway)
  file <- c()
  if(length(tmp_pathway) > 0){
    file <- tmp_pathway
  }
}


if(length(file) > 0){
  #Load excalibur functions
  source(paste(path_Excalibur_function, "stat_Framework.R", sep = ""))
  source(paste(path_Excalibur_function, "parallel_result_combination.R", sep = ""))
  source(paste(path_Excalibur_function, "over_representation_analysis.R", sep = ""))

  if(version != "optimal" && version != "all"){
    name_test <- stat_Framework(get_test = TRUE,
                                version = version)
    if(length(name_test) != length(version)){
      stop("this version is currently not available, please use version == 0 (optimal)")
    }
  }

  for(i in 1:length(file)){
    input_dir <- paste(path_to_results, file[i], "/", sep = "")
    output_dir <- paste(path_store_results, file[i], "/", sep = "")
    dir.create(output_dir)
    parallel_result_combination(path_to_results = input_dir,
                                version = version,
                                Binary = TRUE,
                                pvalue_threshold = 0.05,
                                level = level,
                                path_to_table = path_Excalibur_Data,
                                path_store_results = output_dir)
  }

}

if(level == "pathway"){
  print(path_to_results)
  ###Retrieve global info
  file <- list.files(path_to_results)
  file <- file[grep("global_info_experiment", file)]
  print(file)
  if(length(file)>1){
      a <- readRDS(paste(path_to_results, file[2], sep = ""))
      b <- readRDS(paste(path_to_results, file[1], sep = ""))
      global_info <- rbind(a,b)
  }else{
      global_info <- readRDS(paste(path_to_results, file, sep = ""))
  }

  ###Retrieve global results
  file <- list.files(path_store_results)
  name_experiment <- file[-grep("general", file)]
  tmp_id_name_exp <- grep("gene", name_experiment)
  tmp_id_name_exp <- c(tmp_id_name_exp, grep("pathway", name_experiment))
  print(tmp_id_name_exp)
  if(length(tmp_id_name_exp) > 0){
      name_experiment <- name_experiment[tmp_id_name_exp]
  }else{
      stop("No folder of interest")
  }
  
  folder <- paste(path_store_results, name_experiment, "/", sep = "")

  for(i in 1:length(folder)){
      print(paste(folder[i], "general_info.xlsx", sep = ""))
      tmp_data <- read_xlsx(paste(folder[i], "general_info.xlsx", sep = ""))
      global_info$Nbr_region_analyzed[which(global_info$name_of_experiment == name_experiment[i])] <- dim(tmp_data)[1]
    
      tmp_list_file <- list.files(folder[i])
      #sig region
      tmp_idx <- which(tmp_list_file == "significant_general_info.xlsx")
      if(length(tmp_idx)>0){
          tmp_data <- read_xlsx(paste(folder[i], "significant_general_info.xlsx", sep = ""))
          global_info$Nbr_significant_region[which(global_info$name_of_experiment == name_experiment[i])] <- dim(tmp_data)[1]
      }else{
          global_info$Nbr_significant_region[which(global_info$name_of_experiment == name_experiment[i])] <- 0
      }
      #not analyzed region
      tmp_idx <- which(tmp_list_file == "general_info_not_analyzed.xlsx")
      if(length(tmp_idx)>0){
          tmp_data <- read_xlsx(paste(folder[i], "general_info_not_analyzed.xlsx", sep = ""))
          global_info$Nbr_region_not_analyzed[which(global_info$name_of_experiment == name_experiment[i])] <- dim(tmp_data)[1]
      }else{
          global_info$Nbr_region_not_analyzed[which(global_info$name_of_experiment == name_experiment[i])] <- 0
      }
  }

  global_info$Proportion_Nbr_significant_region_VS_Nbr_region_analyzed <- global_info$Nbr_significant_region / global_info$Nbr_region_analyzed

  global_info$Proportion_Nbr_region_analyzed_VS_Nbr_region <- global_info$Nbr_region_analyzed / global_info$Nbr_region

  
  write_xlsx(global_info, paste(path_general, "global_info_analysis.xlsx", sep=""))
}