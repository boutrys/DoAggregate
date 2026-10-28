#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop("4 arguments must be supplied ()",
   call.=FALSE)
}


####### Parameter of main function
path_to_global <- args[1]
path_to_results <- args[2]
path_output <- args[3]
RLIBRARIES <- args[4]


# Let's be sure Excalibur use correct R libraries
setwd("/tmp/")
.libPaths(RLIBRARIES,FALSE)
print(.libPaths())

library(writexl)
library(readxl)


###Retrieve global info
file <- list.files(path_to_global)
file <- file[grep("global_info_experiment", file)]
if(length(file)>1){
    a <- readRDS(paste(path_to_global, file[2], sep = ""))
    b <- readRDS(paste(path_to_global, file[1], sep = ""))
    global_info <- rbind(a,b)
}else{
    global_info <- readRDS(paste(path_to_global, file, sep = ""))
}

###Retrieve global results
file <- list.files(path_to_results)
name_experiment <- file[-grep("general", file)]
tmp_id_name_exp <- grep("gene", name_experiment)
tmp_id_name_exp <- c(tmp_id_name_exp, grep("over_represented", name_experiment))
if(length(tmp_id_name_exp) > 0){
    name_experiment <- name_experiment[tmp_id_name_exp]
}else{
    stop("No folder of interest")
}

folder <- paste(path_to_results, name_experiment, "/", sep = "")

for(i in 1:length(folder)){
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

write_xlsx(global_info, paste(path_output, "global_info_analysis.xlsx", sep=""))