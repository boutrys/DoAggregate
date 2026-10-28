#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=8) {
  stop("8 arguments must be supplied ()",
   call.=FALSE)
}


####### Parameter of main function
path_to_data <- args[1]
path_data_annotate <- args[2]
path_Excalibur_function <- args[3]
path_tmp_results <- args[4] 
path_results <- args[5]
Aggregate_over_representation <- args[6]
path_Excalibur_Data <- args[7]
RLIBRARIES <- args[8]
MSigDB <- TRUE


data_patient_annotated <- readRDS(paste(path_to_data, "GENE_SET_data_patient_annotated.rds", sep = ""))
data_control_annotated <- readRDS(paste(path_to_data, "GENE_SET_data_control_annotated.rds", sep = ""))
print(data_patient_annotated)


# Let's be sure Excalibur use correct R libraries
setwd("/tmp/")
.libPaths(RLIBRARIES,FALSE)
print(.libPaths())

library(clusterProfiler)
library(org.Hs.eg.db)
library(GOSemSim)
library(writexl)
library(tidyverse)
library(ReactomePA)
library(DOSE)
library(ComplexHeatmap)
require(circlize)
library(Cairo)


gene_list <- unique(c(data_patient_annotated$gene_symbol,
                      data_control_annotated$gene_symbol))
print(gene_list)


if(Aggregate_over_representation == 0){
  source(paste(path_Excalibur_function, "gene_to_pathway.R", sep = ""))
  res <- gene_to_pathway(gene_list = gene_list,
                          output_path = path_results,
                          min_inputGene_inPathway = 8,
                          max_inputGene_inPathway = 8,
                          RatioInputGeneVSGeneInPathway = 0,
                          path_to_database = path_Excalibur_Data)
  pathway_to_test <- res
  saveRDS(pathway_to_test, paste(path_tmp_results, "pathway_to_test.rds", sep = ""))
}else if(Aggregate_over_representation == 1){
  source(paste(path_Excalibur_function, "over_representation_analysis.R", sep = ""))
  if (length(gene_list) > 1) {
    res <- over_representation_analysis(geneList = gene_list,
                              gene_annot = "SYMBOL",
                              go_of_interest = c("BP", "MF", "CC"),
                              max_GO_similarity_accepted = 0.5,
                              cutoff = 0.05,
                              max_gene_per_item = 500,
                              max_item_plot = 20,
                              coocurrence = TRUE,
                              max_coocurrence = 100,
                              MSigDB = TRUE,
                              path_to_table = path_Excalibur_Data,
                              path_to_store_results = paste(path_results, "over_representation_analysis", sep = ""))
    over_represented_pathway <- res
  }else{
    over_represented_pathway <- c()
  }
  saveRDS(over_represented_pathway, paste(path_tmp_results, "pathway_to_test.rds", sep = ""))
}




