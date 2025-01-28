#' Main Function for Genetic Data Analysis Pipeline
#'
#' @name main
#' @description
#' This comprehensive function orchestrates the entire workflow for genetic data analysis. It imports data,
#' applies pre-treatment steps, generates genotype matrices, conducts data preparation, and carries out
#' over-representation analysis. It is designed to handle both patient and control data, apply filtering based
#' on specified criteria, and execute advanced genetic analyses including PCA and hypothesis testing.
#'
#' @param PATH_PATIENTS File path to the patient data.
#' @param PATH_CONTROLS File path to the control data.
#' @param PATH_GVCF File path to genomic VCF data.
#' @param mandatory_column_mapping A list mapping column names in the input data to required format.
#' @param path_results Directory path to store results.
#' @param importation_method Integer indicating the method of data importation (0 for default, 1 for GVCF).
#' @param RLIBRARIES Paths to R libraries used in the analysis.
#' @param list_genes File path to a list of genes of interest, if any.
#' @param rare_maf_threshold Threshold for defining rare minor allele frequency.
#' @param var_covariate Integer indicating whether to use covariates in the analysis.
#' @param var_covariate_pca Boolean indicating whether to perform PCA on covariates.
#' @param min_var Minimum number of variants to include for analysis.
#' @param max_var Maximum number of variants to consider in the analysis.
#' @param do_data_preparation Boolean indicating whether to perform data preparation.
#' @param DO_over_representation Boolean indicating whether to perform over-representation analysis.
#' @param version String indicating the version of the analysis method.
#' @param filtering_criteria List of criteria for filtering the data.
#' @param filtering_level List of levels corresponding to each filtering criterion.
#' @param filtering_direction List indicating the direction of each filtering criterion.
#' @param coverage_matrix_path Optional path to a pre-existing coverage matrix.
#' @param mapqv_threshold Mapping quality threshold.
#' @param min_coverage Minimum coverage threshold.
#' @param remove_individual Boolean indicating whether to remove individuals based on poor data quality.
#' @param remove_position Boolean indicating whether to remove positions based on coverage.
#' @param remove_private Boolean indicating whether to exclude private variants.
#' @param pvalue_threshold P-value threshold for statistical tests.
#' @param Binary Boolean indicating whether the analysis is binary.
#' @param path_Excalibur_function Path to the Excalibur function library.
#' @return No return value, side effects include writing results to files and possibly printing to console.
#' @import RCurl
#' @import writexl
#' @import readxl
#' @import data.table
#' @import MASS
#' @import ggdendro
#' @import tidyverse
#' @import KATSP
#' @import minqa
#' @import CompQuadForm
#' @import SKAT
#' @import AssotesteR
#' @import CATT
#' @import iGasso
#' @import WGScan
#' @import REBET
#' @import RVtests
#' @import Cairo
#' @import vcfR
#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @import GOSemSim
#' @import ReactomePA
#' @import DOSE
#' @import ComplexHeatmap
#' @import circlize
#' @import readr
#' @import dplyr
#' @import grid
#' @importFrom magrittr %>%
#' @export
main <- function(PATH_PATIENTS = "",
                 PATH_CONTROLS = "",
                 PATH_GVCF = "",
                 mandatory_column_mapping = list(sample = "Indiv",
                                                 chr = "CHROM",
                                                 pos = "POS",
                                                 reference = "REF",
                                                 alternative = "ALT",
                                                 zygosity = "",
                                                 gene_symbol = "Gene.ensGene",
                                                 gnomad_wes_af = "gnomad40_AF",
                                                 gt_GT = "gt_GT",
                                                 gt_DP = "gt_DP"),
                 path_results = "",
                 importation_method = 0,
                 RLIBRARIES = .libPaths(),
                 list_genes = 0,
                 rare_maf_threshold = 0.01,
                 var_covariate = 0,
                 var_covariate_pca = 1,
                 min_var = 2,
                 max_var = 10000,
                 do_data_preparation = 1,
                 DO_over_representation = 1,
                 version= "optimal",
                 filtering_criteria = 0,
                 filtering_level = 0,
                 filtering_direction = 0,
                 coverage_matrix_path = 0,
                 mapqv_threshold = 20,
                 min_coverage = 20,
                 remove_individual = 1,
                 remove_position = 1,
                 remove_private = 0,
                 pvalue_threshold = 0.05,
                 Binary = TRUE,
                 path_Excalibur_function = "/Users/sboutry/Documents/excalibur/Rpackage/DoAggregate/R"){

  ### Let's be sure Excalibur use correct R libraries
  #.libPaths(RLIBRARIES,FALSE)




  ###Load excalibur functions
  source(paste(path_Excalibur_function, "checkFormat.R", sep = "/"))
  source(paste(path_Excalibur_function, "gvcf_pretreatment.R", sep = "/"))
  source(paste(path_Excalibur_function, "genotypeMatrix.R", sep = "/"))
  source(paste(path_Excalibur_function, "genotype_covered_matrix.R", sep = "/"))
  source(paste(path_Excalibur_function, "data_preparation.R", sep = "/"))
  source(paste(path_Excalibur_function, "filtering_option.R", sep = "/"))
  source(paste(path_Excalibur_function, "hypotheses_generation.R", sep = "/"))
  source(paste(path_Excalibur_function, "region_selection.R", sep = "/"))
  source(paste(path_Excalibur_function, "analysis.R", sep = "/"))
  if(DO_over_representation == 1){
    source(paste(path_Excalibur_function, "over_representation_analysis.R", sep = "/"))
  }
  source(paste(path_Excalibur_function, "stat_Framework.R", sep = "/"))
  source(paste(path_Excalibur_function, "ADATest.R", sep = "/"))
  source(paste(path_Excalibur_function, "single_region_analysis.R", sep = "/"))
  source(paste(path_Excalibur_function, "result_combination.R", sep = "/"))


  ###Checking input paramters
  if(importation_method == 0){
    importation_method <- "default"
  }else if(importation_method == 1){
    #using GVCF
    importation_method <- "GVCF"
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


  ###DATA Importation and preparation
  ### Prepare main folder to store results
  dir.create(path_results)
  dir.create(paste(path_results, "general", sep = "/"))
  path_general_results <- paste(path_results, "/general/", sep = "")

  print("start reading files")
  if(importation_method == "default"){
    data_patient <- read_tsv(PATH_PATIENTS, col_types = cols(chr = col_character(), reference = col_character(), alternative = col_character(), zygosity = col_character()))
    data_control <- read_tsv(PATH_CONTROLS, col_types = cols(chr = col_character(), reference = col_character(), alternative = col_character(), zygosity = col_character()))
    results_checkFormat <- checkFormat(patient = data_patient,
                                       control = data_control)
    data_patient <- results_checkFormat[[1]]
    data_control <- results_checkFormat[[2]]
  }else if(importation_method == "GVCF"){
    res <- gvcf_pretreatment(path_data = PATH_GVCF,
                             path_output = path_general_results,
                             cases = PATH_PATIENTS,
                             controls = PATH_CONTROLS,
                             mandatory_column_mapping =  mandatory_column_mapping,
                             MAF = 0.01,
                             het_var = c("0/1", "0|1"),
                             homo_var = c("1/1", "1|1"),
                             homo_ref = c("0/0", "0|0"),
                             cleaning = FALSE,
                             remove_individual = TRUE,
                             remove_position = TRUE,
                             min_coverage = 10)
    data_patient <- res[[1]]
    data_control <- res[[2]]
  }else{
    stop("still to do, e.g. vcf files and so on")
  }
  print("end reading files")

  if(list_genes != 0){
    gene_list <- read_tsv(list_genes, col_names = FALSE)
    colnames(gene_list) <- "region"
    data_patient <- data_patient[which(data_patient$gene_symbol %in% gene_list$region),]
    data_control <- data_control[which(data_control$gene_symbol %in% gene_list$region),]
  }



  start_preparation <- Sys.time()
  ########################### DATA PREPARATION #############################
  if(do_data_preparation){
    result_data_preparation <- data_preparation(data_patient = data_patient,
                                                data_control = data_control,
                                                mapqv_threshold = mapqv_threshold,
                                                min_coverage = min_coverage,
                                                remove_individual = remove_individual,
                                                remove_position = remove_position,
                                                path_to_store_plots = path_general_results,
                                                coverage_matrix_path = coverage_matrix_path)
    data_patient_annotated <- result_data_preparation[[1]]
    data_control_annotated <- result_data_preparation[[2]]
    genotype_matrix <- result_data_preparation[[3]]
    missing_annotation <- result_data_preparation[[4]]
    variant_target <- result_data_preparation[[5]]

    if(coverage_matrix_path != 0){
      print("start trying to write data_patient_annotated")
      clean_coverage_matrix <- result_data_preparation[[6]]
      rejected_variant <- result_data_preparation[[7]]
      rejected_data_patient <- result_data_preparation[[8]]
      rejected_data_control <- result_data_preparation[[9]]
      #Potential bug
      #write_tsv(clean_coverage_matrix, file = paste(path_general_results, "clean_coverage_matrix.tsv", sep = ""))
      #write_tsv(rejected_variant, file = paste(path_general_results, "rejected_variants.tsv", sep = ""))
      #write_tsv(rejected_data_patient, file = paste(path_general_results, "rejected_data_patient.tsv", sep = ""))
      #write_tsv(rejected_data_control, file = paste(path_general_results, "rejected_data_control.tsv", sep = ""))

    }
    #Potential bug

    #err <- try(write_tsv(data_patient_annotated, file = paste(path_general_results, "data_patient_annotated.tsv", sep = "")),
    #            silent = TRUE)

    #if(length(grep("error", err)) > 1){
    #  print("MY_ERROR : cant write data_patient_annotated because of out of memory")
    #}
    print("end writing data_control_annotated")
    #write_tsv(data_control_annotated, file = paste(path_general_results, "data_control_annotated.tsv", sep = ""))
    #write_tsv(missing_annotation, file = paste(path_general_results, "missing_annotation.tsv", sep = ""))
    #write_tsv(variant_target, file = paste(path_general_results, "variant_list.tsv", sep = ""))
  }else{
    data_patient_annotated <- data_patient
    data_control_annotated <- data_control
    Nbr_patient <- length(unique(data_patient_annotated$sample))
    Nbr_control <- length(unique(data_control_annotated$sample))
    result_geno_matrix <- genotypeMatrix(data_patient = data_patient_annotated,
                                         data_control = data_control_annotated,
                                         Nbr_patient_init = Nbr_patient,
                                         Nbr_control_init = Nbr_control)
    genotype_matrix <- result_geno_matrix[[1]]
    patient_name <- unique(data_patient_annotated$patient)
    control_name <- unique(data_control_annotated$patient)
    genotype_matrix <- data.frame(genotype_matrix, row.names = c(patient_name,control_name))
    list_variant <- result_geno_matrix[[4]]
    colnames(genotype_matrix) <- list_variant$id
  }
  stop_preparation <- Sys.time()
  tot_preparation <- stop_preparation - start_preparation


  ###Over representation analysis
  if(DO_over_representation == 1){
    gene_list <- unique(c(data_patient_annotated$gene_symbol,
                          data_control_annotated$gene_symbol))
    if (length(gene_list) > 1) {
      res <- over_representation_analysis(geneList = gene_list,
                                          go_of_interest = c("BP", "MF", "CC"),
                                          max_GO_similarity_accepted = 0.5,
                                          cutoff = 0.05,
                                          max_gene_per_item = 20,
                                          max_item_plot = 20,
                                          coocurrence = TRUE,
                                          MSigDB = FALSE,
                                          max_coocurrence = 100,
                                          path_to_store_results = path_general_results)
      over_represented_pathway <- res
    }else{
      over_represented_pathway <- c()
    }
  }


  ###Perform all experiments
  if(DO_over_representation == 1){
    KindOfRegion <- c("gene", "over_represented_pathway")
  }else{
    KindOfRegion <- c("gene")
  }


  ### Launch analysis
  analysis(data_patient_annotated = data_patient_annotated,
           data_control_annotated = data_control_annotated,
           genotype_matrix = genotype_matrix,
           KindOfRegion = KindOfRegion,
           filtering_criteria = filtering_criteria,
           filtering_level = filtering_level,
           filtering_direction = filtering_direction,
           min_var = min_var,
           max_var = max_var,
           remove_private = remove_private,
           over_represented_pathway = over_represented_pathway,
           version = version,
           var_covariate = var_covariate,
           var_covariate_pca = var_covariate_pca,
           rare_maf_threshold = rare_maf_threshold,
           pvalue_threshold = pvalue_threshold,
           Binary = Binary,
           path_results = path_results,
           path_general_results = path_general_results)

}#end of function

