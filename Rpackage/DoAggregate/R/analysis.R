#' Run General Analysis on Several Experiments
#'
#' @description This function performs various genetic analyses across multiple experiments, incorporating filtering,
#' hypothesis generation, and statistical testing. It allows for customization of input data, filtering criteria, analysis types, and output storage.
#'
#' @param data_patient_annotated A data frame containing annotated patient data.
#' @param data_control_annotated A data frame containing annotated control data.
#' @param genotype_matrix A matrix representing the genotype data.
#' @param KindOfRegion A vector indicating the regions of interest for the analysis.
#' @param filtering_criteria A list specifying filtering criteria for the data.
#' @param filtering_level A list specifying the level of filtering to apply.
#' @param filtering_direction A list indicating the filtering direction.
#' @param min_var Integer. The minimum number of variants required for a region to be analyzed. Default is 2.
#' @param max_var Integer. The maximum number of variants allowed for a region to be analyzed. Default is 10000.
#' @param remove_private Logical. If TRUE, removes private variants from the analysis. Default is FALSE.
#' @param over_represented_pathway A list or vector indicating pathways that are over-represented.
#' @param version Character. Specifies the version of the statistical framework to use. Default is "optimal".
#' @param var_covariate Numeric. A variable used for covariate adjustment.
#' @param var_covariate_pca Numeric. A variable used for PCA-based covariate adjustment.
#' @param rare_maf_threshold Numeric. The minor allele frequency threshold for defining rare variants. Default is 0.01.
#' @param pvalue_threshold Numeric. The p-value threshold for statistical significance. Default is 0.05.
#' @param Binary Logical. If TRUE, applies a binary framework for analysis. Default is TRUE.
#' @param path_results Character. Path to save detailed results of each experiment.
#' @param path_general_results Character. Path to save the general results summary.
#'
#' @return A list containing the results of the analyses, including general experiment information and plots.
#'
#' @examples
#' # Example usage:
#' @export
analysis <- function(data_patient_annotated = c(),
                     data_control_annotated = c(),
                     genotype_matrix = c(),
                     KindOfRegion = c(),
                     filtering_criteria = list(),
                     filtering_level = list(),
                     filtering_direction = list(),
                     min_var = 2,
                     max_var = 10000,
                     remove_private = FALSE,
                     over_represented_pathway = c(),
                     version = "optimal",
                     var_covariate = 0,
                     var_covariate_pca = 1,
                     rare_maf_threshold = 0.01,
                     pvalue_threshold = 0.05,
                     Binary = TRUE,
                     path_results = "",
                     path_general_results = ""){
  #TODO in input parameter change the last one tmp_folder_output
  #Then integrate function into main.R

  start <- Sys.time()
  max_var <- as.integer(max_var)
  patient_list <- unique(data_patient_annotated$sample)
  control_list <- unique(data_control_annotated$sample)
  Nbr_patient <- length(patient_list)
  Nbr_control <- length(control_list)


  ########################### HYPOTHERES GENERATION BASED ON FILTERING #############################
  res_filtering <- filtering_option(KindOfRegion = KindOfRegion,
                                    filtering_criteria = filtering_criteria,
                                    filtering_level = filtering_level)
  #Prepare SUBfolder to store results
  vec_hypotheses <- res_filtering[[1]]
  Nbr_exp <- res_filtering[[2]]
  KindOfRegion <- res_filtering[[3]]
  name_of_experience <- vec_hypotheses$genetic_region

  #Keep general info about each experiments
  global_info_experiment <- data.frame(matrix(NA, ncol = 15, nrow = Nbr_exp))
  colnames(global_info_experiment) <- c("name_of_experiment",
                                        "Nbr_region",
                                        "Nbr_region_analyzed",
                                        "Nbr_significant_region",
                                        "Nbr_region_not_analyzed",
                                        "Proportion_Nbr_significant_region_VS_Nbr_region_analyzed",
                                        "Proportion_Nbr_region_analyzed_VS_Nbr_region",
                                        "Nbr_position_unique",
                                        "Nbr_position_unique_patient",
                                        "Nbr_position_unique_control",
                                        "Nbr_variant_unique",
                                        "Nbr_variant_unique_patient",
                                        "Nbr_variant_unique_control",
                                        "Total_variant_patient",
                                        "Total_variant_control")
  global_info_experiment$name_of_experiment <- name_of_experience

  #Computational time
  name_time_data <- c(paste("tot_single_analysis", seq(1:Nbr_exp), sep = "_"),
                      "data_preparation",
                      "tot_gene_set_analysis",
                      "total_time")
  time_data <- data.frame(matrix(NA, ncol = length(name_time_data), nrow = 1))
  colnames(time_data)<- name_time_data


  #Going though all the different experiments
  for(z in 1:Nbr_exp){
    ###Create directory to store results
    dir.create(paste(path_results,vec_hypotheses$genetic_region[z], sep = "/"))
    output_folder_directory <- paste(path_results,vec_hypotheses$genetic_region[z], sep = "/")


    ########################### HYPOTHERES GENERATION #############################
    results_hypotheses <- hypotheses_generation(data_patient_annotated = data_patient_annotated,
                                                data_control_annotated = data_control_annotated,
                                                filtering_criteria = vec_hypotheses$filtering_criteria[z],
                                                filtering_level = vec_hypotheses$filtering_level[z],
                                                filtering_direction = filtering_direction,
                                                level_of_analysis = KindOfRegion[z])
    filter_annotate_data_patient <- results_hypotheses[[1]]
    filter_annotate_data_control <- results_hypotheses[[2]]
    region_of_interest <- results_hypotheses[[3]]


    #Only analyze over-represented region
    if(name_of_experience[z] == "over_represented_pathway"){
      tmp_table_gene <- over_represented_pathway[[2]]
      region_of_interest <- over_represented_pathway[[1]]
      if(length(region_of_interest) < 1){
        region_of_interest <- NA
        tmp_table_gene <- NA
      }
    }else if(name_of_experience[z] == "sig_aggregate_gene"){
      stop("still not developped")
      #target_region <- readRDS(paste(tmp_path_over_represented, "target_region.rds", sep = ""))
      #tmp_region <- target_region
      #if(length(tmp_region) > 1){
      #  region_of_interest <- tmp_region
      #}else{
      #  region_of_interest <- NA
      #}
    }


    global_info_experiment$Nbr_position_unique[z] <- length(unique(c(filter_annotate_data_patient$position,
                                                                     filter_annotate_data_control$position)))
    global_info_experiment$Nbr_position_unique_patient[z] <- length(unique(filter_annotate_data_patient$position))
    global_info_experiment$Nbr_position_unique_control[z] <- length(unique(filter_annotate_data_control$position))
    global_info_experiment$Nbr_variant_unique[z] <- length(unique(c(filter_annotate_data_patient$variant_id,
                                                                    filter_annotate_data_control$variant_id)))
    global_info_experiment$Nbr_variant_unique_patient[z] <- length(unique(filter_annotate_data_patient$variant_id))
    global_info_experiment$Nbr_variant_unique_control[z] <- length(unique(filter_annotate_data_control$variant_id))
    global_info_experiment$Total_variant_patient[z] <- dim(filter_annotate_data_patient)[1]
    global_info_experiment$Total_variant_control[z] <- dim(filter_annotate_data_control)[1]


    ########################### ANALYSIS #############################
    Nbr_regions <- length(which(!is.na(region_of_interest)))
    if(Nbr_regions > 0){
      if(name_of_experience[z] == "over_represented_pathway" || name_of_experience[z] == "over_represented_GO"){
        annotation_name <- c("Nbr_test_applied",
                             "name_genetic_region",
                             "gene_involved",
                             "Nbr_patient_with_mutation",
                             "Nbr_control_with_mutation",
                             "Proportion_patients_with_mutation",
                             "Proportion_controls_with_mutation",
                             "Total_nbr_variant",
                             "Total_nbr_variant_in_patients",
                             "Total_nbr_variant_in_controls",
                             "Nbr_unique_variant",
                             "Nbr_unique_variant_in_patient",
                             "Nbr_unique_variant_in_controls",
                             "Nbr_unique_variant_in_both",
                             "Total_nbr_heterozygous_variant_patient",
                             "Total_nbr_homozygous_variant_patient",
                             "Total_nbr_heterozygous_variant_control",
                             "Total_nbr_homozygous_variant_control",
                             "Proportion_of_reference_homozygous",
                             "Proportion_of_heterozygous",
                             "Proportion_of_alternatif_homozygous")
        General_info <- data.frame(matrix(data=NA,nrow = Nbr_regions, ncol = length(annotation_name)))
        colnames(General_info) <- annotation_name
        General_info$gene_involved <- tmp_table_gene
      }else{
        annotation_name <- c("Nbr_test_applied",
                             "name_genetic_region",
                             "Nbr_patient_with_mutation",
                             "Nbr_control_with_mutation",
                             "Proportion_patients_with_mutation",
                             "Proportion_controls_with_mutation",
                             "Total_nbr_variant",
                             "Total_nbr_variant_in_patients",
                             "Total_nbr_variant_in_controls",
                             "Nbr_unique_variant",
                             "Nbr_unique_variant_in_patient",
                             "Nbr_unique_variant_in_controls",
                             "Nbr_unique_variant_in_both",
                             "Total_nbr_heterozygous_variant_patient",
                             "Total_nbr_homozygous_variant_patient",
                             "Total_nbr_heterozygous_variant_control",
                             "Total_nbr_homozygous_variant_control",
                             "Proportion_of_reference_homozygous",
                             "Proportion_of_heterozygous",
                             "Proportion_of_alternatif_homozygous")
        General_info <- data.frame(matrix(data=NA,nrow = Nbr_regions, ncol = length(annotation_name)))
        colnames(General_info) <- annotation_name
      }
    }else{
      Nbr_regions <- 0
    }


    #initiate dataframe to collect results
    name_test <- stat_Framework(Binary = Binary,
                                get_test = TRUE,
                                version = version)
    if(version != "optimal" && version != "all"){
      if(length(name_test) != length(version)){
        stop("this version is currently not available, please use version == 0 (optimal)")
      }
    }
    Nbr_test <- length(name_test)
    summary_pvalue <- data.frame(
      name = name_test
    )
    summary_pvalue_excell <- data.frame()
    summary_pvalue_excell <- summary_pvalue
    summary_time_excell <- summary_pvalue


    #  ITERATION ACCROSS REGIONS OF INTEREST  #
    time_count <- 0
    if(Nbr_regions > 0){
      TrueRegion <- c()
      for (loop in 1:Nbr_regions) {
        #region selection
        time_1 <- Sys.time()
        if(name_of_experience[z] == "over_represented_pathway"){
          tmp_geneID <- tmp_table_gene[loop]
          results_region_selection <- region_selection(data_patient = filter_annotate_data_patient,
                                                       data_control = filter_annotate_data_control,
                                                       genotype_matrix = genotype_matrix,
                                                       region_type = KindOfRegion[z],
                                                       region_name = region_of_interest[loop],
                                                       geneID = tmp_geneID)
          focus_data_patient <- results_region_selection[[1]]
          focus_data_control <- results_region_selection[[2]]
          focus_genotype_matrix <- results_region_selection[[3]]
        }else{
          results_region_selection <- region_selection(data_patient = filter_annotate_data_patient,
                                                       data_control = filter_annotate_data_control,
                                                       genotype_matrix = genotype_matrix,
                                                       region_type = KindOfRegion[z],
                                                       region_name = region_of_interest[loop])
          focus_data_patient <- results_region_selection[[1]]
          focus_data_control <- results_region_selection[[2]]
          focus_genotype_matrix <- results_region_selection[[3]]
        }

        time_2 <- Sys.time()
        dif_time <- time_2 - time_1
        time_count <- time_count + dif_time

        ##If private variant are removed from analysis
        if(remove_private){
          indiv_per_var <- colSums(focus_genotype_matrix)
          focus_genotype_matrix <- focus_genotype_matrix[,which(indiv_per_var > 1)]
          focus_data_patient <- focus_data_patient[which(focus_data_patient$variant_id %in% colnames(focus_genotype_matrix)),]
          focus_data_control <- focus_data_control[which(focus_data_control$variant_id %in% colnames(focus_genotype_matrix)),]
        }

        if(dim(focus_genotype_matrix)[2] >= min_var){
          ################## Annotate results ####################
          Nbr_total_mutation <- dim(focus_genotype_matrix)[1] * dim(focus_genotype_matrix)[2]
          General_info$name_genetic_region[loop] <- region_of_interest[loop]
          General_info$Nbr_patient_with_mutation[loop] <- length(unique(focus_data_patient$sample))
          General_info$Nbr_control_with_mutation[loop] <- length(unique(focus_data_control$sample))
          General_info$Proportion_patients_with_mutation[loop] <- length(unique(focus_data_patient$sample)) / Nbr_patient
          General_info$Proportion_controls_with_mutation[loop] <- length(unique(focus_data_control$sample)) / Nbr_control
          General_info$Total_nbr_variant_in_patients[loop] <- dim(focus_data_patient)[1]
          General_info$Total_nbr_variant_in_controls[loop] <- dim(focus_data_control)[1]
          General_info$Total_nbr_variant[loop] <- General_info$Total_nbr_variant_in_patients[loop] + General_info$Total_nbr_variant_in_controls[loop]
          General_info$Nbr_unique_variant[loop] <- dim(focus_genotype_matrix)[2]
          General_info$Nbr_unique_variant_in_patient[loop] <- length(unique(focus_data_patient$variant_id))
          General_info$Nbr_unique_variant_in_controls[loop] <- length(unique(focus_data_control$variant_id))
          General_info$Nbr_unique_variant_in_both[loop] <- length(which(duplicated(c(unique(focus_data_patient$variant_id),unique(focus_data_control$variant_id)))))
          General_info$Total_nbr_heterozygous_variant_patient[loop] <- length(which(focus_data_patient$zygosity == "Heterozygous"))
          General_info$Total_nbr_homozygous_variant_patient[loop] <- length(which(focus_data_patient$zygosity == "Homozygous"))
          General_info$Total_nbr_heterozygous_variant_control[loop] <- length(which(focus_data_control$zygosity == "Heterozygous"))
          General_info$Total_nbr_homozygous_variant_control[loop] <- length(which(focus_data_control$zygosity == "Homozygous"))
          General_info$Proportion_of_reference_homozygous[loop] <- length(focus_genotype_matrix[focus_genotype_matrix == 0]) / Nbr_total_mutation
          General_info$Proportion_of_heterozygous[loop] <- length(focus_genotype_matrix[focus_genotype_matrix == 1]) / Nbr_total_mutation
          General_info$Proportion_of_alternatif_homozygous[loop] <- length(focus_genotype_matrix[focus_genotype_matrix == 2]) / Nbr_total_mutation


          if(dim(focus_genotype_matrix)[2] < max_var){

            ###Perform the region analysis
            #TODO change parallel_region_analysis function to output what needed
            res_single_region_analysis <- single_region_analysis(focus_data_patient = focus_data_patient,
                                                                 focus_data_control = focus_data_control,
                                                                 Nbr_patient = Nbr_patient,
                                                                 Nbr_control = Nbr_control,
                                                                 Nbr_test = Nbr_test,
                                                                 Binary = Binary,
                                                                 focus_genotype_matrix = focus_genotype_matrix,
                                                                 max_var = max_var,
                                                                 var_covariate = var_covariate,
                                                                 var_covariate_pca = var_covariate_pca,
                                                                 rare_maf_threshold = rare_maf_threshold,
                                                                 pvalue_threshold = 0.05,
                                                                 version = version,
                                                                 region_of_interest = region_of_interest[loop])
            time_2 <- Sys.time()
            dif_time <- time_2 - time_1
            pvalue <- res_single_region_analysis[[1]]
            res_time <- res_single_region_analysis[[3]]
            #collecting results
            General_info$Nbr_test_applied[which(General_info$name_genetic_region == region_of_interest[loop])] <- Nbr_test - length(which(is.na(pvalue)))
            tmp_pvalue <- data.frame(pvalue)
            colnames(tmp_pvalue) <- region_of_interest[loop]
            summary_pvalue_excell <-  cbind(summary_pvalue_excell,tmp_pvalue)
            tmp_time <- data.frame(res_time)
            colnames(tmp_time) <- region_of_interest[loop]
            summary_time_excell <- cbind(summary_time_excell, tmp_time)
            TrueRegion <- c(TrueRegion, region_of_interest[loop])



          }
          time_1 <- Sys.time()
          time_2 <- Sys.time()
          dif_time <- time_2 - time_1

        }#end of IF variant left in region
      }#end of for going through all regions of interest


      ###Results combination
      nbr_sig_region <- result_combination(General_info = General_info,
                                           TrueRegion = TrueRegion,
                                           summary_pvalue = summary_pvalue,
                                           summary_pvalue_excell = summary_pvalue_excell,
                                           Nbr_test = Nbr_test,
                                           version = version,
                                           Binary = Binary,
                                           pvalue_threshold = pvalue_threshold,
                                           level = name_of_experience[z],
                                           filter_annotate_data_patient = filter_annotate_data_patient,
                                           filter_annotate_data_control = filter_annotate_data_control,
                                           genotype_matrix = genotype_matrix,
                                           path_store_results = output_folder_directory)


      ###Global info collection
      global_info_experiment$Nbr_region[z] <- length(region_of_interest)
      global_info_experiment$Nbr_region_analyzed[z] <- length(TrueRegion)
      global_info_experiment$Nbr_significant_region[z] <- nbr_sig_region
      global_info_experiment$Nbr_region_not_analyzed[z] <- global_info_experiment$Nbr_region[z] - global_info_experiment$Nbr_region_analyzed[z]
      global_info_experiment$Proportion_Nbr_significant_region_VS_Nbr_region_analyzed[z] <- global_info_experiment$Nbr_significant_region[z] / global_info_experiment$Nbr_region_analyzed[z]
      global_info_experiment$Proportion_Nbr_region_analyzed_VS_Nbr_region[z] <- global_info_experiment$Nbr_region_analyzed[z] / global_info_experiment$Nbr_region[z]

    }#end of if region to be analyzed for experiment

  }#end of for going through all z experiment


  ###save global results
  write_xlsx(global_info_experiment, paste(path_general_results, "global_info_analysis.xlsx", sep = "/"))



}#end of function
