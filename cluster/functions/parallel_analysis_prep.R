############################################# Function run general analysis on several experiments #############################################
# The function allow perform all the desired analysis 
# Input : Nbr_exp: Number of experiment to be performed
#         path_to_save_time_plot: (optional) a string indicating where the plot must be saved 
#                            if no path indicated, then results won't be stored
#         class_to_plot: (optional) a string with 4 different possible values in order to indicate which class of test are of interest
#                           default value = "all", otherwise could be "burden", "variance-component" or "omnibus"
# Output : plot from ggplot2
parallel_analysis_prep <- function(tmp_folder_input = "",
                          KindOfRegion = c(),
                          filtering_criteria = list(),
                          filtering_level = list(),
                          filtering_direction = list(),
                          max_var = 100000,
                          remove_private = FALSE,
                          tmp_path_pathway = "",
                          tmp_folder_output = ""){
  
  start <- Sys.time()
  #Import data needed
  data_patient_annotated <- readRDS(paste(tmp_folder_input, "data_patient_annotated.rds", sep = ""))
  data_control_annotated <- readRDS(paste(tmp_folder_input, "data_control_annotated.rds", sep = ""))
  genotype_matrix <- readRDS(paste(tmp_folder_input, "clean_genotype_matrix.rds", sep = ""))

  max_var <- as.integer(max_var)
  patient_list <- unique(data_patient_annotated$sample)
  control_list <- unique(data_control_annotated$sample)
  Nbr_patient <- length(patient_list)
  Nbr_control <- length(control_list)
  
  
  ########################### HYPOTHERES GENERATION BASED ON FILTERING #############################

  res_filtering <- parallel_filtering_option(KindOfRegion = KindOfRegion,
                                    filtering_criteria = filtering_criteria,
                                    filtering_level = filtering_level)
  #Prepare SUBfolder to store results
  vec_hypotheses <- res_filtering[[1]]
  Nbr_exp <- res_filtering[[2]]
  KindOfRegion <- res_filtering[[3]]

  folder_name_experience <- paste(vec_hypotheses$genetic_region, "/", sep = "")
  output_folder_directory <- paste(tmp_folder_output,folder_name_experience, sep = "")
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
  
  #Going though all the different experiments
  for(z in 1:Nbr_exp){
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
    if(name_of_experience[z] == "over_represented_pathway" || name_of_experience[z] == "pathway"){
      print(tmp_path_pathway)
      data_tmp <- readRDS(paste(tmp_path_pathway, "pathway_to_test.rds", sep = ""))
      over_represented_pathway <- data_tmp[[1]]
      tmp_table_gene <- data_tmp[[2]]
      region_of_interest <- over_represented_pathway
      if(length(region_of_interest) < 1){
        region_of_interest <- NA
        tmp_table_gene <- NA
      }
    }else if(name_of_experience[z] == "sig_aggregate_gene"){
      target_region <- readRDS(paste(tmp_path_pathway, "target_region.rds", sep = ""))
      tmp_region <- target_region
      if(length(tmp_region) > 1){
        region_of_interest <- tmp_region
      }else{
        region_of_interest <- NA
      }
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
    
      
    ########################### SINGLE ANALYSIS #############################
    Nbr_regions <- length(which(!is.na(region_of_interest)))
    if(Nbr_regions > 0){
      dir.create(paste(tmp_folder_output,vec_hypotheses$genetic_region[z], sep = ""))
      path_to_save_plot <- output_folder_directory[z]
      tmp_save <- paste(path_to_save_plot, "data_to_analyze/", sep = "")
      dir.create(tmp_save)
      if(name_of_experience[z] == "over_represented_pathway" || name_of_experience[z] == "pathway"){
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
    
    



    #  ITERATION ACCROSS REGIONS OF INTEREST  #
    time_count <- 0
    if(Nbr_regions > 0){
          for (loop in 1:Nbr_regions) {
      #region selection
      time_1 <- Sys.time()
      if(name_of_experience[z] == "over_represented_pathway" || name_of_experience[z] == "pathway"){
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
      }else if(name_of_experience[z] == "pathway"){
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

      if(dim(focus_genotype_matrix)[2] != 0){
        if(dim(focus_genotype_matrix)[2] < max_var){
            time_1 <- Sys.time()
            ###save tmp files 
            res <- list(focus_data_patient, 
                        focus_data_control, 
                        focus_genotype_matrix,
                        region_of_interest[loop],
                        Nbr_patient, 
                        Nbr_control)
            saveRDS(res, paste(tmp_save, loop,"_data.rds", sep = ""))
            time_2 <- Sys.time()
            dif_time <- time_2 - time_1

        }
        time_1 <- Sys.time()
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
        time_2 <- Sys.time()
        dif_time <- time_2 - time_1

      }#end of IF variant left in region
    }#end of for going through all regions of interest
    saveRDS(General_info, paste(path_to_save_plot,"general_info_singleAnalysis.rds",sep=""))   
    global_info_experiment$Nbr_region[z] <- length(region_of_interest)
    }#end of if region to be analyzed for experiment

  }#end of for going through all z experiment

  #save global results

  if(name_of_experience[z] == "pathway" | name_of_experience[z] == "over_represented_pathway"){
    saveRDS(global_info_experiment,paste(tmp_folder_output,"global_info_experiment_pathway.rds",sep = ""))
    stop <- Sys.time()
    tot <- stop - start
    saveRDS(tot, paste(tmp_folder_output, "time_analysis_preparation_pathway.rds", sep = ""))
  }else{
    saveRDS(global_info_experiment,paste(tmp_folder_output,"global_info_experiment.rds",sep = ""))
    stop <- Sys.time()
    tot <- stop - start
    saveRDS(tot, paste(tmp_folder_output, "time_analysis_preparation.rds", sep = ""))
  }


}#end of function