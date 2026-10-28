############################################# Function run general analysis on several experiments #############################################
# The function allow perform all the desired analysis 
# Input : Nbr_exp: Number of experiment to be performed
#         path_to_save_time_plot: (optional) a string indicating where the plot must be saved 
#                            if no path indicated, then results won't be stored
#         class_to_plot: (optional) a string with 4 different possible values in order to indicate which class of test are of interest
#                           default value = "all", otherwise could be "burden", "variance-component" or "omnibus"
# Output : plot from ggplot2
main_function <- function(path_to_patient = "",
                          path_to_control = "",
                          KindOfRegion = list(),
                          know_target_region = 0,
                          filtering_criteria = list(),
                          filtering_level = list(),
                          filtering_direction = list(),
                          min_var = 10,
                          max_var = 1000,
                          rare_maf_threshold = 0.01,
                          var_covariate = 0,
                          var_covariate_pca = TRUE,
                          importation_method = "default",
                          mapqv_threshold = 20,
                          min_coverage = 20,
                          remove_individual = TRUE,
                          remove_position = TRUE,
                          path_dir_results = getwd(),
                          coverage_matrix_path = 0,
                          do_data_preparation = TRUE,
                          path_data_annotate = "",
                          pvalue_threshold = 0.05,
                          DO_over_representation = TRUE,
                          only_over_represented = TRUE,
                          version = "optimal"){
  
  #Data importation
  if(importation_method == "default"){
    data_patient <- read_tsv(path_to_patient, col_types = cols(chr = col_character(), zygosity = col_character()))
    data_control <- read_tsv(path_to_control, col_types = cols(chr = col_character(), zygosity = col_character()))
    results_checkFormat <- checkFormat(patient = data_patient, 
                                       control = data_control)
    data_patient <- results_checkFormat[[1]]
    data_control <- results_checkFormat[[2]]
  }else{
    results_load_input <- load_input(know_target_region)
    path_to_patient <- results_load_input[[1]]
    path_to_control <- results_load_input[[2]]
    data_patient <- read_tsv(path_to_patient, col_types = cols(chr = col_character(), zygosity = col_character()))
    data_control <- read_tsv(path_to_control, col_types = cols(chr = col_character(), zygosity = col_character()))
    results_checkFormat <- checkFormat(patient = data_patient, 
                                       control = data_control)
    data_patient <- results_checkFormat[[1]]
    data_control <- results_checkFormat[[2]]
  }


  if(know_target_region != 0){
    gene_list <- read_tsv(know_target_region, col_names = FALSE)
    colnames(gene_list) <- "region"
    data_patient <- data_patient[which(data_patient$gene_symbol %in% gene_list$region),]
    data_control <- data_control[which(data_control$gene_symbol %in% gene_list$region),]
  }
  
  #Prepare main folder to store results
  dir.create(paste(path_dir_results, "general", sep = ""))
  path_general_results <- paste(path_dir_results, "general/", sep = "")
  
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
                                                coverage_matrix_path = coverage_matrix_path,
                                                path_to_table = path_data_annotate)
    data_patient_annotated <- result_data_preparation[[1]]
    data_control_annotated <- result_data_preparation[[2]]
    clean_genotype_matrix <- result_data_preparation[[3]]
    missing_annotation <- result_data_preparation[[4]]
    variant_target <- result_data_preparation[[5]]
    
    if(coverage_matrix_path != 0){
      clean_coverage_matrix <- result_data_preparation[[6]]
      rejected_variant <- result_data_preparation[[7]]
      rejected_data_patient <- result_data_preparation[[8]]
      rejected_data_control <- result_data_preparation[[9]]

      write_tsv(clean_coverage_matrix, file = paste(path_general_results, "clean_coverage_matrix.tsv", sep = ""))
      write_tsv(rejected_variant, file = paste(path_general_results, "rejected_variants.tsv", sep = ""))
      write_tsv(rejected_data_patient, file = paste(path_general_results, "rejected_data_patient.tsv", sep = ""))
      write_tsv(rejected_data_control, file = paste(path_general_results, "rejected_data_control.tsv", sep = ""))
    
    }

    write_tsv(data_patient_annotated, file = paste(path_general_results, "data_patient_annotated.tsv", sep = ""))
    write_tsv(data_control_annotated, file = paste(path_general_results, "data_control_annotated.tsv", sep = ""))
    write_tsv(missing_annotation, file = paste(path_general_results, "missing_annotation.tsv", sep = ""))
    write_tsv(variant_target, file = paste(path_general_results, "variant_list.tsv", sep = ""))
  }else{
    data_patient_annotated <- data_patient
    data_control_annotated <- data_control
    Nbr_patient <- length(unique(data_patient_annotated$sample))
    Nbr_control <- length(unique(data_control_annotated$sample))
    result_geno_matrix <- genotypeMatrix(data_patient = data_patient_annotated, 
                                         data_control = data_control_annotated, 
                                         Nbr_patient_init = Nbr_patient, 
                                         Nbr_control_init = Nbr_control)
    clean_genotype_matrix <- result_geno_matrix[[1]]
    patient_name <- unique(data_patient_annotated$patient)
    control_name <- unique(data_control_annotated$patient)
    clean_genotype_matrix <- data.frame(clean_genotype_matrix, row.names = c(patient_name,control_name))
    list_variant <- result_geno_matrix[[4]]
    colnames(clean_genotype_matrix) <- list_variant$id
    rm(list_variant)
  }
  stop_preparation <- Sys.time()
  tot_preparation <- stop_preparation - start_preparation
  
  start <- Sys.time()
  patient_list <- unique(data_patient_annotated$sample)
  control_list <- unique(data_control_annotated$sample)
  Nbr_patient <- length(patient_list)
  Nbr_control <- length(control_list)
  
  
  ########################### HYPOTHERES GENERATION BASED ON FILTERING #############################
  res_filtering <- filtering_option(KindOfRegion = KindOfRegion,
                                    filtering_criteria = filtering_criteria,
                                    filtering_level = filtering_level,
                                    only_over_represented = only_over_represented)
  #Prepare SUBfolder to store results
  vec_hypotheses <- res_filtering[[1]]
  Nbr_exp <- res_filtering[[2]]
  KindOfRegion <- res_filtering[[3]]
  folder_name_experience <- paste(vec_hypotheses$genetic_region, "/", sep = "")
  for (i in 1:Nbr_exp) {
    dir.create(paste(path_dir_results,vec_hypotheses$genetic_region[i], sep = ""))
  }
  output_folder_directory <- paste(path_dir_results,folder_name_experience, sep = "")
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
  global_info_experiment$name_of_experiment <- vec_hypotheses$name_of_experience
  name_time_data <- c(paste("tot_single_analysis", seq(1:Nbr_exp), sep = "_"),
                      "data_preparation",
                      "tot_gene_set_analysis",
                      "total_time")
  time_data <- data.frame(matrix(NA, ncol = length(name_time_data), nrow = 1))
  colnames(time_data)<- name_time_data
  
  
  #vector initialization
  list_general_info_singleAnalysis <- list()
  list_general_info_not_analyzed_singleAnalysis <- list()
  list_qvalue_singleAnalysis <- list()
  tot_gene_set_analysis <- NA
  
  ########################### PRE FILTERING STEP OF GENETIC REGION TO BE ANALYZED ###########################
  if (length(unique(data_patient_annotated$gene_symbol)) > 1 && DO_over_representation) {
    res <- gene_set_analysis(data_patient = data_patient_annotated,
                             data_control = data_control_annotated,
                             level_go = c(3,4,5),
                             go_of_interest = c("BP", "MF", "CC"),
                             max_GO_similarity_accepted = 0.5,
                             cutoff = 0.05,
                             path_to_table = path_data_annotate, 
                             initial_path_to_store_results = path_general_results)
    over_represented_pathway <- res[[1]]
    over_represented_GO <- res[[2]]
  }else{
    over_represented_pathway <- c()
    over_represented_GO <- c()
  }

  
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
    if(name_of_experience[z] == "over_represented_pathway"){
      tmp_region <- region_of_interest[which(region_of_interest %in% over_represented_pathway)]
      if(length(tmp_region) > 1){
        region_of_interest <- tmp_region
      }else{
        region_of_interest <- NA
      }
    }else if(name_of_experience[z] == "over_represented_GO"){
      tmp_region <- region_of_interest[which(region_of_interest %in% over_represented_GO)]
      if(length(tmp_region) > 1){
        region_of_interest <- tmp_region
      }else{
        region_of_interest <- NA
      }
    }
    else if(name_of_experience[z] == "sig_aggregate_gene"){
      if(exists("target_region")){
        region_of_interest <- target_region
      }else{
        region_of_interest <- NA
      }
    }
    if(length(region_of_interest) == 0){
      warning(paste("experiment number", z, "not performed because no data"))
    }else{
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
      
      #set correct path to save results
      path_to_save_plot <- output_folder_directory[z]
      
      #ZFA algo global analysis
      #summary_zfa <- zfa_analysis(data_patient = filter_annotate_data_patient,
      #                            data_control = filter_annotate_data_control,
      #                            genotype_matrix = clean_genotype_matrix,
      #                            bin = 256,
      #                            CommonRare_Cutoff = 0.01,
      #                            path_to_results = path_to_save_plot)
      #if(max(summary_zfa$Nbr_genes) > 200){
      #  write_xlsx(summary_zfa[,-8], paste(path_to_save_plot,"zfa_global_result.xlsx",sep=""))
      #}else{
      #  write_xlsx(summary_zfa, paste(path_to_save_plot,"zfa_global_result.xlsx",sep=""))
      #}
      
      
      ########################### SINGLE ANALYSIS #############################
      start_single_analysis <- Sys.time()

patient = filter_annotate_data_patient, 
                                 control = filter_annotate_data_control,
                                 patient_list = patient_list,
                                 control_list = control_list,
                                 KindOfRegion = KindOfRegion[z], 
                                 NameOfRegion = region_of_interest, 
                                 Binary = TRUE, 
                                 genotype_matrix = clean_genotype_matrix,
                                 min_var = min_var,
                                 max_var = max_var,
                                 var_covariate = var_covariate,
                                 var_covariate_pca = var_covariate_pca,
                                 rare_maf_threshold = rare_maf_threshold,
                                 version = version

  
  #usefull info
  Nbr_regions <- length(NameOfRegion)
  Nbr_patient <- length(patient_list)
  Nbr_control <- length(control_list)
  #if(filtered){
  #  General_info <-matrix(data=NA,nrow = Nbr_regions, ncol = 21)
  #  colnames(General_info) <- c(paste(KindOfRegion),"Nbr individuals with minor alleles", "total minor allele count",
  #                              "Proportion patients with mutation", "Proportion controls with mutation", 
  #                              "Total nbr variant", "Total nbr variant in patients", "Total nbr variant in controls", 
  #                              "Nbr unique variant","Nbr unique variant in controls",
  #                              "Total nbr heterozygous variant patient", "Total nbr homozygous variant patient",
  #                              "Total nbr heterozygous variant control", "Total nbr homozygous variant control",
  #                              "Proportion of reference homozygous", "Proportion of heterozygous", "Proportion of alternatif homozygous",
  #                              "cumul diff zygo before", "cumul diff zygo new",
  #                              "List of variants only found in patient", "list of variants used")
  #} 
  
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
  #"List_of_variants_only_found_in_patient", 
  #"List_of_variants_only_found_in_control",
  #"list_of_variants_used"
  
  General_info <- data.frame(matrix(data=NA,nrow = Nbr_regions, ncol = length(annotation_name)))
  colnames(General_info) <- annotation_name
  rejected_region <- General_info

  
  #initiate dataframe to collect results
  name_of_test <- stat_Framework(Binary = Binary,
                                 get_test = TRUE,
                                 version = version)
  summary_pvalue <- data.frame(
    test = name_of_test
  )
  summary_pvalue_excell <- data.frame()
  summary_pvalue_excell <- summary_pvalue
  Nbr_test <- length(name_of_test)
  
  #  ITERATION ACCROSS REGIONS OF INTEREST  #
  blabla <-   c(1:Nbr_regions)
  TrueRegion <- c()
  for (loop in blabla) {
    
    #region selection
    results_region_selection <- region_selection(data_patient = patient, 
                                                 data_control = control,
                                                 genotype_matrix = genotype_matrix,
                                                 region_type = KindOfRegion, 
                                                 region_name = NameOfRegion[loop])
    focus_data_patient <- results_region_selection[[1]] 
    focus_data_control <- results_region_selection[[2]]
    focus_genotype_matrix <- results_region_selection[[3]]
    
    #no variant left in region , discard region
    if(dim(focus_genotype_matrix)[2] != 0){
      ################## Annotate results ####################
      Nbr_total_mutation <- dim(focus_genotype_matrix)[1] * dim(focus_genotype_matrix)[2]
      General_info$name_genetic_region[loop] <- NameOfRegion[loop]
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
      #General_info$List_of_variants_only_found_in_patient[loop] <- paste(focus_data_patient$variant_id[which(!focus_data_patient$variant_id %in% focus_data_control$variant_id)], collapse = "/")
      #General_info$List_of_variants_only_found_in_control[loop] <- paste(focus_data_control$variant_id[which(!focus_data_control$variant_id %in% focus_data_patient$variant_id)], collapse = "/")
      #General_info$list_of_variants_used[loop] <- paste(colnames(focus_genotype_matrix), collapse = "/")
      
      ############ Run analysis for focus region ############
      if(dim(focus_genotype_matrix)[2] < max_var){
        #estimate covariate via PCA
        err <- try(tmp <- prcomp(focus_genotype_matrix, scale. = TRUE),
                   silent = TRUE)
        if(length(err) == 1){
          err <- try(tmp <- prcomp(focus_genotype_matrix),
                     silent = TRUE)
          if(length(err) == 1){
            var_covariate_pca <- FALSE
          }else{
            pca <- tmp
          }
        }else{
          pca <- tmp
        }
        
        if(var_covariate_pca){
          err <- try(tmp_covariate_pca <- pca$x[,1:2],
                     silent = TRUE)
          if(length(err) == 1){
            covariate_pca <- 1
          }else{
            covariate_pca <- tmp_covariate_pca
          }
        }else{
          covariate_pca <- 1
        }
        
        #phenotype vectore construction
        phenotype <- c(rep(1,Nbr_patient), rep(0,Nbr_control))
        
        #set MAF of variants
        tmp_table_maf <- data.frame(matrix(NA, ncol = 2, nrow = dim(focus_genotype_matrix)[2]))
        #to do big bug only focus on patient boss
        tmp_table_maf[,1] <- colnames(focus_genotype_matrix)
        maf_vec <- c(focus_data_patient$gnomad_wes_af, focus_data_control$gnomad_wes_af)
        for (tmp_id_maf in 1:dim(tmp_table_maf)[1]) {
          tmp <- unique(maf_vec[which(c(focus_data_patient$variant_id, focus_data_control$variant_id) == tmp_table_maf[tmp_id_maf,1])])
          if(length(tmp) > 1){
            print(paste("different MAF info for same variant on region", loop, "for variant", tmp_id_maf))
            tmp <- min(tmp, na.rm = TRUE)
          }
          tmp_table_maf[tmp_id_maf,2] <- tmp
        }
        tmp_table_maf <- tmp_table_maf[order(tmp_table_maf$X1),]
        
        #get position
        position <- c()
        for (pos in 1:dim(focus_genotype_matrix)[2]) {
          tmp_pos <- colnames(focus_genotype_matrix)[pos]
          position <- c(position, as.integer(strsplit(tmp_pos, split = "_")[[1]][2]))
        }
        
        #SKAT framework
        skat_results <- stat_Framework(genotype_matrix = focus_genotype_matrix, 
                                       phenotype = phenotype, 
                                       weight_variant = tmp_table_maf[,2],
                                       Binary = Binary, 
                                       covariate = covariate_pca, 
                                       rare_maf_threshold = rare_maf_threshold,
                                       filter.pval = 0.05,
                                       position = position,
                                       get_test = FALSE, 
                                       version = version)
        pvalue <- skat_results[[1]] 
        #minor_allele_info <-skat_results[[2]]
        
        #collecting results
        General_info$Nbr_test_applied[loop] <- Nbr_test - length(which(is.na(pvalue)))
        tmp_pvalue <- data.frame(pvalue)
        colnames(tmp_pvalue) <- NameOfRegion[loop]
        summary_pvalue_excell <-  cbind(summary_pvalue_excell,tmp_pvalue)
        TrueRegion <- c(TrueRegion, General_info$name_genetic_region[loop])

        msg <- sprintf("%d/%d", loop, Nbr_regions)
        print(paste(msg, "region has been analyzed", sep = "  "))
      }#end of IF enough variants to be aggregated
    }#end of IF variant left in region
  }#end of for going through all regions of interest
  
  
  idx_kept <- which(!General_info$name_genetic_region %in% TrueRegion)
  if(length(idx_kept) > 0){
    General_info_not_analyzed <- General_info[idx_kept,]
    General_info <- General_info[-idx_kept,] 
  }else{
    General_info_not_analyzed <- c()
  }
  
  
  ###Multiple testing correction
  end_file <- dim(summary_pvalue_excell)[2]
  if(length(General_info$name_genetic_region) == 0){
    TrueRegion <- c()
  }

  if(end_file > 2){
    dum_pvalue <- summary_pvalue_excell[,2:end_file]
    qvalue <- p.adjust(unlist(dum_pvalue), method = p.adjust.methods[5]) # Benhjamini Hochberg
    qvalue <- data.frame(matrix(data = qvalue,ncol = dim(dum_pvalue)[2],nrow = Nbr_test))
    summary_qvalue_excell <- bind_cols(summary_pvalue,qvalue)
    colnames(summary_qvalue_excell) <- c("test", TrueRegion)
  }else{
    summary_qvalue_excell <- summary_pvalue_excell
  }
  
  if(length(TrueRegion) > 0){
    #results of analysis
    res_of_analysis <- data.frame(matrix(NA, ncol = 5, nrow = length(TrueRegion)))
    colnames(res_of_analysis) <- c("Excalibur_qvalue", "best_qvalue", "best_Stat_test", "Sig_Stat_test", "Nbr_sig_test")
    for (tmp_best in 2:dim(summary_qvalue_excell)[2]) {
      res_of_analysis$Excalibur_qvalue[tmp_best-1] <- summary_qvalue_excell[which(summary_qvalue_excell$test == "Excalibur"),tmp_best]
      res_of_analysis$best_qvalue[tmp_best-1] <- min(summary_qvalue_excell[,tmp_best], na.rm = TRUE)
      tmp_idx_test <- which(summary_qvalue_excell[,tmp_best] == res_of_analysis$best_qvalue[tmp_best-1])
      if(res_of_analysis$best_qvalue[tmp_best-1] <= pvalue_threshold){
        order_test_idx <- base::order(summary_qvalue_excell[,tmp_best])
        tmp_idx_test <- order_test_idx[which(summary_qvalue_excell[order_test_idx, tmp_best] <= pvalue_threshold)]
        res_of_analysis$Nbr_sig_test[tmp_best-1] <- length(tmp_idx_test)
        tmp_name_test <- as.character(summary_qvalue_excell$test[tmp_idx_test])
        tmp_best_test <- summary_qvalue_excell$test[which(summary_qvalue_excell[,tmp_best] == res_of_analysis$best_qvalue[tmp_best-1])]
        if(length(tmp_best_test) == 1){
          res_of_analysis$best_Stat_test[tmp_best-1] <- tmp_best_test
        }else{
          res_of_analysis$best_Stat_test[tmp_best-1] <- paste(tmp_best_test, collapse = "/")
        }
        if(length(tmp_name_test) == 1){
          res_of_analysis$Sig_Stat_test[tmp_best-1] <- tmp_name_test
        }else{
          res_of_analysis$Sig_Stat_test[tmp_best-1] <- paste(tmp_name_test, collapse = "/")
        }
      }else{
        res_of_analysis$Nbr_sig_test[tmp_best-1] <- 0
        tmp_name_test <- as.character(summary_qvalue_excell$Sig_Stat_test_name[tmp_idx_test])
        if(length(tmp_name_test) == 1){
          res_of_analysis$Sig_Stat_test[tmp_best-1] <- tmp_name_test
        }else{
          res_of_analysis$Sig_Stat_test[tmp_best-1] <- paste(tmp_name_test, collapse = "/")
        }
      }
    }
    General_info <- cbind(res_of_analysis, General_info)
  }
  
  
  return(list(General_info,
              General_info_not_analyzed,
              TrueRegion,
              summary_pvalue_excell,
              summary_qvalue_excell))



      results <- list(General_info,
              General_info_not_analyzed,
              TrueRegion,
              summary_pvalue_excell,
              summary_qvalue_excell)) 

      general_info_singleAnalysis <- results[[1]]
      general_info_not_analyzed_singleAnalysis <- results[[2]]
      TrueRegions <- results[[3]]
      pvalue_singleAnalysis <- results[[4]]
      qvalue_singleAnalysis <- results[[5]]
      stop_single_analysis <- Sys.time()
      tot_single_analysis <- stop_single_analysis - start_single_analysis
      time_data[1,z] <- paste(round(as.double(tot_single_analysis), digits = 2), units(tot_single_analysis))
      
      #TO REMOVE
      saveRDS(general_info_singleAnalysis, paste(path_to_save_plot,"general_info_singleAnalysis.rds",sep=""))
      saveRDS(general_info_not_analyzed_singleAnalysis, paste(path_to_save_plot,"general_info_not_analyzed_singleAnalysis.rds",sep=""))
      saveRDS(TrueRegions, paste(path_to_save_plot,"TrueRegions.rds",sep=""))
      saveRDS(pvalue_singleAnalysis, paste(path_to_save_plot,"pvalue_singleAnalysis.rds",sep=""))
      saveRDS(qvalue_singleAnalysis, paste(path_to_save_plot,"qvalue_singleAnalysis.rds",sep=""))

      #summary info of experiment
      global_info_experiment$Nbr_region[z] <- length(region_of_interest)
      global_info_experiment$Nbr_region_analyzed[z] <-length(TrueRegions)
      if(length(general_info_not_analyzed_singleAnalysis) > 0){
        global_info_experiment$Nbr_region_not_analyzed[z] <- dim(general_info_not_analyzed_singleAnalysis)[1]
      }else{
        global_info_experiment$Nbr_region_not_analyzed[z] <- 0
      }
      
      print(paste("We are currently working on experiment number ",z,sep = "  "))
      
      if(length(TrueRegions) > 0){
        #Compute the ranking
        ranking <- as.numeric(as.factor(base::rank(general_info_singleAnalysis$Excalibur_qvalue)))
        ranking <- data.frame(ranking, stringsAsFactors = FALSE)
        general_info_singleAnalysis <- bind_cols(ranking,general_info_singleAnalysis)
        #order by ranking 
        if(dim(general_info_singleAnalysis)[1] > 0){
          general_info_singleAnalysis <- general_info_singleAnalysis[order(general_info_singleAnalysis$ranking),]
        }
        
        
        #Ploting the results
        #plot_results <- plotSingleAnalysis(data_to_plot=qvalue_singleAnalysis)
        
        #Plot the ranking
        #plot_rank_results <- plotRankSingleAnalysis(data_pvalue=qvalue_singleAnalysis, data_rank=rank_list)
        
        
        #only significant results 
        significant_general_info_singleAnalysis <- filter(general_info_singleAnalysis, Excalibur_qvalue <= pvalue_threshold)
        global_info_experiment$Nbr_significant_region[z] <- dim(significant_general_info_singleAnalysis)[1]
        global_info_experiment$Proportion_Nbr_significant_region_VS_Nbr_region_analyzed[z] <- global_info_experiment$Nbr_significant_region[z] / global_info_experiment$Nbr_region_analyzed[z]
        global_info_experiment$Proportion_Nbr_region_analyzed_VS_Nbr_region[z] <- dim(general_info_singleAnalysis)[1] / global_info_experiment$Nbr_region[z]
        
        #Gene set analysis on significant results IFF level of analysis is gene
        if(vec_hypotheses[z,3][[1]] == "gene" && length(significant_general_info_singleAnalysis$name_genetic_region) > 1){
          start_gene_set_analysis <- Sys.time()
          target_region <- significant_general_info_singleAnalysis$name_genetic_region
          #focus only on significant aggregated gene 
          target_data_patient <- filter(filter_annotate_data_patient, gene_symbol %in% target_region)
          target_data_control <- filter(filter_annotate_data_control, gene_symbol %in% target_region)
          res <- gene_set_analysis(data_patient = target_data_patient,
                                   data_control = target_data_control,
                                   level_go = c(3,4,5),
                                   go_of_interest = c("BP", "MF", "CC"),
                                   max_GO_similarity_accepted = 0.5,
                                   cutoff = 0.05,
                                   path_to_table = path_data_annotate, 
                                   initial_path_to_store_results = output_folder_directory[z])
          sig_over_represented_pathway <- res[[1]]
          sig_over_represented_GO <- res[[2]]
          stop_gene_set_analysis <- Sys.time()
          tot_gene_set_analysis <- stop_gene_set_analysis - start_gene_set_analysis
        }
        
        
        #ensure data format
        #name <- data.frame(colnames(qvalue_singleAnalysis[-1]))

        #pvalue_singleAnalysis <- bind_cols(name,as.data.frame(t(pvalue_singleAnalysis[-1])))
        #qvalue_singleAnalysis <- bind_cols(name,as.data.frame(t(qvalue_singleAnalysis[-1])))
        
        #Saving results
        err <- try(write_xlsx(pvalue_singleAnalysis,paste(path_to_save_plot,"pvalue.xlsx",sep="")),
                  silent = TRUE)
        if(err != paste(path_to_save_plot,"pvalue.xlsx",sep="")){
          write_tsv(pvalue_singleAnalysis, file = paste(path_to_save_plot,"pvalue.tsv",sep=""))
        }
        err <- try(write_xlsx(qvalue_singleAnalysis,paste(path_to_save_plot,"qvalue.xlsx",sep="")),
                  silent = TRUE)
        if(err != paste(path_to_save_plot,"qvalue.xlsx",sep="")){
          write_tsv(qvalue_singleAnalysis, file = paste(path_to_save_plot,"qvalue.tsv",sep=""))
        }
        err <- try(write_xlsx(general_info_singleAnalysis,paste(path_to_save_plot,"general_info.xlsx",sep = "")),
                  silent = TRUE)
        if(err != paste(path_to_save_plot,"general_info.xlsx",sep = "")){
          write_tsv(general_info_singleAnalysis, file = paste(path_to_save_plot,"general_info.tsv",sep = ""))
        }

        if(length(general_info_not_analyzed_singleAnalysis)>0){
          err <- try(write_xlsx(general_info_not_analyzed_singleAnalysis,paste(path_to_save_plot,"general_info_not_analyzed.xlsx",sep = "")),
                    silent = TRUE)
          if(err != paste(path_to_save_plot,"general_info_not_analyzed.xlsx",sep = "")){
            write_tsv(general_info_not_analyzed_singleAnalysis, file = paste(path_to_save_plot,"general_info_not_analyzed.tsv",sep = ""))
          }
        }
        if(dim(significant_general_info_singleAnalysis)[1] > 0){
          err <- try(write_xlsx(significant_general_info_singleAnalysis,paste(path_to_save_plot,"significant_general_info.xlsx",sep = "")),
                    silent = TRUE)
          if(err != paste(path_to_save_plot,"significant_general_info.xlsx",sep = "")){
            write_tsv(significant_general_info_singleAnalysis, file = paste(path_to_save_plot,"significant_general_info.tsv",sep = ""))
          }
        }
        
        
      }else{
        #If no region was analyzed 
        global_info_experiment$Nbr_significant_region[z] <- 0
        global_info_experiment$Proportion_Nbr_significant_region_VS_Nbr_region_analyzed[z] <- NA
        global_info_experiment$Proportion_Nbr_region_analyzed_VS_Nbr_region[z] <- NA
        
        if(dim(general_info_not_analyzed_singleAnalysis)[1]>0){
          #save results not analyzed
          err <- try(write_xlsx(general_info_not_analyzed_singleAnalysis,paste(path_to_save_plot,"general_info_not_analyzed.xlsx",sep = "")),
                    silent = TRUE)
          if(err != paste(path_to_save_plot,"general_info_not_analyzed.xlsx",sep = "")){
            write_tsv(general_info_not_analyzed_singleAnalysis, file = paste(path_to_save_plot,"general_info_not_analyzed.tsv",sep = ""))
          }
        }
      }
      
      #keep all data from experiment z in memory
      list_general_info_singleAnalysis[[z]] <- general_info_singleAnalysis
      list_general_info_not_analyzed_singleAnalysis[[z]] <- general_info_not_analyzed_singleAnalysis
      list_qvalue_singleAnalysis[[z]] <- qvalue_singleAnalysis
    }#end of else in the case there is data and region to be analyzed for experiment z
    
  }#end for loop going though all the different experiments
  
  #TO DO assemble/analyse results from differents experiments
  
  #save global results
  write_xlsx(global_info_experiment,paste(path_general_results,"global_info_experiment.xlsx",sep = ""))
  stop <- Sys.time()
  tot <- stop - start
  
  if(is.na(tot_gene_set_analysis)){
    time_data[1,(z+1):dim(time_data)[2]] <- c(paste(round(as.double(tot_preparation), digits = 2), units(tot_preparation)),
                                            "0",
                                            paste(round(as.double(tot), digits = 2), units(tot)))
  }else{
    time_data[1,(z+1):dim(time_data)[2]] <- c(paste(round(as.double(tot_preparation), digits = 2), units(tot_preparation)),
                                            paste(round(as.double(tot_gene_set_analysis), digits = 2), units(tot_gene_set_analysis)),
                                            paste(round(as.double(tot), digits = 2), units(tot)))
  }
  
  write_xlsx(time_data, paste(path_general_results, "time_analysis.xlsx", sep = ""))

  return(time_data)
}