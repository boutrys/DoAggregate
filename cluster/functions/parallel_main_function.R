############################################# Function run general analysis on several experiments #############################################
# The function allow perform all the desired analysis 
# Input : Nbr_exp: Number of experiment to be performed
#         path_to_save_time_plot: (optional) a string indicating where the plot must be saved 
#                            if no path indicated, then results won't be stored
#         class_to_plot: (optional) a string with 4 different possible values in order to indicate which class of test are of interest
#                           default value = "all", otherwise could be "burden", "variance-component" or "omnibus"
# Output : plot from ggplot2
parallel_main_function <- function(path_to_patient = "",
                          path_to_control = "",
                          know_target_region = 0,
                          importation_method = "default",
                          mapqv_threshold = 20,
                          min_coverage = 20,
                          remove_individual = TRUE,
                          remove_position = TRUE,
                          path_dir_results = getwd(),
                          coverage_matrix_path = 0,
                          do_data_preparation = TRUE,
                          path_data_annotate = "",
                          tmp_folder_data_preparation = getwd(),
                          path_gvcf = "",
                          mandatory_column_mapping =  list(sample = "Indiv", 
                                                                    chr = "CHROM", 
                                                                    pos = "POS", 
                                                                    reference = "REF", 
                                                                    alternative = "ALT", 
                                                                    zygosity = "",
                                                                    gene_symbol = "SYMBOL", 
                                                                    gnomad_wes_af = "gnomad40_AF",
                                                                    gt_GT = "gt_GT", 
                                                                    gt_DP = "gt_DP"),
                            MAF = 0.01){
  
  #Data importation
  print("start reading files")
  if(importation_method == "default"){
    data_patient <- read_tsv(path_to_patient, col_types = cols(chr = col_character(), reference = col_character(), alternative = col_character(), zygosity = col_character()))
    data_control <- read_tsv(path_to_control, col_types = cols(chr = col_character(), reference = col_character(), alternative = col_character(), zygosity = col_character()))
    results_checkFormat <- checkFormat(patient = data_patient, 
                                       control = data_control)
    data_patient <- results_checkFormat[[1]]
    data_control <- results_checkFormat[[2]]
  }else if(importation_method == "gvcf"){
    res_gvcf <- gvcf_pretreatment(path_data = path_gvcf,
                                  path_output = path_dir_results,
                                  cases = path_to_patient,
                                  controls = path_to_control,
                                  mandatory_column_mapping =  mandatory_column_mapping,
                                  MAF = MAF,
                                  het_var = c("0/1", "0|1"), 
                                  homo_var = c("1/1", "1|1"),
                                  homo_ref = c("0/0", "0|0"),
                                  cleaning = do_data_preparation,
                                  remove_individual = remove_individual,
                                  remove_position = remove_position,
                                  min_coverage = min_coverage)
    data_patient <- res_gvcf[[1]]
    data_control <- res_gvcf[[2]]
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
  print("end reading files")

  if(know_target_region != 0){
    gene_list <- read_tsv(know_target_region, col_names = FALSE)
    colnames(gene_list) <- "region"
    data_patient <- data_patient[which(data_patient$gene_symbol %in% gene_list$region),]
    data_control <- data_control[which(data_control$gene_symbol %in% gene_list$region),]
  }


  ########################### USE BLACKLIST TO CLEAN INPUT DATA ###########################
  data_patient$id_variant <- paste(data_patient$chr, data_patient$pos, data_patient$reference, data_patient$alternative, sep = "_")
  data_control$id_variant <- paste(data_control$chr, data_control$pos, data_control$reference, data_control$alternative, sep = "_")
  # Optional blacklist of variants to remove (e.g. recurrent artefacts in an in-house cohort).
  # Set the environment variable EXCALIBUR_BLACKLIST to a TSV file with CHROM, POS, REF and ALT columns to use it.
  path_blacklist <- Sys.getenv("EXCALIBUR_BLACKLIST", "")
  if(nzchar(path_blacklist) && file.exists(path_blacklist)){
    blackList <- read_tsv(path_blacklist)
    blackList$CHROM <- paste("chr", blackList$CHROM, sep = "")
    blackList$variant_id <- paste(blackList$CHROM, blackList$POS, blackList$REF, blackList$ALT, sep = "_")
    to_rm_patient <- which(data_patient$id_variant %in% blackList$variant_id)
    to_rm_control <- which(data_control$id_variant %in% blackList$variant_id)
    if(length(to_rm_patient) > 0){
      blacklisted_data_patient <- data_patient[to_rm_patient,]
      data_patient <- data_patient[-to_rm_patient,]
      write_tsv(blacklisted_data_patient, file = paste(path_dir_results, "blacklisted_data_cases.tsv", sep = ""))
    }
    if(length(to_rm_control) > 0){
      blacklisted_data_control <- data_control[to_rm_control,]
      data_control <- data_control[-to_rm_control,]
      write_tsv(blacklisted_data_control, file = paste(path_dir_results, "blacklisted_data_controls.tsv", sep = ""))
    }
  }


  start_preparation <- Sys.time()
  ########################### DATA PREPARATION #############################
  if(do_data_preparation && importation_method == "default"){
    result_data_preparation <- data_preparation(data_patient = data_patient,
                                                data_control = data_control,
                                                mapqv_threshold = mapqv_threshold,
                                                min_coverage = min_coverage,
                                                remove_individual = remove_individual,
                                                remove_position = remove_position,
                                                path_to_store_plots = path_dir_results,
                                                coverage_matrix_path = coverage_matrix_path,
                                                path_to_table = path_data_annotate)
    data_patient_annotated <- result_data_preparation[[1]]
    data_control_annotated <- result_data_preparation[[2]]
    clean_genotype_matrix <- result_data_preparation[[3]]
    missing_annotation <- result_data_preparation[[4]]
    variant_target <- result_data_preparation[[5]]
    
    if(coverage_matrix_path != 0){
      print("start trying to write data_patient_annotated")
      clean_coverage_matrix <- result_data_preparation[[6]]
      rejected_variant <- result_data_preparation[[7]]
      rejected_data_patient <- result_data_preparation[[8]]
      rejected_data_control <- result_data_preparation[[9]]
      #Potential bug
      #write_tsv(clean_coverage_matrix, file = paste(path_dir_results, "clean_coverage_matrix.tsv", sep = ""))
      #write_tsv(rejected_variant, file = paste(path_dir_results, "rejected_variants.tsv", sep = ""))
      #write_tsv(rejected_data_patient, file = paste(path_dir_results, "rejected_data_patient.tsv", sep = ""))
      #write_tsv(rejected_data_control, file = paste(path_dir_results, "rejected_data_control.tsv", sep = ""))
    
    }
    #Potential bug
    
    #err <- try(write_tsv(data_patient_annotated, file = paste(path_dir_results, "data_patient_annotated.tsv", sep = "")),
    #            silent = TRUE)
 
    #if(length(grep("error", err)) > 1){
    #  print("MY_ERROR : cant write data_patient_annotated because of out of memory")
    #}
    print("end writing data_control_annotated")         
    #write_tsv(data_control_annotated, file = paste(path_dir_results, "data_control_annotated.tsv", sep = ""))
    #write_tsv(missing_annotation, file = paste(path_dir_results, "missing_annotation.tsv", sep = ""))
    #write_tsv(variant_target, file = paste(path_dir_results, "variant_list.tsv", sep = ""))
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
  }
  stop_preparation <- Sys.time()
  tot_preparation <- stop_preparation - start_preparation

  #save two time to avoid parallel jobs to work on the same data
  saveRDS(data_patient_annotated, paste(tmp_folder_data_preparation, "data_patient_annotated.rds", sep = ""))
  saveRDS(data_control_annotated, paste(tmp_folder_data_preparation, "data_control_annotated.rds", sep = ""))
  saveRDS(data_patient_annotated, paste(tmp_folder_data_preparation, "GENE_SET_data_patient_annotated.rds", sep = ""))
  saveRDS(data_control_annotated, paste(tmp_folder_data_preparation, "GENE_SET_data_control_annotated.rds", sep = ""))
  saveRDS(clean_genotype_matrix, paste(tmp_folder_data_preparation, "clean_genotype_matrix.rds", sep = ""))
  saveRDS(tot_preparation, paste(tmp_folder_data_preparation, "time_data_preparation.rds", sep = ""))
}#end of function