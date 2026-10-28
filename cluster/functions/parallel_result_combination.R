############################################# Function run general analysis on several experiments #############################################
# The function allow perform all the desired analysis 
# Input : Nbr_exp: Number of experiment to be performed
#         path_to_save_time_plot: (optional) a string indicating where the plot must be saved 
#                            if no path indicated, then results won't be stored
#         class_to_plot: (optional) a string with 4 different possible values in order to indicate which class of test are of interest
#                           default value = "all", otherwise could be "burden", "variance-component" or "omnibus"
# Output : plot from ggplot2
parallel_result_combination <- function(path_to_results = "",
                                        version = "optimal",
                                        Binary = TRUE,
                                        pvalue_threshold = 0.05,
                                        level = "gene",
                                        path_to_table = "",
                                        path_store_results = ""){
  
  #Data importation
  General_info <- readRDS(paste(path_to_results, "general_info_singleAnalysis.rds", sep = ""))
  path_input <- paste(path_to_results, "results/", sep = "")
  path_time <- paste(path_to_results, "time/", sep = "")
  path_data <- paste(path_to_results, "data_to_analyze/", sep = "")
  name_file <- list.files(path_input)
  name_file_time <- list.files(path_time)

  #usefull info
  Nbr_regions <- length(name_file)
  start_single_analysis <- Sys.time()

  #initiate dataframe to collect results
  rejected_region <- General_info
  name_of_test <- stat_Framework(Binary = Binary,
                                  get_test = TRUE,
                                  version = version)
  summary_pvalue <- data.frame(
    test = name_of_test
  )
  summary_pvalue_excell <- data.frame()
  summary_pvalue_excell <- summary_pvalue
  summary_time <- data.frame()
  summary_time <- summary_pvalue
  Nbr_test <- length(name_of_test)
  
  #  ITERATION ACCROSS REGIONS TESTED  #
  iter <-   c(1:Nbr_regions)
  TrueRegion <- c()
  for (loop in iter) {
    res <- readRDS(paste(path_input, name_file[loop], sep = ""))
    pvalue <- res[[1]]
    region_of_interest <- res[[2]]
    res_time <- readRDS(paste(path_time, name_file_time[loop], sep = ""))
    #collecting results
    General_info$Nbr_test_applied[which(General_info$name_genetic_region == region_of_interest)] <- Nbr_test - length(which(is.na(pvalue)))
    tmp_pvalue <- data.frame(pvalue)
    colnames(tmp_pvalue) <- region_of_interest
    summary_pvalue_excell <-  cbind(summary_pvalue_excell,tmp_pvalue)
    tmp_time <- data.frame(res_time)
    colnames(tmp_time) <- region_of_interest
    summary_time <- cbind(summary_time, tmp_time)
    TrueRegion <- c(TrueRegion, region_of_interest)
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
    qvalue_final <- data.frame(matrix(NA, ncol = length(2:end_file), nrow = Nbr_test))
    for(i in 1:dim(summary_pvalue_excell)[1]){
      dum_pvalue <- summary_pvalue_excell[i,2:end_file]
       qvalue <- p.adjust(unlist(dum_pvalue), method = p.adjust.methods[5]) # Benhjamini Hochberg
       
       qvalue_final[i,] <- qvalue
    }
    #dum_pvalue <- summary_pvalue_excell[,2:end_file]
    #qvalue <- p.adjust(unlist(dum_pvalue), method = p.adjust.methods[5]) # Benhjamini Hochberg
    #qvalue <- data.frame(matrix(data = qvalue,ncol = dim(dum_pvalue)[2],nrow = Nbr_test))
    #summary_qvalue_excell <- bind_cols(summary_pvalue,qvalue)
    summary_qvalue_excell <- bind_cols(summary_pvalue,qvalue_final)
    colnames(summary_qvalue_excell) <- c("test", TrueRegion)
  }else{
    summary_qvalue_excell <- summary_pvalue_excell
  }

  if(length(TrueRegion) > 0){
    #results of analysis
    res_of_analysis <- data.frame(matrix(NA, ncol = 5, nrow = length(TrueRegion)))
    colnames(res_of_analysis) <- c("Excalibur_qvalue", "best_qvalue", "best_Stat_test", "Sig_Stat_test", "Nbr_sig_test")
    for (tmp_best in 2:dim(summary_qvalue_excell)[2]) {
      tmp_idx <- which(General_info$name_genetic_region == colnames(summary_qvalue_excell)[tmp_best])
      res_of_analysis$Excalibur_qvalue[tmp_idx] <- summary_qvalue_excell[which(summary_qvalue_excell$test == "Excalibur"),tmp_best]
      res_of_analysis$best_qvalue[tmp_idx] <- min(summary_qvalue_excell[,tmp_best], na.rm = TRUE)
      tmp_idx_test <- which(summary_qvalue_excell[,tmp_best] == res_of_analysis$best_qvalue[tmp_idx])
      if(res_of_analysis$best_qvalue[tmp_idx] <= pvalue_threshold){
        order_test_idx <- base::order(summary_qvalue_excell[,tmp_best])
        tmp_idx_test <- order_test_idx[which(summary_qvalue_excell[order_test_idx, tmp_best] <= pvalue_threshold)]
        res_of_analysis$Nbr_sig_test[tmp_idx] <- length(tmp_idx_test)
        tmp_name_test <- as.character(summary_qvalue_excell$test[tmp_idx_test])
        tmp_best_test <- summary_qvalue_excell$test[which(summary_qvalue_excell[,tmp_best] == res_of_analysis$best_qvalue[tmp_idx])]
        if(length(tmp_best_test) == 1){
          res_of_analysis$best_Stat_test[tmp_idx] <- tmp_best_test
        }else{
          res_of_analysis$best_Stat_test[tmp_idx] <- paste(tmp_best_test, collapse = "/")
        }
        if(length(tmp_name_test) == 1){
          res_of_analysis$Sig_Stat_test[tmp_idx] <- tmp_name_test
        }else{
          res_of_analysis$Sig_Stat_test[tmp_idx] <- paste(tmp_name_test, collapse = "/")
        }
      }else{
        res_of_analysis$Nbr_sig_test[tmp_idx] <- 0
        tmp_name_test <- as.character(summary_qvalue_excell$Sig_Stat_test_name[tmp_idx_test])
        if(length(tmp_name_test) == 1){
          res_of_analysis$Sig_Stat_test[tmp_idx] <- tmp_name_test
        }else{
          res_of_analysis$Sig_Stat_test[tmp_idx] <- paste(tmp_name_test, collapse = "/")
        }
      }
    }
    General_info <- cbind(res_of_analysis, General_info)
  }


  if(length(TrueRegion) > 0){
    #Compute the ranking
    ranking <- as.numeric(as.factor(base::rank(General_info$Excalibur_qvalue)))
    ranking <- data.frame(ranking, stringsAsFactors = FALSE)
    General_info <- bind_cols(ranking,General_info)
    #order by ranking 
    if(dim(General_info)[1] > 0){
      General_info <- General_info[order(General_info$ranking),]
    }
        
    #only significant results 
    significant_General_info <- filter(General_info, Excalibur_qvalue <= pvalue_threshold)
        
  
    #Saving results
    #err <- try(write_xlsx(summary_pvalue_excell,paste(path_store_results,"pvalue.xlsx",sep="")),
    #              silent = TRUE)
    #if(err != paste(path_store_results,"pvalue.xlsx",sep="")){
    #  write_tsv(summary_pvalue_excell, file = paste(path_store_results,"pvalue.tsv",sep=""))
    #}
    #err <- try(write_xlsx(summary_time,paste(path_store_results,"time.xlsx",sep="")),
    #              silent = TRUE)
    #if(err != paste(path_store_results,"time.xlsx",sep="")){
    #  write_tsv(summary_time, file = paste(path_store_results,"time.tsv",sep=""))
    #}
    #err <- try(write_xlsx(summary_qvalue_excell,paste(path_store_results,"qvalue.xlsx",sep="")),
    #              silent = TRUE)
    #if(err != paste(path_store_results,"qvalue.xlsx",sep="")){
    #  write_tsv(summary_qvalue_excell, file = paste(path_store_results,"qvalue.tsv",sep=""))
    #}
    err <- try(write_xlsx(General_info,paste(path_store_results,"general_info.xlsx",sep = "")),
                  silent = TRUE)
    if(err != paste(path_store_results,"general_info.xlsx",sep = "")){
      write_tsv(General_info, file = paste(path_store_results,"general_info.tsv",sep = ""))
    }

    if(length(General_info_not_analyzed)>0){
      err <- try(write_xlsx(General_info_not_analyzed,paste(path_store_results,"general_info_not_analyzed.xlsx",sep = "")),
                    silent = TRUE)
      if(err != paste(path_store_results,"general_info_not_analyzed.xlsx",sep = "")){
        write_tsv(General_info_not_analyzed, file = paste(path_store_results,"general_info_not_analyzed.tsv",sep = ""))
      }

    }
    
    if(dim(significant_General_info)[1] > 0){
      err <- try(write_xlsx(significant_General_info,paste(path_store_results,"significant_general_info.xlsx",sep = "")),
                    silent = TRUE)
      if(err != paste(path_store_results,"significant_general_info.xlsx",sep = "")){
        write_tsv(significant_General_info, file = paste(path_store_results,"significant_general_info.tsv",sep = ""))
      }
      file_data <- list.files(path_data)
      sig_name <- significant_General_info$name_genetic_region

      #Over representation analysis of significant gene
      if(level == "gene" && length(sig_name) > 1){
        path_over_representation <- paste(path_store_results, "sig_gene_over_representation", sep = "/")
        dir.create(path_over_representation)
        res_ov <- over_representation_analysis(geneList = sig_name,
                                              go_of_interest = c("BP", "MF", "CC"),
                                              max_GO_similarity_accepted = 0.5,
                                              cutoff = 0.05,
                                              max_gene_per_item = 500,
                                              max_item_plot = 20,
                                              coocurrence = TRUE,
                                              max_coocurrence = 100,
                                              MSigDB = TRUE,
                                              path_to_table = path_Excalibur_Data,
                                              path_to_store_results = path_over_representation)
      }

      #Save focus cases significant results
      data_focus_patient <- significant_General_info[which(significant_General_info$Total_nbr_variant_in_patients > significant_General_info$Total_nbr_variant_in_controls),]
      if(dim(data_focus_patient)[1] > 0){
        err <- try(write_xlsx(data_focus_patient,paste(path_store_results,"Cases_significant_results.xlsx",sep = "")),
                    silent = TRUE)
        if(err != paste(path_store_results,"Cases_significant_results.xlsx",sep = "")){
          write_tsv(data_focus_patient, file = paste(path_store_results,"Cases_significant_results.tsv",sep = ""))
        }

        sig_region_patient <- data_focus_patient$name_genetic_region
        #Over representation analysis of significant gene in cases
        if(level == "gene" && dim(data_focus_patient)[1] > 1){
          path_over_representation <- paste(path_store_results, "Cases_sig_gene_over_representation", sep = "/")
          dir.create(path_over_representation)
          res_ov <- over_representation_analysis(geneList = sig_region_patient,
                                                go_of_interest = c("BP", "MF", "CC"),
                                                max_GO_similarity_accepted = 0.5,
                                                cutoff = 0.05,
                                                max_gene_per_item = 500,
                                                max_item_plot = 20,
                                                coocurrence = TRUE,
                                                max_coocurrence = 100,
                                                MSigDB = TRUE,
                                                path_to_table = path_Excalibur_Data,
                                                path_to_store_results = path_over_representation)
        }

        #Save raw data of significant regions focus on patient
        sig_reg_dir <- paste(path_store_results, "data_significant_region/", sep ="")
        dir.create(sig_reg_dir)
        for(k in 1:length(file_data)){
          tmp <- readRDS(paste(path_data, file_data[k], sep = ""))
          tmp_name <- tmp[[4]]
          tmp_idx <- which(sig_region_patient == tmp_name)
          if(length(tmp_idx) > 0){
            if(length(grep("/", tmp_name)) > 0){
              tmp_name <- strsplit(tmp_name, split = "/")[[1]][1]
            }
            tmp_save_data <- rbind(tmp[[1]], tmp[[2]])
            tmp_name <- substr(tmp_name, 1, 20)
            tmp_name <- paste(tmp_name,"_data.xlsx",sep = "")
            tmp_name <- str_remove(tmp_name, pattern = "'")
            tmp_name <- str_remove(tmp_name, pattern = "'")
            err <- try(write_xlsx(tmp_save_data,paste(sig_reg_dir, tmp_name,sep = "")),
                      silent = TRUE)
            if(err != paste(sig_reg_dir, tmp_name,sep = "")){
            write_tsv(tmp_save_data, file = paste(sig_reg_dir, tmp_name,sep = ""))
            }#end of saving file
          }#end of if we found the region
        }#end of while going through all rds files
      }

      
    }#end of significant region saving
  }else{
    if(dim(General_info_not_analyzed)[1]>0){
      #save results not analyzed
      err <- try(write_xlsx(General_info_not_analyzed,paste(path_store_results,"general_info_not_analyzed.xlsx",sep = "")),
                    silent = TRUE)
      if(err != paste(path_store_results,"general_info_not_analyzed.xlsx",sep = "")){
        write_tsv(General_info_not_analyzed, file = paste(path_store_results,"general_info_not_analyzed.tsv",sep = ""))
      }
    }
  }

}#end of function