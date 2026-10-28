#' Perform Comprehensive Analysis Across Experiments
#'
#' @name result_combination
#' @description This function performs a variety of statistical tests across multiple genetic regions or experiments,
#' aggregates the results, and optionally saves significant findings and summary plots.
#' It applies multiple testing correction and returns various statistics such as p-values and q-values.
#'
#' @param General_info A data frame containing general information about each genetic region or experiment.
#' @param TrueRegion A vector of region names that are specifically analyzed.
#' @param summary_pvalue A data frame summarizing p-values for each test and region.
#' @param summary_pvalue_excell A data frame similar to summary_pvalue but formatted for Excel output.
#' @param Nbr_test The number of tests performed.
#' @param version Character string denoting the version of analysis, default is "optimal".
#' @param Binary Logical indicating whether the analysis is binary.
#' @param pvalue_threshold The threshold for determining significance.
#' @param level The level at which analysis is performed (e.g., 'gene', 'pathway').
#' @param filter_annotate_data_patient Filtered patient data used for annotations.
#' @param filter_annotate_data_control Filtered control data used for annotations.
#' @param genotype_matrix A matrix or data frame of genotype information.
#' @param path_store_results The path where results should be stored. If not provided, results are not saved.
#' @return Depending on the function's processing, it may return the number of significant names detected or simply the performance metric (time taken).
#'
#' @importFrom stats p.adjust
#' @export
#' @importFrom utils globalVariables
utils::globalVariables(c("Excalibur_qvalue"))


result_combination <- function(General_info = c(),
                               TrueRegion = c(),
                               summary_pvalue = c(),
                               summary_pvalue_excell = c(),
                               Nbr_test = 0,
                               version = "optimal",
                               Binary = TRUE,
                               pvalue_threshold = 0.05,
                               level = "gene",
                               filter_annotate_data_patient = c(),
                               filter_annotate_data_control = c(),
                               genotype_matrix = c(),
                               path_store_results = ""){
  start_single_analysis <- Sys.time()


  #todo check bellow because i removed code here no TRUEREGION variable


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
    #err <- try(write_xlsx(summary_pvalue_excell,paste(path_store_results,"pvalue.xlsx",sep="/")),
    #              silent = TRUE)
    #if(err != paste(path_store_results,"pvalue.xlsx",sep="/")){
    #  write_tsv(summary_pvalue_excell, file = paste(path_store_results,"pvalue.tsv",sep="/"))
    #}
    #err <- try(write_xlsx(summary_time,paste(path_store_results,"time.xlsx",sep="/")),
    #              silent = TRUE)
    #if(err != paste(path_store_results,"time.xlsx",sep="/")){
    #  write_tsv(summary_time, file = paste(path_store_results,"time.tsv",sep="/"))
    #}
    #err <- try(write_xlsx(summary_qvalue_excell,paste(path_store_results,"qvalue.xlsx",sep="/")),
    #              silent = TRUE)
    #if(err != paste(path_store_results,"qvalue.xlsx",sep="/")){
    #  write_tsv(summary_qvalue_excell, file = paste(path_store_results,"qvalue.tsv",sep="/"))
    #}
    err <- try(write_xlsx(General_info,paste(path_store_results,"general_info.xlsx",sep = "/")),
               silent = TRUE)
    if(err != paste(path_store_results,"general_info.xlsx",sep = "/")){
      write_tsv(General_info, file = paste(path_store_results,"general_info.tsv",sep = "/"))
    }

    if(length(General_info_not_analyzed)>0){
      err <- try(write_xlsx(General_info_not_analyzed,paste(path_store_results,"general_info_not_analyzed.xlsx",sep = "/")),
                 silent = TRUE)
      if(err != paste(path_store_results,"general_info_not_analyzed.xlsx",sep = "/")){
        write_tsv(General_info_not_analyzed, file = paste(path_store_results,"general_info_not_analyzed.tsv",sep = "/"))
      }

    }

    if(dim(significant_General_info)[1] > 0){
      err <- try(write_xlsx(significant_General_info,paste(path_store_results,"significant_general_info.xlsx",sep = "/")),
                 silent = TRUE)
      if(err != paste(path_store_results,"significant_general_info.xlsx",sep = "/")){
        write_tsv(significant_General_info, file = paste(path_store_results,"significant_general_info.tsv",sep = "/"))
      }
      sig_name <- significant_General_info$name_genetic_region

      #Over representation analysis of significant gene
      if(level == "gene" && length(sig_name) > 1){
        path_over_representation <- paste(path_store_results, "sig_gene_over_representation", sep = "/")
        dir.create(path_over_representation)
        res_ov <- over_representation_analysis(geneList = sig_name,
                                               go_of_interest = c("BP", "MF", "CC"),
                                               max_GO_similarity_accepted = 0.5,
                                               cutoff = 0.05,
                                               max_gene_per_item = 20,
                                               max_item_plot = 20,
                                               coocurrence = TRUE,
                                               MSigDB = FALSE,
                                               max_coocurrence = 100,
                                               path_to_store_results = path_over_representation)
      }

      #Save raw data of significant regions
      sig_reg_dir <- paste(path_store_results, "data_significant_region", sep = "/")
      dir.create(sig_reg_dir)
      for(k in 1:length(sig_name)){
        tmp_name <- sig_name[k]
        tmp_idx <- which(sig_name == tmp_name)
        if(length(tmp_idx) > 0){
          if(level %in% c("GO", "pathway", "over_represented_pathway", "over_represented_GO")){
            tmp_geneID <- strsplit(significant_General_info$gene_involved[k], split = "/")[[1]][1]
            results_region_selection <- region_selection(data_patient = filter_annotate_data_patient,
                                                         data_control = filter_annotate_data_control,
                                                         genotype_matrix = genotype_matrix,
                                                         region_type = level,
                                                         region_name = tmp_name,
                                                         geneID = tmp_geneID)
          }else{
            results_region_selection <- region_selection(data_patient = filter_annotate_data_patient,
                                                         data_control = filter_annotate_data_control,
                                                         genotype_matrix = genotype_matrix,
                                                         region_type = level,
                                                         region_name = tmp_name)
          }

          focus_data_patient <- results_region_selection[[1]]
          focus_data_control <- results_region_selection[[2]]
          tmp_save_data <- rbind(focus_data_patient, focus_data_control)
          if(length(grep("/", tmp_name)) > 0){
            tmp_name <- strsplit(tmp_name, split = "/")[[1]][1]
          }
          tmp_name <- substr(tmp_name, 1, 20)
          tmp_name <- gsub(pattern = " ", replacement = "_", tmp_name)
          err <- try(write_xlsx(tmp_save_data,paste(sig_reg_dir, "/", tmp_name, "_data.xlsx", sep = "")),
                     silent = TRUE)
          if(err != paste(sig_reg_dir, "/", tmp_name, "_data.xlsx", sep = "")){
            write_tsv(tmp_save_data, file = paste(sig_reg_dir, "/", tmp_name, "_data.xlsx", sep = ""))
          }#end of saving file
        }#end of if we found the region
      }#end of while going through all rds files
    }#end of significant region saving
  }else{
    if(dim(General_info_not_analyzed)[1]>0){
      #save results not analyzed
      err <- try(write_xlsx(General_info_not_analyzed,paste(path_store_results,"general_info_not_analyzed.xlsx",sep = "/")),
                 silent = TRUE)
      if(err != paste(path_store_results,"general_info_not_analyzed.xlsx",sep = "/")){
        write_tsv(General_info_not_analyzed, file = paste(path_store_results,"general_info_not_analyzed.tsv",sep = "/"))
      }
    }
  }


  if(exists("sig_name")){
    return(length(sig_name))
  }else{
    return(0)
  }


}#end of function
