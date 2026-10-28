#' Run Single Case/Control Analysis on a Specified Region
#'
#' @name single_region_analysis
#' @description This function conducts a statistical analysis on a specific genetic region using case and control data.
#' It processes the data based on the provided genetic matrix and performs statistical tests to determine significant variants.
#'
#' @param focus_data_patient Data frame of patient data including at least patient, chromosome, position, reference, alternative, zygosity, and gene.
#' @param focus_data_control Data frame of control data including similar columns as patient data.
#' @param Nbr_patient Number of patients.
#' @param Nbr_control Number of controls.
#' @param Nbr_test Number of tests to be performed.
#' @param Binary Logical indicating if the analysis is binary.
#' @param focus_genotype_matrix Genotype matrix focused on the specific region.
#' @param max_var Maximum number of variants to be considered in the PCA.
#' @param var_covariate Variance to be used as covariate.
#' @param var_covariate_pca Logical indicating whether PCA should be used for covariate adjustment.
#' @param rare_maf_threshold Threshold for rare minor allele frequency (MAF).
#' @param pvalue_threshold Threshold for significance in p-values.
#' @param version Specifies the version of the analysis to be used, default is 'optimal'.
#' @param region_id Identifier for the region being analyzed.
#' @param region_of_interest Name or identifier of the region of interest.
#' @param time_path Path to save timing data.
#' @return A list containing p-values for the tests, the region of interest, and computational time.
#'
#' @importFrom stats prcomp
#' @export

single_region_analysis <- function(focus_data_patient = c(),
                                   focus_data_control = c(),
                                   Nbr_patient = 0,
                                   Nbr_control = 0,
                                   Nbr_test = 0,
                                   Binary = TRUE,
                                   focus_genotype_matrix = c(),
                                   max_var = max_var,
                                   var_covariate = 0,
                                   var_covariate_pca = FALSE,
                                   rare_maf_threshold = 0.01,
                                   pvalue_threshold = 0.05,
                                   version = "optimal",
                                   region_id = 0,
                                   region_of_interest = "",
                                   time_path = ""){
  ############ Run analysis for focus region ############
  if((dim(focus_genotype_matrix)[2] < max_var) && (dim(focus_genotype_matrix)[2] > 1)){
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

    tmp_table_maf[,1] <- colnames(focus_genotype_matrix)
    maf_vec <- c(focus_data_patient$gnomad_wes_af, focus_data_control$gnomad_wes_af)
    for (tmp_id_maf in 1:dim(tmp_table_maf)[1]) {
      tmp <- unique(maf_vec[which(c(focus_data_patient$variant_id, focus_data_control$variant_id) == tmp_table_maf[tmp_id_maf,1])])
      if(length(tmp) > 1){
        print(paste("different MAF info for same variant on region for variant", tmp_id_maf))
        tmp <- min(tmp, na.rm = TRUE)
      }
      tmp_table_maf[tmp_id_maf,2] <- tmp
    }
    tmp_table_maf <- tmp_table_maf[order(tmp_table_maf$X1),]
    tmp_table_maf[,2] <- as.double(tmp_table_maf[,2])

    #get position
    position <- c()
    for (pos in 1:dim(focus_genotype_matrix)[2]) {
      tmp_pos <- colnames(focus_genotype_matrix)[pos]
      position <- c(position, as.integer(strsplit(tmp_pos, split = "_")[[1]][2]))
    }

    #STAT framework
    excalibur_results <- stat_Framework(genotype_matrix = focus_genotype_matrix,
                                        phenotype = phenotype,
                                        weight_variant = tmp_table_maf[,2],
                                        Binary = Binary,
                                        covariate = covariate_pca,
                                        rare_maf_threshold = rare_maf_threshold,
                                        filter.pval = 0.05,
                                        position = position,
                                        get_test = FALSE,
                                        version = version)
    pvalue <- excalibur_results[[1]]
    computational_time <- excalibur_results[[2]]
    output <- list(pvalue, region_of_interest, computational_time)
    return(output)
    #saveRDS(computational_time, paste(time_path, region_id, "_computational_time.rds", sep = ""))
  }else{
    #stop(paste("Region id ", region_id, " has more variant than the max, see variable max_var (currently fixed at : ", max_var, sep =""))
  }


}#end of function
