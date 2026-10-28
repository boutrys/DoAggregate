############################################# Function to run a single case/control analysis on a given kind of region #############################################
# The function allow to run an analysis through the statistical framework using a single case/control on a given kind of region
# Input : patient: data frame with at least column such as
#                                   patient, chromosome, position, reference, altenative, zygosity, gene
#         control: data frame with at least column such as
#                                   patient, chromosome, position, reference, altenative, zygosity, gene
#         KindOfRegion: [chr] must be equal to "chr", "pos" or "gene_symbol" (for the moment)
#         NameOfRegion: [list] name of the targeted region
# Output : data frame where column are 12 ranking 
#          data frame where each column is a region and each raw is a pvalue from a particular test
parallel_region_analysis <- function(focus_data_patient = c(),
                            focus_data_control = c(),
                            Nbr_patient = 0,
                            Nbr_control = 0, 
                            Binary = TRUE, 
                            focus_genotype_matrix = c(),
                            max_var = max_var,
                            var_covariate = 0,
                            var_covariate_pca = FALSE,
                            rare_maf_threshold = 0.01,
                            pvalue_threshold = 0.05,
                            version = "optimal",
                            path_output = getwd(),
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
    output <- list(pvalue, region_of_interest)
    saveRDS(output, paste(path_output, region_id, "_output.rds", sep = ""))
    saveRDS(computational_time, paste(time_path, region_id, "_computational_time.rds", sep = ""))
  }else{
    #stop(paste("Region id ", region_id, " has more variant than the max, see variable max_var (currently fixed at : ", max_var, sep =""))
  }


}#end of function