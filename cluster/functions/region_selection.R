############################################# Function to select the genetic region of interest #############################################
# The function allow to filter a data frame to retained only information within the desired region of interest
# Input : data: data frame with at least column such as
#                                   patient, chromosome, position, reference, altenative, zygosity, gene
#         region_type: [chr] must be equal to "chr", "pos" or "gene_symbol" (for the moment)
#         region_name: [list] of name of the region of interest, if no value then all region of the region type will be used (all genes, all chromosomes, etc)
# Output : data frame where column are patient, chromosome, position, reference, altenative, zygosity
region_selection <- function(data_patient = c(), 
                             data_control = c(),
                             genotype_matrix = c(),
                             region_type = c(), 
                             region_name = c(),
                             geneID = c()){
  if(missing(region_type)){
    msg<-sprintf("No region type given")
    warning(msg,call.=TRUE)
  }
  if(missing(region_name)){
    msg<-sprintf("No pregion name given")
    warning(msg,call.=TRUE)
  }
  
  if(length(region_name)==0){
    return(list(data_patient, data_control, genotype_matrix))
    
  }else {
    
    if(region_type == "GO" || region_type == "pathway"){
      geneID <- strsplit(geneID, split = "/")[[1]]
      region_name <- geneID 
    }
    focus_data_patient <- filter(data_patient, gene_symbol %in% unlist(region_name))
    focus_data_control <- filter(data_control, gene_symbol %in% unlist(region_name))
    }
    
    #remove column of genotype matrix not in the region of interest
    if(is.null(focus_data_patient$variant_id)){
      focus_data_patient$variant_id <- paste(focus_data_patient$chr, 
                                             focus_data_patient$pos,
                                             focus_data_patient$reference,
                                             focus_data_patient$alternative,
                                             sep = "_")
      focus_data_control$variant_id <- paste(focus_data_control$chr, 
                                             focus_data_control$pos,
                                             focus_data_control$reference,
                                             focus_data_control$alternative,
                                             sep = "_")
    }
    kept_id <- which(colnames(genotype_matrix) %in% unique(c(focus_data_patient$variant_id,focus_data_control$variant_id)))
    if(length(kept_id) == 1){
      focus_genotype_matrix <- data.frame(genotype_matrix[, kept_id], stringsAsFactors = FALSE)
      rownames(focus_genotype_matrix) <- rownames(genotype_matrix)
      colnames(focus_genotype_matrix) <- colnames(genotype_matrix)[kept_id]
    }else if(length(kept_id) == 0){
      #stop("region_selection problem, no variant left in region")
      focus_genotype_matrix <- data.frame(matrix(NA, ncol = 0, nrow = 0))
    }else{
      focus_genotype_matrix <- genotype_matrix[, kept_id]
    }
    
    #focus_genotype_matrix <- focus_genotype_matrix[which(row.names(focus_genotype_matrix) %in% c(unique(focus_data_patient$patient), unique(focus_data_control$patient))),]
    return(list(focus_data_patient, focus_data_control, focus_genotype_matrix))
  
  
}#end of function
#######################################################################################################################################
