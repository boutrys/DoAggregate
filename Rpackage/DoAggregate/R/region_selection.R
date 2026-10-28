#' Select Genetic Region of Interest
#'
#' @name region_selection
#' @description This function filters a dataset to retain only information within a specified region of interest,
#' such as a chromosome, position, or gene symbol. It is primarily used to narrow down the data
#' for focused genetic analysis in patient and control datasets along with genotype matrices.
#'
#' @param data_patient A data frame containing patient genetic data.
#' @param data_control A data frame containing control genetic data.
#' @param genotype_matrix A matrix or data frame of genotype information.
#' @param region_type A character string indicating the type of region to filter by; options include "chr", "pos", or "gene_symbol".
#' @param region_name A list of names specifying the regions of interest. If left empty, all regions of the specified type are included.
#' @param geneID Optional; used if region type pertains to pathways or GO terms to specify gene IDs.
#' @return A list containing the filtered patient data, control data, and genotype matrix.
#' @export
#' @importFrom utils globalVariables
utils::globalVariables(c("gene_symbol"))


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

    if(region_type == "GO" || region_type == "pathway" || region_type == "over_represented_GO" || region_type == "over_represented_pathway"){
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
