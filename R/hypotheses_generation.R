#' Generate Hypotheses from Annotated Patient and Control Data
#'
#' This function filters annotated datasets of patients and controls based on specified criteria and then generates hypotheses
#' at different levels of genetic analysis such as gene, pathway, or gene ontology (GO) terms. The function allows flexible
#' filtering with dynamic criteria and supports analysis at multiple biological levels, facilitating targeted hypothesis generation.
#'
#' @param data_patient_annotated A data frame with annotated data for patients.
#' @param data_control_annotated A data frame with annotated data for controls.
#' @param filtering_criteria A list containing the criteria names to filter the datasets.
#' @param filtering_level A list containing the threshold levels for each filtering criterion.
#' @param filtering_direction A list indicating the direction of filtering ('+' for above or '-' for below the threshold).
#' @param level_of_analysis A vector indicating the levels of analysis ('gene', 'pathway', 'GO') to generate hypotheses.
#' @return A list containing the filtered patient data, filtered control data, and a vector of regions of interest based on the level of analysis.
#' @export
############################################# Function for hypotheses generation #############################################
hypotheses_generation <- function(data_patient_annotated = c(),
                                  data_control_annotated = c(),
                                  filtering_criteria = list(),
                                  filtering_level = list(),
                                  filtering_direction = list(),
                                  level_of_analysis = c("gene", "pathway", "GO")){
  if(!filtering_criteria == "Default"){
    ############## Hypotheses based on filtering criterias ##############
    #patient dataset
    if(!is.na(filtering_criteria)){
      for (i in 1:length(filtering_criteria[[1]])) {
        tmp_col <- which(colnames(data_patient_annotated) == filtering_criteria[[1]][i])
        if(length(tmp_col) > 0){
          if(filtering_direction[[1]][i] == "-"){
            data_patient_annotated <- data_patient_annotated[which(data_patient_annotated[,tmp_col] <= filtering_level[[1]][i]),]
          }else{
            data_patient_annotated <- data_patient_annotated[which(data_patient_annotated[,tmp_col] >= filtering_level[[1]][i]),]
          }
        }else{
          stop(paste(filtering_criteria[[1]][i], "missing column annotation for filtering criteria in patient dataset", sep = " "))
        }
      }#end of for going through all criteria
    }#end of if filtering criteria

    #control dataset
    if(!is.na(filtering_criteria)){
      for (i in 1:length(filtering_criteria[[1]])) {
        tmp_col <- which(colnames(data_control_annotated) == filtering_criteria[[1]][i])
        if(length(tmp_col) > 0){
          if(filtering_direction[[1]][i] == "-"){
            data_control_annotated <- data_control_annotated[which(data_control_annotated[,tmp_col] <= filtering_level[[1]][i]),]
          }else{
            data_control_annotated<- data_control_annotated[which(data_control_annotated[,tmp_col] >= filtering_level[[1]][i]),]
          }
        }else{
          stop(paste(filtering_criteria[[1]][i], "missing column annotation for filtering criteria in patient dataset", sep = " "))
        }
      }#end of for going through all criteria
    }#end of if filtering criteria
  }


  ############## Hypotheses based on level of genetic region ##############
  region_of_interest <- c()

  #gene level
  if(level_of_analysis[[1]] == "gene"){
    region_of_interest <- unique(data_patient_annotated$gene_symbol)
  }
  #old way to analyze all pathway/GO terms present in data, thanks to fullAnnotator using our own tables downloaded from consensusPathDB and transformed with in-house scripts
  if(FALSE){
    #pathway level
    if(level_of_analysis[[1]] == "pathway"){
      Nbr_pathway <- length(unique(data_patient_annotated$Pathway_list))
      for(i in 1:Nbr_pathway){
        tmp_pathway_list <- unique(strsplit(data_patient_annotated$Pathway_list[i], split = ";")[[1]][-1])
        region_of_interest <- c(region_of_interest, tmp_pathway_list)
      }
      region_of_interest <- unique(region_of_interest)
    }
    #GO level
    if(level_of_analysis[[1]] == "GO"){
      Nbr_GO <- length(unique(data_patient_annotated$GO_list))
      for(i in 1:Nbr_GO){
        tmp_GO_list <- unique(strsplit(data_patient_annotated$GO_list[i], split = ";")[[1]][-1])
        tmp_GO_list <- tmp_GO_list[-which(tmp_GO_list == "")]
        tmp_GO_list <- tmp_GO_list[which(!is.na(tmp_GO_list))]
        region_of_interest <- c(region_of_interest, tmp_GO_list)
      }
      region_of_interest <- unique(region_of_interest)

    }
  }


  #keep only desired column
  #data_patient_annotated <- dplyr::select(data_patient_annotated, patient, pathology, chr, ..)
  #data_control_annotated <- dplyr::select(data_control_annotated, patient, pathology, chr, ..)

  return(list(data_patient_annotated,
              data_control_annotated,
              region_of_interest))
}
