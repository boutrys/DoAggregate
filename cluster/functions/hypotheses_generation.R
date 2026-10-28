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
  
  
  #keep only desired column
  #data_patient_annotated <- dplyr::select(data_patient_annotated, patient, pathology, chr, ..)
  #data_control_annotated <- dplyr::select(data_control_annotated, patient, pathology, chr, ..)
  
  return(list(data_patient_annotated, 
              data_control_annotated, 
              region_of_interest))
}