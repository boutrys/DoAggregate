############################################# load_input #############################################
# The function allow load all needed input
load_input <- function(know_target_region = FALSE){
  
  #retrieve path to data of patient file 
  path_patient_data <- dlg_open(title = "select the file containing the data of the patient group")$res
  if(length(path_patient_data)<1){
    stop("no data for patient seleted, restart please!")
  }
  #retrieve path to data of control file
  path_control_data <- dlg_open(title = "select the file containing the data of the control group")$res
  if(length(path_control_data)<1){
    stop("no data for patient seleted, restart please!")
  }
  
  #retrieve path to data of control file
  if(know_target_region == TRUE){
    path_list_regions <- dlg_open(title = "select the file containing the data of list of regions of interest")$res
    if(length(path_list_regions)<1){
      stop("no data for patient seleted, restart please!")
    }
    return(list(path_patient_data, path_control_data, path_list_regions))
  }else{
    return(list(path_patient_data, path_control_data))
  }
}
