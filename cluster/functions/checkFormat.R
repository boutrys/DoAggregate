############################################# Function to check format of data  #############################################
# The function allow to ensure that input data have desired format
# Input : patient: data frame with at least column such as
#                                   patient, chromosome, position, reference, altenative, zygosity, gene_symbol
#         control: data frame with at least column such as
#                                   patient, chromosome, position, reference, altenative, zygosity, gene_symbol
# Output : the two input data frame with desired format of column
checkFormat <- function(patient = c(), control = c()){
  #file must have at least a patient column in character, 
  # chr in character, pos in interger, reference&alternative in character, zygosity in integer
  patient <- data.frame(patient, stringsAsFactors = FALSE)
  control <- data.frame(control, stringsAsFactors = FALSE)
  
  #filter non annotated data
  patient <- patient %>% filter(!is.na(sample)) %>%
    filter(!is.na(gene_symbol)) %>% filter(!is.na(chr)) %>% 
    filter(!is.na(pos)) %>% filter(!is.na(reference)) %>% 
    filter(!is.na(alternative))  %>% filter(!is.na(zygosity))
  control <- control %>% filter(!is.na(sample)) %>% 
    filter(!is.na(gene_symbol)) %>% filter(!is.na(chr)) %>% 
    filter(!is.na(pos)) %>% filter(!is.na(reference)) %>% 
    filter(!is.na(alternative))  %>% filter(!is.na(zygosity))
  
  patient$chr <- as.character(patient$chr, as.is = TRUE)
  control$chr <- as.character(control$chr, as.is = TRUE)
  
  id_patient_X <- which(patient$chr == "23")
  id_patient_Y <- which(patient$chr == "24")
  id_control_X <- which(control$chr == "23")
  id_control_Y <- which(control$chr == "24")
  if(length(id_patient_X) > 0){patient[id_patient_X,]$chr <- "X"}
  if(length(id_patient_Y) > 0){patient[id_patient_Y,]$chr <- "Y"}
  if(length(id_control_X) > 0){control[id_control_X,]$chr <- "X"}
  if(length(id_control_Y) > 0){control[id_control_Y,]$chr <- "Y"}
  
  return(list(patient,control))
}