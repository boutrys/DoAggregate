#' Construct Genotype Matrix
#'
#' @name genotypeMatrix
#' @description This function constructs a genotype matrix from patient and control data. The function expects data frames
#' containing genetic variant information such as chromosome position, reference, and alternate alleles,
#' as well as zygosity. The output is a genotype matrix where columns represent variants and rows represent individuals.
#' Additional output includes lists of unique variants specific to patients and controls.
#'
#'
#' @param data_patient A data frame containing the columns: sample, chr, pos, reference, alternative, zygosity.
#' @param data_control A data frame containing the columns: sample, chr, pos, reference, alternative, zygosity.
#' @param Nbr_patient_init The initial number of patients expected (used for matrix dimensioning).
#' @param Nbr_control_init The initial number of controls expected (used for matrix dimensioning).
#' @return A list containing the genotype matrix and three lists of variant identifiers: those specific to patients,
#'         those specific to controls, and a combined list from both sets.
#' @export
#' @importFrom utils globalVariables
utils::globalVariables(c("sample", "chr", "pos", "reference", "alternative", "zygosity"))


genotypeMatrix <- function(data_patient=c(), data_control=c(), Nbr_patient_init=0, Nbr_control_init=0){
  data_patient <- data.frame(data_patient) #avoid problem using bind_rows
  data_control <- data.frame(data_control) #avoid problem using bind_rows

  data_patient <- data_patient %>% dplyr::select(sample,chr,pos,reference,alternative,zygosity) %>%
    filter(!is.na(sample)) %>% filter(!is.na(chr)) %>% filter(!is.na(pos)) %>%
    filter(!is.na(reference)) %>% filter(!is.na(alternative))  %>% filter(!is.na(zygosity))
  data_control <- data_control %>% dplyr::select(sample,chr,pos,reference,alternative,zygosity) %>%
    filter(!is.na(sample)) %>% filter(!is.na(chr)) %>% filter(!is.na(pos)) %>%
    filter(!is.na(reference)) %>% filter(!is.na(alternative))  %>% filter(!is.na(zygosity))

  total <- Nbr_patient_init + Nbr_control_init

  if((dim(data_patient)[1]==0)&(dim(data_control)[1]==0)){
    newG <- matrix(0,nrow = total,ncol = 1)
    return(list(newG, c(), c()))
  }else{
    Nbr_patient_region <- length(unique(data_patient$sample))
    Nbr_control_region <- length(unique(data_control$sample))


    id_hete_patient <- which(data_patient$zygosity=="Heterozygous")
    id_homo_patient <- which(data_patient$zygosity=="Homozygous")
    id_hete_control <- which(data_control$zygosity=="Heterozygous")
    id_homo_control <- which(data_control$zygosity=="Homozygous")
    if(length(id_hete_patient) > 0){data_patient[id_hete_patient,]$zygosity <- 1}
    if(length(id_homo_patient) > 0){data_patient[id_homo_patient,]$zygosity <- 2}
    if(length(id_hete_control) > 0){data_control[id_hete_control,]$zygosity <- 1}
    if(length(id_homo_control) > 0){data_control[id_homo_control,]$zygosity <- 2}

    #merging the patient and the control
    newdata <- bind_rows(data_patient,data_control)

    #take the unique variant present within the data
    variant_unique_patient <- data.frame() #initialize empty data frame
    variant_unique_control <- data.frame() #initialize empty data frame
    variant_unique_patient <- distinct(data_patient,chr,pos,reference,alternative)
    variant_unique_control <- distinct(data_control,chr,pos,reference,alternative)
    variant_unique <- data.frame() #initialize empty data frame
    variant_unique <- distinct(newdata,chr,pos,reference,alternative)
    variant_unique$position <- paste(variant_unique$chr, variant_unique$pos, sep = "_")
    variant_unique_intersection <- dplyr::intersect(variant_unique_patient,variant_unique_control)
    variant_unique_specificTopatient <- dplyr::setdiff(variant_unique_patient,variant_unique_intersection)
    variant_unique_specificTocontrol <- dplyr::setdiff(variant_unique_control,variant_unique_intersection)

    #preprocessing to hasten computation
    list_variant_unique_specificTopatient <- variant_unique_specificTopatient %>% dplyr::mutate(id = paste(chr,pos,reference,alternative, sep = "_"))
    list_variant_unique_specificTocontrol <- variant_unique_specificTocontrol %>% dplyr::mutate(id = paste(chr,pos,reference,alternative, sep = "_"))
    variant_table <- variant_unique %>% dplyr::mutate(id = paste(chr,pos,reference,alternative, sep = "_"))

    Nbr_variant <- dim(variant_unique)[1]
    unique_patient <- unique(newdata$sample) #list of the patient

    G = matrix(0,nrow = total,ncol = Nbr_variant) #initialization of genotype matrix with all zero entry
    for(i in 1:length(unique_patient)){
      name <- unique_patient[i]
      variant_patient_hete <- filter(newdata, sample == name, zygosity == "1") %>% dplyr::mutate(id = paste(chr,pos,reference,alternative, sep = "_"))
      variant_patient_homo <- filter(newdata, sample == name, zygosity == "2") %>% dplyr::mutate(id = paste(chr,pos,reference,alternative, sep = "_"))
      inter_data <- c()
      inter_data_hete <- which(variant_table$id %in% variant_patient_hete$id)
      inter_data_homo <- which(variant_table$id %in% variant_patient_homo$id)
      G[i,inter_data_hete] <-  1
      G[i,inter_data_homo] <-  2
    }
    #adjust genotype matrix format
    if((Nbr_patient_init == Nbr_patient_region)&(Nbr_control_init == Nbr_control_region)){
      newG <- G
    } else if(Nbr_patient_region == 0 | Nbr_control_region == 0){
      if(Nbr_patient_region == 0){
        newG = matrix(0,nrow = total,ncol = Nbr_variant)
        index_control_begin <- Nbr_patient_init + 1
        index_control_end <- index_control_begin + Nbr_control_region - 1
        newG[index_control_begin:index_control_end,] <- G[(Nbr_patient_region+1):(Nbr_patient_region + Nbr_control_region),]
      } else if(Nbr_control_region == 0){
        newG = matrix(0,nrow = total,ncol = Nbr_variant)
        newG[1:Nbr_patient_region,] <- G[1:Nbr_patient_region,]
      } else {
        msg<-sprintf("PROBLEM FOR GENOTYPE MATRIX COMPUTATION, ONE MISSING CASE")
        warning(msg,call.=TRUE)
        }
    } else{
      newG = matrix(0,nrow = total,ncol = Nbr_variant)
      newG[1:Nbr_patient_region,] <- G[1:Nbr_patient_region,] #assign the value to the patient present in this region
      index_control_begin <- Nbr_patient_init + 1
      index_control_end <- index_control_begin + Nbr_control_region - 1
      newG[index_control_begin:index_control_end,] <- G[(Nbr_patient_region+1):(Nbr_patient_region + Nbr_control_region),] #assign the value to the controls present in this region
    }

    return(list(newG, list_variant_unique_specificTopatient$id, list_variant_unique_specificTocontrol$id, variant_table))
  }#end of else
}
