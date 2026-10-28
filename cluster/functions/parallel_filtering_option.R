############################################# Function to filter take different filtering option into account for several experiment #############################################
# The function allow to prepare the different filtering criteria and option to be considered for each analysis
# Input : 
#         
# Output : 
parallel_filtering_option <- function(KindOfRegion = list(),
                             filtering_criteria = list(),
                             filtering_level = list()){
  if(length(KindOfRegion) > 1){
    KindOfRegion <- paste(KindOfRegion, collapse="/")
    if(KindOfRegion == "gene/pathway"){
    name_of_experience <- c("gene", "pathway")
    KindOfRegion <- c("gene", "pathway")
    }else if(KindOfRegion == "over_represented_pathway"){
    name_of_experience <- c("over_represented_pathway")
    KindOfRegion <- c("pathway") 
    }
  }else if(KindOfRegion == "gene"){
    name_of_experience <- c("gene")
    KindOfRegion <- c("gene")
  }else if(KindOfRegion == "pathway"){
    name_of_experience <- c("pathway")
    KindOfRegion <- c("pathway")
  }else if(KindOfRegion == "over_represented_pathway"){
    name_of_experience <- c("over_represented_pathway")
    KindOfRegion <- c("pathway") 
  }else if(KindOfRegion == "sig_aggregate_gene"){
    name_of_experience <- c("sig_aggregate_gene")
    KindOfRegion <- c("gene")
  }

  if(length(filtering_level) > 0){
    Nbr_exp <- length(name_of_experience) * length(filtering_level)
    vec_hypotheses <- data.frame(matrix(NA,ncol = 3, nrow = Nbr_exp, byrow = TRUE))
    colnames(vec_hypotheses) <- c("filtering_criteria", "filtering_level", "genetic_region")
    vec_hypotheses$genetic_region <- name_of_experience
    vec_hypotheses$filtering_criteria <- filtering_criteria
    tmp_level <- c()
    for (i in 1:length(filtering_level)) {
      tmp_level <- c(tmp_level, rep(filtering_level[i], length(KindOfRegion)))
    }
    vec_hypotheses$filtering_level <- tmp_level
    vec_hypotheses$name_of_experience <- NA
    for (i in 1:Nbr_exp) {
      tmp <- unlist(vec_hypotheses$filtering_criteria[i])
      tmp_level <- unlist(vec_hypotheses$filtering_level[i])
      tmp_name <- c()
      if(length(tmp) == length(tmp_level)){
        for (j in 1:length(tmp)) {
          tmp_name <- paste(tmp_name, paste(tmp[j], tmp_level[j], sep = "_"), sep = "")
        }
        vec_hypotheses$name_of_experience[i] <- paste(tmp_name, vec_hypotheses$genetic_region[i], sep = "_")
      }else{
        stop("Number of filtering criteria and number of filtering level do not match!")
      }
    }
    KindOfRegion <- rep(KindOfRegion, length(filtering_level))
  }else{
    Nbr_exp <- length(name_of_experience)
    vec_hypotheses <- data.frame(matrix(NA,ncol = 3, nrow = Nbr_exp, byrow = TRUE))
    colnames(vec_hypotheses) <- c("filtering_criteria", "filtering_level", "genetic_region")
    vec_hypotheses$genetic_region <- name_of_experience
    vec_hypotheses$filtering_criteria <- "Default"
    vec_hypotheses$filtering_level <- "Default"
  }

  return(list(vec_hypotheses, Nbr_exp, KindOfRegion))
}#end of function