#' Run General Analysis on Several Experiments
#'
#' @name data_preparation
#' @description Prepares data for analysis by annotating VCF files with VEP, creating a BED file, and building genotype and coverage matrices. The function includes steps for harmonizing coverage data and cleaning input data based on thresholds and user preferences.
#'
#' @param data_patient A data frame containing patient data, with required columns:
#' \itemize{
#'   \item \code{sample}: character
#'   \item \code{chr}: character
#'   \item \code{pos}: integer
#'   \item \code{reference}: character
#'   \item \code{alternative}: character
#'   \item \code{zygosity}: integer
#' }
#' @param data_control A data frame containing control data, with the same required columns as \code{data_patient}.
#' @param mapqv_threshold Minimum mapping quality value for variants (default: 20).
#' @param min_coverage Minimum coverage threshold for filtering variants (default: 20).
#' @param remove_individual Logical. Whether to remove individuals not meeting coverage criteria (default: \code{TRUE}).
#' @param remove_position Logical. Whether to remove positions not meeting coverage criteria (default: \code{TRUE}).
#' @param path_to_store_plots A string specifying the directory to store generated plots and results (default: "").
#' @param coverage_matrix_path Path to the coverage matrix file, or 0 if no coverage matrix is provided (default: 0).
#'
#' @return A list containing processed data and results:
#' \itemize{
#'   \item \code{data_patient_annotated}: Annotated patient data.
#'   \item \code{data_control_annotated}: Annotated control data.
#'   \item \code{clean_genotype_matrix}: Cleaned genotype matrix.
#'   \item \code{missing_annotation}: Placeholder for missing annotations (currently \code{NA}).
#'   \item \code{variant_target}: Filtered variant information.
#'   \item \code{clean_coverage_matrix}: Cleaned coverage matrix (if applicable).
#'   \item \code{rejected_variant}: Variants removed due to coverage thresholds.
#'   \item \code{rejected_data_patient}: Rejected patient data due to coverage thresholds.
#'   \item \code{rejected_data_control}: Rejected control data due to coverage thresholds.
#' }
#'
#' @import stringr
#' @export
#' @importFrom utils globalVariables
utils::globalVariables(c("sample", "chr", "pos", "reference", "alternative", "zygosity", "position"))

data_preparation <- function(data_patient = c(),
                             data_control = c(),
                             mapqv_threshold = 20,
                             min_coverage = 20,
                             remove_individual = TRUE,
                             remove_position = TRUE,
                             path_to_store_plots = "",
                             coverage_matrix_path = 0){

  print("data_preparation start data preparation R")
  patient_name = unique(data_patient$sample)
  control_name = unique(data_control$sample)
  Nbr_patient <- length(patient_name)
  Nbr_control <- length(control_name)
  #genotype matrix
  result_geno_matrix <- genotypeMatrix(data_patient = data_patient,
                                       data_control = data_control,
                                       Nbr_patient_init = Nbr_patient,
                                       Nbr_control_init = Nbr_control)
  matrix_genotype <- data.frame(result_geno_matrix[[1]], row.names = c(patient_name,control_name))
  #list_variant_patient <- result_geno_matrix[[2]]  #useless, we don't use them for the moment
  #list_variant_control <- result_geno_matrix[[3]]
  list_variant <- result_geno_matrix[[4]]

  colnames(matrix_genotype) <- list_variant$id
  position_unique <- dplyr::select(list_variant[which(!duplicated(list_variant$position)),], c(chr,pos))


  ###coverage harmonization
  if(coverage_matrix_path != 0){
    print("data_preparation coverage matrix harmonization")

    matrix_coverage <- read_tsv(coverage_matrix_path, col_types = cols(chr = col_character()))

    tmp_chr <- grep(matrix_coverage$chr, pattern="chr")
    if(length(tmp_chr > 0)){
      matrix_coverage$chr <- str_remove_all(matrix_coverage$chr, pattern = "chr")
    }


    individu_id <- which(colnames(matrix_coverage) %in% c(patient_name, control_name))
    matrix_coverage <- matrix_coverage[,c(1,2,individu_id)]
    tmp_individu_cov <- colnames(matrix_coverage)[-c(1:2)]
    if(length(tmp_individu_cov) == 0){
      individu_not_cov <- c(patient_name, control_name)
    }else{
      tmp_id <- which(!tmp_individu_cov %in% c(patient_name, control_name))
      if(length(tmp_id) > 0){
        individu_not_cov <- tmp_individu_cov[tmp_id]
      }else{
        individu_not_cov <- 0
      }
    }

    tmp_position_cov <- paste(matrix_coverage$chr, matrix_coverage$pos, sep = "_")
    id_pos_notcov <- which(!list_variant$position %in% tmp_position_cov)
    if(length(id_pos_notcov) > 0){
      position_Notcov <- list_variant$position[id_pos_notcov]
    }else{
      position_Notcov <- 0
    }

    if(individu_not_cov != 0 | position_Notcov != 0){
      bug_path <- paste(path_to_store_plots, "Problems/", sep="")
      dir.create(bug_path)
      if(individu_not_cov != 0 && position_Notcov == 0){
        write_tsv(data.frame(sample = individu_not_cov), paste(bug_path, "individu_not_covered.tsv", sep=""), col_names = FALSE)
        write_tsv(matrix_coverage[,c(1,2)], paste(bug_path, "bed_to_do.tsv", sep=""))
        stop("Missing individuals in coverage matrix, Add missing individuals to the coverage matrix please")
      }else if(individu_not_cov != 0 && position_Notcov != 0){
        stop("Missing individuals AND positions in coverage matrix, redo the coverage matrix please")
      }else if(individu_not_cov == 0 && position_Notcov != 0){
        write_tsv(data.frame(sample = c(patient_name, control_name)), paste(bug_path, "individu_to_do.tsv", sep=""), col_names = FALSE)
        write_tsv(list_variant[id_pos_notcov, c(1,2)], paste(bug_path, "bed_to_do.tsv", sep=""))
        stop("Missing positions in coverage matrix, Add missing positions to the coverage matrix please")
      }
    }

    #remove non necessary info in input coverage matrix
    row_to_rm <- which(!tmp_position_cov %in% list_variant$position)
    if(length(row_to_rm) > 0){
      matrix_coverage <- matrix_coverage[-row_to_rm,]
    }

    #clean genotype matrix given coverage matrix
    clean_result <- genotype_covered_matrix(coverage_matrix = matrix_coverage,
                                            genotype_matrix = matrix_genotype,
                                            variant = list_variant,
                                            patient_name = patient_name,
                                            control_name = control_name,
                                            min_coverage = min_coverage,
                                            remove_individual = remove_individual,
                                            remove_position = remove_position,
                                            path_store_result = path_to_store_plots)
    clean_genotype_matrix <- clean_result[[1]]
    clean_coverage_matrix <- clean_result[[2]]
    variant_target <- clean_result[[3]]

    rejected_variant <- clean_result[[4]]
    if(remove_individual){
      rejected_individu <- clean_result[[5]]
    }
    #data management
    if(dim(rejected_individu)[1] > 0){
      data_patient <- filter(data_patient, ! sample %in% rejected_individu$sample)
      data_control <- filter(data_control, ! sample %in% rejected_individu$sample)
    }
    data_patient$position <- paste(data_patient$chr, data_patient$pos, sep = "_")
    clean_data_patient <- filter(data_patient, position %in% variant_target$position)
    rejected_data_patient <- filter(data_patient, !position %in% variant_target$position)
    data_control$position <- paste(data_control$chr, data_control$pos, sep = "_")
    clean_data_control <- filter(data_control, position %in% variant_target$position)
    rejected_data_control <- filter(data_control, !position %in% variant_target$position)

  }else{#No cleaning using coverage matrix
    print("data_preparation No cleaning to be done")

    clean_genotype_matrix <- matrix_genotype
    variant_target <- list_variant

    data_patient$position <- paste(data_patient$chr, data_patient$pos, sep = "_")
    clean_data_patient <- filter(data_patient, position %in% variant_target$position)

    data_control$position <- paste(data_control$chr, data_control$pos, sep = "_")
    clean_data_control <- filter(data_control, position %in% variant_target$position)
  }


  data_patient_annotated <- clean_data_patient
  data_control_annotated <- clean_data_control
  missing_annotation <- NA


  #Cohort analysis (Patient -> dendogram/PCA/... clustering patient)
  #similarity_cohort_matrix <- similarity_cohort(clean_genotype_matrix = clean_genotype_matrix,
  #                                              path_to_store_plots = path_to_store_plots)
  #writexl::write_xlsx(as.data.frame(similarity_cohort_matrix), paste(path_to_store_plots, "cohort_similarity_matrix.xlsx", sep = ""))

  #variant_id annotation
  data_patient_annotated$variant_id <- paste(data_patient_annotated$position, data_patient_annotated$reference, data_patient_annotated$alternative, sep = "_")
  data_control_annotated$variant_id <- paste(data_control_annotated$position, data_control_annotated$reference, data_control_annotated$alternative, sep = "_")

  if(coverage_matrix_path != 0){
    return(list(data_patient_annotated, data_control_annotated,
                clean_genotype_matrix, missing_annotation, variant_target,clean_coverage_matrix, rejected_variant,
                rejected_data_patient, rejected_data_control))
  }else{
    return(list(data_patient_annotated, data_control_annotated,
                clean_genotype_matrix, missing_annotation, variant_target))
  }

}
