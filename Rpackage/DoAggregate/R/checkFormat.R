#' Check Format of Data
#'
#' @name checkFormat
#' @description Ensures that the input data frames have the desired format.
#' It checks for the presence of required columns, removes rows with missing values in key columns,
#' and standardizes chromosome values (e.g., converting "23" and "24" to "X" and "Y").
#'
#' @param patient A data frame containing at least the following columns:
#' \itemize{
#'   \item \code{patient}: character
#'   \item \code{chromosome}: character
#'   \item \code{position}: integer
#'   \item \code{reference}: character
#'   \item \code{alternative}: character
#'   \item \code{zygosity}: integer
#'   \item \code{gene_symbol}: character
#' }
#' @param control A data frame with the same required columns as \code{patient}.
#'
#' @return A list containing the two cleaned and formatted data frames:
#' \itemize{
#'   \item \code{patient}: Cleaned and formatted patient data frame.
#'   \item \code{control}: Cleaned and formatted control data frame.
#' }
#'
#'
#' @export
#' @importFrom utils globalVariables
utils::globalVariables(c("gene_symbol", "chr", "pos", "reference", "alternative", "zygosity", "sample"))

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
