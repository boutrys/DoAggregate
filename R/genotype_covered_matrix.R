#' Generate a Cleaned Genotype and Coverage Matrix
#'
#' @name genotype_covered_matrix
#' @description This function processes genotype and coverage matrices to clean and filter the data based on variant coverage, removing low-quality individuals and positions. It also calculates data reliability and outputs diagnostic plots and cleaned datasets.
#'
#' @param coverage_matrix A data frame representing the coverage matrix, with positions in rows and individuals in columns.
#' @param genotype_matrix A data frame representing the genotype matrix, with positions in rows and individuals in columns.
#' @param variant A data frame containing variant information, including columns `id` and `position`.
#' @param patient_name A vector of patient sample names corresponding to columns in the matrices.
#' @param control_name A vector of control sample names corresponding to columns in the matrices.
#' @param min_coverage An integer specifying the minimum coverage threshold for a position to be considered reliable. Default is 20.
#' @param remove_individual A logical value indicating whether to remove individuals with a high number of poorly covered positions. Default is TRUE.
#' @param remove_position A logical value indicating whether to remove positions with a high number of poorly covered individuals. Default is TRUE.
#' @param path_store_result A character string specifying the directory path to store result files (e.g., plots, diagnostics).
#' @return A list containing cleaned genotype and coverage matrices, filtered variants, rejected variants, and (if `remove_individual` is TRUE) information about rejected individuals.
#'
#' @import ggplot2
#' @importFrom stats sd
#' @export
#' @importFrom utils globalVariables
utils::globalVariables(c("Nbr_badly_covered_position_per_individual", "group"))

genotype_covered_matrix <- function(coverage_matrix = c(),
                                    genotype_matrix = c(),
                                    variant = c(),
                                    patient_name = c(),
                                    control_name = c(),
                                    min_coverage = 20,
                                    remove_individual = TRUE,
                                    remove_position = TRUE,
                                    path_store_result = ""){
  #pre treatment and fusion of the 2 coverage matrices
  if(length(which(colnames(genotype_matrix) %in% variant$id)) == length(variant$position)){
    table_data_reliability <- data.frame(matrix(NA, ncol = 3, nrow = 5))
    colnames(table_data_reliability) <- c("Info", "Before_cleaning", "after_cleaning")
    table_data_reliability$Info <- c("Nbr_patient",
                                     "Nbr_control",
                                     "Nbr_unique_variant",
                                     "Total_nbr_variant",
                                     "reliability")
    total_nbr_var <- length(which(genotype_matrix > 0))
    table_data_reliability$Before_cleaning[1:4] <- c(length(patient_name), length(control_name), dim(genotype_matrix)[2], total_nbr_var)

    coverage_matrix$position <- paste(coverage_matrix$chr,coverage_matrix$pos,sep = "_")
    position_unique <- variant$position

    focus_coverage_matrix_restricted <- data.frame(t(coverage_matrix[which(coverage_matrix$position %in% position_unique),which(colnames(coverage_matrix) %in% c(patient_name,control_name))]))
    colnames(focus_coverage_matrix_restricted) <- coverage_matrix$position[which(coverage_matrix$position %in% position_unique)]
    rm(coverage_matrix)
    focus_coverage_matrix_restricted <- focus_coverage_matrix_restricted[,order(colnames(focus_coverage_matrix_restricted))]
    genotype_matrix <- genotype_matrix[,order(colnames(genotype_matrix))]
    focus_coverage_matrix_restricted <- focus_coverage_matrix_restricted[match(rownames(genotype_matrix), rownames(focus_coverage_matrix_restricted)),]
    #remove possible undesired positions and variants
    rm_pos_idx <- setdiff(position_unique, colnames(focus_coverage_matrix_restricted))
    if(length(rm_pos_idx) > 0){
      genotype_matrix <- genotype_matrix[,-which(colnames(genotype_matrix) %in% variant$id[which(position_unique %in% rm_pos_idx)])]
    }
    rm_var_idx <- setdiff(position_unique, colnames(focus_coverage_matrix_restricted))
    if(length(rm_var_idx) > 0){
      variant <- variant[-which(variant$id %in% variant$id[which(position_unique %in% rm_var_idx)]),]
      position_unique <- variant$position
    }

    duplicate_position_restricted <- position_unique[which(duplicated(position_unique))]
    duplicate_position_restricted <- duplicate_position_restricted[order(duplicate_position_restricted)]
    temp_table <- focus_coverage_matrix_restricted[,which(colnames(focus_coverage_matrix_restricted) %in% duplicate_position_restricted)]
    if(length(which(duplicated(duplicate_position_restricted))) > 0){
      tmp_pos <- duplicate_position_restricted[which(duplicated(duplicate_position_restricted))]
      to_add <- data.frame(matrix(NA, ncol = length(tmp_pos), nrow = dim(temp_table)[1]))
      colnames(to_add) <- tmp_pos
      for(i in 1:length(tmp_pos)){
        to_add[,i] <- data.frame(focus_coverage_matrix_restricted[,which(colnames(focus_coverage_matrix_restricted) == tmp_pos[i])])
      }
      rm(tmp_pos)
      temp_table <- bind_cols(temp_table, to_add)
      rm(to_add)
    }
    new_focus_coverage_matrix_restricted <- cbind(focus_coverage_matrix_restricted,temp_table)
    rm(temp_table)
    new_focus_coverage_matrix_restricted <- new_focus_coverage_matrix_restricted[,order(colnames(new_focus_coverage_matrix_restricted))]
    rm(duplicate_position_restricted)
    rownames(new_focus_coverage_matrix_restricted) <- rownames(focus_coverage_matrix_restricted)
    rm(focus_coverage_matrix_restricted)

    reliability <- length(which(new_focus_coverage_matrix_restricted > min_coverage)) / (dim(new_focus_coverage_matrix_restricted)[1]*dim(new_focus_coverage_matrix_restricted)[2])
    table_data_reliability$Before_cleaning[5] <- reliability

    #removing individuals
    if(remove_individual){
      bad_element_per_position <- rep(0, dim(new_focus_coverage_matrix_restricted)[2])
      bad_element_per_patient <- rep(0, dim(new_focus_coverage_matrix_restricted)[1])
      for (k in 1:dim(new_focus_coverage_matrix_restricted)[2]) {
        idx_individu <- which(new_focus_coverage_matrix_restricted[,k] < min_coverage)
        bad_element_per_patient[idx_individu] <- bad_element_per_patient[idx_individu] + 1
        bad_element_per_position[k] <- length(idx_individu)
      }
      median_bad_elem_per_patient <- stats::median(bad_element_per_patient)
      sd_bad_elem_per_patient <- sd(bad_element_per_patient)
      table_plot <- data.frame(matrix(NA, ncol = 2, nrow = length(c(patient_name,control_name))))
      colnames(table_plot) <- c("Nbr_badly_covered_position_per_individual", "group")
      table_plot$Nbr_badly_covered_position_per_individual <- bad_element_per_patient
      table_plot$group <- rep("control", dim(genotype_matrix)[1])
      table_plot$group[1:length(patient_name)] <- "patient"
      g <- ggplot(table_plot, aes(Nbr_badly_covered_position_per_individual)) +
        geom_histogram(aes(fill=group)) +
        geom_vline(xintercept = median_bad_elem_per_patient) +
        geom_vline(xintercept = median_bad_elem_per_patient + sd_bad_elem_per_patient, linetype = "dashed")
      options(bitmapType="cairo")
      ggsave(g, filename = paste(path_store_result, "bad_pos_per_individuals_before_removal.png", sep = ""), width = 7, height = 7, device = "png")
      idx_individu_rm <- which(bad_element_per_patient > median_bad_elem_per_patient + sd_bad_elem_per_patient)
      if(length(idx_individu_rm) > 0){
        new_focus_coverage_matrix_restricted <- new_focus_coverage_matrix_restricted[-idx_individu_rm,]
        genotype_matrix <- genotype_matrix[-idx_individu_rm,]
      }
      idx_pos_to_rm_bc_patient <- which(colSums(genotype_matrix) == 0)
      if(length(idx_pos_to_rm_bc_patient) > 0){
        genotype_matrix <- genotype_matrix[,-idx_pos_to_rm_bc_patient]
        new_focus_coverage_matrix_restricted <- new_focus_coverage_matrix_restricted[,-idx_pos_to_rm_bc_patient]
      }

      rejected_individu <- data.frame(matrix(NA, ncol = 2, nrow = length(idx_individu_rm)))
      colnames(rejected_individu) <- c("sample", "nbr_badly_covered_position")
      rejected_individu$sample <- c(patient_name, control_name)[idx_individu_rm]
      rejected_individu$nbr_badly_covered_position <- bad_element_per_patient[idx_individu_rm]
      patient_name <- patient_name[patient_name %in% rownames(new_focus_coverage_matrix_restricted)]
      control_name <- control_name[control_name %in% rownames(new_focus_coverage_matrix_restricted)]
    }

    #Update value OR first computation
    bad_element_per_position <- rep(0, dim(new_focus_coverage_matrix_restricted)[2])
    bad_element_per_patient <- rep(0, dim(new_focus_coverage_matrix_restricted)[1])
    for (k in 1:dim(new_focus_coverage_matrix_restricted)[2]) {
      idx_individu <- which(new_focus_coverage_matrix_restricted[,k] < min_coverage)
      bad_element_per_patient[idx_individu] <- bad_element_per_patient[idx_individu] + 1
      bad_element_per_position[k] <- length(idx_individu)
    }
    table_plot <- data.frame(matrix(NA, ncol = 2, nrow = length(c(patient_name,control_name))))
    colnames(table_plot) <- c("Nbr_badly_covered_position_per_individual", "group")
    table_plot$Nbr_badly_covered_position_per_individual <- bad_element_per_patient
    table_plot$group <- rep("control", dim(genotype_matrix)[1])
    table_plot$group[1:length(patient_name)] <- "patient"
    g <- ggplot(table_plot, aes(Nbr_badly_covered_position_per_individual)) + geom_histogram(aes(fill=group))
    if(remove_individual){
      ggsave(g, filename = paste(path_store_result, "bad_pos_per_individuals_after_individu_removal.png", sep = ""), width = 7, height = 7)
    }else{
      ggsave(g, filename = paste(path_store_result, "bad_pos_per_individuals.png", sep = ""), width = 7, height = 7)
    }

    #removing positions
    if(remove_position){
      median_bad_elem_per_position <- stats::median(bad_element_per_position)
      sd_bad_elem_per_position <- sd(bad_element_per_position)
      table_plot <- data.frame(bad_element_per_position)
      g <- ggplot(table_plot, aes(bad_element_per_position)) +
        geom_histogram() +
        geom_vline(xintercept = median_bad_elem_per_position) +
        geom_vline(xintercept = median_bad_elem_per_position + sd_bad_elem_per_position, linetype = "dashed")
      ggsave(g, filename = paste(path_store_result, "bad_individual_per_position.png", sep = ""), width = 7, height = 7)
      idx_pos_to_remove <- which(bad_element_per_position > median_bad_elem_per_position + sd_bad_elem_per_position)
      if(length(idx_pos_to_remove) > 0){
        reduced_coverage_matrix <- new_focus_coverage_matrix_restricted[,-idx_pos_to_remove]
        genotype_matrix_new <- genotype_matrix[,-idx_pos_to_remove]
      }else{
        reduced_coverage_matrix <- new_focus_coverage_matrix_restricted
        genotype_matrix_new <- genotype_matrix
      }
      rm(new_focus_coverage_matrix_restricted)
      rm(genotype_matrix)
      id_no_var_left <- which(rowSums(genotype_matrix_new) == 0)
      if(length(id_no_var_left) > 0){
        genotype_matrix_new <- genotype_matrix_new[-id_no_var_left,]
        id_rm <- which(rownames(reduced_coverage_matrix) %in% names(id_no_var_left))
        reduced_coverage_matrix <- reduced_coverage_matrix[-id_rm,]
        patient_name <- which(!patient_name %in% names(id_no_var_left))
        control_name <- which(!control_name %in% names(id_no_var_left))
        }
    }
    if(!remove_position){
      reduced_coverage_matrix <- new_focus_coverage_matrix_restricted
      genotype_matrix_new <- genotype_matrix
      rm(new_focus_coverage_matrix_restricted)
      rm(genotype_matrix)
    }

    #Update value OR first computation
    bad_element_per_position <- rep(0, dim(reduced_coverage_matrix)[2])
    bad_element_per_patient <- rep(0, dim(reduced_coverage_matrix)[1])
    for (k in 1:dim(reduced_coverage_matrix)[2]) {
      idx_individu <- which(reduced_coverage_matrix[,k] < min_coverage)
      bad_element_per_patient[idx_individu] <- bad_element_per_patient[idx_individu] + 1
      bad_element_per_position[k] <- length(idx_individu)
    }
    table_plot <- data.frame(bad_element_per_position)
    g <- ggplot(table_plot, aes(bad_element_per_position)) +
      geom_histogram()
    if(remove_position){
      ggsave(g, filename = paste(path_store_result, "bad_individual_per_position_after_position_removal.png", sep = ""), width = 7, height = 7)
      table_plot <- data.frame(matrix(NA, ncol = 2, nrow = length(c(patient_name,control_name))))
      colnames(table_plot) <- c("Nbr_badly_covered_position_per_individual", "group")
      table_plot$Nbr_badly_covered_position_per_individual <- bad_element_per_patient
      table_plot$group <- rep("control", dim(genotype_matrix_new)[1])
      table_plot$group[1:length(patient_name)] <- "patient"
      g <- ggplot(table_plot, aes(Nbr_badly_covered_position_per_individual)) + geom_histogram(aes(fill=group))
      ggsave(g, filename = paste(path_store_result, "bad_pos_per_individuals_final.png", sep = ""), width = 7, height = 7)
    }else{
      ggsave(g, filename = paste(path_store_result, "bad_individual_per_position.png", sep = ""), width = 7, height = 7)
    }

    reduced_variant <- variant %>% filter(variant$position %in% colnames(reduced_coverage_matrix))
    rejected_variant <- variant %>% filter(!variant$position %in% colnames(reduced_coverage_matrix))

    total_nbr_var <- length(which(genotype_matrix_new > 0))
    reliability <- length(which(reduced_coverage_matrix > min_coverage)) / (dim(reduced_coverage_matrix)[1]*dim(reduced_coverage_matrix)[2])
    table_data_reliability$after_cleaning <- c(length(patient_name), length(control_name), dim(genotype_matrix_new)[2], total_nbr_var, reliability)
    write_xlsx(table_data_reliability, paste(path_store_result, "impact_of_cleaning.xlsx", sep = ""))

    if(remove_individual){
      return(list(genotype_matrix_new,
                  reduced_coverage_matrix,
                  reduced_variant,
                  rejected_variant,
                  rejected_individu))
    }else{
      return(list(genotype_matrix_new,
                  reduced_coverage_matrix,
                  reduced_variant,
                  rejected_variant))
    }

  }else{
    stop("Problem not solved yet, contact the stupid Simon Boutry")
  }
}
