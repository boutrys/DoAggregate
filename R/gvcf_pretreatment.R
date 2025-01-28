#' Pre-treat VCF Data
#'
#' @name gvcf_pretreatment
#' @description This function preprocesses VCF (and gVCF) data for variant analysis. It reads the gVCF file, renames columns according
#' to a provided mapping, and filters variants based on the minor allele frequency (MAF). It also identifies and renames genes,
#' performs data cleaning based on coverage and allele frequency, and splits the data into cases and controls.
#'
#' @param path_data The file path to the input gVCF data.
#' @param path_output The output directory for saving processed files and reports.
#' @param cases The file path to the list of case samples.
#' @param controls The file path to the list of control samples.
#' @param mandatory_column_mapping A list mapping gVCF column names to required names for processing.
#' @param MAF The minor allele frequency threshold for variant filtering.
#' @param het_var A vector of genotype codes representing heterozygous variants.
#' @param homo_var A vector of genotype codes representing homozygous variants.
#' @param homo_ref A vector of genotype codes representing homozygous reference.
#' @param cleaning A boolean flag to perform data cleaning based on coverage and variant quality.
#' @param remove_individual A boolean flag to remove individuals based on poor coverage across many positions.
#' @param remove_position A boolean flag to remove poorly covered positions across many individuals.
#' @param min_coverage The minimum coverage required for a position to be considered reliable.
#' @return A list containing two data frames: one for case data and one for control data after preprocessing.
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @importFrom stats sd median
#' @export
#' @importFrom utils globalVariables
utils::globalVariables(c("gnomad_wes_af", "id_variant", "gt_DP", "Nbr_badly_covered_position_per_individual", "group", "zygosity"))

gvcf_pretreatment <- function(path_data = "",
                              path_output = "",
                              cases = "",
                              controls = "",
                              mandatory_column_mapping =  list(sample = "Indiv",
                                                               chr = "CHROM",
                                                               pos = "POS",
                                                               reference = "REF",
                                                               alternative = "ALT",
                                                               zygosity = "",
                                                               gene_symbol = "SYMBOL",
                                                               gnomad_wes_af = "gnomad40_AF",
                                                               gt_GT = "gt_GT",
                                                               gt_DP = "gt_DP"),
                              MAF = 0.01,
                              het_var = c("0/1", "0|1"),
                              homo_var = c("1/1", "1|1"),
                              homo_ref = c("0/0", "0|0"),
                              cleaning = TRUE,
                              remove_individual = TRUE,
                              remove_position = TRUE,
                              min_coverage = 10){

  ###Read input data
  data <- read.vcfR(path_data)
  transformed <- vcfR2tidy(data, single_frame = TRUE)


  ###Pre process data
  rm(data)
  colnames(transformed$dat)[which(colnames(transformed$dat) == mandatory_column_mapping$sample)] <- names(mandatory_column_mapping)[1]
  colnames(transformed$dat)[which(colnames(transformed$dat) == mandatory_column_mapping$chr)] <- names(mandatory_column_mapping)[2]
  colnames(transformed$dat)[which(colnames(transformed$dat) == mandatory_column_mapping$pos)] <- names(mandatory_column_mapping)[3]
  if(length(which(colnames(transformed$dat) == "REF...4")) > 0){
    colnames(transformed$dat)[which(colnames(transformed$dat) == "REF...4")] <- "REF"
  }
  colnames(transformed$dat)[which(colnames(transformed$dat) == mandatory_column_mapping$reference)] <- names(mandatory_column_mapping)[4]
  colnames(transformed$dat)[which(colnames(transformed$dat) == mandatory_column_mapping$alternative)] <- names(mandatory_column_mapping)[5]
  if(mandatory_column_mapping$zygosity == ""){
    transformed$dat$zygosity <- NA
    transformed$dat$zygosity[which(transformed$dat$gt_GT %in% het_var)] <- "Heterozygous"
    transformed$dat$zygosity[which(transformed$dat$gt_GT %in% homo_var)] <- "Homozygous"
    transformed$dat$zygosity[which(transformed$dat$gt_GT %in% homo_ref)] <- "Reference"
  }else{
    colnames(transformed$dat)[which(colnames(transformed$dat) == mandatory_column_mapping$zygosity)] <- names(mandatory_column_mapping)[6]
  }
  colnames(transformed$dat)[which(colnames(transformed$dat) == mandatory_column_mapping$gene_symbol)] <- names(mandatory_column_mapping)[7]
  colnames(transformed$dat)[which(colnames(transformed$dat) == mandatory_column_mapping$gnomad_wes_af)] <- names(mandatory_column_mapping)[8]
  colnames(transformed$dat)[which(colnames(transformed$dat) == mandatory_column_mapping$gt_GT)] <- names(mandatory_column_mapping)[9]
  colnames(transformed$dat)[which(colnames(transformed$dat) == mandatory_column_mapping$gt_DP)] <- names(mandatory_column_mapping)[10]

  transformed$dat$gnomad_wes_af <- as.double(transformed$dat$gnomad_wes_af)


  ### Filtering of data
  #Using MAF
  data_filtered <- transformed$dat %>%
    filter(gnomad_wes_af <= MAF | is.na(gnomad_wes_af))
  rm(transformed)


  ###Check gene symbol
  geneNameToChange <- grep(pattern = "\\", data_filtered$gene_symbol, fixed = TRUE)
  list_to_change <- unique(data_filtered$gene_symbol[geneNameToChange])

  for (i in 1:length(list_to_change)) {
    tmp_id <- which(data_filtered$gene_symbol[geneNameToChange] == list_to_change[i])
    tmp_new_name <- strsplit(data_filtered$gene_symbol[geneNameToChange[tmp_id[1]]], split = "\\", fixed = TRUE)[[1]][1]
    data_filtered$gene_symbol[geneNameToChange[tmp_id]] <- tmp_new_name
  }


  ###Load list of cases and controls
  cases <- read_tsv(cases, col_names = FALSE)
  controls <- read_tsv(controls, col_names = FALSE)


  ### Data cleaning
  if(cleaning){
    table_data_reliability <- data.frame(matrix(NA, ncol = 3, nrow = 5))
    colnames(table_data_reliability) <- c("Info", "Before_cleaning", "after_cleaning")
    table_data_reliability$Info <- c("Nbr_patient",
                                     "Nbr_control",
                                     "Nbr_unique_variant",
                                     "Total_nbr_variant",
                                     "reliability")
    total_nbr_var <- length(which(data_filtered$zygosity %in% c("Heterozygous", "Homozygous")))
    patient_name <- cases$X1[which(cases$X1 %in% data_filtered$sample)]
    control_name <- controls$X1[which(controls$X1 %in% data_filtered$sample)]
    data_filtered$id_variant <- paste(data_filtered$chr, data_filtered$pos, data_filtered$reference, data_filtered$alternative, sep = "_")
    list_variant <- unique(data_filtered$id_variant)
    nbr_variant <- length(list_variant)
    table_data_reliability$Before_cleaning[1:4] <- c(length(patient_name), length(control_name), nbr_variant, total_nbr_var)


    ### Build coverage
    data_coverage <- data_filtered[,c("sample", "id_variant" ,"gt_DP")] %>%
      pivot_wider(
        names_from = id_variant,
        values_from = gt_DP,
        values_fill = list(gt_DP = NA_integer_)
      )
    #reliability of the data given coverage
    reliability <- length(which(data_coverage[,-1] >= min_coverage)) / (dim(data_coverage[,-1])[1]*dim(data_coverage[,-1])[2])
    table_data_reliability$Before_cleaning[5] <- reliability


    ###remove individuals that have badly covered position more than median + sd of badly covered position
    bad_position_per_individu <- c()
    nbr_individu <- length(patient_name) + length(control_name)
    for (i in 1:nbr_individu) {
      bad_position_per_individu <- c(bad_position_per_individu, length(which(data_coverage[i,-1] < min_coverage)))
    }
    if(remove_individual){
      table_plot <- data.frame(matrix(NA, ncol = 2, nrow = length(c(patient_name,control_name))))
      colnames(table_plot) <- c("Nbr_badly_covered_position_per_individual", "group")
      table_plot$Nbr_badly_covered_position_per_individual <- bad_position_per_individu
      table_plot$group[which(data_coverage$sample %in% cases$X1)] <- "cases"
      table_plot$group[which(data_coverage$sample %in% controls$X1)] <- "controls"
      median_bad_elem_per_individu <- median(bad_position_per_individu)
      sd_bad_elem_per_individu <- sd(bad_position_per_individu)
      g <- ggplot(table_plot, aes(Nbr_badly_covered_position_per_individual)) +
        geom_histogram(aes(fill=group)) +
        geom_vline(xintercept = median_bad_elem_per_individu) +
        geom_vline(xintercept = median_bad_elem_per_individu + sd_bad_elem_per_individu, linetype = "dashed")
      options(bitmapType="cairo")
      ggsave(g, filename = paste(path_output, "bad_pos_per_individuals_before_removal.png", sep = "/"), width = 7, height = 7, device = "png")

      id_rm_individu <- which(bad_position_per_individu > (median_bad_elem_per_individu + sd_bad_elem_per_individu))
      if(length(id_rm_individu) > 0){
        #Update values
        data_coverage <- data_coverage[-id_rm_individu,]
        data_filtered <- data_filtered[-which(data_filtered$sample %in% data_coverage$sample[id_rm_individu]),]
        remove_case <- patient_name[which(!patient_name %in% unique(data_filtered$sample))]
        if(length(remove_case) > 0){
          write_xlsx(data.frame(remove_cases = remove_case), paste(path_output, "removed_cases.xlsx", sep = "/"))
        }
        remove_control <- control_name[which(!control_name %in% data_filtered$sample)]
        if(length(remove_control) > 0){
          write_xlsx(data.frame(remove_controls = remove_control), paste(path_output, "removed_controls.xlsx", sep = "/"))
        }
        patient_name <- cases$X1[which(cases$X1 %in% data_filtered$sample)]
        control_name <- controls$X1[which(controls$X1 %in% data_filtered$sample)]
        nbr_individu <- length(patient_name) + length(control_name)
        bad_position_per_individu <- c()
        for (i in 1:nbr_individu) {
          bad_position_per_individu <- c(bad_position_per_individu, length(which(data_coverage[i,-1] < min_coverage)))
        }
      }

      table_plot <- data.frame(matrix(NA, ncol = 2, nrow = length(c(patient_name,control_name))))
      colnames(table_plot) <- c("Nbr_badly_covered_position_per_individual", "group")
      table_plot$Nbr_badly_covered_position_per_individual <- bad_position_per_individu
      table_plot$group[which(data_coverage$sample %in% cases$X1)] <- "cases"
      table_plot$group[which(data_coverage$sample %in% controls$X1)] <- "controls"
      g <- ggplot(table_plot, aes(Nbr_badly_covered_position_per_individual)) + geom_histogram(aes(fill=group))
      ggsave(g, filename = paste(path_output, "bad_pos_per_individuals_after_individu_removal.png", sep = "/"), width = 7, height = 7)


    }else{
      table_plot <- data.frame(matrix(NA, ncol = 2, nrow = length(c(patient_name,control_name))))
      colnames(table_plot) <- c("Nbr_badly_covered_position_per_individual", "group")
      table_plot$Nbr_badly_covered_position_per_individual <- bad_position_per_individu
      table_plot$group[which(data_coverage$sample %in% cases$X1)] <- "cases"
      table_plot$group[which(data_coverage$sample %in% controls$X1)] <- "controls"
      g <- ggplot(table_plot, aes(Nbr_badly_covered_position_per_individual)) + geom_histogram(aes(fill=group))
      ggsave(g, filename = paste(path_output, "bad_pos_per_individuals.png", sep = "/"), width = 7, height = 7)
    }


    ###remove position for which the number of badly covered individual is  more than median + sd of badly covered individual
    bad_individu_per_position <- c()
    end <- nbr_variant + 1
    for (i in 2:end) {
      bad_individu_per_position <- c(bad_individu_per_position, length(which(data_coverage[,i] < min_coverage)))
    }
    if(remove_position){
      median_bad_elem_per_position <- median(bad_individu_per_position)
      sd_bad_elem_per_position <- sd(bad_individu_per_position)
      table_plot <- data.frame(bad_individu_per_position)
      g <- ggplot(table_plot, aes(bad_individu_per_position)) +
        geom_histogram() +
        geom_vline(xintercept = median_bad_elem_per_position) +
        geom_vline(xintercept = median_bad_elem_per_position + sd_bad_elem_per_position, linetype = "dashed")
      ggsave(g, filename = paste(path_output, "bad_individual_per_position.png", sep = "/"), width = 7, height = 7)

      id_rm_position <- which(bad_individu_per_position > (median_bad_elem_per_position + sd_bad_elem_per_position))

      if(length(id_rm_position) > 0){
        #Update values
        coverage_id_rm_position <- id_rm_position + 1
        data_coverage <- data_coverage[,-coverage_id_rm_position]
        nbr_variant <- nbr_variant - length(id_rm_position)
        list_variant <- list_variant[-id_rm_position]
        data_filtered <- data_filtered[-which(!data_filtered$id_variant %in%  list_variant),]
        bad_individu_per_position <- c()
        end <- nbr_variant + 1
        for (i in 2:end) {
          bad_individu_per_position <- c(bad_individu_per_position, length(which(data_coverage[,i] < min_coverage)))
        }
        bad_position_per_individu <- c()
        for (i in 1:nbr_individu) {
          bad_position_per_individu <- c(bad_position_per_individu, length(which(data_coverage[i,-1] < min_coverage)))
        }
      }
      table_plot <- data.frame(bad_individu_per_position)
      g <- ggplot(table_plot, aes(bad_individu_per_position)) +
        geom_histogram()
      ggsave(g, filename = paste(path_output, "bad_individual_per_position_after_position_removal.png", sep = "/"), width = 7, height = 7)
      table_plot <- data.frame(matrix(NA, ncol = 2, nrow = length(c(patient_name,control_name))))
      colnames(table_plot) <- c("Nbr_badly_covered_position_per_individual", "group")
      table_plot$Nbr_badly_covered_position_per_individual <- bad_position_per_individu
      table_plot$group[which(data_coverage$sample %in% cases$X1)] <- "cases"
      table_plot$group[which(data_coverage$sample %in% controls$X1)] <- "controls"
      g <- ggplot(table_plot, aes(Nbr_badly_covered_position_per_individual)) + geom_histogram(aes(fill=group))
      ggsave(g, filename = paste(path_output, "bad_pos_per_individuals_final.png", sep = "/"), width = 7, height = 7)

    }else{
      table_plot <- data.frame(bad_individu_per_position)
      g <- ggplot(table_plot, aes(bad_individu_per_position)) +
        geom_histogram()
      ggsave(g, filename = paste(path_output, "bad_individual_per_position.png", sep = "/"), width = 7, height = 7)
    }

    ###Update reliability of data after cleaning
    total_nbr_var <- length(which(data_filtered$zygosity %in% c("Heterozygous", "Homozygous")))
    reliability <- length(which(data_coverage[,-1] >= min_coverage)) / (dim(data_coverage[,-1])[1]*dim(data_coverage[,-1])[2])
    table_data_reliability$after_cleaning <- c(length(patient_name), length(control_name), nbr_variant, total_nbr_var, reliability)
    write_xlsx(table_data_reliability, paste(path_output, "impact_of_cleaning.xlsx", sep = "/"))

  }else{
    #No cleaning done on data
  }


  ###Split data cases and controls
  patient_data <- data_filtered[which(data_filtered$sample %in% cases$X1),]
  control_data <- data_filtered[which(data_filtered$sample %in% controls$X1),]
  nbr_cases <- length(unique(patient_data$sample))
  nbr_controls <- length(unique(control_data$sample))
  nbr_individu <- nbr_cases + nbr_controls


  ###Keep only variant that are Heterozygous or Homozygous alternative
  data_filtered <- data_filtered[which(data_filtered$zygosity %in% c("Heterozygous", "Homozygous")),]
  patient_data <- patient_data[which(patient_data$zygosity %in% c("Heterozygous", "Homozygous")),]
  control_data <- control_data[which(control_data$zygosity %in% c("Heterozygous", "Homozygous")),]


  ###Save the cleaned data into separated tsv for cases and controls
  write_tsv(patient_data, paste(path_output, "patient_data.tsv", sep = "/"))
  write_tsv(control_data, paste(path_output, "control_data.tsv", sep = "/"))


  ###Output of the function
  return(list(patient_data,
              control_data))


}

