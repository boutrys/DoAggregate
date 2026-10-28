rm(list = ls())
start <- Sys.time()


###Input parameter
path_data <- "data/big_example/inputs/annovar.vcf.hg38_multianno.reduced.vcf.gz"
path_output <- "data/big_example/inputs"
MAF <- 0.01
het_var <- c("0/1", "0|1")
homo_var <- c("1/1", "1|1")
homo_ref <- c("0/0", "0|0")


#what we want
column_to_have <- c("sample", "chr", "pos", "reference", "alternative", "zygosity", "gene_symbol")


###libraries
library(vcfR)
library(tidyverse)


###Read input data
data <- read.vcfR(path_data)

transformed <- vcfR2tidy(data, single_frame = TRUE) 

dim(transformed$dat)

###Pre process data
rm(data)
transformed$dat$gnomad40_AF <- as.double(transformed$dat$gnomad40_AF)
colnames(transformed$dat)[which(colnames(transformed$dat) == "REF...4")] <- "REF"


### Filtering of data
data_filtered <- transformed$dat %>%
  filter(gnomad40_AF <= MAF | is.na(gnomad40_AF))
rm(transformed)
dim(data_filtered)
data_filtered$zygosity <- NA
nbr_individu <- length(unique(data_filtered$Indiv))

data_filtered$zygosity[which(data_filtered$gt_GT %in% het_var)] <- "Heterozygous"
data_filtered$zygosity[which(data_filtered$gt_GT %in% homo_var)] <- "Homozygous"
data_filtered$zygosity[which(data_filtered$gt_GT %in% homo_ref)] <- "Reference"

variant_NA <- unique(data_filtered$ID[which(is.na(data_filtered$zygosity))])

data_wellCovered <- data_filtered[which(!data_filtered$ID %in% variant_NA),]
dim(data_wellCovered)
rm(data_filtered)

colnames(data_wellCovered)[which(colnames(data_wellCovered) == "CHROM")] <- "chr"
colnames(data_wellCovered)[which(colnames(data_wellCovered) == "POS")] <- "pos"
colnames(data_wellCovered)[which(colnames(data_wellCovered) == "REF")] <- "reference"
colnames(data_wellCovered)[which(colnames(data_wellCovered) == "ALT")] <- "alternative"
colnames(data_wellCovered)[which(colnames(data_wellCovered) == "Gene.ensGene")] <- "gene_symbol"
colnames(data_wellCovered)[which(colnames(data_wellCovered) == "Indiv")] <- "sample"
colnames(data_wellCovered)[which(colnames(data_wellCovered) == "gnomad40_AF")] <- "gnomad_wes_af"

geneNameToChange <- grep(pattern = "\\", data_wellCovered$gene_symbol, fixed = TRUE)
list_to_change <- unique(data_wellCovered$gene_symbol[geneNameToChange])

for (i in 1:length(list_to_change)) {
  tmp_id <- which(data_wellCovered$gene_symbol[geneNameToChange] == list_to_change[i])
  tmp_new_name <- strsplit(data_wellCovered$gene_symbol[geneNameToChange[tmp_id[1]]], split = "\\", fixed = TRUE)[[1]][1]
  data_wellCovered$gene_symbol[geneNameToChange[tmp_id]] <- tmp_new_name
}


#remove reference and keep only variants
data_wellCovered <- data_wellCovered[which(data_wellCovered$zygosity %in% c("Heterozygous", "Homozygous")),]
dim(data_wellCovered)


#remove badly covered variant
data_wellCovered <- data_wellCovered[which(data_wellCovered$gt_DP >= 10),]

#Split data cases and controls
sample_list <- unique(data_wellCovered$sample)
cases <- sample_list[grep(pattern = "eNA", sample_list)]
controls <- sample_list[grep(pattern = "eHG", sample_list)]

patient_data <- data_wellCovered[which(data_wellCovered$sample %in% cases),]
control_data <- data_wellCovered[which(data_wellCovered$sample %in% controls),]


write_tsv(patient_data, paste(path_output, "patient_data.tsv", sep = "/"))
write_tsv(control_data, paste(path_output, "control_data.tsv", sep = "/"))
stop <- Sys.time()
tot <- stop - start
print(tot)