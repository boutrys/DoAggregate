###Create a smaller public dataset

###Input parameter
path_data <- "data/big_example/inputs/annovar.vcf.hg38_multianno.reduced.vcf.gz"
path_output <- "data/small_example/small_public_test_data.vcf"
gene_to_keep <- c("AKR1E2", "ADK", "GANC", "GNB2", "GIGYF1", "ZNF845", "FAM41C", "OR4F16", "CELA3B")


###libraries
library(vcfR)
library(tidyverse)


###Read input data
data <- read.vcfR(path_data)


###Data filtering
id_to_keep <- c()
for (i in 1:dim(data@fix)[1]) {
  tmp_id <- grep(pattern = "Gene.ensGene", strsplit(split = ";", data@fix[i,8])$INFO)
  tmp <- strsplit(split = "=", strsplit(split = ";", data@fix[i,8])$INFO[tmp_id])[[1]][2] 
  if(length(grep(pattern = "\\", tmp, fixed = TRUE)) == 1){
    tmp <- strsplit(split = "\\", tmp, fixed = TRUE)
  }
  if(length(which(gene_to_keep %in% tmp)) > 0){
    id_to_keep <- c(id_to_keep, i)
  }
}

data@fix <- data@fix[id_to_keep,]
data@gt <- data@gt[id_to_keep,]

write.vcf(data, file = path_output, mask = FALSE, APPEND = FALSE)
