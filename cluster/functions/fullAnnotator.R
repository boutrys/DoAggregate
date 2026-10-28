############################################# ANNOTATOR #############################################
# The function allow to fully annotated with VEP the VCF, and annotate on GO/PATHWAYS/DO term
fullAnnotator <- function(data = c(), path_to_table = ""){
  
  go_table <- read.table(paste(path_to_table, "gene_go.txt.gz", sep = ""), stringsAsFactors = FALSE)
  pathway_table <- read.table(paste(path_to_table, "gene_pathway.txt.gz", sep = ""), stringsAsFactors = FALSE)
  
  #VEP annotation tool
  #STIL TO DO
  geneList <- unique(data$gene_symbol)

  ########################### CONSENSUS PATH DB  ###################
  empty_col <- c(matrix(NA, ncol = 1, nrow = dim(data)[1]))
  data$Pathway_list <- empty_col
  data$Nbr_pathway_involved <- empty_col
  gene_not_annotated_pathway <- c()
  pathway_table <- filter(pathway_table, gene_name %in% geneList)
  count_missing_pathway <- 0
  for (i in 1:length(geneList)) {
    good_idx_pathway <- which(pathway_table$gene_name == geneList[i])
    idx <- which(data$gene_symbol == geneList[i])
    if(length(good_idx_pathway) > 0){
      data$Pathway_list[idx] <- as.character(pathway_table$Pathway_list[good_idx_pathway])
      data$Nbr_pathway_involved[idx] <- pathway_table$Nbr_pathway_involved[good_idx_pathway]
    }else{
      gene_not_annotated_pathway <- paste(gene_not_annotated_pathway, geneList[i], sep = ";")
      count_missing_pathway <- count_missing_pathway + length(idx)
    }
  }
  
  
  ########################## GO ANNOTATION  ########################
  tmp_table <- data.frame(matrix(NA, ncol = dim(go_table)[2]-4, nrow = dim(data)[1]))
  gene_not_annotated_GO <- c()
  go_table <- filter(go_table, gene_name %in% geneList)
  go_table[,16:45] <- data.frame(lapply(go_table[,16:45], as.character), stringsAsFactors=FALSE)
  count_missing_GO <- 0
  for (j in 1:length(geneList)) {
    good_idx_goterm <- which(go_table$gene_name == geneList[j])
    idx <- which(data$gene_symbol == geneList[j])
    if(length(good_idx_goterm) > 0){
      tmp <- go_table[good_idx_goterm,-c(46:49)]
      tmp_table[idx,] <- tmp
    }else{
      gene_not_annotated_GO <- paste(gene_not_annotated_GO, geneList[j], sep = ";")
      count_missing_GO <- count_missing_GO + length(idx)
    }
  }
  colnames(tmp_table) <- colnames(go_table)[-c(46:49)]
  data <- bind_cols(data,tmp_table)
  data$GO_list <- apply(data[,grep(pattern = "ID_level_", x = colnames(data))], 1, paste, collapse = ";")
  
  gene_not_annotated_pathway <- substring(gene_not_annotated_pathway, 2)
  gene_not_annotated_GO <- substring(gene_not_annotated_GO, 2)
  missing_annotation <- data.frame(matrix(NA, ncol = 4, nrow = 1))
  colnames(missing_annotation) <- c("Nbr_variant_missing_pathway", "Nbr_variant_missing_GO",
                                    "Missing_pathway", "Missing_GO")
  missing_annotation$Nbr_variant_missing_pathway <- count_missing_pathway
  missing_annotation$Nbr_variant_missing_GO <- count_missing_GO
  if(length(gene_not_annotated_pathway) > 0){
    missing_annotation$Missing_pathway <- gene_not_annotated_pathway
  }else{
    missing_annotation$Missing_pathway <- NA
  }
  if(length(gene_not_annotated_GO) > 0){
    missing_annotation$Missing_GO <- gene_not_annotated_GO
  }else{
    missing_annotation$Missing_GO <- NA
  }

  return(list(data, 
              missing_annotation))
}