#!/usr/bin/Rscript
#/*****************************************************************************************
#  *
#  * Database mapping- Copyright (C) <2024-2026> <Ecole Polythecnique Federal de Lausanne (EPFL)>
#  * 	
#  * List of the contributors to the development of Database mapping: see LICENSE file.
#* Description and complete License: see LICENSE file.
#* 	
#  * This program (Database mapping) is free software: 
#  * you can redistribute it and/or modify it under the terms of the 
#* GNU General Public License as published by the Free Software Foundation, 
#* either version 3 of the License, or (at your option) any later version.
#* 
#  * This program is distributed in the hope that it will be useful,
#* but WITHOUT ANY WARRANTY; without even the implied warranty of
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#* GNU General Public License for more details.
#* 
#  * You should have received a copy of the GNU General Public License
#* along with this program (see COPYING file).  If not, 
#* see <http://www.gnu.org/licenses/>.
#* 
#  *****************************************************************************************/
#  
#  /**
#  *
#  * @author Simon Boutry
#*
#  */
gene_to_pathway <- function(gene_list = c(),
                            output_path = "",
                            min_inputGene_inPathway = 2,
                            max_inputGene_inPathway = 100,
                            RatioInputGeneVSGeneInPathway = 0.8,
                            path_to_database = ""){
  #Maybe add cutoff minimum RatioInputGeneVSGeneInPathway or adapt DoAggregate and for each significant pathway (group of genes) make over representation (change over representation to use our new database as pathway)
  
  
  ###Create new folder to store results 
  #dir.create(output_path)
  
  ###Load database (Carefull with names of tables)
  database <- read_tsv(paste(path_to_database, "Final_database.tsv", sep = "/"))
  table_geneToANpathway <- read_tsv(paste(path_to_database, "Final_table_geneToANpathway.tsv", sep = "/"))
  
  
  ###Gene not found in database from input gene list
  input_gene_notFound <- gene_list[which(!gene_list %in% table_geneToANpathway$Gene)]
  write_xlsx(data.frame(Gene = input_gene_notFound), paste(output_path, "gene_not_found_in_database.xlsx", sep = "/"))
  gene_list <- gene_list[-which(!gene_list %in% table_geneToANpathway$Gene)]

  
  ###Filter table_geneToANpathway based on gene_list using min and max gene per pathway
  table_geneToANpathway <- table_geneToANpathway[which(table_geneToANpathway$Gene %in% gene_list),]
  #remove too little number of input gene mapping to same pathway
  AN_to_keep <- as.double(names(which(table(table_geneToANpathway$AN) >= min_inputGene_inPathway)))
  table_geneToANpathway <- table_geneToANpathway[which(table_geneToANpathway$AN %in% AN_to_keep),]
  #remove too big number of input gene mapping to same pathway
  AN_to_keep <- as.double(names(which(table(table_geneToANpathway$AN) <= max_inputGene_inPathway)))
  table_geneToANpathway <- table_geneToANpathway[which(table_geneToANpathway$AN %in% AN_to_keep),]
  
  
  ###Filter database based on AN to keep
  database <- database[which(database$AN %in% AN_to_keep),]
  
  
  ###Add annotation based on input gene list
  new_annotation <- c("InputGene", "nbrInputGene", "RatioInputGeneVSGeneInPathway")
  new_data <- data.frame(matrix(NA, nrow = dim(database)[1], ncol = length(new_annotation)))  
  colnames(new_data) <- new_annotation
  new_data <- cbind(new_data, database)
  for (i in 1:dim(new_data)[1]) {
    tmp_gene <- table_geneToANpathway$Gene[which(table_geneToANpathway$AN == new_data$AN[i])]
    new_data$InputGene[i] <- paste(tmp_gene, collapse = "/")
    new_data$nbrInputGene[i] <- length(tmp_gene)
    new_data$RatioInputGeneVSGeneInPathway[i] <- length(tmp_gene) / new_data$Nbr_gene[i]
  }
  
  
  ###Filter pathways based on RatioInputGeneVSGeneInPathway
  to_keep <- which(new_data$RatioInputGeneVSGeneInPathway > RatioInputGeneVSGeneInPathway)
  if(length(to_keep) > 0){
    new_data <- new_data[to_keep,]
    
    
    ###Create pathway ID annotation, same pathway ID for same InputGene and retrieve necessary info for DoAggregate
    pathway <- data.frame(matrix(NA, ncol = 1, nrow = dim(new_data)[1]))
    colnames(pathway) <- "pathway"
    new_data <- cbind(pathway, new_data)
    #duplicated_InputGene <- names(which(table(new_data$InputGene) > 1))
    unique_InputGene <- unique(new_data$InputGene)
    row_to_keep <- c()
    for (i in 1:length(unique_InputGene)) {
      tmp_id <- which(new_data$InputGene == unique_InputGene[i])
      new_data$pathway[tmp_id] <- i
      if(length(tmp_id) > 1){
        tmp_id <- tmp_id[which(new_data$RatioInputGeneVSGeneInPathway[tmp_id] == max(new_data$RatioInputGeneVSGeneInPathway[tmp_id]))]
        tmp_id <- tmp_id[1] #if all have the maximum RatioInputGeneVSGeneInPathway then just take the first one
        row_to_keep <- c(row_to_keep, tmp_id)
      }else{
        row_to_keep <- c(row_to_keep, tmp_id)
      }
    }
    
    
    ###Keep output for aggregation at pathway level
    output_for_DoAggregate <- list(c(new_data$name[row_to_keep]),
                                   c(new_data$InputGene[row_to_keep]))
    
    
    ###pathway ID annotation and order by database based on them and based on RatioInputGeneVSGeneInPathway
    new_data <- new_data[order(new_data$pathway, new_data$RatioInputGeneVSGeneInPathway), ]
    
    
    ###Remove AN from database to save in general folder
    new_data <- new_data[,-which(colnames(new_data) == "AN")]
    write_xlsx(new_data, paste(output_path, "pathway_database.xlsx", sep = "/"))
    
    
    ###output from function
    return(output_for_DoAggregate)
    
  }else{
    #No pathway passes the threshold
    new_data <- NULL
    return(NULL)
  }
  
  
  
  
  
}#end of function