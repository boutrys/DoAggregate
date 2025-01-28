#' Over Representation Analysis for Gene Lists
#'
#' @name over_representation_analysis
#' @description This function performs over-representation analysis on provided gene lists using various databases
#' like MSigDB, GO, and others. It supports extensive customization for the analysis, including
#' background gene set adjustment, similarity cutoffs, and the ability to handle co-occurrence matrices.
#' The function is designed to identify significant pathways and gene ontology terms that are
#' over-represented in the provided gene list relative to a background set.
#'
#' @param geneList Vector of gene symbols to analyze.
#' @param gene_annot Annotation type of the genes provided (default: "SYMBOL").
#' @param background_gene Background gene set for normalization in the analysis.
#' @param go_of_interest Vector specifying which gene ontology categories to consider ("BP", "MF", "CC").
#' @param max_GO_similarity_accepted Maximum similarity threshold for considering GO terms as distinct.
#' @param cutoff Significance cutoff for enrichment results.
#' @param max_gene_per_item Maximum number of genes allowed per GO/pathway term before filtering.
#' @param max_item_plot Maximum number of terms to display in plots.
#' @param coocurrence Whether to perform and plot gene co-occurrence analysis.
#' @param max_coocurrence Maximum number of top terms to include in co-occurrence analysis.
#' @param MSigDB Boolean to decide whether to use the MSigDB database for analysis.
#' @param path_to_table Path to the database tables needed for the analysis.
#' @param path_to_store_results Directory path where the results and plots will be saved.
#' @return A list of data frames containing the over-represented pathways and GO terms along with their details.
#'
#' @import grDevices
#' @export
#' @importFrom utils globalVariables
utils::globalVariables(c("p.adjust"))

#NB
#Plots for group only show the best pvalue of the group, while all others are also significant. For Go term GO_Plot, group, only the most significant of the group is plot.
#Potential problem not handle yet : some pathway might have name including "/". But we also use "/" to seprate name of pathway within group of pathway with same geneID!
#Output : for the moment over_represented_GO is the term that have been found in one of the GO database at least, and over_represented_pathway are the term found in other database than GO
#coocurrence matrix and plot are by default no performed, be carefull, if gene list is big, it is computationaly intensive
#coocurrence plot only display the top 100 most present genes
over_representation_analysis <- function(geneList = c(),
                                         gene_annot = "SYMBOL",
                                         background_gene = c(),
                                         go_of_interest = c("BP", "MF", "CC"),
                                         max_GO_similarity_accepted = 0.5,
                                         cutoff = 0.05,
                                         max_gene_per_item = 500,
                                         max_item_plot = 20,
                                         coocurrence = TRUE,
                                         max_coocurrence = 100,
                                         MSigDB = TRUE,
                                         path_to_table = "",
                                         path_to_store_results = ""){
  start <- Sys.time()

  dir.create(path_to_store_results)
  options(bitmapType="cairo")


  ##############################################################################################
  ############################################### GLOBAL GENE SET ANALYSIS   ###################
  ##############################################################################################
  if(MSigDB){

    ###Load MSigDB database
    pathway_overrepresentation_table <- read_tsv(paste(path_to_table, "pathway_overRepresentation.tsv", sep = "/"))
    database_pathway <- read_tsv(paste(path_to_table, "Final_database.tsv", sep = "/"))


    ###check if there is more than one gene in input and more than one that map a pathway
    geneList <- unique(geneList)
    if(gene_annot != "SYMBOL"){
      geneList <- bitr(geneList, fromType=gene_annot, toType="SYMBOL", OrgDb="org.Hs.eg.db")$SYMBOL
      if(length(background_gene) > 0){
        background_gene <- bitr(background_gene, fromType=gene_annot, toType="SYMBOL", OrgDb="org.Hs.eg.db")$SYMBOL
      }
    }
    id_gene_not_found <- which(!geneList %in% pathway_overrepresentation_table$gene_name)
    if(length(id_gene_not_found) > 0){
      write_xlsx(data.frame(Gene = geneList[id_gene_not_found]), paste(path_to_store_results, "gene_not_found_in_database.xlsx", sep = "/"))
      geneList <- geneList[-id_gene_not_found]
    }


    if(length(geneList) > 1){
      ###Over representation analysis
      if(length(background_gene) == 0){
        ego <- clusterProfiler::enricher(gene = geneList,
                                         TERM2GENE = pathway_overrepresentation_table,
                                         pAdjustMethod = "BH",
                                         qvalueCutoff = cutoff)
      }else{
        ego <- clusterProfiler::enricher(gene = geneList,
                                         universe = background_gene,
                                         TERM2GENE = pathway_overrepresentation_table,
                                         pAdjustMethod = "BH",
                                         qvalueCutoff = cutoff)
      }


      if(!is.null(ego)){
        #Keep only significant results
        ego_result <- ego@result
        short_ego <- filter(ego_result , ego_result$p.adjust < cutoff)


        #Remove pathways with more than max_gene_per_item included (big pathway filter)
        id_too_big_rm <- which(short_ego$Count > max_gene_per_item)
        if(length(id_too_big_rm) > 0){
          remove_too_big <-  short_ego[id_too_big_rm,]
          to_save_remove_too_big <- remove_too_big[,-c(1:2)]
          colnames(to_save_remove_too_big)[c(1:2,6)] <- c("GeneRatioFromInputGene", "Ratio_TotalGene_TotalGeneDatabase", "InputGene")
          to_save_remove_too_big <- cbind(to_save_remove_too_big, database_pathway[which(database_pathway$name %in% rownames(to_save_remove_too_big)),])
          writexl::write_xlsx(to_save_remove_too_big, path = paste(path_to_store_results, "pathway_removed_too_big.xlsx", sep = "/"))
          short_ego <- short_ego[-id_too_big_rm,]
        }

        #Remove terms with only one gene from the gene list
        short_ego <- short_ego[which(short_ego$Count > 1),]


        #Simplify results
        if(dim(short_ego)[1] > 0){
          #Remove redundant geneID list from results, keep only most significant one
          duplicate_geneID <- duplicated(short_ego$geneID)
          to_rm <- which(duplicate_geneID)
          name_too_long <- FALSE
          if(length(to_rm) > 0){
            removed_duplicate <- short_ego[to_rm,]
            short_ego <- short_ego[which(!duplicate_geneID),]
            to_save_removed_duplicate <- removed_duplicate[,-c(1:2)]
            colnames(to_save_removed_duplicate)[c(1:2,6)] <- c("GeneRatioFromInputGene", "Ratio_TotalGene_TotalGeneDatabase", "InputGene")
            to_save_removed_duplicate <- cbind(to_save_removed_duplicate, database_pathway[which(database_pathway$name %in% rownames(to_save_removed_duplicate)),])
            writexl::write_xlsx(to_save_removed_duplicate, path = paste(path_to_store_results, "pathway_removed_duplicated_gene_list.xlsx", sep = "/"))


            #Create group for pathway terms with same geneID list
            dup_gene_list <- unique(removed_duplicate$geneID)
            nbr_dup_gene_list <- length(dup_gene_list)
            short_ego_name_shortened <- short_ego
            for (n in 1:nbr_dup_gene_list) {
              new_name <- paste(gsub(pattern = "/", replacement = "-", x = short_ego$Description[which(short_ego$geneID == dup_gene_list[n])]),
                                paste(gsub(pattern = "/", replacement = "-", x = removed_duplicate$Description[which(removed_duplicate$geneID == dup_gene_list[n])]), collapse = "-"),
                                sep = "-")
              short_ego$Description[which(short_ego$geneID == dup_gene_list[n])] <- new_name
              short_ego_name_shortened$Description[which(short_ego_name_shortened$geneID == dup_gene_list[n])] <- new_name
              if(nchar(new_name) > 50){
                new_name <- substr(new_name, 1, 50)
                short_ego_name_shortened$Description[which(short_ego_name_shortened$geneID == dup_gene_list[n])] <- new_name
                name_too_long <- TRUE
              }
            }
          }


          #Plotting the results
          ego_simplified <- ego
          ego_simplified@result <- short_ego
          if(name_too_long){
            ego_simplified@result <- short_ego_name_shortened
          }
          current_item_plot <- dim(ego_simplified@result)[1]
          if(current_item_plot > max_item_plot){
            current_item_plot <- max_item_plot
          }
          Cairo(paste(path_to_store_results, "Dot_plot_pathway_overrepresentation.png", sep = "/"),
                width = 1000,
                height = 1000)
          print(enrichplot::dotplot(ego_simplified, showCategory = current_item_plot))
          dev.off()

          tmp_width <- 600 + length(unique(strsplit(paste(short_ego$geneID, collapse = "/"), split = "/")[[1]])) * 10
          tmp_height <- 600 + current_item_plot * 10
          if(tmp_width > 2000){
            tmp_width <- 2000
          }
          if(tmp_height > 2000){
            tmp_height <- 2000
          }
          Cairo(paste(path_to_store_results, "Plot_Gene_occurence_pathway.png", sep = "/"),
                width = tmp_width,
                height = tmp_height)
          print(enrichplot::heatplot(ego_simplified, showCategory = current_item_plot))
          dev.off()

          Cairo(paste(path_to_store_results, "Network_overrepresentation.png", sep = "/"),
                width = 1000,
                height = 1000)
          p <- cnetplot(ego_simplified,
                        showCategory = current_item_plot,
                        categorySize="pvalue",
                        colorEdge = TRUE,
                        color_category='firebrick',
                        color_gene='steelblue')
          print(p)
          dev.off()

          Cairo(paste(path_to_store_results, "Network_Circle_overrepresentation.png", sep = "/"),
                width = 1000,
                height = 1000)
          p <- cnetplot(ego_simplified,
                        showCategory = current_item_plot,
                        categorySize="pvalue",
                        circular = TRUE,
                        colorEdge = TRUE,
                        color_category='firebrick',
                        color_gene='steelblue')
          print(p)
          dev.off()


          #TODO combine info of our database and remove useless columns from short_ego to save interesting info
          #transform results if group are forme
          short_ego$NbrPathwayCollapsed <- 0
          to_check <- grep(pattern = "-", short_ego$Description)
          for (ww in to_check) {
            short_ego$NbrPathwayCollapsed[ww] <- length(strsplit(split = "-", short_ego$Description[ww])[[1]])
          }

          #Save results
          to_save_short_ego <- short_ego %>%
            left_join(database_pathway, by=c("ID"="name"))
          to_save_short_ego <- to_save_short_ego[,-c(4,9)]
          to_save_short_ego <- to_save_short_ego[,c(1,8,2:7,18,19,9:14,16,17,15)]
          colnames(to_save_short_ego)[c(1, 3, 9)] <- c("name_pathway", "collapsed_pathway", "Total_Nbr_gene")
          writexl::write_xlsx(to_save_short_ego, path = paste(path_to_store_results, "Summary_results_over_representation.xlsx", sep = "/"))

        }#IF dim(short_ego)[1] > 0


        #Gene never over represented in original input
        id_rm <- which(!geneList %in% unique(strsplit(paste(short_ego$geneID, collapse = "/"), split = "/")[[1]]))
        if(length(id_rm) > 0){
          name_no_OR <- geneList[id_rm]
          writexl::write_xlsx(data.frame(Gene = name_no_OR), path = paste(path_to_store_results, "Gene_not_over_represented.xlsx", sep = "/"))
        }


        if(dim(short_ego)[1] > 0){
          if(coocurrence){
            start_matrix <- Sys.time()
            print("Starting gene co-occurence analysis after over-representation analysis")
            #restrict analysis on top most significant GO/pathway
            if(dim(short_ego)[1] > max_coocurrence){
              short_ego <- short_ego[c(1:max_coocurrence),]
            }


            #By gene based on input geneList
            nbr_gene <- length(geneList)
            geneBygene_matrix <- matrix(0, ncol = nbr_gene, nrow = nbr_gene)
            colnames(geneBygene_matrix) <- geneList
            rownames(geneBygene_matrix) <- geneList
            for (i in 1:dim(short_ego)[1]) {
              tmp_list_gene <- strsplit(short_ego$geneID[i], split = "/")[[1]]
              for (j in 1:length(tmp_list_gene)) {
                tmp_gene <- tmp_list_gene[j]
                geneBygene_matrix[which(rownames(geneBygene_matrix) == tmp_gene),which(colnames(geneBygene_matrix) %in% tmp_list_gene)] <- geneBygene_matrix[which(rownames(geneBygene_matrix) == tmp_gene),which(colnames(geneBygene_matrix) %in% tmp_list_gene)] + 1
              }
            }
            id_rm <- which(colSums(geneBygene_matrix) == 0)
            if(length(id_rm) > 0){
              geneBygene_matrix <- geneBygene_matrix[-id_rm, -id_rm]
            }
            geneBygene_data <- data.frame(Gene = rownames(geneBygene_matrix))
            geneBygene_data <- cbind(geneBygene_data, as.data.frame(geneBygene_matrix))
            writexl::write_xlsx(geneBygene_data, path = paste(path_to_store_results, "Cooccurence_genes.xlsx", sep = "/"))
            end_matrix <- Sys.time()


            start_plot <- Sys.time()
            if(dim(geneBygene_matrix)[2] > 100){
              max_id <- colSums(geneBygene_matrix)
              condition_plot <- min(max_id[order(max_id, decreasing = TRUE)][1:100])
              to_keep_plot <- which(max_id >= condition_plot)[1:100]
              geneBygene_matrix_reduce <- geneBygene_matrix[to_keep_plot,to_keep_plot]
            }else{
              geneBygene_matrix_reduce <- geneBygene_matrix
            }

            #ploting coocurrence of genes
            name_of_plot <- paste(path_to_store_results, "Coocurrence_genes.png", sep = "/")
            Cairo(name_of_plot,
                  width = 1500, height = 1200)

            h1 <- Heatmap(geneBygene_matrix_reduce,
                          heatmap_legend_param = list(
                            title = "CoOcurrence"),
                          column_km = 1,
                          cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                            if(geneBygene_matrix_reduce[i, j] >= 1)
                              grid.text(sprintf("%d", geneBygene_matrix_reduce[i, j]), x, y, gp = gpar(fontsize = 10))
                          }
            )
            draw(h1)
            dev.off()
            print("End gene co-occurence analysis after over-representation analysis")
            end_plot <- Sys.time()
          }#enf of If coocurrence
        }


        #Returning results
        Description <- c(short_ego$ID)
        geneID <- c(short_ego$geneID)
        over_represented_pathway <- list(Description, geneID)

        stop <- Sys.time()
        print(stop - start)
        return(over_represented_pathway)


        ###End of if significantly over represented pathways
      }else{
        return(list(NA, NA))
      }
    }else{
      print("less than two genes found in pathway database so no over representation performed")
      return(list(NA, NA))
    }






    ###End of if using MSigDB pathways

  }else{

    ###Use clusterProfiler accessible databases

    print("Start to map gene symbol to ENTREZID")
    ids_geneList <- bitr(geneList, fromType=gene_annot, toType="ENTREZID", OrgDb="org.Hs.eg.db")
    if(length(background_gene) > 0){
      background_gene <- bitr(background_gene, fromType=gene_annot, toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
    }
    print("End to map gene symbol to ENTREZID")

    if(dim(ids_geneList)[1] > 1){
      if(length(unique(ids_geneList[,1])) < length(unique(geneList))){
        gene_not_handled <- setdiff(geneList, unique(ids_geneList[,1]))
        writexl::write_xlsx(data.frame(gene_not_handled), path = paste(path_to_store_results, "list_gene_not_handled_for_analysis.xlsx", sep = "/"))
        geneList <- geneList[-which(geneList %in% gene_not_handled)]
      }

      list_ego <- list()

      #GO OVER-REPRESENTATION TEST (on unique gene list : see below to use original dataframe)
      path_GO_over_representation <- paste(path_to_store_results, "GO_over_representation", sep = "/")
      dir.create(path_GO_over_representation)
      short_ego_track <- list()
      for (j in 1:length(go_of_interest)) {
        current_folder <- paste(path_GO_over_representation, go_of_interest[j], sep = "/")
        current_folder <- paste(current_folder, "/", sep = "")
        if(length(background_gene) == 0){
          ego <- clusterProfiler::enrichGO(gene = base::unique(ids_geneList[,2]),
                                           OrgDb = org.Hs.eg.db,
                                           ont = go_of_interest[j],
                                           pAdjustMethod = "BH",
                                           qvalueCutoff = cutoff,
                                           readable = TRUE)
        }else{
          ego <- clusterProfiler::enrichGO(gene = base::unique(ids_geneList[,2]),
                                           OrgDb = org.Hs.eg.db,
                                           universe = background_gene,
                                           ont = go_of_interest[j],
                                           pAdjustMethod = "BH",
                                           qvalueCutoff = cutoff,
                                           readable = TRUE)
        }
        list_ego[[j]] <- ego
        if(!is.null(ego)){
          ego_result <- ego@result
          dir.create(current_folder)

          #Keep only significant results
          short_ego <- filter(ego_result , ego_result$p.adjust < cutoff)

          #Remove GO term with more than max_gene_per_item included (big GO term filter)
          id_too_big_rm <- which(short_ego$Count > max_gene_per_item)
          if(length(id_too_big_rm) > 0){
            remove_too_big <-  short_ego[id_too_big_rm,]
            writexl::write_xlsx(remove_too_big, path = paste(current_folder, "/", go_of_interest[j], "_GO_removed_too_big.xlsx", sep = ""))
            short_ego <- short_ego[-id_too_big_rm,]
          }

          #Remove terms with only one gene from the gene list
          short_ego <- short_ego[which(short_ego$Count > 1),]


          #Simplify results
          if(dim(short_ego)[1] > 0){
            #Simplify the GO terms, by removing redundant terms based on their semantic similarities
            hsGO <- GOSemSim::godata('org.Hs.eg.db', ont=go_of_interest[j])

            similarity_matrix  <- mgoSim(short_ego$ID,
                                         short_ego$ID,
                                         semData = hsGO,
                                         combine = NULL)
            for (l in 1:dim(similarity_matrix)[1]) {
              for (k in 1:dim(similarity_matrix)[2]) {
                if(l == k){
                  similarity_matrix[l,k] <- -1
                }
              }
            }
            over_similarity <- base::colSums(similarity_matrix > max_GO_similarity_accepted)
            tmp_over_similarity <- over_similarity
            tmp_similarity_matrix <- similarity_matrix
            current_max_sim <- max(tmp_similarity_matrix)
            idx_to_remove <- c()
            nbr_loop <- 0
            while (current_max_sim > max_GO_similarity_accepted) {
              nbr_loop <- nbr_loop + 1
              current_idx <- which(tmp_over_similarity == max(tmp_over_similarity))
              if(length(current_idx) > 1){
                current_idx <- current_idx[length(current_idx)]
              }
              tmp_idx <- which(short_ego$ID == colnames(tmp_similarity_matrix)[current_idx])
              tmp_similarity_matrix <- tmp_similarity_matrix[-current_idx,-current_idx]
              if(is.null(dim(tmp_similarity_matrix))){
                tmp_over_similarity <- 0
              }else{
                tmp_over_similarity <- base::colSums(tmp_similarity_matrix > max_GO_similarity_accepted)
              }
              idx_to_remove <- c(idx_to_remove, tmp_idx)
              current_max_sim <- max(tmp_similarity_matrix)
            }
            if(length(idx_to_remove) > 0){
              removed_from_short_ego <- short_ego[idx_to_remove,]
              short_ego <- short_ego[-idx_to_remove,]
              writexl::write_xlsx(removed_from_short_ego, path = paste(current_folder, "/", go_of_interest[j], "_GO_removed_similar_representation.xlsx", sep = ""))
            }


            #Remove redundant geneID list from results, keep only most significant one
            duplicate_geneID <- duplicated(short_ego$geneID)
            to_rm <- which(duplicate_geneID)
            name_too_long <- FALSE
            if(length(to_rm) > 0){
              removed_duplicate <- short_ego[to_rm,]
              short_ego <- short_ego[which(!duplicate_geneID),]
              writexl::write_xlsx(removed_duplicate, path = paste(current_folder, "/", go_of_interest[j], "_GO_removed_duplicated_gene_list.xlsx", sep = ""))

              #Create group for GO terms with same geneID list
              dup_gene_list <- unique(removed_duplicate$geneID)
              nbr_dup_gene_list <- length(dup_gene_list)
              short_ego_name_shortened <- short_ego
              for (n in 1:nbr_dup_gene_list) {
                new_name <- paste(gsub(pattern = "/", replacement = "-", x = short_ego$Description[which(short_ego$geneID == dup_gene_list[n])]),
                                  paste(gsub(pattern = "/", replacement = "-", x = removed_duplicate$Description[which(removed_duplicate$geneID == dup_gene_list[n])]), collapse = "_"),
                                  sep = "_")
                short_ego$Description[which(short_ego$geneID == dup_gene_list[n])] <- new_name
                short_ego_name_shortened$Description[which(short_ego_name_shortened$geneID == dup_gene_list[n])] <- new_name
                if(nchar(new_name) > 50){
                  new_name <- substr(new_name, 1, 50)
                  short_ego_name_shortened$Description[which(short_ego_name_shortened$geneID == dup_gene_list[n])] <- new_name
                  name_too_long <- TRUE
                }
              }
            }


            #Plotting the results
            ego_simplified <- ego
            ego_simplified@result <- short_ego
            current_item_plot <- dim(ego_simplified@result)[1]
            if(current_item_plot > max_item_plot){
              current_item_plot <- max_item_plot
            }
            #TODO why bug when only one term to plot
            if(dim(short_ego)[1] > 1){
              Cairo(paste(current_folder, "GO_Plot_", go_of_interest[j], "_GO_overrepresentation.png", sep = ""),
                    width = 1000,
                    height = 1000)
              print(enrichplot::goplot(ego_simplified, showCategory = current_item_plot))
              dev.off()
            }

            if(name_too_long){
              ego_simplified@result <- short_ego_name_shortened
            }

            Cairo(paste(current_folder, "Dot_plot_", go_of_interest[j], "_GO_overrepresentation.png", sep = ""),
                  width = 1000,
                  height = 1000)
            print(enrichplot::dotplot(ego_simplified, showCategory = current_item_plot))
            dev.off()

            tmp_width <- 600 + length(unique(strsplit(paste(short_ego$geneID, collapse = "/"), split = "/")[[1]])) * 10
            tmp_height <- 600 + current_item_plot * 10
            if(tmp_width > 2000){
              tmp_width <- 2000
            }
            if(tmp_height > 2000){
              tmp_height <- 2000
            }
            Cairo(paste(current_folder, "Plot_Gene_occurence_GO_", go_of_interest[j], ".png", sep = ""),
                  width = tmp_width,
                  height = tmp_height)
            print(enrichplot::heatplot(ego_simplified, showCategory = current_item_plot))
            dev.off()


            if(length(to_rm) > 0){
              #Create group for GO terms with same geneID list
              for (n in 1:nbr_dup_gene_list) {
                new_name <- paste(gsub(pattern = "/", replacement = "-", x = short_ego$ID[which(short_ego$geneID == dup_gene_list[n])]),
                                  paste(gsub(pattern = "/", replacement = "-", x = removed_duplicate$ID[which(removed_duplicate$geneID == dup_gene_list[n])]), collapse = "_"),
                                  sep = "_")
                short_ego$ID[which(short_ego$geneID == dup_gene_list[n])] <- new_name
              }
            }

            #Add database information
            short_ego$Database <- paste("GO", go_of_interest[j], sep = "-")

            #Save results
            writexl::write_xlsx(short_ego, path = paste(current_folder, "/", go_of_interest[j], "_GO_over_representation.xlsx", sep = ""))

          }#IF dim(short_ego)[1] > 0

          short_ego_track[[j]] <- short_ego
        }else{
          short_ego_track[[j]] <- NULL
        }

      }#For loop for each BP, CC, MF GO terms


      ###Others Databas for pathway analysis
      pathway_analysis <- paste(path_to_store_results, "pathway_over_representation", sep = "/")
      dir.create(pathway_analysis)
      short_pathway_track <- list()


      ###KEGG analysis
      KEGG_path <- paste(pathway_analysis, "KEGG", sep = "/")
      if(length(background_gene) == 0){
        ego <- clusterProfiler::enrichKEGG(gene = base::unique(ids_geneList[,2]),
                                           organism = "hsa",
                                           keyType = "kegg",
                                           pAdjustMethod = "BH",
                                           qvalueCutoff = cutoff,
                                           use_internal_data = FALSE)
      }else{
        ego <- clusterProfiler::enrichKEGG(gene = base::unique(ids_geneList[,2]),
                                           organism = "hsa",
                                           keyType = "kegg",
                                           universe = background_gene,
                                           pAdjustMethod = "BH",
                                           qvalueCutoff = cutoff,
                                           use_internal_data = FALSE)
      }

      list_ego[[4]] <- ego
      if(!is.null(ego)){
        ego <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID')
        ego_result <- ego@result
        dir.create(KEGG_path)

        #Keep only significant results
        short_ego <- filter(ego_result , ego_result$p.adjust < cutoff)


        #Remove pathways with more than max_gene_per_item included (big pathway filter)
        id_too_big_rm <- which(short_ego$Count > max_gene_per_item)
        if(length(id_too_big_rm) > 0){
          remove_too_big <-  short_ego[id_too_big_rm,]
          writexl::write_xlsx(remove_too_big, path = paste(KEGG_path, "KEGG_removed_too_big.xlsx", sep = "/"))
          short_ego <- short_ego[-id_too_big_rm,]
        }

        #Remove terms with only one gene from the gene list
        short_ego <- short_ego[which(short_ego$Count > 1),]


        #Simplify results
        if(dim(short_ego)[1] > 0){
          #Remove redundant geneID list from results, keep only most significant one
          duplicate_geneID <- duplicated(short_ego$geneID)
          to_rm <- which(duplicate_geneID)
          name_too_long <- FALSE
          if(length(to_rm) > 0){
            removed_duplicate <- short_ego[to_rm,]
            short_ego <- short_ego[which(!duplicate_geneID),]
            writexl::write_xlsx(removed_duplicate, path = paste(KEGG_path, "KEGG_removed_duplicated_gene_list.xlsx", sep = "/"))


            #Create group for pathway terms with same geneID list
            dup_gene_list <- unique(removed_duplicate$geneID)
            nbr_dup_gene_list <- length(dup_gene_list)
            short_ego_name_shortened <- short_ego
            for (n in 1:nbr_dup_gene_list) {
              new_name <- paste(gsub(pattern = "/", replacement = "-", x = short_ego$Description[which(short_ego$geneID == dup_gene_list[n])]),
                                paste(gsub(pattern = "/", replacement = "-", x = removed_duplicate$Description[which(removed_duplicate$geneID == dup_gene_list[n])]), collapse = "_"),
                                sep = "_")
              short_ego$Description[which(short_ego$geneID == dup_gene_list[n])] <- new_name
              short_ego_name_shortened$Description[which(short_ego_name_shortened$geneID == dup_gene_list[n])] <- new_name
              if(nchar(new_name) > 50){
                new_name <- substr(new_name, 1, 50)
                short_ego_name_shortened$Description[which(short_ego_name_shortened$geneID == dup_gene_list[n])] <- new_name
                name_too_long <- TRUE
              }
            }
          }


          #Plotting the results
          ego_simplified <- ego
          ego_simplified@result <- short_ego
          if(name_too_long){
            ego_simplified@result <- short_ego_name_shortened
          }
          current_item_plot <- dim(ego_simplified@result)[1]
          if(current_item_plot > max_item_plot){
            current_item_plot <- max_item_plot
          }
          Cairo(paste(KEGG_path, "Dot_plot_KEGG_overrepresentation.png", sep = "/"),
                width = 1000,
                height = 1000)
          print(enrichplot::dotplot(ego_simplified, showCategory = current_item_plot))
          dev.off()

          tmp_width <- 600 + length(unique(strsplit(paste(short_ego$geneID, collapse = "/"), split = "/")[[1]])) * 10
          tmp_height <- 600 + current_item_plot * 10
          if(tmp_width > 2000){
            tmp_width <- 2000
          }
          if(tmp_height > 2000){
            tmp_height <- 2000
          }
          Cairo(paste(KEGG_path, "Plot_Gene_occurence_KEGG_pathway.png", sep = "/"),
                width = tmp_width,
                height = tmp_height)
          print(enrichplot::heatplot(ego_simplified, showCategory = current_item_plot))
          dev.off()

          Cairo(paste(KEGG_path, "Network_KEGG_overrepresentation.png", sep = "/"),
                width = 1000,
                height = 1000)
          p <- cnetplot(ego_simplified,
                        showCategory = current_item_plot,
                        categorySize="pvalue",
                        colorEdge = TRUE,
                        color_category='firebrick',
                        color_gene='steelblue')
          print(p)
          dev.off()

          Cairo(paste(KEGG_path, "Network_Circle_KEGG_overrepresentation.png", sep = "/"),
                width = 1000,
                height = 1000)
          p <- cnetplot(ego_simplified,
                        showCategory = current_item_plot,
                        categorySize="pvalue",
                        circular = TRUE,
                        colorEdge = TRUE,
                        color_category='firebrick',
                        color_gene='steelblue')
          print(p)
          dev.off()


          if(length(to_rm) > 0){
            #Create group for GO terms with same geneID list
            for (n in 1:nbr_dup_gene_list) {
              new_name <- paste(gsub(pattern = "/", replacement = "-", x = short_ego$ID[which(short_ego$geneID == dup_gene_list[n])]),
                                paste(gsub(pattern = "/", replacement = "-", x = removed_duplicate$ID[which(removed_duplicate$geneID == dup_gene_list[n])]), collapse = "_"),
                                sep = "_")
              short_ego$ID[which(short_ego$geneID == dup_gene_list[n])] <- new_name
            }
          }


          #transform results if group are formed
          short_ego$Database <- "KEGG"

          #Save results
          writexl::write_xlsx(short_ego, path = paste(KEGG_path, "KEGG_over_representation.xlsx", sep = "/"))

        }#IF dim(short_ego)[1] > 0

        short_pathway_track[[1]] <- short_ego
      }else{
        short_pathway_track[[1]] <- NULL
      }


      ###WikiPathways analysis
      WikiPathways_path <- paste(pathway_analysis, "WikiPathways", sep = "/")
      if(length(background_gene) == 0){
        ego <- clusterProfiler::enrichWP(gene = base::unique(ids_geneList[,2]),
                                         organism = "Homo sapiens",
                                         pAdjustMethod = "BH",
                                         qvalueCutoff = cutoff)
      }else{
        ego <- clusterProfiler::enrichWP(gene = base::unique(ids_geneList[,2]),
                                         organism = "Homo sapiens",
                                         universe = background_gene,
                                         pAdjustMethod = "BH",
                                         qvalueCutoff = cutoff)
      }


      list_ego[[5]] <- ego
      if(!is.null(ego)){
        ego <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID')
        ego_result <- ego@result
        dir.create(WikiPathways_path)

        #Keep only significant results
        short_ego <- filter(ego_result , ego_result$p.adjust < cutoff)


        #Remove pathways with more than max_gene_per_item included (big pathway filter)
        id_too_big_rm <- which(short_ego$Count > max_gene_per_item)
        if(length(id_too_big_rm) > 0){
          remove_too_big <-  short_ego[id_too_big_rm,]
          writexl::write_xlsx(remove_too_big, path = paste(WikiPathways_path, "WikiPathways_removed_too_big.xlsx", sep = "/"))
          short_ego <- short_ego[-id_too_big_rm,]
        }

        #Remove terms with only one gene from the gene list
        short_ego <- short_ego[which(short_ego$Count > 1),]

        #Simplify results
        if(dim(short_ego)[1] > 0 && length(which(short_ego$Count > 1)) > 0){
          #Remove redundant geneID list from results, keep only most significant one
          duplicate_geneID <- duplicated(short_ego$geneID)
          to_rm <- which(duplicate_geneID)
          name_too_long <- FALSE
          if(length(to_rm) > 0){
            removed_duplicate <- short_ego[to_rm,]
            short_ego <- short_ego[which(!duplicate_geneID),]
            writexl::write_xlsx(removed_duplicate, path = paste(WikiPathways_path, "WikiPathways_removed_duplicated_gene_list.xlsx", sep = "/"))

            #Create group for pathway terms with same geneID list
            dup_gene_list <- unique(removed_duplicate$geneID)
            nbr_dup_gene_list <- length(dup_gene_list)
            short_ego_name_shortened <- short_ego
            for (n in 1:nbr_dup_gene_list) {
              new_name <- paste(gsub(pattern = "/", replacement = "-", x = short_ego$Description[which(short_ego$geneID == dup_gene_list[n])]),
                                paste(gsub(pattern = "/", replacement = "-", x = removed_duplicate$Description[which(removed_duplicate$geneID == dup_gene_list[n])]), collapse = "_"),
                                sep = "_")
              short_ego$Description[which(short_ego$geneID == dup_gene_list[n])] <- new_name
              short_ego_name_shortened$Description[which(short_ego_name_shortened$geneID == dup_gene_list[n])] <- new_name
              if(nchar(new_name) > 50){
                new_name <- substr(new_name, 1, 50)
                short_ego_name_shortened$Description[which(short_ego_name_shortened$geneID == dup_gene_list[n])] <- new_name
                name_too_long <- TRUE
              }
            }
          }


          #Plotting the results
          ego_simplified <- ego
          ego_simplified@result <- short_ego
          if(name_too_long){
            ego_simplified@result <- short_ego_name_shortened
          }
          current_item_plot <- dim(ego_simplified@result)[1]
          if(current_item_plot > max_item_plot){
            current_item_plot <- max_item_plot
          }
          Cairo(paste(WikiPathways_path, "Dot_plot_WikiPathways_overrepresentation.png", sep = "/"),
                width = 1000,
                height = 1000)
          print(enrichplot::dotplot(ego_simplified, showCategory = current_item_plot))
          dev.off()

          tmp_width <- 600 + length(unique(strsplit(paste(short_ego$geneID, collapse = "/"), split = "/")[[1]])) * 10
          tmp_height <- 600 + current_item_plot * 10
          if(tmp_width > 2000){
            tmp_width <- 2000
          }
          if(tmp_height > 2000){
            tmp_height <- 2000
          }
          Cairo(paste(WikiPathways_path, "Plot_Gene_occurence_WikiPathways_pathway.png", sep = "/"),
                width = tmp_width,
                height = tmp_height)
          print(enrichplot::heatplot(ego_simplified, showCategory = current_item_plot))
          dev.off()

          Cairo(paste(WikiPathways_path, "Network_WikiPathways_overrepresentation.png", sep = "/"),
                width = 1000,
                height = 1000)
          p <- cnetplot(ego_simplified,
                        showCategory = current_item_plot,
                        categorySize="pvalue",
                        colorEdge = TRUE,
                        color_category='firebrick',
                        color_gene='steelblue')
          print(p)
          dev.off()

          Cairo(paste(WikiPathways_path, "Network_Circle_WikiPathways_overrepresentation.png", sep = "/"),
                width = 1000,
                height = 1000)
          p <- cnetplot(ego_simplified,
                        showCategory = current_item_plot,
                        categorySize="pvalue",
                        circular = TRUE,
                        colorEdge = TRUE,
                        color_category='firebrick',
                        color_gene='steelblue')
          print(p)
          dev.off()


          if(length(to_rm) > 0){
            #Create group for terms with same geneID list
            for (n in 1:nbr_dup_gene_list) {
              new_name <- paste(gsub(pattern = "/", replacement = "-", x = short_ego$ID[which(short_ego$geneID == dup_gene_list[n])]),
                                paste(gsub(pattern = "/", replacement = "-", x = removed_duplicate$ID[which(removed_duplicate$geneID == dup_gene_list[n])]), collapse = "_"),
                                sep = "_")
              short_ego$ID[which(short_ego$geneID == dup_gene_list[n])] <- new_name
            }
          }

          short_ego$Database <- "WikiPathways"


          #Save results
          writexl::write_xlsx(short_ego, path = paste(WikiPathways_path, "WikiPathways_over_representation.xlsx", sep = "/"))

        }#IF dim(short_ego)[1] > 0

        short_pathway_track[[2]] <- short_ego
      }else{
        short_pathway_track[[2]] <- NULL
      }


      ###Reactome analysis
      Reactome_path <- paste(pathway_analysis, "Reactome", sep = "/")
      if(length(background_gene) == 0){
        ego <- ReactomePA::enrichPathway(gene = base::unique(ids_geneList[,2]),
                                         organism = "human",
                                         pAdjustMethod = "BH",
                                         qvalueCutoff = cutoff,
                                         readable = TRUE)
      }else{
        ego <- ReactomePA::enrichPathway(gene = base::unique(ids_geneList[,2]),
                                         organism = "human",
                                         universe = background_gene,
                                         pAdjustMethod = "BH",
                                         qvalueCutoff = cutoff,
                                         readable = TRUE)
      }

      list_ego[[6]] <- ego
      if(!is.null(ego)){
        ego_result <- ego@result
        dir.create(Reactome_path)


        #Keep only significant results
        short_ego <- filter(ego_result , ego_result$p.adjust < cutoff)


        #Remove pathways with more than max_gene_per_item included (big pathway filter)
        id_too_big_rm <- which(short_ego$Count > max_gene_per_item)
        if(length(id_too_big_rm) > 0){
          remove_too_big <-  short_ego[id_too_big_rm,]
          writexl::write_xlsx(remove_too_big, path = paste(Reactome_path, "Reactome_removed_too_big.xlsx", sep = "/"))
          short_ego <- short_ego[-id_too_big_rm,]
        }

        #Remove terms with only one gene from the gene list
        short_ego <- short_ego[which(short_ego$Count > 1),]


        #Simplify results
        if(dim(short_ego)[1] > 0 && length(which(short_ego$Count > 1)) > 0){
          #Remove redundant geneID list from results, keep only most significant one
          duplicate_geneID <- duplicated(short_ego$geneID)
          to_rm <- which(duplicate_geneID)
          name_too_long <- FALSE
          if(length(to_rm) > 0){
            removed_duplicate <- short_ego[to_rm,]
            short_ego <- short_ego[which(!duplicate_geneID),]
            writexl::write_xlsx(removed_duplicate, path = paste(Reactome_path, "Reactome_removed_duplicated_gene_list.xlsx", sep = "/"))

            #Create group for pathway terms with same geneID list
            dup_gene_list <- unique(removed_duplicate$geneID)
            nbr_dup_gene_list <- length(dup_gene_list)
            short_ego_name_shortened <- short_ego
            for (n in 1:nbr_dup_gene_list) {
              new_name <- paste(gsub(pattern = "/", replacement = "-", x = short_ego$Description[which(short_ego$geneID == dup_gene_list[n])]),
                                paste(gsub(pattern = "/", replacement = "-", x = removed_duplicate$Description[which(removed_duplicate$geneID == dup_gene_list[n])]), collapse = "_"),
                                sep = "_")
              short_ego$Description[which(short_ego$geneID == dup_gene_list[n])] <- new_name
              short_ego_name_shortened$Description[which(short_ego_name_shortened$geneID == dup_gene_list[n])] <- new_name
              if(nchar(new_name) > 50){
                new_name <- substr(new_name, 1, 50)
                short_ego_name_shortened$Description[which(short_ego_name_shortened$geneID == dup_gene_list[n])] <- new_name
                name_too_long <- TRUE
              }
            }
          }


          #Plotting the results
          ego_simplified <- ego
          ego_simplified@result <- short_ego
          if(name_too_long){
            ego_simplified@result <- short_ego_name_shortened
          }
          current_item_plot <- dim(ego_simplified@result)[1]
          if(current_item_plot > max_item_plot){
            current_item_plot <- max_item_plot
          }
          Cairo(paste(Reactome_path, "Dot_plot_Reactome_overrepresentation.png", sep = "/"),
                width = 1000,
                height = 1000)
          print(enrichplot::dotplot(ego_simplified, showCategory = current_item_plot))
          dev.off()

          tmp_width <- 600 + length(unique(strsplit(paste(short_ego$geneID, collapse = "/"), split = "/")[[1]])) * 10
          tmp_height <- 600 + current_item_plot * 10
          if(tmp_width > 2000){
            tmp_width <- 2000
          }
          if(tmp_height > 2000){
            tmp_height <- 2000
          }
          Cairo(paste(Reactome_path, "Plot_Gene_occurence_Reactome_pathway.png", sep = "/"),
                width = tmp_width,
                height = tmp_height)
          print(enrichplot::heatplot(ego_simplified, showCategory = current_item_plot))
          dev.off()


          Cairo(paste(Reactome_path, "Network_Reactome_overrepresentation.png", sep = "/"),
                width = 1000,
                height = 1000)
          p <- cnetplot(ego_simplified,
                        showCategory = current_item_plot,
                        categorySize="pvalue",
                        colorEdge = TRUE,
                        color_category='firebrick',
                        color_gene='steelblue')
          print(p)
          dev.off()

          Cairo(paste(Reactome_path, "Network_Circle_Reactome_overrepresentation.png", sep = "/"),
                width = 1000,
                height = 1000)
          p <- cnetplot(ego_simplified,
                        showCategory = current_item_plot,
                        categorySize="pvalue",
                        circular = TRUE,
                        colorEdge = TRUE,
                        color_category='firebrick',
                        color_gene='steelblue')
          print(p)
          dev.off()



          if(length(to_rm) > 0){
            #Create group for terms with same geneID list
            for (n in 1:nbr_dup_gene_list) {
              new_name <- paste(gsub(pattern = "/", replacement = "-", x = short_ego$ID[which(short_ego$geneID == dup_gene_list[n])]),
                                paste(gsub(pattern = "/", replacement = "-", x = removed_duplicate$ID[which(removed_duplicate$geneID == dup_gene_list[n])]), collapse = "_"),
                                sep = "_")
              short_ego$ID[which(short_ego$geneID == dup_gene_list[n])] <- new_name
            }
          }
          short_ego$Database <- "Reactome"


          #Save results
          writexl::write_xlsx(short_ego, path = paste(Reactome_path, "Reactome_over_representation.xlsx", sep = "/"))

        }#IF dim(short_ego)[1] > 0

        short_pathway_track[[3]] <- short_ego
      }else{
        short_pathway_track[[3]] <- NULL
      }


      ###DiseaseOntology analysis
      DiseaseOntology_path <- paste(pathway_analysis, "DiseaseOntology", sep = "/")
      if(length(background_gene) == 0){
        ego <- DOSE::enrichDO(gene = base::unique(ids_geneList[,2]),
                              ont = "DO",
                              pAdjustMethod = "BH",
                              qvalueCutoff = cutoff,
                              readable = TRUE)
      }else{
        ego <- DOSE::enrichDO(gene = base::unique(ids_geneList[,2]),
                              ont = "DO",
                              universe = background_gene,
                              pAdjustMethod = "BH",
                              qvalueCutoff = cutoff,
                              readable = TRUE)
      }
      list_ego[[7]] <- ego
      if(!is.null(ego)){
        ego_result <- ego@result
        dir.create(DiseaseOntology_path)


        #Keep only significant results
        short_ego <- filter(ego_result , ego_result$p.adjust < cutoff)


        #Remove pathways with more than max_gene_per_item included (big pathway filter)
        id_too_big_rm <- which(short_ego$Count > max_gene_per_item)
        if(length(id_too_big_rm) > 0){
          remove_too_big <-  short_ego[id_too_big_rm,]
          writexl::write_xlsx(remove_too_big, path = paste(DiseaseOntology_path, "DiseaseOntology_removed_too_big.xlsx", sep = "/"))
          short_ego <- short_ego[-id_too_big_rm,]
        }

        #Remove terms with only one gene from the gene list
        short_ego <- short_ego[which(short_ego$Count > 1),]


        #Simplify results
        if(dim(short_ego)[1] > 0 && length(which(short_ego$Count > 1)) > 0){
          #Remove redundant geneID list from results, keep only most significant one
          duplicate_geneID <- duplicated(short_ego$geneID)
          to_rm <- which(duplicate_geneID)
          name_too_long <- FALSE
          if(length(to_rm) > 0){
            removed_duplicate <- short_ego[to_rm,]
            short_ego <- short_ego[which(!duplicate_geneID),]
            writexl::write_xlsx(removed_duplicate, path = paste(DiseaseOntology_path, "DiseaseOntology_removed_duplicated_gene_list.xlsx", sep = "/"))

            #Create group for pathway terms with same geneID list
            dup_gene_list <- unique(removed_duplicate$geneID)
            nbr_dup_gene_list <- length(dup_gene_list)
            short_ego_name_shortened <- short_ego
            for (n in 1:nbr_dup_gene_list) {
              new_name <- paste(gsub(pattern = "/", replacement = "-", x = short_ego$Description[which(short_ego$geneID == dup_gene_list[n])]),
                                paste(gsub(pattern = "/", replacement = "-", x = removed_duplicate$Description[which(removed_duplicate$geneID == dup_gene_list[n])]), collapse = "_"),
                                sep = "_")
              short_ego$Description[which(short_ego$geneID == dup_gene_list[n])] <- new_name
              short_ego_name_shortened$Description[which(short_ego_name_shortened$geneID == dup_gene_list[n])] <- new_name
              if(nchar(new_name) > 50){
                new_name <- substr(new_name, 1, 50)
                short_ego_name_shortened$Description[which(short_ego_name_shortened$geneID == dup_gene_list[n])] <- new_name
                name_too_long <- TRUE
              }
            }
          }


          #Plotting the results
          ego_simplified <- ego
          ego_simplified@result <- short_ego
          if(name_too_long){
            ego_simplified@result <- short_ego_name_shortened
          }
          current_item_plot <- dim(ego_simplified@result)[1]
          if(current_item_plot > max_item_plot){
            current_item_plot <- max_item_plot
          }
          Cairo(paste(DiseaseOntology_path, "Dot_plot_DiseaseOntology_overrepresentation.png", sep = "/"),
                width = 1000,
                height = 1000)
          print(enrichplot::dotplot(ego_simplified, showCategory = current_item_plot))
          dev.off()

          tmp_width <- 600 + length(unique(strsplit(paste(short_ego$geneID, collapse = "/"), split = "/")[[1]])) * 10
          tmp_height <- 600 + current_item_plot * 10
          if(tmp_width > 2000){
            tmp_width <- 2000
          }
          if(tmp_height > 2000){
            tmp_height <- 2000
          }
          Cairo(paste(DiseaseOntology_path, "Plot_Gene_occurence_DiseaseOntology_pathway.png", sep = "/"),
                width = tmp_width,
                height = tmp_height)
          print(enrichplot::heatplot(ego_simplified, showCategory = current_item_plot))
          dev.off()

          Cairo(paste(DiseaseOntology_path, "Network_DiseaseOntology_overrepresentation.png", sep = "/"),
                width = 1000,
                height = 1000)
          p <- cnetplot(ego_simplified,
                        showCategory = current_item_plot,
                        categorySize="pvalue",
                        colorEdge = TRUE,
                        color_category='firebrick',
                        color_gene='steelblue')
          print(p)
          dev.off()

          Cairo(paste(DiseaseOntology_path, "Network_Circle_DiseaseOntology_overrepresentation.png", sep = "/"),
                width = 1000,
                height = 1000)
          p <- cnetplot(ego_simplified,
                        showCategory = current_item_plot,
                        categorySize="pvalue",
                        circular = TRUE,
                        colorEdge = TRUE,
                        color_category='firebrick',
                        color_gene='steelblue')
          print(p)
          dev.off()


          if(length(to_rm) > 0){
            #Create group for terms with same geneID list
            for (n in 1:nbr_dup_gene_list) {
              new_name <- paste(gsub(pattern = "/", replacement = "-", x = short_ego$ID[which(short_ego$geneID == dup_gene_list[n])]),
                                paste(gsub(pattern = "/", replacement = "-", x = removed_duplicate$ID[which(removed_duplicate$geneID == dup_gene_list[n])]), collapse = "_"),
                                sep = "_")
              short_ego$ID[which(short_ego$geneID == dup_gene_list[n])] <- new_name
            }
          }
          short_ego$Database <- "DiseaseOntology"


          #Save results
          writexl::write_xlsx(short_ego, path = paste(DiseaseOntology_path, "DiseaseOntology_over_representation.xlsx", sep = "/"))

        }#IF dim(short_ego)[1] > 0

        short_pathway_track[[4]] <- short_ego
      }else{
        short_pathway_track[[4]] <- NULL
      }


      ###NetCanGen analysis : Network of cancer gene database
      NetCanGen_path <- paste(pathway_analysis, "NetCanGen", sep = "/")
      ego <- DOSE::enrichNCG(gene = base::unique(ids_geneList[,2]),
                             pAdjustMethod = "BH",
                             qvalueCutoff = cutoff,
                             readable = TRUE)
      list_ego[[8]] <- ego
      if(!is.null(ego)){
        ego_result <- ego@result
        dir.create(NetCanGen_path)


        #Keep only significant results
        short_ego <- filter(ego_result , ego_result$p.adjust < cutoff)


        #Remove pathways with more than max_gene_per_item included (big pathway filter)
        id_too_big_rm <- which(short_ego$Count > max_gene_per_item)
        if(length(id_too_big_rm) > 0){
          remove_too_big <-  short_ego[id_too_big_rm,]
          writexl::write_xlsx(remove_too_big, path = paste(NetCanGen_path, "NetCanGen_removed_too_big.xlsx", sep = "/"))
          short_ego <- short_ego[-id_too_big_rm,]
        }

        #Remove terms with only one gene from the gene list
        short_ego <- short_ego[which(short_ego$Count > 1),]


        #Simplify results
        if(dim(short_ego)[1] > 0 && length(which(short_ego$Count > 1)) > 0){
          #Remove redundant geneID list from results, keep only most significant one
          duplicate_geneID <- duplicated(short_ego$geneID)
          to_rm <- which(duplicate_geneID)
          name_too_long <- FALSE
          if(length(to_rm) > 0){
            removed_duplicate <- short_ego[to_rm,]
            short_ego <- short_ego[which(!duplicate_geneID),]
            writexl::write_xlsx(removed_duplicate, path = paste(NetCanGen_path, "NetCanGen_removed_duplicated_gene_list.xlsx", sep = "/"))

            #Create group for pathway terms with same geneID list
            dup_gene_list <- unique(removed_duplicate$geneID)
            nbr_dup_gene_list <- length(dup_gene_list)
            short_ego_name_shortened <- short_ego
            for (n in 1:nbr_dup_gene_list) {
              new_name <- paste(gsub(pattern = "/", replacement = "-", x = short_ego$Description[which(short_ego$geneID == dup_gene_list[n])]),
                                paste(gsub(pattern = "/", replacement = "-", x = removed_duplicate$Description[which(removed_duplicate$geneID == dup_gene_list[n])]), collapse = "_"),
                                sep = "_")
              short_ego$Description[which(short_ego$geneID == dup_gene_list[n])] <- new_name
              short_ego_name_shortened$Description[which(short_ego_name_shortened$geneID == dup_gene_list[n])] <- new_name
              if(nchar(new_name) > 50){
                new_name <- substr(new_name, 1, 50)
                short_ego_name_shortened$Description[which(short_ego_name_shortened$geneID == dup_gene_list[n])] <- new_name
                name_too_long <- TRUE
              }
            }
          }


          #Plotting the results
          ego_simplified <- ego
          ego_simplified@result <- short_ego
          if(name_too_long){
            ego_simplified@result <- short_ego_name_shortened
          }
          current_item_plot <- dim(ego_simplified@result)[1]
          if(current_item_plot > max_item_plot){
            current_item_plot <- max_item_plot
          }
          Cairo(paste(NetCanGen_path, "Dot_plot_NetCanGen_overrepresentation.png", sep = "/"),
                width = 1000,
                height = 1000)
          print(enrichplot::dotplot(ego_simplified, showCategory = current_item_plot))
          dev.off()

          tmp_width <- 600 + length(unique(strsplit(paste(short_ego$geneID, collapse = "/"), split = "/")[[1]])) * 10
          tmp_height <- 600 + current_item_plot * 10
          if(tmp_width > 2000){
            tmp_width <- 2000
          }
          if(tmp_height > 2000){
            tmp_height <- 2000
          }
          Cairo(paste(NetCanGen_path, "Plot_Gene_occurence_NetCanGen_pathway.png", sep = "/"),
                width = tmp_width,
                height = tmp_height)
          print(enrichplot::heatplot(ego_simplified, showCategory = current_item_plot))
          dev.off()


          Cairo(paste(NetCanGen_path, "Network_NetCanGen_overrepresentation.png", sep = "/"),
                width = 1000,
                height = 1000)
          p <- cnetplot(ego_simplified,
                        showCategory = current_item_plot,
                        categorySize="pvalue",
                        colorEdge = TRUE,
                        color_category='firebrick',
                        color_gene='steelblue')
          print(p)
          dev.off()

          Cairo(paste(NetCanGen_path, "Network_Circle_NetCanGen_overrepresentation.png", sep = "/"),
                width = 1000,
                height = 1000)
          p <- cnetplot(ego_simplified,
                        showCategory = current_item_plot,
                        categorySize="pvalue",
                        circular = TRUE,
                        colorEdge = TRUE,
                        color_category='firebrick',
                        color_gene='steelblue')
          print(p)
          dev.off()


          if(length(to_rm) > 0){
            #Create group for terms with same geneID list
            for (n in 1:nbr_dup_gene_list) {
              new_name <- paste(gsub(pattern = "/", replacement = "-", x = short_ego$ID[which(short_ego$geneID == dup_gene_list[n])]),
                                paste(gsub(pattern = "/", replacement = "-", x = removed_duplicate$ID[which(removed_duplicate$geneID == dup_gene_list[n])]), collapse = "_"),
                                sep = "_")
              short_ego$ID[which(short_ego$geneID == dup_gene_list[n])] <- new_name
            }
          }
          short_ego$Database <- "NetCanGen"


          #Save results
          writexl::write_xlsx(short_ego, path = paste(NetCanGen_path, "NetCanGen_over_representation.xlsx", sep = "/"))

        }#IF dim(short_ego)[1] > 0

        short_pathway_track[[5]] <- short_ego
      }else{
        short_pathway_track[[5]] <- NULL
      }


      ###DisGeNET analysis : Disease Gene Network database
      DisGeNET_path <- paste(pathway_analysis, "DisGeNET", sep = "/")
      if(length(background_gene) == 0){
        ego <- DOSE::enrichDGN(gene = base::unique(ids_geneList[,2]),
                               pAdjustMethod = "BH",
                               qvalueCutoff = cutoff,
                               readable = TRUE)
      }else{
        ego <- DOSE::enrichDGN(gene = base::unique(ids_geneList[,2]),
                               pAdjustMethod = "BH",
                               universe = background_gene,
                               qvalueCutoff = cutoff,
                               readable = TRUE)
      }

      list_ego[[9]] <- ego
      if(!is.null(ego)){
        ego_result <- ego@result


        #Keep only significant results
        short_ego <- filter(ego_result , ego_result$p.adjust < cutoff)

        dir.create(DisGeNET_path)
        #Remove pathways with more than max_gene_per_item included (big pathway filter)
        id_too_big_rm <- which(short_ego$Count > max_gene_per_item)
        if(length(id_too_big_rm) > 0){
          remove_too_big <-  short_ego[id_too_big_rm,]
          writexl::write_xlsx(remove_too_big, path = paste(DisGeNET_path, "DisGeNET_removed_too_big.xlsx", sep = "/"))
          short_ego <- short_ego[-id_too_big_rm,]
        }

        #Remove terms with only one gene from the gene list
        short_ego <- short_ego[which(short_ego$Count > 1),]


        #Simplify results
        if(dim(short_ego)[1] > 0 && length(which(short_ego$Count > 1)) > 0){
          #Remove redundant geneID list from results, keep only most significant one
          duplicate_geneID <- duplicated(short_ego$geneID)
          to_rm <- which(duplicate_geneID)
          name_too_long <- FALSE
          if(length(to_rm) > 0){
            removed_duplicate <- short_ego[to_rm,]
            short_ego <- short_ego[which(!duplicate_geneID),]
            writexl::write_xlsx(removed_duplicate, path = paste(DisGeNET_path, "DisGeNET_removed_duplicated_gene_list.xlsx", sep = "/"))

            #Create group for pathway terms with same geneID list
            dup_gene_list <- unique(removed_duplicate$geneID)
            nbr_dup_gene_list <- length(dup_gene_list)
            short_ego_name_shortened <- short_ego
            for (n in 1:nbr_dup_gene_list) {
              new_name <- paste(gsub(pattern = "/", replacement = "-", x = short_ego$Description[which(short_ego$geneID == dup_gene_list[n])]),
                                paste(gsub(pattern = "/", replacement = "-", x = removed_duplicate$Description[which(removed_duplicate$geneID == dup_gene_list[n])]), collapse = "_"),
                                sep = "_")
              short_ego$Description[which(short_ego$geneID == dup_gene_list[n])] <- new_name
              short_ego_name_shortened$Description[which(short_ego_name_shortened$geneID == dup_gene_list[n])] <- new_name
              if(nchar(new_name) > 50){
                new_name <- substr(new_name, 1, 50)
                short_ego_name_shortened$Description[which(short_ego_name_shortened$geneID == dup_gene_list[n])] <- new_name
                name_too_long <- TRUE
              }
            }
          }


          #Plotting the results
          ego_simplified <- ego
          ego_simplified@result <- short_ego
          if(name_too_long){
            ego_simplified@result <- short_ego_name_shortened
          }
          current_item_plot <- dim(ego_simplified@result)[1]
          if(current_item_plot > max_item_plot){
            current_item_plot <- max_item_plot
          }
          Cairo(paste(DisGeNET_path, "Dot_plot_DisGeNET_overrepresentation.png", sep = "/"),
                width = 1000,
                height = 1000)
          print(enrichplot::dotplot(ego_simplified, showCategory = current_item_plot))
          dev.off()

          tmp_width <- 600 + length(unique(strsplit(paste(short_ego$geneID, collapse = "/"), split = "/")[[1]])) * 10
          tmp_height <- 600 + current_item_plot * 10
          if(tmp_width > 2000){
            tmp_width <- 2000
          }
          if(tmp_height > 2000){
            tmp_height <- 2000
          }
          Cairo(paste(DisGeNET_path, "Plot_Gene_occurence_DisGeNET_pathway.png", sep = "/"),
                width = tmp_width,
                height = tmp_height)
          print(enrichplot::heatplot(ego_simplified, showCategory = current_item_plot))
          dev.off()


          Cairo(paste(DisGeNET_path, "Network_DisGeNET_overrepresentation.png", sep = "/"),
                width = 1000,
                height = 1000)
          p <- cnetplot(ego_simplified,
                        showCategory = current_item_plot,
                        categorySize="pvalue",
                        colorEdge = TRUE,
                        color_category='firebrick',
                        color_gene='steelblue')
          print(p)
          dev.off()

          Cairo(paste(DisGeNET_path, "Network_Circle_DisGeNET_overrepresentation.png", sep = "/"),
                width = 1000,
                height = 1000)
          p <- cnetplot(ego_simplified,
                        showCategory = current_item_plot,
                        categorySize="pvalue",
                        circular = TRUE,
                        colorEdge = TRUE,
                        color_category='firebrick',
                        color_gene='steelblue')
          print(p)
          dev.off()


          if(length(to_rm) > 0){
            #Create group for terms with same geneID list
            for (n in 1:nbr_dup_gene_list) {
              new_name <- paste(gsub(pattern = "/", replacement = "-", x = short_ego$ID[which(short_ego$geneID == dup_gene_list[n])]),
                                paste(gsub(pattern = "/", replacement = "-", x = removed_duplicate$ID[which(removed_duplicate$geneID == dup_gene_list[n])]), collapse = "_"),
                                sep = "_")
              short_ego$ID[which(short_ego$geneID == dup_gene_list[n])] <- new_name
            }
          }
          short_ego$Database <- "DisGeNET"


          #Save results
          writexl::write_xlsx(short_ego, path = paste(DisGeNET_path, "DisGeNET_over_representation.xlsx", sep = "/"))

        }#IF dim(short_ego)[1] > 0

        short_pathway_track[[6]] <- short_ego
      }else{
        short_pathway_track[[6]] <- NULL
      }


      ###Results combination
      nbr_result <- length(short_ego_track) + length(short_pathway_track)
      name_column <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count", "Database")
      combined_results <- data.frame(matrix(NA, ncol = length(name_column), nrow = 0))
      colnames(combined_results) <- name_column
      count <- 1
      for (i in 1:nbr_result) {
        if(i <= length(short_ego_track)){
          #Retrieving GO term results
          tmp <- short_ego_track[[i]]
          if(!is.null(tmp)){
            #if(dim(tmp)[2] < 10 && dim(tmp)[1] > 0){
            #  tmp$Database <- paste("GO", go_of_interest[i], sep = "_")
            #}
            if(dim(tmp)[1] > 0){
              tmp <- tmp[,which(colnames(tmp) %in% name_column)]
              combined_results <- rbind(combined_results, tmp)
            }
          }
        }else{
          #Retrieving Pathway results
          tmp <- short_pathway_track[[count]]
          if(!is.null(tmp)){
            #if(dim(tmp)[2] < 10 && dim(tmp)[1] > 0){
            #  tmp$Database <- pathway_database[count]
            #}
            if(dim(tmp)[1] > 0){
              tmp <- tmp[,which(colnames(tmp) %in% name_column)]
              combined_results <- rbind(combined_results, tmp)
            }
          }

          count <- count + 1
        }
      }


      #reorder combined resulta to lowest p.adjust to biggest p.adjust
      combined_results <- combined_results %>%
        arrange(p.adjust)

      #Remove redundant geneID list from combined results, keep only the first one
      duplicate_geneID <- duplicated(combined_results$geneID)
      to_rm <- which(duplicate_geneID)
      name_too_long <- FALSE
      if(length(to_rm) > 0){
        removed_duplicate <- combined_results[to_rm,]
        combined_results <- combined_results[which(!duplicate_geneID),]
        writexl::write_xlsx(removed_duplicate, path = paste(path_to_store_results, "Combined_removed_duplicated_gene_list.xlsx", sep = "/"))

        #Create group for pathway terms with same geneID list
        dup_gene_list <- unique(removed_duplicate$geneID)
        nbr_dup_gene_list <- length(dup_gene_list)
        combined_results_name_shortened <- combined_results
        for (n in 1:nbr_dup_gene_list) {
          new_name <- paste(gsub(pattern = "/", replacement = "-", x = combined_results$Description[which(combined_results$geneID == dup_gene_list[n])]),
                            paste(gsub(pattern = "/", replacement = "-", x = removed_duplicate$Description[which(removed_duplicate$geneID == dup_gene_list[n])]), collapse = "_"),
                            sep = "_")
          combined_results$Description[which(combined_results$geneID == dup_gene_list[n])] <- new_name
          combined_results_name_shortened$Description[which(combined_results_name_shortened$geneID == dup_gene_list[n])] <- new_name
          if(nchar(new_name) > 50){
            new_name <- substr(new_name, 1, 50)
            combined_results_name_shortened$Description[which(combined_results_name_shortened$geneID == dup_gene_list[n])] <- new_name
            name_too_long <- TRUE
          }
        }
      }


      #Plotting the results
      if(is.null(ego)){
        ego <- list_ego[[length(list_ego)]]
      }
      ego_simplified <- ego
      ego_simplified@result <- combined_results
      ego_simplified@result$ID <- rownames(ego_simplified@result)
      if(dim(ego_simplified@result)[1] > 0){
        if(name_too_long){
          ego_simplified@result <- combined_results_name_shortened
        }
        current_item_plot <- dim(ego_simplified@result)[1]
        if(current_item_plot > max_item_plot){
          current_item_plot <- max_item_plot
        }
        Cairo(paste(path_to_store_results, "Dot_plot_combined_overrepresentation.png", sep = "/"),
              width = 1000,
              height = 1000)
        print(enrichplot::dotplot(ego_simplified, showCategory = current_item_plot))
        dev.off()


        if(length(to_rm) > 0){
          #Create group for terms with same geneID list
          for (n in 1:nbr_dup_gene_list) {
            new_name <- paste(gsub(pattern = "/", replacement = "-", x = combined_results$ID[which(combined_results$geneID == dup_gene_list[n])]),
                              paste(gsub(pattern = "/", replacement = "-", x = removed_duplicate$ID[which(removed_duplicate$geneID == dup_gene_list[n])]), collapse = "_"),
                              sep = "_")
            combined_results$ID[which(combined_results$geneID == dup_gene_list[n])] <- new_name
            #Database
            new_name <- paste(gsub(pattern = "/", replacement = "-", x = combined_results$Database[which(combined_results$geneID == dup_gene_list[n])]),
                              paste(gsub(pattern = "/", replacement = "-", x = removed_duplicate$Database[which(removed_duplicate$geneID == dup_gene_list[n])]), collapse = "_"),
                              sep = "_")
            combined_results$Database[which(combined_results$geneID == dup_gene_list[n])] <- new_name
          }
        }

      }


      #Gene never over represented in original input
      id_rm <- which(!geneList %in% unique(strsplit(paste(combined_results$geneID, collapse = "/"), split = "/")[[1]]))
      if(length(id_rm) > 0){
        name_no_OR <- geneList[id_rm]
        writexl::write_xlsx(data.frame(Gene = name_no_OR), path = paste(path_to_store_results, "Gene_not_over_represented.xlsx", sep = "/"))
      }


      ###Summary
      ORA_summary_col <- c("geneList", "GeneRatio", "BestPvalueAdjusted", "Nbr_ID", "ID", "Description", "Database", "Nbr_Database")
      ORA_summary <- data.frame(matrix(NA, ncol = length(ORA_summary_col), nrow = dim(combined_results)[1]))
      colnames(ORA_summary) <- ORA_summary_col
      ORA_summary$geneList <- combined_results$geneID
      ORA_summary$GeneRatio <- combined_results$GeneRatio
      ORA_summary$BestPvalueAdjusted <- combined_results$qvalue
      ORA_summary$ID <- combined_results$ID
      ORA_summary$Description <- combined_results$Description
      ORA_summary$Database <- combined_results$Database
      if(dim(combined_results)[1] > 0){
        for (i in 1:dim(ORA_summary)[1]) {
          tmp_ID <- ORA_summary$ID[i]
          tmp_ID <- length(strsplit(tmp_ID, split = "_")[[1]])
          ORA_summary$Nbr_ID[i] <- tmp_ID
          tmp_Database <- ORA_summary$Database[i]
          tmp_Database <- length(strsplit(tmp_Database, split = "_")[[1]])
          ORA_summary$Nbr_Database[i] <- tmp_Database
        }
      }


      if(dim(ORA_summary)[1] > 0){
        writexl::write_xlsx(ORA_summary, path = paste(path_to_store_results, "Summary_results_over_representation.xlsx", sep = "/"))
        if(coocurrence){
          start_matrix <- Sys.time()
          print("Starting gene co-occurence analysis after over-representation analysis")
          #restrict analysis on top most significant GO/pathway
          if(dim(ORA_summary)[1] > max_coocurrence){
            ORA_summary <- ORA_summary[c(1:max_coocurrence),]
          }


          if(gene_annot != "SYMBOL"){
            geneList <- bitr(geneList, fromType=gene_annot, toType="SYMBOL", OrgDb="org.Hs.eg.db")$SYMBOL
          }
          #By gene based on input geneList
          nbr_gene <- length(geneList)
          geneBygene_matrix <- matrix(0, ncol = nbr_gene, nrow = nbr_gene)
          colnames(geneBygene_matrix) <- geneList
          rownames(geneBygene_matrix) <- geneList
          for (i in 1:dim(ORA_summary)[1]) {
            tmp_list_gene <- strsplit(ORA_summary$geneList[i], split = "/")[[1]]
            for (j in 1:length(tmp_list_gene)) {
              tmp_gene <- tmp_list_gene[j]
              geneBygene_matrix[which(rownames(geneBygene_matrix) == tmp_gene),which(colnames(geneBygene_matrix) %in% tmp_list_gene)] <- geneBygene_matrix[which(rownames(geneBygene_matrix) == tmp_gene),which(colnames(geneBygene_matrix) %in% tmp_list_gene)] + 1
            }
          }
          id_rm <- which(colSums(geneBygene_matrix) == 0)
          if(length(id_rm) > 0){
            geneBygene_matrix <- geneBygene_matrix[-id_rm, -id_rm]
          }
          geneBygene_data <- data.frame(Gene = rownames(geneBygene_matrix))
          geneBygene_data <- cbind(geneBygene_data, as.data.frame(geneBygene_matrix))
          writexl::write_xlsx(geneBygene_data, path = paste(path_to_store_results, "Cooccurence_genes.xlsx", sep = "/"))
          end_matrix <- Sys.time()


          start_plot <- Sys.time()
          if(dim(geneBygene_matrix)[2] > 100){
            max_id <- colSums(geneBygene_matrix)
            condition_plot <- min(max_id[order(max_id, decreasing = TRUE)][1:100])
            to_keep_plot <- which(max_id >= condition_plot)[1:100]
            geneBygene_matrix_reduce <- geneBygene_matrix[to_keep_plot,to_keep_plot]
          }else{
            geneBygene_matrix_reduce <- geneBygene_matrix
          }

          #ploting coocurrence of genes
          name_of_plot <- paste(path_to_store_results, "Coocurrence_genes.png", sep = "/")
          Cairo(name_of_plot,
                width = 1500, height = 1200)

          h1 <- Heatmap(geneBygene_matrix_reduce,
                        heatmap_legend_param = list(
                          title = "CoOcurrence"),
                        column_km = 1,
                        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                          if(geneBygene_matrix_reduce[i, j] >= 1)
                            grid.text(sprintf("%d", geneBygene_matrix_reduce[i, j]), x, y, gp = gpar(fontsize = 10))
                        }
          )
          draw(h1)
          dev.off()
          print("End gene co-occurence analysis after over-representation analysis")
          end_plot <- Sys.time()
        }#enf of If coocurrence
      }


      #If there is term with same ID but coming from different databases and have different gene list
      id_dup <- which(duplicated(combined_results$Description))
      if(length(id_dup) > 0){
        for(i in 1:length(id_dup)){
          tmp_id <- which(combined_results$Description == combined_results$Description[id_dup[i]])
          for(j in 1:length(tmp_id)){
            combined_results$Description[tmp_id[j]] <- paste(combined_results$Database[tmp_id[j]],combined_results$Description[tmp_id[j]], sep = "")
          }

          print(combined_results$Description[tmp_id])
        }
      }


      #Returning results
      #Pathways
      Description <- c(combined_results$Description)
      geneID <- c(combined_results$geneID)
      Database <- c(combined_results$Database)
      over_represented_pathway <- list(Description, geneID, Database)

      stop <- Sys.time()
      print(stop - start)
      return(over_represented_pathway)
    }else{
      print("only one gene left for over representation analysis, therefore step is not performed")
      return(list(NA, NA, NA))
    }

  }### End of if clusterProfiler databases




}#end of function

