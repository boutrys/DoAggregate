############################################# Function to adjust the genotype matrix via coverage matrix #############################################
# The function can adapt the Genotype matrix given the genotype matrix

coverage_matrix_creation <- function(list_patient = c(),
                                     target_position = c(),
                                     path_to_bam = "", 
                                     mapqv_threshold = 20){
  #using php script to go into highlander database 
  if(path_to_bam == "exomes_1000g" | path_to_bam == "exomes_haplotype_caller"){
    target_position$posEnd <- target_position$pos
    current_time_bed <- unclass(Sys.time())
    temporary_bed_file <- as.character(as.integer(current_time_bed))
    temporary_file <- tempfile(pattern = temporary_bed_file, tmpdir = tempdir(), fileext = ".txt")
    write_tsv(target_position,temporary_file, col_names = FALSE) 
    
    list_patient <- base::paste(list_patient, collapse = "+")
    current_time <- unclass(Sys.time())
    output <- paste(as.character(as.integer(current_time)), ".tsv", sep = "")
    query_list <- list()
    query_list[[1]] <- upload_file(temporary_file)
    names(query_list)[1] <- "bed"
    query_list[[2]] <- list_patient
    names(query_list)[2] <- "patients"
    query_list[[3]] <- path_to_bam
    names(query_list)[3] <- "analysis"
    query_list[[4]] <- mapqv_threshold
    names(query_list)[4] <- "mapqv"
    query_list[[5]] <- output
    names(query_list)[5] <- "output"
    name_url <- paste(Sys.getenv("EXCALIBUR_COVERAGE_SERVER", "http://localhost"), "/php/coveragematrix.php", sep = "")
    
    p <- POST(name_url, body = query_list)
    #p <- POST(name_url, body = query_list, verbose(data_out = FALSE, data_in = TRUE, ssl = TRUE), encode = "multipart")
    resp <- FALSE
    while(!resp){
      resp <- url.exists(url = paste(Sys.getenv("EXCALIBUR_COVERAGE_SERVER", "http://localhost"), "/bam/bamout/", output, sep = ""))
      Sys.sleep(60)
    }
    answer <- GET(url = paste(Sys.getenv("EXCALIBUR_COVERAGE_SERVER", "http://localhost"), "/bam/bamout/", output, sep = ""))
    coverage_matrix <- httr::content(answer, as = "parsed")
    return(coverage_matrix)
  }else{
    #still to do, if bam are in a folder
    stop("no implementation yet for bam stored in folder")
  }
}
