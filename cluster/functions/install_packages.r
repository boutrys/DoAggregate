#!/usr/bin/Rscript
### If new version R, script to reinstall all needed packages
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("2 arguments must be supplied ()",
   call.=FALSE)
}


####### Parameter of main function
path_Rlibrary <- args[1]
do_install <- args[2]


.libPaths(path_Rlibrary, FALSE)
print(.libPaths())
print(getRversion())


if(do_install == 1){
    # Set the CRAN mirror URL
    cran_mirror <- "https://cloud.r-project.org"
    # Set the mirror option
    options(repos = c(CRAN = cran_mirror))
    

    #CRAN package to install via install.packages("NameOfPackage")
    install.packages("RCurl", lib=path_Rlibrary)
    install.packages("writexl", lib=path_Rlibrary)
    install.packages("readxl", lib=path_Rlibrary)
    install.packages("data.table", lib=path_Rlibrary)
    install.packages("ggdendro", lib=path_Rlibrary)
    install.packages("tidyverse", lib=path_Rlibrary)
    install.packages("circlize", lib=path_Rlibrary)
    install.packages("minqa", lib=path_Rlibrary)
    install.packages("CompQuadForm", lib=path_Rlibrary)
    install.packages("SKAT", lib=path_Rlibrary)
    install.packages("CATT", lib=path_Rlibrary)
    install.packages("iGasso", lib=path_Rlibrary)
    install.packages("WGScan", lib=path_Rlibrary)
    install.packages("Cairo", lib=path_Rlibrary)
    install.packages("vcfR", lib=path_Rlibrary)


    #Bioconductor packages
    if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", lib=path_Rlibrary)
    BiocManager::install("REBET")
    BiocManager::install("GOSemSim")
    BiocManager::install("org.Hs.eg.db")
    BiocManager::install("DOSE")
    BiocManager::install("biomaRt")
    BiocManager::install("ReactomePA")
    BiocManager::install("enrichplot")
    BiocManager::install("ComplexHeatmap")


    ##necessary for install_github
    install.packages("devtools", lib=path_Rlibrary)
    library(devtools)

    #Github packages
    install_github("baolinwu/KATSP")
    install_github('gastonstat/AssotesteR')
    remotes::install_github("YuLab-SMU/clusterProfiler") 


    #Download tar.gz file at https://cran.r-project.org/src/contrib/Archive/RVtests/
    #install.packages("C:/Users/boutrys/OneDrive - UCL/DRAFT_ARTICLES/Statistical_framework_ensemble/To_publish_Plos_computational_biology/Response_to_reviewers/SIMULATION/functions/library_to_install/RVtests_1.2.tar.gz", repos=NULL, type="source")
    install.packages("${path_Rlibrary}/RVtests_1.2.tar.gz", repos=NULL, type="source", lib=path_Rlibrary)

}








