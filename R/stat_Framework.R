#' Statistical Framework for Genetic Association Testing
#'
#' @name stat_Framework
#' @description This function provides a comprehensive statistical analysis for genetic association studies.
#' It can perform a variety of statistical tests based on the input genotype and phenotype data,
#' including burden tests, variance-component tests, and omnibus tests. The function allows
#' specifying whether to perform a binary or continuous analysis and supports adjustments for
#' covariates and minor allele frequency (MAF) thresholds.
#'
#' @param genotype_matrix A matrix of genotype data where rows represent samples and columns represent genetic variants.
#' @param phenotype A vector of phenotypes corresponding to each sample in the genotype matrix.
#' @param weight_variant A vector of weights for the variants, used in weighted statistical tests.
#' @param Binary Logical; if TRUE, performs tests appropriate for binary traits, otherwise for continuous traits.
#' @param covariate Covariates to be included in the model, if any.
#' @param rare_maf_threshold The threshold for defining rare variants based on minor allele frequency.
#' @param filter.pval The p-value threshold for filtering results.
#' @param position A vector indicating the position of each variant, used in some types of tests.
#' @param get_test Logical; if TRUE, returns the names of available tests instead of performing them.
#' @param version Specifies the version of the tests to run ('all', 'optimal', or a specific subset).
#' @param max_var_per_test Maximum number of variants allowed per test to ensure computational feasibility.
#' @return A list containing three elements: p-values from the tests, computational times for each test, and any errors encountered during testing.
#'
#' @importFrom stats p.adjust p.adjust.methods
#' @export


#in further developpement, we should split burden, skat and omnibus tests into different function,
#or put an additionnal argument to allow to compute only for a specific class of test
#put additional argument for the model, here we always suppose y~1
stat_Framework <- function(genotype_matrix = c(),
                           phenotype = c(),
                           weight_variant = c(),
                           Binary = TRUE,
                           covariate = 1,
                           rare_maf_threshold = 0.01,
                           filter.pval = 0.05,
                           position = c(),
                           get_test = FALSE,
                           version = "optimal",
                           max_var_per_test = 300){


  #Not running any test, just return number of test and their names
  if(get_test){
    if(Binary){
      if(version == "all"){
        name_of_test <- c("Excalibur_baseline",
                          "robust_burden",
                          "robust_SKAT",
                          "robust_SKATO",
                          "p_linear_liu_burden",
                          "p_linear_weighted_liumod_burden",
                          "p_linear_davies_skat",
                          "p_linear_weighted_davies_skat",
                          "pBin_linear_IBS_SKAT_MA",
                          "pBin_weighted_IBS_SKAT_MA",
                          "pBin_2wayIX_SKAT_MA",
                          "pBin_weighted_quadratic_SKAT_ERA",
                          "pBin_linear_SKATO_MA",
                          "pBin_linear_weighted_SKATO_MA",
                          "p_ascore",
                          "p_ascore_ord",
                          "p_assu",
                          "p_assu_ord",
                          "p_asum",
                          "p_asum_ord",
                          "p_bst",
                          "p_calpha",
                          "p_calpha_asymptopic",
                          "p_carv_hard",
                          "p_carv_variable",
                          "p_carv_stepup",
                          "p_cast_fisher",
                          "p_cast_chisq",
                          "p_cmat",
                          "p_cmc",
                          "p_rbt",
                          "p_rvt1",
                          "p_rvt2",
                          "p_rwas",
                          "p_score",
                          "p_score_asymptopic",
                          "p_ssu",
                          "p_ssu_asymptopic",
                          "p_ssuw",
                          "p_ssuw_asymptopic",
                          "p_ttest_asymptopic",
                          "p_vt",
                          "p_wss",
                          "p_wst",
                          "p_wst_asymptopic",
                          "p_catt",
                          "p_KAT",
                          "p_SKATplus",
                          "p_wgscan_region",
                          "p_rebet",
                          "p_pcr",
                          "p_rr",
                          "p_spls",
                          "p_t1p",
                          "p_t5p",
                          "p_wep",
                          "p_score_vtp",
                          "p_wod01",
                          "p_wod05",
                          "p_ada")
        return(name_of_test)
      }else if(version == "optimal"){
        name_of_test <- c("Excalibur",
                          "robust_burden",
                          "robust_SKAT",
                          "robust_SKATO",
                          "p_ascore",
                          "p_ascore_ord",
                          "p_assu",
                          "p_assu_ord",
                          "p_asum",
                          "p_asum_ord",
                          "p_bst",
                          "p_calpha",
                          "p_cast_chisq",
                          "p_cmat",
                          "p_cmc",
                          "p_rbt",
                          "p_rvt1",
                          "p_rvt2",
                          "p_rwas",
                          "p_score",
                          "p_score_asymptopic",
                          "p_ssu",
                          "p_ssu_asymptopic",
                          "p_ssuw",
                          "p_ssuw_asymptopic",
                          "p_ttest_asymptopic",
                          "p_vt",
                          "p_wss",
                          "p_wst",
                          "p_KAT",
                          "p_rebet",
                          "p_rr",
                          "p_spls",
                          "p_t1p",
                          "p_t5p",
                          "p_wep",
                          "p_score_vtp")
        return(name_of_test)
      }else{
        #if a user wants to perform using specific methods
        complete_list_test <- c("Excalibur",
                          "robust_burden",
                          "robust_SKAT",
                          "robust_SKATO",
                          "p_linear_liu_burden",
                          "p_linear_weighted_liumod_burden",
                          "p_linear_davies_skat",
                          "p_linear_weighted_davies_skat",
                          "pBin_linear_IBS_SKAT_MA",
                          "pBin_weighted_IBS_SKAT_MA",
                          "pBin_2wayIX_SKAT_MA",
                          "pBin_weighted_quadratic_SKAT_ERA",
                          "pBin_linear_SKATO_MA",
                          "pBin_linear_weighted_SKATO_MA",
                          "p_ascore",
                          "p_ascore_ord",
                          "p_assu",
                          "p_assu_ord",
                          "p_asum",
                          "p_asum_ord",
                          "p_bst",
                          "p_calpha",
                          "p_calpha_asymptopic",
                          "p_carv_hard",
                          "p_carv_variable",
                          "p_carv_stepup",
                          "p_cast_fisher",
                          "p_cast_chisq",
                          "p_cmat",
                          "p_cmc",
                          "p_rbt",
                          "p_rvt1",
                          "p_rvt2",
                          "p_rwas",
                          "p_score",
                          "p_score_asymptopic",
                          "p_ssu",
                          "p_ssu_asymptopic",
                          "p_ssuw",
                          "p_ssuw_asymptopic",
                          "p_ttest_asymptopic",
                          "p_vt",
                          "p_wss",
                          "p_wst",
                          "p_wst_asymptopic",
                          "p_catt",
                          "p_KAT",
                          "p_SKATplus",
                          "p_wgscan_region",
                          "p_rebet",
                          "p_pcr",
                          "p_rr",
                          "p_spls",
                          "p_t1p",
                          "p_t5p",
                          "p_wep",
                          "p_score_vtp",
                          "p_wod01",
                          "p_wod05",
                          "p_ada")
        name_of_test <- complete_list_test[which(complete_list_test %in% version)]
        if(length(name_of_test) == length(version) && length(version) > 0){
          return(name_of_test)
        }else{
          stop("error in the name of test selected, no test which such name found")
        }
      }
    }else{
      name_of_test <- c("Excalibur",
                        "p_1_linear_davies_burden",
                        "p_2_linear_weighted_liu_burden",
                        "p_3_linear_weighted_liumod_skat",
                        "p_4_weighted_quadratic_liumod_skat",
                        "p_5_linear_liu_skat",
                        "p_6_linear_IBS_liu_skat",
                        "p_7_linear_skato",
                        "p_8_linear_weighted_skato")
      return(name_of_test)
    }

  }else{
    #preprocess and initialization
    genotype_matrix <- as.matrix(genotype_matrix)

    errors <- data.frame()
    err_robust_burden <- 0
    err_robust_SKAT <- 0
    err_robust_SKATO <- 0
    err_p_linear_liu_burden <- 0
    err_p_linear_weighted_liumod_burden <- 0
    err_p_linear_davies_skat <- 0
    err_p_linear_weighted_davies_skat <- 0
    err_pBin_linear_IBS_SKAT_MA <- 0
    err_pBin_weighted_IBS_SKAT_MA <- 0
    err_pBin_2wayIX_SKAT_MA <- 0
    err_pBin_weighted_quadratic_SKAT_ERA <- 0
    err_pBin_linear_SKATO_MA <- 0
    err_pBin_linear_weighted_SKATO_MA <- 0
    err_p_ascore <- 0
    err_p_ascore_ord <- 0
    err_p_assu <- 0
    err_p_assu_ord <- 0
    err_p_asum <- 0
    err_p_asum_ord <- 0
    err_p_bst <- 0
    err_p_calpha <- 0
    err_p_calpha_asymptopic <- 0
    err_p_carv_hard <- 0
    err_p_carv_variable <- 0
    err_p_carv_stepup <- 0
    err_p_cast_fisher <- 0
    err_p_cast_chisq <- 0
    err_p_cmat <- 0
    err_p_cmc <- 0
    err_p_rbt <- 0
    err_p_rvt1 <- 0
    err_p_rvt2 <- 0
    err_p_rwas <- 0
    err_p_score <- 0
    err_p_score_asymptopic <- 0
    err_p_ssu <- 0
    err_p_ssu_asymptopic <- 0
    err_p_ssuw <- 0
    err_p_ssuw_asymptopic <- 0
    err_p_ttest_asymptopic <- 0
    err_p_vt <- 0
    err_p_wss <- 0
    err_p_wst <- 0
    err_p_wst_asymptopic <- 0
    err_p_catt <- 0
    err_p_KAT <- 0
    err_p_SKATplus <- 0
    err_p_wgscan_region <- 0
    err_p_rebet <- 0
    err_p_pcr <- 0
    err_p_rr <- 0
    err_p_spls <- 0
    err_p_t1p <- 0
    err_p_t5p <- 0
    err_p_wep <- 0
    err_p_score_vtp <- 0
    err_p_wod01 <- 0
    err_p_wod05 <- 0
    err_p_ada <- 0


    if(Binary){

      ###Do the null model
      #CATT preprocess data
      catt_matrix <- matrix(NA, ncol = dim(genotype_matrix)[2], nrow = 2)
      catt_matrix[1,] <- colSums(genotype_matrix[which(phenotype == 0),])
      catt_matrix[2,] <- colSums(genotype_matrix[which(phenotype == 1),])

      if(length(covariate) == 1){
        obj <-SKAT_Null_Model_MomentAdjust(phenotype ~ 1, data = data.frame(genotype_matrix), type.Resampling = "bootstrap")
        obj_wgscan <- WGScan.prelim(Y = phenotype, out_type="D")
      }else{
        obj <-SKAT_Null_Model_MomentAdjust(phenotype ~ covariate, data = data.frame(genotype_matrix), type.Resampling = "bootstrap")
        obj_wgscan <- WGScan.prelim(Y = phenotype, X = covariate, out_type="D")
      }

      ########################  Perform all tests  ########################
      tmp_pvalue <- c()
      tmp_computational_time <- c()
      tmp_error <- c()

      if(version == "all" | version == "optimal" | length(which(version == "robust_burden"))> 0){
        start_time <- Sys.time()
        if(length(which(is.na(weight_variant))) > 0){
          err <- try(tmp <-  SKAT::SKATBinary_Robust(genotype_matrix, obj, method = "Burden"),
                  silent = TRUE)
        }else{
          err <- try(tmp <-  SKAT::SKATBinary_Robust(genotype_matrix, obj, method = "Burden", weights = weight_variant),
                  silent = TRUE)
        }
        if(length(err) == 1){
          robust_burden <- NA
          err_robust_burden <- err_robust_burden + 1
        }else{
          robust_burden <- tmp$p.value
        }
        end_time <- Sys.time()
        time_robust_burden <- as.double(difftime(end_time, start_time, units = "secs"))
        tmp_pvalue <- c(tmp_pvalue, robust_burden)
        tmp_computational_time <- c(tmp_computational_time, time_robust_burden)
        tmp_error <- c(tmp_error, err_robust_burden)
      }


      if(version == "all" | version == "optimal" | length(which(version == "robust_SKAT"))> 0){
        start_time <- Sys.time()
        if(length(which(is.na(weight_variant))) > 0){
          err <- try(tmp <-  SKAT::SKATBinary_Robust(genotype_matrix, obj, method = "SKAT"),
                    silent = TRUE)
        }else{
          err <- try(tmp <-  SKAT::SKATBinary_Robust(genotype_matrix, obj, method = "SKAT", weights = weight_variant),
                  silent = TRUE)
        }
        if(length(err) == 1){
          robust_SKAT <- NA
          err_robust_SKAT <- err_robust_SKAT + 1
        }else{
          robust_SKAT <- tmp$p.value
        }
        end_time <- Sys.time()
        time_robust_SKAT <- as.double(difftime(end_time, start_time, units = "secs"))
        tmp_pvalue <- c(tmp_pvalue, robust_SKAT)
        tmp_computational_time <- c(tmp_computational_time, time_robust_SKAT)
        tmp_error <- c(tmp_error, err_robust_SKAT)
      }


      if(version == "all" | version == "optimal" | length(which(version == "robust_SKATO"))> 0){
        start_time <- Sys.time()
        if(length(which(is.na(weight_variant))) > 0){
          err <- try(tmp <-  SKAT::SKATBinary_Robust(genotype_matrix, obj, method = "SKATO"),
                    silent = TRUE)
        }else{
          err <- try(tmp <-  SKAT::SKATBinary_Robust(genotype_matrix, obj, method = "SKATO", weights = weight_variant),
                  silent = TRUE)
        }
        if(length(err) == 1){
          robust_SKATO <- NA
          err_robust_SKATO <- err_robust_SKATO + 1
        }else{
          robust_SKATO <- tmp$p.value
        }
        end_time <- Sys.time()
        time_robust_SKATO <- as.double(difftime(end_time, start_time, units = "secs"))
        tmp_pvalue <- c(tmp_pvalue, robust_SKATO)
        tmp_computational_time <- c(tmp_computational_time, time_robust_SKATO)
        tmp_error <- c(tmp_error, err_robust_SKATO)
      }


      if(version == "all" | length(which(version == "p_linear_liu_burden"))> 0){
        start.time_p_linear_liu_burden <- Sys.time()
        err <- try(tmp <-  SKAT::SKAT(genotype_matrix, obj, kernel = "linear", method = "liu", r.corr = 1, weights = weight_variant),
                  silent = TRUE)
        if(length(err) == 1){
          p_linear_liu_burden <- NA
          err_p_linear_liu_burden <- err_p_linear_liu_burden + 1
        }else{
          p_linear_liu_burden <- tmp$p.value
        }
        end.time_p_linear_liu_burden <- Sys.time()
        time_p_linear_liu_burden <- as.double(difftime(end.time_p_linear_liu_burden,start.time_p_linear_liu_burden, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_linear_liu_burden)
        tmp_computational_time <- c(tmp_computational_time, time_p_linear_liu_burden)
        tmp_error <- c(tmp_error, err_p_linear_liu_burden)
      }


      if(version == "all" | length(which(version == "p_linear_weighted_liumod_burden"))> 0){
        start.time_p_linear_weighted_liumod_burden <- Sys.time()
        err <- try(tmp <-  SKAT::SKAT(genotype_matrix, obj, kernel = "linear.weighted", method = "liu.mod", r.corr = 1, weights = weight_variant),
                  silent = TRUE)
        if(length(err) == 1){
          p_linear_weighted_liumod_burden <- NA
          err_p_linear_weighted_liumod_burden <- err_p_linear_weighted_liumod_burden + 1
        }else{
          p_linear_weighted_liumod_burden <- tmp$p.value
        }
        end.time_p_linear_weighted_liumod_burden <- Sys.time()
        time_p_linear_weighted_liumod_burden <- as.double(difftime(end.time_p_linear_weighted_liumod_burden,start.time_p_linear_weighted_liumod_burden, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_linear_weighted_liumod_burden)
        tmp_computational_time <- c(tmp_computational_time, time_p_linear_weighted_liumod_burden)
        tmp_error <- c(tmp_error, err_p_linear_weighted_liumod_burden)
      }


      if(version == "all" | length(which(version == "p_linear_davies_skat"))> 0){
        #Linear kernel for variance-component test
        start.time_p_linear_davies_skat <- Sys.time()
        err <- try(tmp <-  SKAT::SKAT(genotype_matrix, obj, kernel = "linear", method = "davies", r.corr = 0, weights = weight_variant),
                  silent = TRUE)
        if(length(err) == 1){
          p_linear_davies_skat <- NA
          err_p_linear_davies_skat <- err_p_linear_davies_skat + 1
        }else{
          p_linear_davies_skat <- tmp$p.value
        }
        end.time_p_linear_davies_skat <- Sys.time()
        time_p_linear_davies_skat <- as.double(difftime(end.time_p_linear_davies_skat, start.time_p_linear_davies_skat, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_linear_davies_skat)
        tmp_computational_time <- c(tmp_computational_time, time_p_linear_davies_skat)
        tmp_error <- c(tmp_error, err_p_linear_davies_skat)
      }


      if(version == "all" | length(which(version == "p_linear_weighted_davies_skat"))> 0){
        start.time_p_linear_weighted_davies_skat <- Sys.time()
        err <- try(tmp <-  SKAT::SKAT(genotype_matrix, obj, kernel = "linear.weighted", method = "davies", r.corr = 0, weights = weight_variant),
                  silent = TRUE)
        if(length(err) == 1){
          p_linear_weighted_davies_skat <- NA
          err_p_linear_weighted_davies_skat <- err_p_linear_weighted_davies_skat + 1
        }else{
          p_linear_weighted_davies_skat <- tmp$p.value
        }
        end.time_p_linear_weighted_davies_skat <- Sys.time()
        time_p_linear_weighted_davies_skat <- as.double(difftime(end.time_p_linear_weighted_davies_skat, start.time_p_linear_weighted_davies_skat, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_linear_weighted_davies_skat)
        tmp_computational_time <- c(tmp_computational_time, time_p_linear_weighted_davies_skat)
        tmp_error <- c(tmp_error, err_p_linear_weighted_davies_skat)
      }


      if(version == "all" | length(which(version == "pBin_linear_IBS_SKAT_MA"))> 0){
        start.time_pBin_linear_IBS_SKAT_MA <- Sys.time()
        err <- try(tmp <- SKAT::SKATBinary(genotype_matrix, obj, kernel = "IBS", method = "SKAT", method.bin = "MA", weights = weight_variant),
                  silent = TRUE)
        if(length(err) == 1){
          pBin_linear_IBS_SKAT_MA <- NA
          err_pBin_linear_IBS_SKAT_MA <- err_pBin_linear_IBS_SKAT_MA + 1
        }else{
          pBin_linear_IBS_SKAT_MA <- tmp$p.value
        }
        end.time_pBin_linear_IBS_SKAT_MA <- Sys.time()
        time_pBin_linear_IBS_SKAT_MA <- as.double(difftime(end.time_pBin_linear_IBS_SKAT_MA, start.time_pBin_linear_IBS_SKAT_MA, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, pBin_linear_IBS_SKAT_MA)
        tmp_computational_time <- c(tmp_computational_time, time_pBin_linear_IBS_SKAT_MA)
        tmp_error <- c(tmp_error, err_pBin_linear_IBS_SKAT_MA)
      }


      if(version == "all" | length(which(version == "pBin_weighted_IBS_SKAT_MA"))> 0){
        start.time_pBin_weighted_IBS_SKAT_MA <- Sys.time()
        err <- try(tmp <- SKAT::SKATBinary(genotype_matrix, obj, kernel = "IBS.weighted", method = "SKAT", method.bin = "MA", weights = weight_variant),
                  silent = TRUE)
        if(length(err) == 1){
          pBin_weighted_IBS_SKAT_MA <- NA
          err_pBin_weighted_IBS_SKAT_MA <- err_pBin_weighted_IBS_SKAT_MA + 1
        }else{
          pBin_weighted_IBS_SKAT_MA <- tmp$p.value
        }
        end.time_pBin_weighted_IBS_SKAT_MA <- Sys.time()
        time_pBin_weighted_IBS_SKAT_MA <- as.double(difftime(end.time_pBin_weighted_IBS_SKAT_MA, start.time_pBin_weighted_IBS_SKAT_MA, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, pBin_weighted_IBS_SKAT_MA)
        tmp_computational_time <- c(tmp_computational_time, time_pBin_weighted_IBS_SKAT_MA)
        tmp_error <- c(tmp_error, err_pBin_weighted_IBS_SKAT_MA)
      }


      if(version == "all" | length(which(version == "pBin_2wayIX_SKAT_MA"))> 0){
        start.time_pBin_2wayIX_SKAT_MA <- Sys.time()
        err <- try(tmp <- SKAT::SKATBinary(genotype_matrix, obj, kernel = "2wayIX", method = "SKAT", method.bin = "MA", weights = weight_variant),
                  silent = TRUE)
        if(length(err) == 1){
          pBin_2wayIX_SKAT_MA <- NA
          err_pBin_2wayIX_SKAT_MA <- err_pBin_2wayIX_SKAT_MA + 1
        }else{
          pBin_2wayIX_SKAT_MA <- tmp$p.value
        }
        end.time_pBin_2wayIX_SKAT_MA <- Sys.time()
        time_pBin_2wayIX_SKAT_MA <- as.double(difftime(end.time_pBin_2wayIX_SKAT_MA, start.time_pBin_2wayIX_SKAT_MA, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, pBin_2wayIX_SKAT_MA)
        tmp_computational_time <- c(tmp_computational_time, time_pBin_2wayIX_SKAT_MA)
        tmp_error <- c(tmp_error, err_pBin_2wayIX_SKAT_MA)
      }


      if(version == "all" | length(which(version == "pBin_weighted_quadratic_SKAT_ERA"))> 0){
        start.time_pBin_weighted_quadratic_SKAT_ERA <- Sys.time()
        err <- try(tmp <- SKAT::SKATBinary(genotype_matrix, obj, kernel = "quadratic", method = "SKAT", method.bin = "ER.A", weights = weight_variant),
                  silent = TRUE)
        if(length(err) == 1){
          pBin_weighted_quadratic_SKAT_ERA <- NA
          err_pBin_weighted_quadratic_SKAT_ERA <- err_pBin_weighted_quadratic_SKAT_ERA + 1
        }else{
          pBin_weighted_quadratic_SKAT_ERA <- tmp$p.value
        }
        end.time_pBin_weighted_quadratic_SKAT_ERA <- Sys.time()
        time_pBin_weighted_quadratic_SKAT_ERA <- as.double(difftime(end.time_pBin_weighted_quadratic_SKAT_ERA, start.time_pBin_weighted_quadratic_SKAT_ERA, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, pBin_weighted_quadratic_SKAT_ERA)
        tmp_computational_time <- c(tmp_computational_time, time_pBin_weighted_quadratic_SKAT_ERA)
        tmp_error <- c(tmp_error, err_pBin_weighted_quadratic_SKAT_ERA)
      }


      if(version == "all" | length(which(version == "pBin_linear_SKATO_MA"))> 0){
        start.time_pBin_linear_SKATO_MA <- Sys.time()
        err <- try(tmp <- SKAT::SKATBinary(genotype_matrix, obj, kernel = "linear", method = "SKATO", method.bin = "MA", weights = weight_variant),
                  silent = TRUE)
        if(length(err) == 1){
          pBin_linear_SKATO_MA <- NA
          err_pBin_linear_SKATO_MA <- err_pBin_linear_SKATO_MA + 1
        }else{
          pBin_linear_SKATO_MA <- tmp$p.value
        }
        end.time_pBin_linear_SKATO_MA <- Sys.time()
        time_pBin_linear_SKATO_MA <- as.double(difftime(end.time_pBin_linear_SKATO_MA, start.time_pBin_linear_SKATO_MA, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, pBin_linear_SKATO_MA)
        tmp_computational_time <- c(tmp_computational_time, time_pBin_linear_SKATO_MA)
        tmp_error <- c(tmp_error, err_pBin_linear_SKATO_MA)
      }


      if(version == "all" | length(which(version == "pBin_linear_weighted_SKATO_MA"))> 0){
        start.time_pBin_linear_weighted_SKATO_MA <- Sys.time()
        err <- try(tmp <- SKAT::SKATBinary(genotype_matrix, obj, kernel = "linear.weighted", method = "SKATO", method.bin = "MA", weights = weight_variant),
                  silent = TRUE)
        if(length(err) == 1){
          pBin_linear_weighted_SKATO_MA <- NA
          err_pBin_linear_weighted_SKATO_MA <- err_pBin_linear_weighted_SKATO_MA + 1
        }else{
          pBin_linear_weighted_SKATO_MA <- tmp$p.value
        }
        end.time_pBin_linear_weighted_SKATO_MA <- Sys.time()
        time_pBin_linear_weighted_SKATO_MA <- as.double(difftime(end.time_pBin_linear_weighted_SKATO_MA, start.time_pBin_linear_weighted_SKATO_MA, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, pBin_linear_weighted_SKATO_MA)
        tmp_computational_time <- c(tmp_computational_time, time_pBin_linear_weighted_SKATO_MA)
        tmp_error <- c(tmp_error, err_pBin_linear_weighted_SKATO_MA)
      }


      if(version == "all" | version == "optimal" | length(which(version == "p_ascore"))> 0){
        #Adaptive Score test  by Hand and Pan (2010)
        start <- Sys.time()
        if(dim(genotype_matrix)[2] < 1000){
          err <- try(tmp <- AssotesteR::ASCORE(phenotype, genotype_matrix),
                    silent = TRUE)
        }else{
          err <- "too much variant"
        }
        if(length(err) == 1){
          p_ascore <- NA
          err_p_ascore <- err_p_ascore + 1
        }else{
          p_ascore <- tmp$perm.pval
        }
        stop <- Sys.time()
        time_p_ascore <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_ascore)
        tmp_computational_time <- c(tmp_computational_time, time_p_ascore)
        tmp_error <- c(tmp_error, err_p_ascore)
      }


      if(version == "all" | version == "optimal" | length(which(version == "p_ascore_ord"))> 0){
        #Ordered Adaptive Score test  by Hand and Pan (2010)
        start <- Sys.time()
        if(dim(genotype_matrix)[2] < max_var_per_test){
          err <- try(tmp <- AssotesteR::ASCORE.Ord(phenotype, genotype_matrix),
                    silent = TRUE)
        }else{
          err <- "too much variant"
        }
        if(length(err) == 1){
          p_ascore_ord <- NA
          err_p_ascore_ord <- err_p_ascore_ord + 1
        }else{
          p_ascore_ord <- tmp$perm.pval
        }
        stop <- Sys.time()
        time_p_ascore_ord <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_ascore_ord)
        tmp_computational_time <- c(tmp_computational_time, time_p_ascore_ord)
        tmp_error <- c(tmp_error, err_p_ascore_ord)
      }


      if(version == "all" | version == "optimal" | length(which(version == "p_assu"))> 0){
        #adaptive SSU test by Han and Pan (2010)
        start <- Sys.time()
        if(dim(genotype_matrix)[2] < max_var_per_test){
          err <- try(tmp <- AssotesteR::ASSU(phenotype, genotype_matrix),
                      silent = TRUE)
        }else{
          err <- "too much variant"
        }
        if(length(err) == 1){
          p_assu <- NA
          err_p_assu <- err_p_assu + 1
        }else{
          p_assu <- tmp$perm.pval
        }
        stop <- Sys.time()
        time_p_assu <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_assu)
        tmp_computational_time <- c(tmp_computational_time, time_p_assu)
        tmp_error <- c(tmp_error, err_p_assu)
      }


      if(version == "all" | version == "optimal" | length(which(version == "p_assu_ord"))> 0){
        #Ordered adaptive SSU test by Han and Pan (2010)
        start <- Sys.time()
        if(dim(genotype_matrix)[2] < max_var_per_test){
          err <- try(tmp <- AssotesteR::ASSU.Ord(phenotype, genotype_matrix),
                    silent = TRUE)
        }else{
          err <- "too much variant"
        }
        if(length(err) == 1){
          p_assu_ord <- NA
          err_p_assu_ord <- err_p_assu_ord + 1
        }else{
          p_assu_ord <- tmp$perm.pval
        }
        stop <- Sys.time()
        time_p_assu_ord <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_assu_ord)
        tmp_computational_time <- c(tmp_computational_time, time_p_assu_ord)
        tmp_error <- c(tmp_error, err_p_assu_ord)
      }


      if(version == "all" | version == "optimal" | length(which(version == "p_asum"))> 0){
        #The adaptive Adaptive Sum test by Han and Pan (2010)
        start <- Sys.time()
        if(dim(genotype_matrix)[2] < max_var_per_test){
          err <- try(tmp <- AssotesteR::ASUM(phenotype, genotype_matrix),
                    silent = TRUE)
        }else{
          err <- "too much variant"
        }
        if(length(err) == 1){
          p_asum <- NA
          err_p_asum <- err_p_asum + 1
        }else{
          p_asum <- tmp$perm.pval
        }
        stop <- Sys.time()
        time_p_asum <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_asum)
        tmp_computational_time <- c(tmp_computational_time, time_p_asum)
        tmp_error <- c(tmp_error, err_p_asum)
      }


      if(version == "all" | version == "optimal" | length(which(version == "p_asum_ord"))> 0){
        #Ordered adaptive Adaptive Sum test by Han and Pan (2010)
        start <- Sys.time()
        if(dim(genotype_matrix)[2] < max_var_per_test){
          err <- try(tmp <- AssotesteR::ASUM.Ord(phenotype, genotype_matrix),
                    silent = TRUE)
        }else{
          err <- "too much variant"
        }
        if(length(err) == 1){
          p_asum_ord <- NA
          err_p_asum_ord <- err_p_asum_ord + 1
        }else{
          p_asum_ord <- tmp$perm.pval
        }
        stop <- Sys.time()
        time_p_asum_ord <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_asum_ord)
        tmp_computational_time <- c(tmp_computational_time, time_p_asum_ord)
        tmp_error <- c(tmp_error, err_p_asum_ord)
      }


      if(version == "all" | version == "optimal" | length(which(version == "p_bst"))> 0){
        #Bayesian Score Test by Goeman et al (2005)
        start <- Sys.time()
        err <- try(tmp <- AssotesteR::BST(phenotype, genotype_matrix),
                  silent = TRUE)
        if(length(err) == 1){
          p_bst <- NA
          err_p_bst <- err_p_bst + 1
        }else{
          p_bst <- tmp$perm.pval
        }
        stop <- Sys.time()
        time_p_bst <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_bst)
        tmp_computational_time <- c(tmp_computational_time, time_p_bst)
        tmp_error <- c(tmp_error, err_p_bst)
      }


      if(version == "all" | version == "optimal" | length(which(version %in% c("p_calpha", "p_calpha_asymptopic")))> 0){
        #C-alpha Score Test by Neale et al (2011)
        start <- Sys.time()
        err <- try(tmp <- AssotesteR::CALPHA(phenotype, genotype_matrix, perm = 100),
                  silent = TRUE)
        if(length(err) == 1){
          p_calpha <- NA
          err_p_calpha <- err_p_calpha + 1
          p_calpha_asymptopic <- NA
          err_p_calpha_asymptopic <- err_p_calpha_asymptopic + 1
        }else{
          p_calpha <- tmp$perm.pval
          p_calpha_asymptopic <- tmp$asym.pval
        }
        stop <- Sys.time()
        time_p_calpha <- as.double(difftime(stop, start, units = "secs")) / 2
        time_p_calpha_asymptopic <- time_p_calpha

        if(version == "optimal" | length(which(version == c("p_calpha")))>0){
          tmp_pvalue <- c(tmp_pvalue, p_calpha)
          tmp_computational_time <- c(tmp_computational_time, time_p_calpha)
          tmp_error <- c(tmp_error, err_p_calpha)
        }else if(version == "all"){
          tmp_pvalue <- c(tmp_pvalue, p_calpha, p_calpha_asymptopic)
          tmp_computational_time <- c(tmp_computational_time, time_p_calpha, time_p_calpha_asymptopic)
          tmp_error <- c(tmp_error, err_p_calpha, err_p_calpha_asymptopic)
        }else if(length(which(version == c("p_calpha_asymptopic")))>0){
          tmp_pvalue <- c(tmp_pvalue, p_calpha_asymptopic)
          tmp_computational_time <- c(tmp_computational_time, time_p_calpha_asymptopic)
          tmp_error <- c(tmp_error, err_p_calpha_asymptopic)
        }
      }


      if(version == "all" | length(which(version == "p_carv_hard"))> 0){
        #Comprehrensive Approach to Analyzing Rare Variants by Hoffmann et al (2010)
        start <- Sys.time()
        err <- try(tmp <- AssotesteR::CARV(phenotype, genotype_matrix, waf = TRUE, signs = TRUE, approach = "hard", maf = rare_maf_threshold),
                  silent = TRUE)
        if(length(err) == 1){
          p_carv_hard <- NA
          err_p_carv_hard <- err_p_carv_hard + 1
        }else{
          p_carv_hard <- tmp$perm.pval
        }
        stop <- Sys.time()
        time_p_carv_hard <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_carv_hard)
        tmp_computational_time <- c(tmp_computational_time, time_p_carv_hard)
        tmp_error <- c(tmp_error, err_p_carv_hard)
      }


      if(version == "all" | length(which(version == "p_carv_variable"))> 0){
        start <- Sys.time()
        err <- try(tmp <- AssotesteR::CARV(phenotype, genotype_matrix, waf = TRUE, signs = TRUE, approach = "variable", maf = rare_maf_threshold),
                  silent = TRUE)
        if(length(err) == 1){
          p_carv_variable <- NA
          err_p_carv_variable <- err_p_carv_variable + 1
        }else{
          p_carv_variable <- tmp$perm.pval
        }
        stop <- Sys.time()
        time_p_carv_variable <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_carv_variable)
        tmp_computational_time <- c(tmp_computational_time, time_p_carv_variable)
        tmp_error <- c(tmp_error, err_p_carv_variable)
      }


      if(version == "all" | length(which(version == "p_carv_stepup"))> 0){
        start <- Sys.time()
        err <- try(tmp <- AssotesteR::CARV(phenotype, genotype_matrix, waf = TRUE, signs = TRUE, approach = "stepup", maf = rare_maf_threshold),
                  silent = TRUE)
        if(length(err) == 1){
          p_carv_stepup <- NA
          err_p_carv_stepup <- err_p_carv_stepup + 1
        }else{
          p_carv_stepup <- tmp$perm.pval
        }
        stop <- Sys.time()
        time_p_carv_stepup <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_carv_stepup)
        tmp_computational_time <- c(tmp_computational_time, time_p_carv_stepup)
        tmp_error <- c(tmp_error, err_p_carv_stepup)
      }


      if(version == "all" | length(which(version == "p_cast_fisher"))> 0){
        #Cohort Allelic Sums Test by S. Morgenthaler et al (2007)
        start <- Sys.time()
        err <- try(tmp <- AssotesteR::CAST(phenotype, genotype_matrix, maf = rare_maf_threshold, test = "fisher"),
                  silent = TRUE)
        if(length(err) == 1){
          p_cast_fisher <- NA
          err_p_cast_fisher <- err_p_cast_fisher + 1
        }else{
          p_cast_fisher <- tmp$asym.pval
        }
        stop <- Sys.time()
        time_p_cast_fisher <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_cast_fisher)
        tmp_computational_time <- c(tmp_computational_time, time_p_cast_fisher)
        tmp_error <- c(tmp_error, err_p_cast_fisher)
      }


      if(version == "all" | version == "optimal" | length(which(version == "p_cast_chisq"))> 0){
        start <- Sys.time()
        err <- try(tmp <- AssotesteR::CAST(phenotype, genotype_matrix, maf = rare_maf_threshold, test = "chisq"),
                  silent = TRUE)
        if(length(err) == 1){
          p_cast_chisq <- NA
          err_p_cast_chisq <- err_p_cast_chisq + 1
        }else{
          p_cast_chisq <- tmp$asym.pval
        }
        stop <- Sys.time()
        time_p_cast_chisq <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_cast_chisq)
        tmp_computational_time <- c(tmp_computational_time, time_p_cast_chisq)
        tmp_error <- c(tmp_error, err_p_cast_chisq)
      }


      if(version == "all" | version == "optimal" | length(which(version == "p_cmat"))> 0){
        #Cumulative Minor Allele Test by Zawistowski et al (2010)
        start <- Sys.time()
        err <- try(tmp <- AssotesteR::CMAT(phenotype, genotype_matrix, maf = rare_maf_threshold, weights = weight_variant),
                  silent = TRUE)
        if(length(err) == 1){
          p_cmat <- NA
          err_p_cmat <- err_p_cmat + 1
        }else{
          p_cmat <- tmp$perm.pval
        }
        stop <- Sys.time()
        time_p_cmat <- as.double(difftime(stop, start, units = "secs"))
        tmp_pvalue <- c(tmp_pvalue, p_cmat)
        tmp_computational_time <- c(tmp_computational_time, time_p_cmat)
        tmp_error <- c(tmp_error, err_p_cmat)
      }


      if(version == "all" | version == "optimal" | length(which(version == "p_cmc"))> 0){
        #Combined Multivariate and Collapsing Method by Li and Leal (2008)
        start <- Sys.time()
        err <- try(tmp <- AssotesteR::CMC(phenotype, genotype_matrix, maf = rare_maf_threshold),
                  silent = TRUE)
        if(length(err) == 1){
          p_cmc <- NA
          err_p_cmc <- err_p_cmc + 1
        }else{
          p_cmc <- tmp$perm.pval
        }
        stop <- Sys.time()
        time_p_cmc <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_cmc)
        tmp_computational_time <- c(tmp_computational_time, time_p_cmc)
        tmp_error <- c(tmp_error, err_p_cmc)
      }


      if(version == "all" | version == "optimal" | length(which(version == "p_rbt"))> 0){
        #Replication Based Test by Ionita-Laza et al (2011)
        start <- Sys.time()
        err <- try(tmp <- AssotesteR::RBT(phenotype, genotype_matrix),
                  silent = TRUE)
        if(length(err) == 1){
          p_rbt <- NA
          err_p_rbt <- err_p_rbt + 1
        }else{
          p_rbt <- tmp$perm.pval
        }
        stop <- Sys.time()
        time_p_rbt <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_rbt)
        tmp_computational_time <- c(tmp_computational_time, time_p_rbt)
        tmp_error <- c(tmp_error, err_p_rbt)
      }


      if(version == "all" | version == "optimal" | length(which(version == "p_rvt1"))> 0){
        #Rare Variant Test 1 for dichotomous traits by Morris and Zeggini (2010)
        start <- Sys.time()
        err <- try(tmp <- AssotesteR::RVT1(phenotype, genotype_matrix, maf = rare_maf_threshold),
                  silent = TRUE)
        if(length(err) == 1){
          p_rvt1 <- NA
          err_p_rvt1 <- err_p_rvt1 + 1
        }else{
          p_rvt1 <- tmp$asym.pval
        }
        stop <- Sys.time()
        time_p_rvt1 <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_rvt1)
        tmp_computational_time <- c(tmp_computational_time, time_p_rvt1)
        tmp_error <- c(tmp_error, err_p_rvt1)
      }


      if(version == "all" | version == "optimal" | length(which(version == "p_rvt2"))> 0){
        #Rare Variant Test 2 for dichotomous traits by Morris and Zeggini (2010)
        start <- Sys.time()
        err <- try(tmp <- AssotesteR::RVT2(phenotype, genotype_matrix, maf = rare_maf_threshold),
                  silent = TRUE)
        if(length(err) == 1){
          p_rvt2 <- NA
          err_p_rvt2 <- err_p_rvt2 + 1
        }else{
          p_rvt2 <- tmp$asym.pval
        }
        stop <- Sys.time()
        time_p_rvt2 <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_rvt2)
        tmp_computational_time <- c(tmp_computational_time, time_p_rvt2)
        tmp_error <- c(tmp_error, err_p_rvt2)
      }


      if(version == "all" | version == "optimal" | length(which(version == "p_rwas"))> 0){
        #Rare-Variant Weighted Aggregate Statistic by Sul et al (2011)
        start <- Sys.time()
        err <- try(tmp <- AssotesteR::RWAS(phenotype, genotype_matrix, maf = rare_maf_threshold, perm = 100),
                  silent = TRUE)
        if(length(err) == 1){
          p_rwas <- NA
          err_p_rwas <- err_p_rwas + 1
        }else{
          p_rwas <- tmp$asym.pval
        }
        stop <- Sys.time()
        time_p_rwas <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_rwas)
        tmp_computational_time <- c(tmp_computational_time, time_p_rwas)
        tmp_error <- c(tmp_error, err_p_rwas)
      }


      if(version == "all" | version == "optimal" | length(which(version %in% c("p_score", "p_score_asymptopic")))> 0){
        #Score Test (from Logistic Regression) by Chapman J et al (2008)
        start <- Sys.time()
        if(dim(genotype_matrix)[2] < max_var_per_test){
          err <- try(tmp <- AssotesteR::SCORE(phenotype, genotype_matrix),
                    silent = TRUE)
        }else{
          err <- "too much variant"
        }
        if(length(err) == 1){
          p_score <- NA
          err_p_score <- err_p_score + 1
          p_score_asymptopic <- NA
          err_p_score_asymptopic <- err_p_score_asymptopic + 1
        }else{
          p_score <- tmp$perm.pval
          p_score_asymptopic <- tmp$asym.pval
        }
        stop <- Sys.time()
        time_p_score <- as.double(difftime(stop, start, units = "secs")) / 2
        time_p_score_asymptopic <- time_p_score

        if(version == "all" | version == "optimal"){
          tmp_pvalue <- c(tmp_pvalue, p_score, p_score_asymptopic)
          tmp_computational_time <- c(tmp_computational_time, time_p_score, time_p_score_asymptopic)
          tmp_error <- c(tmp_error, err_p_score, err_p_score_asymptopic)
        }else if(length(which(version == "p_score"))> 0){
          tmp_pvalue <- c(tmp_pvalue, p_score)
          tmp_computational_time <- c(tmp_computational_time, time_p_score)
          tmp_error <- c(tmp_error, err_p_score)
        }else if(length(which(version == "p_score_asymptopic"))> 0){
          tmp_pvalue <- c(tmp_pvalue, p_score_asymptopic)
          tmp_computational_time <- c(tmp_computational_time, time_p_score_asymptopic)
          tmp_error <- c(tmp_error, err_p_score_asymptopic)
        }
      }


      if(version == "all" | version == "optimal" | length(which(version %in% c("p_ssu", "p_ssu_asymptopic")))> 0){
        #Sum of Squared Score U Statistic by Pan (2009)
        start <- Sys.time()
        if(dim(genotype_matrix)[2] < max_var_per_test){
          err <- try(tmp <- AssotesteR::SSU(phenotype, genotype_matrix),
                    silent = TRUE)
        }else{
          err <- "too much variant"
        }
        if(length(err) == 1){
          p_ssu <- NA
          err_p_ssu <- err_p_ssu + 1
          p_ssu_asymptopic <- NA
          err_p_ssu_asymptopic <- err_p_ssu_asymptopic + 1
        }else{
          p_ssu <- tmp$perm.pval
          p_ssu_asymptopic <- tmp$asym.pval
        }
        stop <- Sys.time()
        time_p_ssu <- as.double(difftime(stop, start, units = "secs")) / 2
        time_p_ssu_asymptopic <- time_p_ssu

        if(version == "all" | version == "optimal"){
          tmp_pvalue <- c(tmp_pvalue, p_ssu, p_ssu_asymptopic)
          tmp_computational_time <- c(tmp_computational_time, time_p_ssu, time_p_ssu_asymptopic)
          tmp_error <- c(tmp_error, err_p_ssu, err_p_ssu_asymptopic)
        }else if(length(which(version == "p_ssu"))> 0){
          tmp_pvalue <- c(tmp_pvalue, p_ssu)
          tmp_computational_time <- c(tmp_computational_time, time_p_ssu)
          tmp_error <- c(tmp_error, err_p_ssu)
        }else if(length(which(version == "p_ssu_asymptopic"))> 0){
          tmp_pvalue <- c(tmp_pvalue, p_ssu_asymptopic)
          tmp_computational_time <- c(tmp_computational_time, time_p_ssu_asymptopic)
          tmp_error <- c(tmp_error, err_p_ssu_asymptopic)
        }
      }


      if(version == "all" | version == "optimal" | length(which(version %in% c("p_ssuw", "p_ssuw_asymptopic")))> 0){
        #Weighted Sum of Squared Score U Statistic by Pan (2009)
        start <- Sys.time()
        err <- try(tmp <- AssotesteR::SSUW(phenotype, genotype_matrix),
                  silent = TRUE)
        if(length(err) == 1){
          p_ssuw <- NA
          err_p_ssuw <- err_p_ssuw + 1
          p_ssuw_asymptopic <- NA
          err_p_ssuw_asymptopic <- err_p_ssuw_asymptopic + 1
        }else{
          p_ssuw <- tmp$perm.pval
          p_ssuw_asymptopic <- tmp$asym.pval
        }
        stop <- Sys.time()
        time_p_ssuw <- as.double(difftime(stop, start, units = "secs")) / 2
        time_p_ssuw_asymptopic <- time_p_ssuw

        if(version == "all" | version == "optimal"){
          tmp_pvalue <- c(tmp_pvalue, p_ssuw, p_ssuw_asymptopic)
          tmp_computational_time <- c(tmp_computational_time, time_p_ssuw, time_p_ssuw_asymptopic)
          tmp_error <- c(tmp_error, err_p_ssuw, err_p_ssuw_asymptopic)
        }else if(length(which(version == "p_ssuw"))> 0){
          tmp_pvalue <- c(tmp_pvalue, p_ssuw)
          tmp_computational_time <- c(tmp_computational_time, time_p_ssuw)
          tmp_error <- c(tmp_error, err_p_ssuw)
        }else if(length(which(version == "p_ssuw_asymptopic"))> 0){
          tmp_pvalue <- c(tmp_pvalue, p_ssuw_asymptopic)
          tmp_computational_time <- c(tmp_computational_time, time_p_ssuw_asymptopic)
          tmp_error <- c(tmp_error, err_p_ssuw_asymptopic)
        }
      }


      if(version == "all" | version == "optimal" | length(which(version == "p_ttest_asymptopic"))> 0){
        #Hotelling T2 Test by Xiong et al (2002)
        start <- Sys.time()
        err <- try(tmp <- AssotesteR::TTEST(phenotype, genotype_matrix),
                  silent = TRUE)
        if(length(err) == 1){
          p_ttest_asymptopic <- NA
          err_p_ttest_asymptopic <- err_p_ttest_asymptopic + 1
        }else{
          p_ttest_asymptopic <- tmp$asym.pval
        }
        stop <- Sys.time()
        time_p_ttest_asymptopic <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_ttest_asymptopic)
        tmp_computational_time <- c(tmp_computational_time, time_p_ttest_asymptopic)
        tmp_error <- c(tmp_error, err_p_ttest_asymptopic)
      }


      if(version == "all" | version == "optimal" | length(which(version == "p_vt"))> 0){
        #Variable Threshold by Price et al (2010)
        start <- Sys.time()
        err <- try(tmp <- AssotesteR::VT(phenotype, genotype_matrix, maf = rare_maf_threshold),
                  silent = TRUE)
        if(length(err) == 1){
          p_vt <- NA
          err_p_vt <- err_p_vt + 1
        }else{
          p_vt <- tmp$perm.pval
        }
        stop <- Sys.time()
        time_p_vt <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_vt)
        tmp_computational_time <- c(tmp_computational_time, time_p_vt)
        tmp_error <- c(tmp_error, err_p_vt)
      }


      if(version == "all" | version == "optimal" | length(which(version == "p_wss"))> 0){
        #Weighted Sum Statistic by Madsen and Browning (2009)
        start <- Sys.time()
        err <- try(tmp <- AssotesteR::WSS(phenotype, genotype_matrix),
                  silent = TRUE)
        if(length(err) == 1){
          p_wss <- NA
          err_p_wss <- err_p_wss + 1
        }else{
          p_wss <- tmp$perm.pval
        }
        stop <- Sys.time()
        time_p_wss <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_wss)
        tmp_computational_time <- c(tmp_computational_time, time_p_wss)
        tmp_error <- c(tmp_error, err_p_wss)
      }


      if(version == "all" | version == "optimal" | length(which(version %in% c("p_wst", "p_wst_asymptopic")))> 0){
        #Weighted Score Test by Wang and Elston (2007)
        start <- Sys.time()
        err <- try(tmp <- AssotesteR::WST(phenotype, genotype_matrix),
                  silent = TRUE)
        if(length(err) == 1){
          p_wst <- NA
          err_p_wst <- err_p_wst + 1
          p_wst_asymptopic <- NA
          err_p_wst_asymptopic <- err_p_wst_asymptopic + 1
        }else{
          p_wst <- tmp$perm.pval
          p_wst_asymptopic <- tmp$asym.pval
        }
        stop <- Sys.time()
        time_p_wst <- as.double(difftime(stop, start, units = "secs")) / 2
        time_p_wst_asymptopic <- time_p_wst

        if(version == "all"){
          tmp_pvalue <- c(tmp_pvalue, p_wst, p_wst_asymptopic)
          tmp_computational_time <- c(tmp_computational_time, time_p_wst, time_p_wst_asymptopic)
          tmp_error <- c(tmp_error, err_p_wst, err_p_wst_asymptopic)
        }else if(length(which(version == "p_wst"))> 0 | version == "optimal"){
          tmp_pvalue <- c(tmp_pvalue, p_wst)
          tmp_computational_time <- c(tmp_computational_time, time_p_wst)
          tmp_error <- c(tmp_error, err_p_wst)
        }else if(length(which(version == "p_wst_asymptopic"))> 0){
          tmp_pvalue <- c(tmp_pvalue, p_wst_asymptopic)
          tmp_computational_time <- c(tmp_computational_time, time_p_wst_asymptopic)
          tmp_error <- c(tmp_error, err_p_wst_asymptopic)
        }
      }


      if(version == "all" | length(which(version == "p_catt"))> 0){
        #CATT by Zhicheng Du et al (2017)
        start <- Sys.time()
        err <- try(tmp <- CATT::CATT(table = catt_matrix),
                  silent = TRUE)
        if(length(err) == 1){
          p_catt <- NA
          err_p_catt <- err_p_catt + 1
        }else{
          p_catt <- tmp$p.value
        }
        stop <- Sys.time()
        time_p_catt <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_catt)
        tmp_computational_time <- c(tmp_computational_time, time_p_catt)
        tmp_error <- c(tmp_error, err_p_catt)
      }


      if(version == "all" | version == "optimal" | length(which(version == "p_KAT"))> 0){
        #Conditional Inference for the Kernel Association Test by Wang, K. (2016)
        start <- Sys.time()
        err <- try(tmp <- iGasso::KAT.coin(y = phenotype,
                                          G = genotype_matrix,
                                          X = covariate,
                                          out_type = "D"),
                  silent = TRUE)
        if(length(err) == 1 | tmp$p.value < 0){
          p_KAT <- NA
          err_p_KAT <- err_p_KAT + 1
        }else{
          p_KAT <- tmp$p.value
        }
        stop <- Sys.time()
        time_p_KAT <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_KAT)
        tmp_computational_time <- c(tmp_computational_time, time_p_KAT)
        tmp_error <- c(tmp_error, err_p_KAT)
      }


      if(version == "all" | length(which(version == "p_SKATplus"))> 0){
        #enhanced power over SKAT SKATplus by Wang, K. (2016)
        start <- Sys.time()
        err <- try(tmp <- iGasso::SKATplus(y = phenotype,
                                          G = genotype_matrix,
                                          X = covariate,
                                          out_type = "D",
                                          tau = 1),
                  silent = TRUE)
        if(length(err) == 1){
          p_SKATplus <- NA
          err_p_SKATplus <- err_p_SKATplus + 1
        }else{
          p_SKATplus <- tmp$p.value
        }
        stop <- Sys.time()
        time_p_SKATplus <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_SKATplus)
        tmp_computational_time <- c(tmp_computational_time, time_p_SKATplus)
        tmp_error <- c(tmp_error, err_p_SKATplus)
      }


      if(version == "all" | length(which(version == "p_wgscan_region"))> 0){
        #WGS-scan score type statistics by Zihuai et al (2019)
        start <- Sys.time()
        err <- try(tmp <- WGScan::WGScan.Region(result.prelim = obj_wgscan,
                                                G = genotype_matrix,
                                                pos = position,
                                                MAF.threshold = rare_maf_threshold),
                  silent = TRUE)
        if(length(err) == 1){
          p_wgscan_region <- NA
          err_p_wgscan_region <- err_p_wgscan_region + 1
        }else{
          p_wgscan_region <- tmp$p.value
        }
        stop <- Sys.time()
        time_p_wgscan_region <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_wgscan_region)
        tmp_computational_time <- c(tmp_computational_time, time_p_wgscan_region)
        tmp_error <- c(tmp_error, err_p_wgscan_region)
      }


      if(version == "all" | version == "optimal" | length(which(version == "p_rebet"))> 0){
        #REBET (subREgion-based BurdEn Test) by Bin Zhu et al (2018)
        start <- Sys.time()
        err <- try(tmp <- REBET::rebet(response = phenotype,
                                      genotypes = genotype_matrix,
                                      subRegions = rep("1", dim(genotype_matrix)[2]),
                                      responseType = "binary"),
                  silent = TRUE)
        if(length(err) == 1){
          p_rebet <- NA
          err_p_rebet <- err_p_rebet + 1
        }else{
          p_rebet <- as.double(tmp$Meta$pval)
        }
        stop <- Sys.time()
        time_p_rebet <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_rebet)
        tmp_computational_time <- c(tmp_computational_time, time_p_rebet)
        tmp_error <- c(tmp_error, err_p_rebet)
      }


      if(version == "all" | length(which(version == "p_pcr"))> 0){
        #PCR Principal Components Regression for RV tests by C. Xu et al (2012)
        start <- Sys.time()
        err <- try(tmp <- RVtests::PCR(x = genotype_matrix,
                                      y = phenotype),
                  silent = TRUE)
        if(length(err) == 1){
          p_pcr <- NA
          err_p_pcr <- err_p_pcr + 1
        }else{
          p_pcr <- as.double(tmp$pvalue.empirical)
        }
        stop <- Sys.time()
        time_p_pcr <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_pcr)
        tmp_computational_time <- c(tmp_computational_time, time_p_pcr)
        tmp_error <- c(tmp_error, err_p_pcr)
      }


      if(version == "all" | version == "optimal" | length(which(version == "p_rr"))> 0){
        #Ridge Regression for RV Tests by C. Xu et al (2012)
        start <- Sys.time()
        err <- try(tmp <- RVtests::RR(x = genotype_matrix,
                                      y = phenotype,
                                      z = covariate,
                                      weights = as.character(weight_variant)),
                  silent = TRUE)
        if(length(err) == 1){
          p_rr <- NA
          err_p_rr <- err_p_rr + 1
        }else{
          p_rr <- as.double(tmp$pvalue.empirical)
        }
        stop <- Sys.time()
        time_p_rr <- as.double(difftime(stop, start, units = "secs"))
        tmp_pvalue <- c(tmp_pvalue, p_rr)
        tmp_computational_time <- c(tmp_computational_time, time_p_rr)
        tmp_error <- c(tmp_error, err_p_rr)
      }


      if(version == "all" | version == "optimal" | length(which(version == "p_spls"))> 0){
        #Sparse PLS for RV Tests by C. Xu et al (2012)
        start <- Sys.time()
        err <- try(tmp <- RVtests::SPLS(x = genotype_matrix,
                                        y = phenotype,
                                        npermutation = 100),
                  silent = TRUE)
        if(length(err) == 1){
          p_spls <- NA
          err_p_spls <- err_p_spls + 1
        }else{
          p_spls <- as.double(tmp$pvalue.empirical)
        }
        stop <- Sys.time()
        time_p_spls <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_spls)
        tmp_computational_time <- c(tmp_computational_time, time_p_spls)
        tmp_error <- c(tmp_error, err_p_spls)
      }


      if(version == "all" | version == "optimal" | length(which(version %in% c("p_t1p", "p_t5p", "p_wep", "p_score_vtp", "p_wod01", "p_wod05")))> 0){
        #SVT and WOD for RV Tests by C. Xu et al (2012)
        start <- Sys.time()
        err <- try(tmp <- RVtests::VTWOD(x = genotype_matrix,
                                        y = phenotype),
                  silent = TRUE)
        if(length(err) == 1){
          p_t1p <- NA
          err_p_t1p <- err_p_t1p + 1
          p_t5p <- NA
          err_p_t5p <- err_p_t5p + 1
          p_wep <- NA
          err_p_wep <- err_p_wep + 1
          p_score_vtp <- NA
          err_p_score_vtp <- err_p_score_vtp + 1
          p_wod01 <- NA
          err_p_wod01 <- err_p_wod01 + 1
          p_wod05 <- NA
          err_p_wod05 <- err_p_wod05 + 1
        }else{
          p_t1p <- as.double(tmp$pvalue.empirical[2])
          p_t5p <- as.double(tmp$pvalue.empirical[4])
          p_wep <- as.double(tmp$pvalue.empirical[6])
          p_score_vtp <- as.double(tmp$pvalue.empirical[8])
          p_wod01 <- as.double(tmp$pvalue.empirical[9])
          p_wod05 <- as.double(tmp$pvalue.empirical[10])
        }
        stop <- Sys.time()
        time_p_t1p <- as.double(difftime(stop, start, units = "secs")) / 6
        time_p_t5p <- time_p_t1p
        time_p_wep <- time_p_t1p
        time_p_score_vtp <- time_p_t1p
        time_p_wod01 <- time_p_t1p
        time_p_wod05 <- time_p_t1p

        if(version == "all"){
          tmp_pvalue <- c(tmp_pvalue, p_t1p, p_t5p, p_wep, p_score_vtp, p_wod01, p_wod05)
          tmp_computational_time <- c(tmp_computational_time, time_p_t1p, time_p_t5p, time_p_wep, time_p_score_vtp, time_p_wod01, time_p_wod05)
          tmp_error <- c(tmp_error, err_p_t1p, err_p_t5p, err_p_wep, err_p_score_vtp, err_p_wod01, err_p_wod05)
        }else if(version == "optimal"){
          tmp_pvalue <- c(tmp_pvalue, p_t1p, p_t5p, p_wep, p_score_vtp)
          tmp_computational_time <- c(tmp_computational_time, time_p_t1p, time_p_t5p, time_p_wep, time_p_score_vtp)
          tmp_error <- c(tmp_error, err_p_t1p, err_p_t5p, err_p_wep, err_p_score_vtp)
        }
        if(length(which(version == "p_t1p"))> 0){
          tmp_pvalue <- c(tmp_pvalue, p_t1p)
          tmp_computational_time <- c(tmp_computational_time, time_p_t1p)
          tmp_error <- c(tmp_error, err_p_t1p)
        }
        if(length(which(version == "p_t5p"))> 0){
          tmp_pvalue <- c(tmp_pvalue, p_t5p)
          tmp_computational_time <- c(tmp_computational_time, time_p_t5p)
          tmp_error <- c(tmp_error, err_p_t5p)
        }
        if(length(which(version == "p_wep"))> 0){
          tmp_pvalue <- c(tmp_pvalue, p_wep)
          tmp_computational_time <- c(tmp_computational_time, time_p_wep)
          tmp_error <- c(tmp_error, err_p_wep)
        }
        if(length(which(version == "p_score_vtp"))> 0){
          tmp_pvalue <- c(tmp_pvalue, p_score_vtp)
          tmp_computational_time <- c(tmp_computational_time, time_p_score_vtp)
          tmp_error <- c(tmp_error, err_p_score_vtp)
        }
        if(length(which(version == "p_wod01"))> 0){
          tmp_pvalue <- c(tmp_pvalue, p_wod01)
          tmp_computational_time <- c(tmp_computational_time, time_p_wod01)
          tmp_error <- c(tmp_error, err_p_wod01)
        }
        if(length(which(version == "p_wod05"))> 0){
          tmp_pvalue <- c(tmp_pvalue, p_wod05)
          tmp_computational_time <- c(tmp_computational_time, time_p_wod05)
          tmp_error <- c(tmp_error, err_p_wod05)
        }
      }


      if(version == "all" | length(which(version == "p_ada"))> 0){
        #ADA Tests by Lin W-Y (2016)
        start <- Sys.time()
        err <- try(tmp <- ADATest(genotype = genotype_matrix,
                                  phenotype = phenotype,
                                  mafThr = rare_maf_threshold),
                  silent = TRUE)
        if(length(err) == 1){
          p_ada <- NA
          err_p_ada <- err_p_ada + 1
        }else{
          p_ada <- tmp$pval
        }
        stop <- Sys.time()
        time_p_ada <- as.double(difftime(stop, start, units = "secs"))

        tmp_pvalue <- c(tmp_pvalue, p_ada)
        tmp_computational_time <- c(tmp_computational_time, time_p_ada)
        tmp_error <- c(tmp_error, err_p_ada)
      }


      ########################  Store results  ########################
      pvalue <- data.frame(pvalue = tmp_pvalue)
      pvalue[which(pvalue[,1] == 0),1] <- NA
      computational_time <- data.frame(computational_time = tmp_computational_time)
      errors <- data.frame(nbr_errors = tmp_error)

      #test_in_excalibur <- stat_Framework(get_test = TRUE, version = "all")[-1]
      #test_in_testing_framework <- stat_Framework(get_test = TRUE, version = "all")[-1]
      #idx_excalibur <- which(test_in_testing_framework %in% test_in_excalibur)
      if(version == "all" | version == "optimal"){
        #Excalibur_baseline or Excalibur
        qvalue <- p.adjust(pvalue[,1], method = p.adjust.methods[4]) # Bonferonni
        value_output <- min(qvalue, na.rm = TRUE)
        if(value_output == Inf){
          #If all test have failed to return a pvalue
          pvalue <- rbind(data.frame(pvalue = NA), pvalue)
          errors <- rbind(data.frame(nbr_errors = 1), errors)
        }else{
          Excalibur <- value_output
          pvalue <- rbind(data.frame(pvalue = Excalibur), pvalue)
          errors <- rbind(data.frame(nbr_errors = 0), errors)
        }
        time_Excalibur <- sum(computational_time[,1])
        computational_time <- rbind(data.frame(computational_time = time_Excalibur), computational_time)
      }


    } else {
      ###Continuous case
      #Do the null model
      if(length(covariate) == 1){
        obj <-SKAT_Null_Model(phenotype~ 1, data = data.frame(genotype_matrix), out_type = "C")
      }else{
        obj <-SKAT_Null_Model(phenotype~ covariate, data = data.frame(genotype_matrix), out_type = "C")
      }

      ########################  Perform all tests  ########################
      #Excalibur
      start_Excalibur <- Sys.time()


      #Linear kernel for Burden test using davies method for p value computation
      start <- Sys.time()
      err <- try(tmp <- SKAT::SKAT(genotype_matrix, obj, kernel = "linear", method = "davies", r.corr = 1),
                 silent = TRUE)
      if(length(err) == 1){
        p_linear_davies_burden <- NA
      }else{
        p_linear_davies_burden <- tmp$p.value
      }
      stop <- Sys.time()
      time_p_linear_davies_burden <- stop - start


      #Linear weighted kernel for Burden test using liu method for p value computation
      start <- Sys.time()
      err <- try(tmp <- SKAT(genotype_matrix, obj, kernel = "linear.weighted", method = "liu", r.corr = 1),
                 silent = TRUE)
      if(length(err) == 1){
        p_linear_weighted_liu_burden <- NA
      }else{
        p_linear_weighted_liu_burden <- tmp$p.value
      }
      stop <- Sys.time()
      time_p_linear_weighted_liu_burden <- stop - start


      #Linear weighted kernel for Variance-Component test using modified liu method for p value computation
      start <- Sys.time()
      err <- try(tmp <- SKAT::SKAT(genotype_matrix, obj, kernel = "linear.weighted", method = "liu.mod", r.corr = 0),
                 silent = TRUE)
      if(length(err) == 1){
        p_linear_weighted_liumod_skat <- NA
      }else{
        p_linear_weighted_liumod_skat <- tmp$p.value
      }
      stop <- Sys.time()
      time_p_linear_weighted_liumod_skat <- stop - start

      #Linear kernel for Variance-Component test using liu method for p value computation
      start <- Sys.time()
      err <- try(tmp <- SKAT::SKAT(genotype_matrix, obj, kernel = "linear", method = "liu", r.corr = 0),
                 silent = TRUE)
      if(length(err) == 1){
        p_linear_liu_skat <- NA
      }else{
        p_linear_liu_skat <- tmp$p.value
      }
      stop <- Sys.time()
      time_p_linear_liu_skat <- stop - start

      #Linear kernel for Omnibus test using optimal method for p value computation
      start <- Sys.time()
      err <- try(tmp <- SKAT::SKAT(genotype_matrix, obj, kernel = "linear", method = "optimal.adj"),
                 silent = TRUE)
      if(length(err) == 1){
        p_linear_skato <- NA
      }else{
        p_linear_skato <- tmp$p.value
      }
      stop <- Sys.time()
      time_p_linear_skato <- stop - start

      #Linear weighted kernel for Omnibus test using optimal method for p value computation
      start <- Sys.time()
      err <- try(tmp <- SKAT::SKAT(genotype_matrix, obj, kernel = "linear.weighted", method = "optimal.adj"),
                 silent = TRUE)
      if(length(err) == 1){
        p_linear_weighted_skato <- NA
      }else{
        p_linear_weighted_skato <- tmp$p.value
      }
      stop <- Sys.time()
      time_p_linear_weighted_skato <- stop - start

      ########################  Store results  ########################
      pvalue <- data.frame()
      pvalue <- data.frame(
        pvalue = c(p_linear_davies_burden,
                   p_linear_weighted_liu_burden,
                   p_linear_weighted_liumod_skat,
                   p_linear_liu_skat,
                   p_linear_skato,
                   p_linear_weighted_skato)
      )

      pvalue[which(pvalue[,1] == 0),1] <- NA

      #Excalibur
      test_in_excalibur <- stat_Framework(Binary = FALSE, get_test = TRUE)[-1]
      test_in_testing_framework <- stat_Framework(Binary = FALSE, get_test = TRUE)[-1]
      idx_excalibur <- which(test_in_testing_framework %in% test_in_excalibur)
      qvalue <- p.adjust(pvalue[idx_excalibur,1], method = p.adjust.methods[5]) # B-H
      value_output <- min(qvalue, na.rm = TRUE)
      Excalibur <- value_output


      pvalue <- rbind(data.frame(pvalue = Excalibur), pvalue)


      computational_time <- data.frame()
      computational_time <- data.frame(
        computational_time = c(time_p_linear_davies_burden,
                               time_p_linear_weighted_liu_burden,
                               time_p_linear_weighted_liumod_skat,
                               time_p_linear_liu_skat,
                               time_p_linear_skato,
                               time_p_linear_weighted_skato)
      )

      time_Excalibur <- sum(computational_time[,1])
      computational_time <- rbind(data.frame(computational_time = time_Excalibur), computational_time)


    }#end of else, continuous case


    pvalue <- pvalue[,1]
    computational_time <- computational_time[,1]
    errors <- errors[,1]
    return(list(pvalue,
                computational_time,
                errors))
  }
}#end of function
