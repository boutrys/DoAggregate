test_that("small test on public dataset on the GVCF", {
  # Test code here
  ###Input parameters
  path_results = "/Users/sboutry/Documents/excalibur/Rpackage/DoAggregate/tests/testthat/debug_27012025"


  load_moduleR = 1
  RLIBRARIES = .libPaths()
  list_genes = 0
  rare_maf_threshold = 0.01
  var_covariate = 0
  var_covariate_pca = 1
  min_var = 2
  max_var = 10000
  do_data_preparation = 1
  DO_over_representation = 1
  version= "optimal"
  filtering_criteria = 0
  filtering_level = 0
  filtering_direction = 0
  coverage_matrix_path = 0
  mapqv_threshold = 20
  min_coverage = 20
  remove_individual = 1
  remove_position = 1
  remove_private = 0
  pvalue_threshold = 0.05
  Binary = TRUE
  path_Excalibur_function = "/Users/sboutry/Documents/excalibur/Rpackage/DoAggregate/R"


  ###Load excalibur functions
  source(paste(path_Excalibur_function, "main.R", sep = "/"))


  ################## other example launching on small public dataset using GVCF file

  PATH_GVCF = "/Users/sboutry/Documents/excalibur/Rpackage/DoAggregate/tests/testthat/small_test/small_public_test_data.vcf"
  case <- "/Users/sboutry/Documents/excalibur/Rpackage/DoAggregate/tests/testthat/small_test/cases.txt"
  control <- "/Users/sboutry/Documents/excalibur/Rpackage/DoAggregate/tests/testthat/small_test/controls.txt"
  mandatory_column_mapping =  list(sample = "Indiv",
                                   chr = "CHROM",
                                   pos = "POS",
                                   reference = "REF",
                                   alternative = "ALT",
                                   zygosity = "",
                                   gene_symbol = "Gene.ensGene",
                                   gnomad_wes_af = "gnomad40_AF",
                                   gt_GT = "gt_GT",
                                   gt_DP = "gt_DP")
  importation_method = 1


  source(paste(path_Excalibur_function, "main.R", sep = "/"))


  start <- Sys.time()
  ###Launching excalibur pipeline
  main(PATH_PATIENTS = case,
       PATH_CONTROLS = control,
       PATH_GVCF = PATH_GVCF,
       mandatory_column_mapping = mandatory_column_mapping,
       path_results = path_results,
       importation_method = importation_method,
       RLIBRARIES = .libPaths(),
       list_genes = 0,
       path_Excalibur_function = path_Excalibur_function)
  stop <- Sys.time()
  print(stop - start)




})
