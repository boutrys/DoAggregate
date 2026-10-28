test_that("DoAggregate runs on the small public GVCF test dataset", {
  skip_on_cran()
  path_results <- file.path(tempdir(), "DoAggregate_small_test")

  mandatory_column_mapping <- list(sample = "Indiv",
                                   chr = "CHROM",
                                   pos = "POS",
                                   reference = "REF",
                                   alternative = "ALT",
                                   zygosity = "",
                                   gene_symbol = "Gene.ensGene",
                                   gnomad_wes_af = "gnomad40_AF",
                                   gt_GT = "gt_GT",
                                   gt_DP = "gt_DP")

  main(PATH_PATIENTS = test_path("small_test", "cases.txt"),
       PATH_CONTROLS = test_path("small_test", "controls.txt"),
       PATH_GVCF = test_path("small_test", "small_public_test_data.vcf"),
       mandatory_column_mapping = mandatory_column_mapping,
       path_results = path_results,
       importation_method = 1,
       RLIBRARIES = .libPaths(),
       list_genes = 0)

  expect_true(file.exists(file.path(path_results, "general", "global_info_analysis.xlsx")))
  expect_true(file.exists(file.path(path_results, "gene", "general_info.xlsx")))
})
