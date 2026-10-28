## Run DoAggregate (serial R package) on the public test datasets of this repository.
## Run from the repository root after installing the package (see Rpackage/DoAggregate/README.md).

library(DoAggregate)

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

### 1) Small example: annotated VCF (nine genes) + lists of case and control samples
start <- Sys.time()
main(PATH_PATIENTS = "data/small_example/cases.txt",
     PATH_CONTROLS = "data/small_example/controls.txt",
     PATH_GVCF = "data/small_example/small_public_test_data.vcf",
     mandatory_column_mapping = mandatory_column_mapping,
     importation_method = 1,
     path_results = "DoAggregate_results_small_example")
print(Sys.time() - start)

### 2) Full exome example: TSV files for cases (26 CEU) and controls (25 GBR)
## Decompress data/big_example/inputs/*.tsv.gz first. This run tests more than 10,000 genes and takes
## much longer in serial mode; the cluster implementation is recommended for this size.
start <- Sys.time()
main(PATH_PATIENTS = "data/big_example/inputs/patient_data.tsv",
     PATH_CONTROLS = "data/big_example/inputs/control_data.tsv",
     importation_method = 0,
     path_results = "DoAggregate_results_big_example")
print(Sys.time() - start)
