# DoAggregate (R package, serial implementation)

DoAggregate runs rare-variant aggregation tests on genes, over-represented pathways and over-represented GO terms,
using the Excalibur ensemble method by default. This package runs the whole pipeline in a single R session. For large
cohorts, the parallel implementation for Slurm clusters (`cluster/` folder of the repository) is much faster.

## Installation

```r
# RVtests is not on CRAN; the source tarball is provided in Rpackage/dependencies/
install.packages("Rpackage/dependencies/RVtests_1.2.tar.gz", repos = NULL, type = "source")

install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "GOSemSim", "ReactomePA", "DOSE",
                       "ComplexHeatmap", "enrichplot"))

install.packages("remotes")
remotes::install_github("boutrys/DoAggregate", subdir = "Rpackage/DoAggregate")
```

The other dependencies (SKAT, AssotesteR, iGasso, WGScan, REBET, CATT, KATSP, vcfR, ...) are installed from CRAN
automatically. A `glmnet` 4.1-8 source tarball is also provided in `Rpackage/dependencies/` for older R versions.

## Usage

```r
library(DoAggregate)

# TSV input: one file for cases and one for controls
main(PATH_PATIENTS = "patient_data.tsv",
     PATH_CONTROLS = "control_data.tsv",
     importation_method = 0,
     path_results = "DoAggregate_results")

# Annotated (G)VCF input: one VCF plus two files listing case and control sample names
test_dir <- system.file("extdata", "small_test", package = "DoAggregate")
main(PATH_PATIENTS = file.path(test_dir, "cases.txt"),
     PATH_CONTROLS = file.path(test_dir, "controls.txt"),
     PATH_GVCF = file.path(test_dir, "small_public_test_data.vcf"),
     importation_method = 1,
     mandatory_column_mapping = list(sample = "Indiv", chr = "CHROM", pos = "POS",
                                     reference = "REF", alternative = "ALT", zygosity = "",
                                     gene_symbol = "Gene.ensGene", gnomad_wes_af = "gnomad40_AF",
                                     gt_GT = "gt_GT", gt_DP = "gt_DP"),
     path_results = "DoAggregate_results")
```

### Mandatory TSV columns

`sample`, `chr`, `pos`, `reference`, `alternative`, `zygosity` (Heterozygous or Homozygous), `gene_symbol` (HGNC)
and `gnomad_wes_af` (population allele frequency). With VCF input, `mandatory_column_mapping` maps these to the VCF
fields, and zygosity is computed from the genotype (`gt_GT`) when no zygosity field is given.

### Main parameters

| Parameter | Default | Meaning |
|---|---|---|
| `rare_maf_threshold` | 0.01 | Maximum allele frequency for a variant to be considered rare |
| `min_var`, `max_var` | 2, 10000 | Minimum and maximum number of unique variants for a region to be tested |
| `list_genes` | 0 | Path to a gene list to restrict the analysis (0 = all genes) |
| `DO_over_representation` | 1 | Run the over-representation analysis and test pathways and GO terms |
| `version` | "optimal" | Excalibur (optimal combination); "all" or a list of test names for other tests |
| `pvalue_threshold` | 0.05 | Significance threshold on Benjamini-Hochberg adjusted p-values |
| `coverage_matrix_path`, `min_coverage`, `mapqv_threshold`, `remove_individual`, `remove_position` | | Optional coverage-based cleaning |

See `?main` and the vignette (`vignette("DoAggregate")`) for all parameters and a description of every output file.

## Tests

`tests/testthat/test_small_public_GVCF.R` runs the pipeline on the small public dataset in `inst/extdata/small_test`
and checks that the main result files are produced (`devtools::test()`).

## License

MIT
