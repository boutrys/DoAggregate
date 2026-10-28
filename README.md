<p align="center"><img src="excalibur.png" alt="Excalibur logo" width="160"></p>

# DoAggregate

**DoAggregate** is a pipeline for rare-variant association testing in case-control cohorts. It takes an annotated list
of variants, cleans the data, runs an over-representation analysis on the genes mutated in cases, and tests each
**gene**, **over-represented pathway** and **over-represented GO term** with aggregation tests. By default it uses
**Excalibur**, an ensemble of 36 aggregation tests, and 59 individual tests are available. Results are written as
annotated Excel files and plots that can be explored down to the variant level.

> Boutry S. *et al.* DoAggregate: a pipeline for aggregation testing for rare-variant on a gene and pathway level.
> bioRxiv (2026).

## Repository layout

| Folder | Content |
|---|---|
| [`Rpackage/`](Rpackage) | **Serial implementation** of DoAggregate as an R package, with installation helpers and an example script |
| [`cluster/`](cluster) | **Parallel implementation** for high-performance computing clusters using the Slurm scheduler |
| [`web_application/`](web_application) | **Web application** (PHP/JavaScript) that submits jobs to the cluster implementation and lets users follow them |
| [`data/`](data) | **Public test data** (1000 Genomes exomes) with the expected DoAggregate results |
| [`article/`](article) | Preprint manuscript, supplementary data, figures, tables and the script that produces them |

## Quick start (serial R package)

```r
# Dependency not available on CRAN/Bioconductor (provided in this repository)
install.packages("Rpackage/dependencies/RVtests_1.2.tar.gz", repos = NULL, type = "source")

# Bioconductor dependencies
install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "GOSemSim", "ReactomePA", "DOSE",
                       "ComplexHeatmap", "enrichplot"))

# DoAggregate
install.packages("remotes")
remotes::install_github("boutrys/DoAggregate", subdir = "Rpackage/DoAggregate")
```

```r
library(DoAggregate)
test_dir <- system.file("extdata", "small_test", package = "DoAggregate")
main(PATH_PATIENTS = file.path(test_dir, "cases.txt"),
     PATH_CONTROLS = file.path(test_dir, "controls.txt"),
     PATH_GVCF = file.path(test_dir, "small_public_test_data.vcf"),
     importation_method = 1,
     path_results = "DoAggregate_results")
```

See [`Rpackage/README.md`](Rpackage/README.md) for the input formats and parameters, and
[`cluster/README.md`](cluster/README.md) to run the parallel version on a Slurm cluster.

## Input

Either two TSV files (cases and controls; one line per variant per sample) or one annotated (G)VCF with two files
listing case and control sample names. Mandatory information per variant: sample, chromosome, position, reference
allele, alternative allele, gene symbol (HGNC), zygosity or genotype, and population allele frequency (e.g. gnomAD).

## Output

Four folders: `general` (cleaning and over-representation results), `gene`, `over_represented_pathway` and
`over_represented_GO`. Each contains annotated Excalibur results for all regions, the significant regions, and the
variant-level data behind each significant region. Every file is described in the Supplementary Data and in
`article/Table_S1_Results_description.xlsx`.

## Citation

If you use DoAggregate, please cite the DoAggregate preprint and the Excalibur article:

- Boutry S, Helaers R, Lenaerts T, Vikkula M. Excalibur: A new ensemble method based on an optimal combination of
  aggregation tests for rare-variant association testing for sequencing data. *PLoS Comput Biol*. 2023;19(9).

## License

MIT, see [LICENSE](LICENSE).
