# Parallel implementation (Slurm cluster)

This version splits the DoAggregate pipeline into Slurm jobs: data preparation, one job per genetic region
(gene, pathway or GO term), then result combination. It is the recommended way to analyse exomes or genomes.

## Contents

- `functions/submit_Excalibur.sh`: entry point; parses the options and submits the pipeline (`parallel_do_Excalibur.sh`).
- `functions/parallel_*.R|sh`: preparation, per-region analysis, result combination and final steps.
- `functions/*.R`: shared pipeline functions (data preparation, over-representation analysis, Excalibur tests, ...).
- `functions/install_packages.r`, `functions/install_r_lib.sh`: install the required R packages on the cluster.
- `functions/VCFtoTSV*.sh`, `functions/combine_gvcf_file*.sh`: helpers to combine GVCFs and convert VCF to TSV.
- `Data/*.tsv.gz`: pathway and gene-set tables used for the over-representation analysis and pathway annotation
  (built from MSigDB collections). Decompress them (`gunzip Data/*.gz`) before use.

## Setup

1. Copy `functions/` and `Data/` to the cluster, e.g. to `${EXCALIBUR}/data/functions/` and `${EXCALIBUR}/data/Data/`.
2. Install the R dependencies with `functions/install_r_lib.sh` (edit `RLIBRARY`, `SCRIPT` and the `module`
   lines for your cluster).
3. Edit the cluster-specific defaults at the top of `submit_Excalibur.sh` (working directory `--chdir`, log
   folders, Slurm partition/QOS, R module and R library path), or pass them as options (`-Q`, `-R`, `-L`, `-M`,
   `-v`, `-w`).
4. Optional: to remove blacklisted variants (e.g. recurrent artefacts of your sequencing platform), set the
   environment variable `EXCALIBUR_BLACKLIST` to a TSV file with `CHROM`, `POS`, `REF` and `ALT` columns.

## Running

```bash
# TSV input (cases and controls)
sbatch functions/submit_Excalibur.sh -E "you@example.org" -A cases.tsv -B controls.tsv -i 0.01 results_folder/

# Annotated VCF input with lists of case and control sample names (public example of this repository)
sbatch functions/submit_Excalibur.sh -E "you@example.org" -I 1 \
  -A data/small_example/cases.txt -B data/small_example/controls.txt \
  -G data/small_example/small_public_test_data.vcf -6 "Gene.ensGene" -i 0.01 results_folder/
```

Run `bash functions/submit_Excalibur.sh -h` for the full list of options (input columns, MAF threshold, minimum
and maximum number of variants per region, coverage-based cleaning, Excalibur version, maximum number of jobs, ...).
`functions/HOW_TO.txt` contains the command lines used for the public example.
