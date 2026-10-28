# Public test data

Whole-exome data of 51 individuals from the 1000 Genomes Project, used in the Results section of the article:
**26 CEU individuals as cases** (`NA*` samples) and **25 GBR individuals as controls** (`HG*` samples). Both groups are
of European ancestry, so this is close to a null comparison that illustrates the pipeline and its outputs.

| Path | Content |
|---|---|
| `big_example/cases.txt`, `controls.txt` | Sample lists |
| `big_example/inputs/annovar.vcf.hg38_multianno.reduced.vcf.gz` | Annotated multi-sample VCF of the 51 exomes (46 MB). Every variant and sample is kept, along with the fields used by DoAggregate (GT, AD, DP, GQ, gene annotation, gnomAD frequencies, CADD, REVEL, ClinVar) |
| `big_example/inputs/patient_data.tsv.gz`, `control_data.tsv.gz` | Rare-variant TSV inputs for cases and controls (decompress before use) |
| `big_example/results/` | DoAggregate results on the full exome dataset |
| `small_example/` | Nine-gene subset as an annotated VCF, with sample lists, cleaning outputs and DoAggregate results |
| `scripts/public_data_preTreatment.R` | Builds the TSV inputs from the annotated multi-sample VCF |
| `scripts/generate_small_test_public_dataset.R` | Builds the nine-gene VCF subset |
| `scripts/reduce_public_vcf.py` | Produces the reduced VCF from the full ANNOVAR output (removes unused annotation columns) |

The full ANNOVAR output (all 227 annotation fields, 131 MB compressed) is larger than GitHub's 100 MB file limit.
It was reduced with `scripts/reduce_public_vcf.py`, which keeps every variant and genotype.
