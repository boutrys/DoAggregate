#!/bin/bash
#SBATCH --chdir /scratch/sboutry/logs/combine_gvcf_file
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --qos=serial
#SBATCH --mem 400G
#SBATCH --time 71:59:00
#SBATCH --mail-user=YOUR_EMAIL
#SBATCH -J combine_gvcf_fileTest
#SBATCH --mail-type=FAIL,END
#SBATCH -o /scratch/sboutry/logs/combine_gvcf_file/%J.out # Standard output
#SBATCH -e /scratch/sboutry/logs/combine_gvcf_file/%J.err # Standard error

# Exit the script if any statement returns a non-true return value.
set -e
set -x


echo "START AT $(date)"

# Set ulimit to increase resource limits
ulimit -c unlimited  # Allow unlimited core dumps
ulimit -s unlimited  # Allow unlimited stack size


##input parameter of the script
#path_input_file=/work/gr-fe/archive/sample_repository/all_exome_gvcfs_hg38/SHCS392exomes ##controls GVCF
#path_input_file=/work/gr-fe/archive/sample_repository/all_exome_gvcfs_hg38/FVH/exomes 
path_input_file=/work/gr-fe/sboutry/NEC/DeepVariant_WES
##patient GVCF
#path_name_individu=/work/gr-fe/sboutry/excalibur/input/11_12_23_gvcf_FVH/patient_name.gvcfs.list
path_name_individu=/work/gr-fe/sboutry/NEC/DeepVariant_WES/DeepVariantWes_patient_name.gvcfs.list
#path_output=/work/gr-fe/sboutry/excalibur/input/11_12_23_gvcf_FVH/patient_data.vcf.gz
path_output=/work/gr-fe/sboutry/NEC/DeepVariant_WES/patient_data.vcf.gz

#Path to database and programs
REF=/work/gr-fe/saadat/Reference_Genome/GRCH38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz
BCFTOOLS=/work/gr-fe/sboutry/tools/bcftools/install/bin/bcftools
export BCFTOOLS_PLUGINS=/work/gr-fe/sboutry/tools/bcftools/bcftools/plugins
TABIX=/work/gr-fe/sboutry/tools/tabix/tabix-0.2.6/tabix
GATK=/work/gr-fe/sboutry/tools/gatk/gatk-4.2.2.0/gatk


###Combine GVCF with GATK
cd ${path_input_file}

${GATK} --java-options "-Xmx180G -XX:ParallelGCThreads=36" CombineGVCFs -R ${REF} --variant ${path_name_individu} -O ${path_output}

echo "FINISH AT $(date)"




