#!/bin/bash
#SBATCH --chdir /scratch/sboutry
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --qos=serial
#SBATCH --mem 50G
#SBATCH --time 71:59:00
#SBATCH --mail-user=YOUR_EMAIL
#SBATCH -J TEST_VCFtoTSV2_patient
#SBATCH --mail-type=FAIL,END
#SBATCH -o /scratch/sboutry/logs/VCFtoTSV/TEST_%J.out # Standard output
#SBATCH -e /scratch/sboutry/logs/VCFtoTSV/TEST_%J.err # Standard error


# Exit the script if any statement returns a non-true return value.
set -e
set -x


#Input parameter
#vcf_file=/scratch/sboutry/testSnpEff.vcf
vcf_file=/work/gr-fe/sboutry/MELIS/NEC_finalMSC_alphaMissense_vep.vcf
# Output TSV file
#final_output_file="/scratch/sboutry/VCFkitTest.tsv"
final_output_file="/scratch/sboutry/NEC_finalMSC_alphaMissense_vep.tsv"

# Activate vep from conda
eval "$(conda shell.bash hook)"
conda activate vcf-kit

vk vcf2tsv long --ANN ${vcf_file} >> ${final_output_file}

conda deactivate