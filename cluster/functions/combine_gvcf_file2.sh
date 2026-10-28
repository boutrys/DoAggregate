#!/bin/bash
#SBATCH --chdir /scratch/sboutry
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --qos=serial
#SBATCH --mem 400G
#SBATCH --time 71:59:00
#SBATCH --mail-user=YOUR_EMAIL
#SBATCH -J combine_gvcf_file
#SBATCH --mail-type=FAIL,END
#SBATCH -o /scratch/sboutry/combineGVCF_%J.out # Standard output
#SBATCH -e /scratch/sboutry/combineGVCF_%J.err # Standard error

# Exit the script if any statement returns a non-true return value.
set -e
set -x


echo "START AT $(date)"

# Set ulimit to increase resource limits
ulimit -c unlimited  # Allow unlimited core dumps
ulimit -s unlimited  # Allow unlimited stack size


#input parameter of the script
#path_input_file=/work/gr-fe/archive/sample_repository/all_exome_gvcfs_hg38/SHCS392exomes #controls GVCF
path_input_file=/work/gr-fe/archive/sample_repository/all_exome_gvcfs_hg38/FVH/exomes #patient GVCF
path_name_individu=/work/gr-fe/sboutry/excalibur/input/11_12_23_gvcf_FVH/patient_name.txt
path_output=/work/gr-fe/sboutry/excalibur/input/11_12_23_gvcf_FVH/patient_data.vcf.gz
tmp_file=/scratch/sboutry

#Path to database and programs
REF=/work/gr-fe/saadat/Reference_Genome/GRCH38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz
BCFTOOLS=/work/gr-fe/sboutry/tools/bcftools/install/bin/bcftools
export BCFTOOLS_PLUGINS=/work/gr-fe/sboutry/tools/bcftools/bcftools/plugins
TABIX=/work/gr-fe/sboutry/tools/tabix/tabix-0.2.6/tabix


#Create a list of all GVCF files to combine
cd ${path_input_file}

#attempt to merge a list of vcf file, do not work when too much VCF files (74)
#$BCFTOOLS merge --file-list ${path_name_individu} -Oz -o ${path_output}/patient_data.vcf.gz

# Initialize the count variable
count=0
#for file in $(cat "$path_name_individu")
while IFS= read -r file
do
    if [ "${count}" == 0 ]
    then
        first_file=$file
    elif [ "${count}" == 1 ]
    then
        echo "Processing file: $first_file and $file"
        $BCFTOOLS merge "$first_file" "$file" -o ${tmp_file}/temporary.vcf.gz
        $TABIX -p vcf ${tmp_file}/temporary.vcf.gz
    else
        echo "Processing file: $file"
        $BCFTOOLS merge ${tmp_file}/temporary.vcf.gz "$file" -o ${tmp_file}/new_temporary.vcf.gz
        $TABIX -p vcf ${tmp_file}/new_temporary.vcf.gz
        rm ${tmp_file}/temporary.vcf.gz
        rm ${tmp_file}/temporary.vcf.gz.tbi
        cp ${tmp_file}/new_temporary.vcf.gz ${tmp_file}/temporary.vcf.gz
        cp ${tmp_file}/new_temporary.vcf.gz.tbi ${tmp_file}/temporary.vcf.gz.tbi
        rm ${tmp_file}/new_temporary.vcf.gz
        rm ${tmp_file}/new_temporary.vcf.gz.tbi
    fi
    
    ((++count))

done < "$path_name_individu"


#Move final results to output
cp ${tmp_file}/temporary.vcf.gz "${path_output}"
cp ${tmp_file}/temporary.vcf.gz.tbi "${path_output}.tbi"
rm ${tmp_file}/temporary.vcf.gz
rm ${tmp_file}/temporary.vcf.gz.tbi


echo "FINISH AT $(date)"




