#!/bin/bash
#SBATCH --chdir /scratch/sboutry
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --qos=serial
#SBATCH --mem 50G
#SBATCH --time 71:59:00
#SBATCH --mail-user=YOUR_EMAIL
#SBATCH -J VCFtoTSV_MElis
#SBATCH --mail-type=FAIL,END
#SBATCH -o /scratch/sboutry/logs/VCFtoTSV/MELIS_%J.out # Standard output
#SBATCH -e /scratch/sboutry/logs/VCFtoTSV/MELIS_%J.err # Standard error


# Exit the script if any statement returns a non-true return value.
set -e
set -x

### The following mandatory variables are set in parallel_do_Excalibur.sh and are passed as environment variables by Slurm
#FOR R SCRIPT
##PATH_INPUT
##PATH_OUTPUT


echo "START AT $(date)"
# Input VCF file
#vcf_file=/work/gr-fe/sboutry/excalibur/input/11_12_23_gvcf_FVH/patient_data.vcf.gz
vcf_file=/work/gr-fe/sboutry/MELIS/NEC_finalMSC_alphaMissense_vep.vcf
COMPRESSED_FILE=0
# Output TSV file
final_output_file="/work/gr-fe/sboutry/MELIS/NEC_finalMSC_alphaMissense_vep.tsv"
nameOfanalysis="MELIS"
VEP=1 #set to 0 if vcf not annotated with VEP
IS_ANNOTATED=1 #set to 0 is vcf is not annotated


#Creating temporary vcf file
output_file="/scratch/sboutry/${nameOfanalysis}tmp_final_data.tsv"
touch /scratch/sboutry/${nameOfanalysis}temporary.vcf
tmp_vcf=/scratch/sboutry/${nameOfanalysis}temporary.vcf
chmod 775 ${tmp_vcf}
chmod 775 ${vcf_file}
if [ "${COMPRESSED_FILE}" == 1 ]
then
	gzip -dc "${vcf_file}" > "${tmp_vcf}"
else
	cp "${vcf_file}" "${tmp_vcf}"
fi


# Find the line with #CHROM and extract the column names
column_names=$(grep '#CHROM' "${tmp_vcf}")

# Extract only the column names after the FORMAT column
sample_list=$(echo "$column_names" | awk -F'\t' '{if(match($0,/\tFORMAT\t/)) {start=RSTART+RLENGTH; print substr($0, start);}}')


# Print column names to the TSV file
tsv_file=/scratch/sboutry/${nameOfanalysis}output.tsv
echo -e "$column_names" > "$tsv_file"

echo -e "$column_names"

# Flag to track whether header has been printed
header_printed=0

# Extract the relevant Genotype data from the VCF file, process sample_list columns, and append to the TSV file
awk -F'\t' -v sample_list="$sample_list" '
BEGIN {
	split(sample_list, samples, "\t")
}
{
if (!/^##/) {
	if (!header_printed) {
	header_printed = 1
	next
	}
	printf "%s", $1
	for (i = 2; i <= 9; i++) {
        printf "\t%s", $i
	}
	for (i = 10; i <= NF; i++) {
	printf "\t%s", $i
	}
	printf "\n"
}
}' "${tmp_vcf}" >> "$tsv_file"


# Output TSV file long format
long_file="/scratch/sboutry/${nameOfanalysis}long.tsv"


# Extract header and sample columns
header=$(head -n 1 "$tsv_file" | cut -f 1-9)
samples=$(head -n 1 "$tsv_file" | cut -f 10-)


echo "starting long tsv"

# make one line per sample and mutation keeping only het/hom variant and hom reference
awk -v header="$header" -v samples="$samples" '
BEGIN {
# Split the header and samples into arrays
split(header, header_arr, "\t")
split(samples, samples_arr, "\t")

# Define valid genotypes
valid_genotypes["0/1"]
valid_genotypes["0|1"]
valid_genotypes["1/1"]
valid_genotypes["1|1"]
valid_genotypes["0/0"]
valid_genotypes["0|0"]
}
NR > 1 {
# Loop through the samples and print the transformed lines for valid genotypes
for (i = 10; i <= NF; i++) {
    split($i, parts, ":")
	if (parts[1] in valid_genotypes) {
	print $1, $2, $3, $4, $5, $6, $7, $8, $9, samples_arr[i-9], parts[1], parts[3]
	}
}
}
' "$tsv_file" > "$long_file"

# Add the header back to the long file
sed -i '1s/^/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\tgenotype\tDP\n/' "$long_file"


#Add zygosity and variant ID columns
awk -v OFS='\t' '
BEGIN {
print "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "sample", "genotype", "DP", "zygosity", "variant_id"
}
NR > 1 {
chrom = $1
pos = $2
ref = $4
alt = $5

# Determine zygosity
if ($11 == "0/0" || $11 == "0|0") {
	zygosity = 0
} else if ($11 == "0/1" || $11 == "0|1") {
	zygosity = 1
} else if ($11 == "1/1" || $11 == "1|1") {
	zygosity = 2
} else {
	zygosity = "NA"
}

# Create variant_id
variant_id = chrom ":" pos ":" ref ":" alt

# Print the transformed line
print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, zygosity, variant_id
}
' "$long_file" > "$output_file"


### Annotation retrieval and processing
if [ $IS_ANNOTATED == "1" ]
then
    base_info=""
    # Use grep to filter the lines, awk to extract the desired part, and a loop to concatenate them
    tmp=$(grep '##INFO=<ID=' "${tmp_vcf}")
    # Extract ID fields and concatenate them into base_info
    base_info=$(echo "$tmp" | grep -o 'ID=[^,]*' | awk -F= '{print $2}' | awk '/ALLELE_END/{exit} {printf "%s|", $0}')
    # Remove the trailing "|"
    base_info=${base_info%|}
    base_info=$(echo "$base_info" | sed 's/NCBoost|//')
    base_info=$(echo "$base_info" | sed 's/AF|//')
    echo "${base_info}"

    # Your vep_column_values
    if [ $VEP == "1" ]
    then
        vep_column_name=$(grep 'Format: ' "${tmp_vcf}")
        vep_column_values=$(echo "$vep_column_name" | awk -F'Format: ' '{print $2}' | sed 's/"//g')
        # Remove the ">" at the end using sed
        vep_column_values=$(echo "$vep_column_values" | sed 's/>$//')
    fi


    # Initialize the output file
    annotation_file=/scratch/sboutry/${nameOfanalysis}_tmp_full_annotation.tsv
    echo "creating and writting in temporary file ${annotation_file}"
    touch ${annotation_file}
    chmod 775 ${annotation_file}
    # Process each line in long_file using awk
    if [ $VEP == "1" ]
    then
        tail -n +2 $output_file | awk -F '\t' -v vep_cols="$vep_column_values" -v base_info="$base_info" -v output_file="$annotation_file" '
        BEGIN {
            # Split vep_column_values into an array
            num_vep_cols = split(vep_cols, vep_columns, "|")
            # Split base_info into an array
            num_base_info = split(base_info, base_info_cols, "|")
        }

        {
            # Extract INFO column and remove ;ALLELE_END
            info_column=$8
            gsub(/;ALLELE_END/, "", info_column)

            # Initialize annotation_name and annotation_value
            annotation_name=""
            annotation_value=""

            # Process INFO column
            num_tmp_values = split(info_column, tmp_values, ";")

            # Iterate through base_info columns
            for (i = 1; i <= num_base_info; i++) {
                name = base_info_cols[i]
                value = ""

                #If annotated with annovar there is two AF column and we want the value of AF from ANNOVAR, so do not take first value
                count = 0
                # Check if the name is present in num_tmp_values
                for (j = 1; j <= num_tmp_values; j++) {
                    split(tmp_values[j], pair, "=")
                    if (pair[1] == name) {
                        if (count == 0 &&  pair[1] == "AF") {
                            count = 1
                        } else {
                            value = pair[2]
                            break
                        }
                        
                    }
                }

                annotation_name = annotation_name name "\t"
                annotation_value = annotation_value value "\t"
            }

            # Check for CSQ
            for (i = 1; i <= num_tmp_values; i++) {

                split(tmp_values[i], pair, "=")
                name = pair[1]
                value = pair[2]

                # Skip if no equal sign is found
                if (index(base_info, name "=") == 0 && value == "") {
                    continue
                }

                if (name == "CSQ") {
                    # Handle CSQ differently
                    split(value, csq_values, "|")
                    for (j = 1; j <= num_vep_cols; j++) {
                        annotation_name = annotation_name vep_columns[j] "\t"
                        annotation_value = annotation_value csq_values[j] "\t"
                    }
                } else if (value == "") {
                    annotation_name = annotation_name name "\t"
                    annotation_value = annotation_value value "\t"
                } 
            }

            # Print the results to the output file
            if (NR == 1) {
                print annotation_name >> output_file
            }
            print annotation_value >> output_file
        }
        '
    else
        tail -n +2 $output_file | awk -F '\t' -v base_info="$base_info" -v output_file="$annotation_file" '
        BEGIN {
            # Split base_info into an array
            num_base_info = split(base_info, base_info_cols, "|")
        }

        {
            # Extract INFO column and remove ;ALLELE_END
            info_column=$8
            gsub(/;ALLELE_END/, "", info_column)

            # Initialize annotation_name and annotation_value
            annotation_name=""
            annotation_value=""

            # Process INFO column
            num_tmp_values = split(info_column, tmp_values, ";")

            # Iterate through base_info columns
            for (i = 1; i <= num_base_info; i++) {
                name = base_info_cols[i]
                value = ""

                #If annotated with annovar there is two AF column and we want the value of AF from ANNOVAR, so do not take first value
                count = 0
                # Check if the name is present in num_tmp_values
                for (j = 1; j <= num_tmp_values; j++) {
                    split(tmp_values[j], pair, "=")
                    if (pair[1] == name) {
                        if (count == 0 &&  pair[1] == "AF") {
                            count = 1
                        } else {
                            value = pair[2]
                            break
                        }
                        
                    }
                }

                annotation_name = annotation_name name "\t"
                annotation_value = annotation_value value "\t"
            }
            
            for (i = 1; i <= num_tmp_values; i++) {

                split(tmp_values[i], pair, "=")
                name = pair[1]
                value = pair[2]

                # Skip if no equal sign is found
                if (index(base_info, name "=") == 0 && value == "") {
                    continue
                }

                if (value == "") {
                    annotation_name = annotation_name name "\t"
                    annotation_value = annotation_value value "\t"
                } 
            }

            # Print the results to the output file
            if (NR == 1) {
                print annotation_name >> output_file
            }
            print annotation_value >> output_file
        }
        '
    fi


    # Check if both files have the same number of lines
    num_lines_output=$(wc -l < "$output_file")
    num_lines_annotation=$(wc -l < "$annotation_file")

    if [ "$num_lines_output" -ne "$num_lines_annotation" ]; then
        echo "Error: The files do not have the same number of lines."
        exit 1
    fi

    # Function to find column number by name
    find_column_number() {
        column_name="$1"
        column_number=$(awk -F'\t' -v name="$column_name" 'NR==1 {for (i=1; i<=NF; i++) if ($i == name) print i; exit}' "$2")
        echo "$column_number"
    }

    # Find column numbers for 'INFO' in $output_file
    info_column_number_output=$(find_column_number 'INFO' "$output_file")

    # Remove 'INFO' column from $output_file
    cut -f1-$((info_column_number_output - 1)),$((info_column_number_output + 1))- "$output_file" > "${output_file}_temp"

    if [ $VEP == "1" ]
    then
        # Find column numbers for 'CSQ' in $annotation_file
        csq_column_number_annotation=$(find_column_number 'CSQ' "$annotation_file")
        # Remove 'CSQ' column from $annotation_file
        cut -f1-$((csq_column_number_annotation - 1)),$((csq_column_number_annotation + 1))- "$annotation_file" > "${annotation_file}_temp"
    else
        cp $annotation_file "${annotation_file}_temp"
    fi

    # Combine files with output columns first
    paste "${output_file}_temp" "${annotation_file}_temp" > "$final_output_file"

    # Remove temporary files
    rm "${output_file}_temp" 
    rm "${annotation_file}_temp"
    rm ${annotation_file}

    echo "Files combined successfully. Output saved to $final_output_file"

else
    ###if no annotations
    cp output_file final_output_file
fi


# removing temporary files
rm ${tmp_vcf}
rm ${tsv_file}
rm ${long_file}
rm ${output_file}
    

echo "END AT $(date)"