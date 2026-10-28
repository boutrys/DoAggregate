#!/bin/bash

# Exit the script if any statement returns a non-true return value.
set -e
set -x
#set -u

### The following mandatory variables are set in submit_Excalibur.sh and are passed as environment variables by Slurm
#FOR SLURM JOBS
###PROJECT
###MAX_JOB
###QUEUE
###RUNNINGDIR
###EMAIL
###THREADS_FINAL
###MEM_FINAL
###MEM_job
###LOG
###ZIP
###ZIPPING

#FOR IMPORTATION AND COVERAG
###PATH_SAMPLES

#FOR R SCRIPT
##PATH_PATIENTS
##PATH_CONTROLS
##path_GVCF
##sample_GVCF
##chr_GVCF
##pos_GVCF
##ref_GVCF
##alt_GVCF
##zyosity_GVCF
##gene_symbol_GVCF
##frequency_GVCF
##GT_GVCF
##DP_GVCF
##know_target_region
##filtering_criteria
##filtering_level
##filtering_direction
##min_var
##max_var
##rare_maf_threshold
##var_covariate
##var_covariate_pca
##importation_method
##mapqv_threshold
##min_coverage
##remove_individual
##remove_position
##INPUT
##RESULTS
##do_data_preparation
##path_data_annotate
##version
##path_Excalibur_function
##path_Excalibur_Data
##DO_over_representation
##RVERSION
##RLIBRARIES
##load_moduleR
##remove_private


### Functions 
## Output starting time, host name and stage name (stages are big chunk of the pipeline)
# First argument is a string describing the stage
function startStage()
{
	host=$(hostname)
	stageStart=$(date +%s)
	echo "### Started on $host at $(date +"%Y-%m-%d %H:%M") - $1"
}

## Output ending time, elpased time, host name and stage name (stages are big chunk of the pipeline)
# First argument is a string describing the stage (use the same than startStage)
function finishStage()
{
    host=$(hostname)
	stageEnd=$(date +%s)
	stageDiff=$(((stageEnd-stageStart)/60))
	echo "### Ended on $host at $(date +"%Y-%m-%d %H:%M") - $1"
	echo "### Elapsed on $host : $stageDiff minutes - $1"
}

## Load modules given after first argument and output starting time, host name and step name (steps are generally one command)
# First argument is a string describing the step
# You can give any number of modules to load starting as second argument
function startStep()
{
	host=$(hostname)
	stepStart=$(date +%s)
	echo "## Started on $host at $(date +"%Y-%m-%d %H:%M") - $1"
    if [ "${load_moduleR}" == 1 ]
    then
        if [ "$#" -gt "1" ]
        then
            module load ${@:2}	
        fi
    else
        if [ "$#" -gt "1" ]
        then
            module use ${@:2}
            #module load gcc/11.3.0
            module load r

        fi
    fi
}

## Purge modules and output ending time, elpased time, host name and step name (steps are generally one command)
# First argument is a string describing the step (use the same than startStep)
function finishStep()
{
	module purge
	host=$(hostname)
	stepEnd=$(date +%s)
	stepDiff=$(((stepEnd-stepStart)/60))
	echo "## Ended on $host at $(date +"%Y-%m-%d %H:%M") - $1"
	echo "## Elapsed on $host : $stepDiff minutes - $1"
}


## Module version
if [ "${load_moduleR}" == 1 ]
then
    export ${RLIBRARIES}
fi


ulimit -c unlimited
umask 002


startStage "Data importation and cleaning and annotation"
#folder to put temporary files 
TMP_FOLDER="${RESULTS}tmp"
mkdir -p ${TMP_FOLDER}
chmod 777 ${TMP_FOLDER} || true


startStep "DATA Importation and preparation" $RVERSION
TMP_FOLDER_data_preparation="${TMP_FOLDER}/data_preparation/"
mkdir -p ${TMP_FOLDER_data_preparation}
chmod 777 ${TMP_FOLDER_data_preparation} || true
#FOLDER_RESULTS="${PROJECTDIR}/Results"
#mkdir -p ${FOLDER_RESULTS}
#chmod 777 ${FOLDER_RESULTS} || true
FOLDER_RESULTS="${RESULTS}"
FOLDER_GENERAL_RESULTS="${FOLDER_RESULTS}general"
mkdir -p ${FOLDER_GENERAL_RESULTS}
chmod 777 ${FOLDER_GENERAL_RESULTS} || true
FOLDER_GENERAL_RESULTS="${FOLDER_GENERAL_RESULTS}/"

R_SCRIPT=${path_Excalibur_function}parallel_script_main.R
#launch main Rscript
COVERAGE_PATH=0
Rscript ${R_SCRIPT} ${PATH_PATIENTS} ${PATH_CONTROLS} ${know_target_region} ${importation_method} ${mapqv_threshold} ${min_coverage} ${remove_individual} ${remove_position} ${FOLDER_GENERAL_RESULTS} ${COVERAGE_PATH} ${do_data_preparation} ${path_data_annotate} ${path_Excalibur_function} ${TMP_FOLDER_data_preparation} ${RLIBRARIES} ${path_GVCF} ${sample_GVCF} ${chr_GVCF} ${pos_GVCF} ${ref_GVCF} ${alt_GVCF} ${zyosity_GVCF} ${gene_symbol_GVCF} ${frequency_GVCF} ${GT_GVCF} ${DP_GVCF} ${rare_maf_threshold}
finishStep "DATA Importation and preparation"

finishStage "Data importation and cleaning and annotation"

startStage "Analysis preparation"
startStep "Analysis preparation" $RVERSION
TMP_FOLDER_analysis_preparation="${TMP_FOLDER}/analysis_preparation/"
mkdir -p ${TMP_FOLDER_analysis_preparation}
chmod 777 ${TMP_FOLDER_analysis_preparation} || true

R_SCRIPT=${path_Excalibur_function}parallel_script_analysis_prep.R
#launch analysis preparation Rscript
KindOfRegion=1 #gene
Rscript ${R_SCRIPT} ${TMP_FOLDER_data_preparation} ${KindOfRegion} ${filtering_criteria} ${filtering_level} ${filtering_direction} ${TMP_FOLDER_analysis_preparation} ${path_Excalibur_function} ${max_var} ${remove_private} ${RLIBRARIES}

SUB_LOG="${LOG}region_logs/"
mkdir -p ${SUB_LOG}
chmod 777 ${SUB_LOG} || true
finishStep "Analysis preparation"
finishStage "Analysis preparation"


TMP_FOLDER_gene_set_analysis="${TMP_FOLDER}/gene_set_analysis/"
mkdir -p ${TMP_FOLDER_gene_set_analysis}
chmod 777 ${TMP_FOLDER_gene_set_analysis} || true
#Variable for Slurm
PROJECT_level="${PROJECT}_pathway"

startStage "Over-representation analysis"
startStep "Gene set analysis" $RVERSION
R_SCRIPT=${path_Excalibur_function}parallel_script_gene_set.R
Rscript ${R_SCRIPT} ${TMP_FOLDER_data_preparation} ${path_data_annotate} ${path_Excalibur_function} ${TMP_FOLDER_gene_set_analysis} ${FOLDER_GENERAL_RESULTS} ${DO_over_representation} ${path_Excalibur_Data} ${RLIBRARIES}
finishStep "Gene set analysis"

finishStage "Over-representation analysis"


startStage "Launching job for Gene analysis"
startStep "Region analysis"
TMP_RESULT_GENE="${TMP_FOLDER_analysis_preparation}gene/results"
mkdir -p ${TMP_RESULT_GENE}
chmod 777 ${TMP_RESULT_GENE} || true
TMP_RESULT_GENE="${TMP_RESULT_GENE}/"
TMP_TIME_GENE="${TMP_FOLDER_analysis_preparation}gene/time"
mkdir -p ${TMP_TIME_GENE}
chmod 777 ${TMP_TIME_GENE} || true
TMP_TIME_GENE="${TMP_TIME_GENE}/"
LIST_JOB=""
LIST_JOB_2=""
count_LIST_JOB=0
PATH_DATA="${TMP_FOLDER_analysis_preparation}gene/data_to_analyze/"
cd ${PATH_DATA}
nbr_region=(*)
nbr_region=${#nbr_region[@]}
files=$(basename "${PATH_DATA}*")
finishStep "Region analysis"
nbr_region_per_job=$(( (nbr_region + MAX_JOB - 1) / MAX_JOB ))
current_region=0
current_batch=1
region_total_looped=0
LIST_id_region=""
for i in ${files}
do
    namefile="${i}"
    IFS="_"
    read -ra id_region <<< "${namefile}"
    IFS=""

    if [ "$LIST_id_region" == "" ]
    then
        LIST_id_region="${id_region}"
    else
        LIST_id_region="${LIST_id_region}/${id_region}"
    fi
    ((current_region+=1))
    ((region_total_looped+=1))
    
    if [ "$current_region" == "$nbr_region_per_job" ] || [ "$region_total_looped" == "$nbr_region" ]
    then
        #Variable for Slurm
        PROJECT_single_region="${PROJECT}_GENE_batch_nbr_${current_batch}"

        startStep "launching ${PROJECT_single_region}"
        
        JOB_GENE_Level=`sbatch --parsable --job-name=${PROJECT_single_region} --qos=${QUEUE} --chdir=${RUNNINGDIR} --mail-user=${EMAIL} --mail-type=FAIL --ntasks=1 --cpus-per-task=1 --mem=${MEM_job} --time=02:00:00 --output=${SUB_LOG}%j.%x --export=ALL,id=${LIST_id_region},PATH_DATA=${PATH_DATA},max_var=${max_var},rare_maf_threshold=${rare_maf_threshold},var_covariate=${var_covariate},var_covariate_pca=${var_covariate_pca},version=${version},path_Excalibur_function=${path_Excalibur_function},RVERSION=${RVERSION},RLIBRARIES=${RLIBRARIES},TMP_RESULT=${TMP_RESULT_GENE},TMP_TIME=${TMP_TIME_GENE},load_moduleR=${load_moduleR} ${path_Excalibur_function}parallel_single_region_analysis.sh`
        
        ((current_batch+=1))
        
        if [ "$LIST_JOB" == "" ]
        then
            LIST_JOB="${LIST_JOB}${JOB_GENE_Level}"
        else
            LIST_JOB="${LIST_JOB},${JOB_GENE_Level}"
        fi
        
        finishStep "launching ${PROJECT_single_region}"
        ((count_LIST_JOB+=1))
        if [ "$count_LIST_JOB" == 100 ]
        then
            LIST_JOB_2="${LIST_JOB}"
            DEPEND="afterok:${LIST_JOB_2}"
            project_wait="wait_job_${PROJECT}"
            JOB_wait=`sbatch --parsable --job-name=${project_wait} --qos=${QUEUE} --chdir=${RUNNINGDIR} --mail-user=${EMAIL} --mail-type=FAIL --ntasks=1 --cpus-per-task=1 --mem=${MEM_job} --time=00:10:00 --output=${LOG}%j.%x --dependency=${DEPEND} --kill-on-invalid-dep=yes --export=ALL ${path_Excalibur_function}wait.sh`
            DEPENDNOTOK="afternotok:${LIST_JOB_2}"
            LIST_JOB="${JOB_wait}"
            count_LIST_JOB=0
        fi

        current_region=0
        LIST_id_region=""
    fi
    
done
finishStage "Launching job for Gene analysis"


if [ "$LIST_JOB" != "" ] || [ "$LIST_JOB_2" != "" ]
then
    ###
    startStage "Results combinations"
    ulimit -s 1000000

    PATH_INPUT=${TMP_FOLDER_analysis_preparation}
    level="gene"
    PATH_OUTPUT="${FOLDER_RESULTS}"
    PATH_GENERAL=${FOLDER_GENERAL_RESULTS}
    DEPEND="afterok:${LIST_JOB}"
    DEPENDNOTOK="afternotok:${LIST_JOB}"

    #Gene level combination
    PROJECT_combi_res="${PROJECT}_combi_results_level_${level}"
    JOB_combi_res=`sbatch --parsable --job-name=${PROJECT_combi_res} --qos=${QUEUE} --chdir=${RUNNINGDIR} --mail-user=${EMAIL} --mail-type=FAIL --ntasks=1 --cpus-per-task=1 --mem=${MEM_job} --time=00:30:00 --output=${LOG}%j.%x --dependency=${DEPEND} --kill-on-invalid-dep=yes --export=ALL,PATH_INPUT=${PATH_INPUT},PATH_OUTPUT=${PATH_OUTPUT},version=${version},level=${level},PATH_GENERAL=${PATH_GENERAL},path_Excalibur_function=${path_Excalibur_function},path_Excalibur_Data=${path_Excalibur_Data},RVERSION=${RVERSION},RLIBRARIES=${RLIBRARIES},load_moduleR=${load_moduleR} ${path_Excalibur_function}parallel_result_combination.sh`
    finishStage "Results combinations"
else
    JOB_combi_res=""
fi



startStep "Pathways aggregation analysis"
startStage "Launching pathway"
#TODO not sure which is best to avoid job dependency problems
LIST_JOB_pathway=""
#LIST_JOB_pathway=${JOB_combi_res}
LIST_JOB_3=""
count_LIST_JOB=0
cd ${TMP_FOLDER_gene_set_analysis}
nbr_files=(*)
nbr_files=${#nbr_files[@]}
if [ "$nbr_files" > 1 ]
then
    startStep "Analysis preparation" $RVERSION
    ###Analysis preparation
    if [ "${DO_over_representation}" == 1 ]
    then
        KindOfRegion=4 #over represented pathway
    else
        KindOfRegion=2 #pathway
    fi
    
    R_SCRIPT=${path_Excalibur_function}parallel_script_analysis_prep.R
    Rscript ${R_SCRIPT} ${TMP_FOLDER_data_preparation} ${KindOfRegion} ${filtering_criteria} ${filtering_level} ${filtering_direction} ${TMP_FOLDER_analysis_preparation} ${path_Excalibur_function} ${max_var} ${remove_private} ${RLIBRARIES} ${TMP_FOLDER_gene_set_analysis}
    echo "1"
    finishStep "Analysis preparation"

    ###Region analysis
    #Pathway over-represented
    FOLDER_pathway_analysis_preparation="${TMP_FOLDER_analysis_preparation}pathway"
    if [ -d "$FOLDER_pathway_analysis_preparation" ]
    then
        startStep "Lauching job for Pathway over-represented"
        TMP_RESULT_pathway="${TMP_FOLDER_analysis_preparation}pathway/results"
        mkdir -p ${TMP_RESULT_pathway}
        chmod 777 ${TMP_RESULT_pathway} || true
        TMP_RESULT_pathway="${TMP_RESULT_pathway}/"
        TMP_TIME_pathway="${TMP_FOLDER_analysis_preparation}pathway/time"
        mkdir -p ${TMP_TIME_pathway}
        chmod 777 ${TMP_TIME_pathway} || true
        TMP_TIME_pathway="${TMP_TIME_pathway}/"
        PATH_DATA="${TMP_FOLDER_analysis_preparation}pathway/data_to_analyze/"
        cd ${PATH_DATA}
        nbr_region=(*)
        nbr_region=${#nbr_region[@]}
        files=$(basename "${PATH_DATA}*")
        nbr_region_per_job=$(( (nbr_region + MAX_JOB - 1) / MAX_JOB ))
        current_region=0
        current_batch=1
        region_total_looped=0
        LIST_id_region=""
        for i in ${files}
        do
            namefile="${i}"
            IFS="_"
            read -ra id_region <<< "${namefile}"
            IFS=""

            if [ "$LIST_id_region" == "" ]
            then
                LIST_id_region="${id_region}"
            else
                LIST_id_region="${LIST_id_region}/${id_region}"
            fi
            
            ((current_region+=1))
            ((region_total_looped+=1))
            
            if [ "$current_region" == "$nbr_region_per_job" ] || [ "$region_total_looped" == "$nbr_region" ]
                then
                #Variable for Slurm
                PROJECT_single_region="${PROJECT_level}_batch_nbr_${current_batch}"

                startStep "lauching ${PROJECT_single_region}"
                JOB_pathway_Level=`sbatch --parsable --job-name=${PROJECT_single_region} --qos=${QUEUE} --chdir=${RUNNINGDIR} --mail-user=${EMAIL} --mail-type=FAIL --ntasks=1 --cpus-per-task=1 --mem=${MEM_job} --time=05:00:00 --output=${SUB_LOG}%j.%x --export=ALL,id=${LIST_id_region},PATH_DATA=${PATH_DATA},max_var=${max_var},rare_maf_threshold=${rare_maf_threshold},var_covariate=${var_covariate},var_covariate_pca=${var_covariate_pca},version=${version},path_Excalibur_function=${path_Excalibur_function},RVERSION=${RVERSION},RLIBRARIES=${RLIBRARIES},TMP_RESULT=${TMP_RESULT_pathway},TMP_TIME=${TMP_TIME_pathway},load_moduleR=${load_moduleR} ${path_Excalibur_function}parallel_single_region_analysis.sh`
                
                ((current_batch+=1))
                
                if [ "$LIST_JOB_pathway" == "" ]
                then
                    LIST_JOB_pathway="${LIST_JOB_pathway}${JOB_pathway_Level}"
                else
                    LIST_JOB_pathway="${LIST_JOB_pathway},${JOB_pathway_Level}"
                fi

                ((count_LIST_JOB+=1))
                if [ "$count_LIST_JOB" == 100 ]
                then
                    LIST_JOB_3="${LIST_JOB_pathway}"
                    DEPEND="afterok:${LIST_JOB_3}"
                    project_wait="wait_job_${PROJECT}"
                    JOB_wait=`sbatch --parsable --job-name=${project_wait} --qos=${QUEUE} --chdir=${RUNNINGDIR} --mail-user=${EMAIL} --mail-type=FAIL --ntasks=1 --cpus-per-task=1 --mem=${MEM_job} --time=00:10:00 --output=${LOG}%j.%x --dependency=${DEPEND} --kill-on-invalid-dep=yes --export=ALL ${path_Excalibur_function}wait.sh`
                    DEPENDNOTOK="afternotok:${LIST_JOB_3}"
                    LIST_JOB_pathway="${JOB_wait}"
                    count_LIST_JOB=0
                fi

                current_region=0
                LIST_id_region=""

                finishStep "done lauching ${PROJECT_single_region}"
            fi
        done
        finishStep "Lauching job for Pathway over-represented"
    fi #end of pathway over-represented
fi
finishStage "Launching pathway"
finishStep "Pathways aggregation analysis"






#TODO Check if bug when no pathway to analyze
if [ "$LIST_JOB_pathway" != "" ] || [ "$LIST_JOB_3" != "" ]
then
    ###
    startStage "Results combinations"
    ulimit -s 1000000
    
    
    #over-represented terms combination
    level="pathway"
    PROJECT_combi_res="${PROJECT}_combi_results_level_${level}"
    if [ "$LIST_JOB_pathway" == "" ]
    then
        DEPEND="afterok:${JOB_combi_res}"
    else
        DEPEND="afterok:${LIST_JOB_pathway},${JOB_combi_res}"
    fi

    JOB_combi_res2=`sbatch --parsable --job-name=${PROJECT_combi_res} --qos=${QUEUE} --chdir=${RUNNINGDIR} --mail-user=${EMAIL} --mail-type=FAIL --ntasks=1 --cpus-per-task=1 --mem=${MEM_job} --time=00:30:00 --output=${LOG}%j.%x --dependency=${DEPEND} --kill-on-invalid-dep=yes --export=ALL,PATH_INPUT=${PATH_INPUT},PATH_OUTPUT=${PATH_OUTPUT},version=${version},level=${level},PATH_GENERAL=${PATH_GENERAL},path_Excalibur_function=${path_Excalibur_function},path_Excalibur_Data=${path_Excalibur_Data},RVERSION=${RVERSION},RLIBRARIES=${RLIBRARIES},load_moduleR=${load_moduleR} ${path_Excalibur_function}parallel_result_combination.sh`
    DEPENDNOTOK="afternotok${LIST_JOB_pathway}:${JOB_combi_res}"
    #sbatch --parsable --job-name=mailres.excalibur.${PROJECT} --qos=${QUEUE} --chdir=/tmp --mail-user=${EMAIL} --mail-type=FAIL --ntasks=1 --cpus-per-task=1 --mem=1g --time=20:00:00 --output=${LOG}%j.%x --dependency=${DEPENDNOTOK} --kill-on-invalid-dep=yes --export=ALL,STATUS=ERROR,TOOL=excalibur,PROJECT=${PROJECT},EMAIL=${EMAIL} ${MAILING}
    finishStage "Results combinations"
else
    JOB_combi_res2=""
fi


if [ "$JOB_combi_res2" != "" ]
then
    startStage "End of analysis"
    startStep "Global results"
    PATH_GLOBAL=${TMP_FOLDER_analysis_preparation}
    PATH_TO_RESULT=${FOLDER_RESULTS}
    PATH_OUTPUT=${FOLDER_GENERAL_RESULTS}
    #DEPEND="afterok:${JOB_combi_res},${JOB_combi_res2}"
    DEPEND="afterok:${JOB_combi_res2}"
    PROJECT_finish="${PROJECT}_finish"

    JOB_finish=`sbatch --parsable --job-name=${PROJECT_finish} --qos=${QUEUE} --chdir=${RUNNINGDIR} --mail-user=${EMAIL} --mail-type=FAIL --ntasks=1 --cpus-per-task=1 --mem=${MEM_job} --time=00:30:00 --output=${LOG}%j.%x --dependency=${DEPEND} --kill-on-invalid-dep=yes --export=ALL,PATH_GLOBAL=${PATH_GLOBAL},PATH_TO_RESULT=${PATH_TO_RESULT},PATH_OUTPUT=${PATH_OUTPUT},path_Excalibur_function=${path_Excalibur_function},RVERSION=${RVERSION},RLIBRARIES=${RLIBRARIES},SUB_LOG=${SUB_LOG},load_moduleR=${load_moduleR} ${path_Excalibur_function}parallel_finish.sh`

    DEPEND="afterok:${JOB_finish}"
    DEPENDNOTOK="afternotok:${JOB_finish}"

    if [ "$ZIP" != "" ]
    then
        MODE=ZIP
        JOB_zip=`sbatch --parsable --job-name=zip.excalibur.${PROJECT} --qos=${QUEUE} --chdir=${RUNNINGDIR} --mail-user=${EMAIL} --mail-type=FAIL,END --ntasks=1 --cpus-per-task=4 --mem=${MEM_job} --time=00:30:00 --output=${LOG}%j.%x --dependency=${DEPEND} --kill-on-invalid-dep=yes --export=ALL,STATUS=DONE,PROJECT=${PROJECT},EMAIL=${EMAIL},FOLDER_RESULTS=${FOLDER_RESULTS},ZIP=${ZIP},TMP_FOLDER=${TMP_FOLDER},load_moduleR=${load_moduleR},SUB_LOG=${SUB_LOG} ${path_Excalibur_function}do_zip.sh`
        ZIP_FILENAME=${ZIP}${PROJECT}.zip
        ZIP_FILENAME=$(basename "$ZIP_FILENAME")
        DEPEND="afterok:${JOB_zip}"
        DEPENDNOTOK="afternotok:${JOB_zip}"
    else
        MODE=LINK
        ZIP_FILENAME=NONE
    fi
    #sbatch --parsable --job-name=mailres.excalibur.${PROJECT} --qos=${QUEUE} --chdir=/tmp --mail-user=${EMAIL} --mail-type=FAIL --ntasks=1 --cpus-per-task=1 --mem=1g --time=20:00:00 --output=${LOG}%j.%x --dependency=${DEPEND} --kill-on-invalid-dep=yes --export=ALL,STATUS=DONE,TOOL=excalibur,PROJECT=${PROJECT},EMAIL=${EMAIL},MODE=${MODE},SAMPLE=${ZIP_FILENAME} ${MAILING}
    #sbatch --parsable --job-name=mailres.excalibur.${PROJECT} --qos=${QUEUE} --chdir=/tmp --mail-user=${EMAIL} --mail-type=FAIL --ntasks=1 --cpus-per-task=1 --mem=1g --time=20:00:00 --output=${LOG}%j.%x --dependency=${DEPENDNOTOK} --kill-on-invalid-dep=yes --export=ALL,STATUS=ERROR,TOOL=excalibur,PROJECT=${PROJECT},EMAIL=${EMAIL} ${MAILING}

    startStep "Global results"
    finishStage "End of analysis"
    
fi





