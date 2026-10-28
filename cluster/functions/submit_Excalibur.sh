#!/bin/bash
#SBATCH --chdir /scratch/sboutry
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --qos=serial
#SBATCH --mem 4G
#SBATCH --time 00:59:00
#SBATCH --mail-user=YOUR_EMAIL
#SBATCH -J submit_excalibur
#SBATCH --mail-type=FAIL
#SBATCH -o /scratch/sboutry/excalibur/logs/submit_log/submit_%J.out # Standard output
#SBATCH -e /scratch/sboutry/excalibur/logs/submit_log/submit_%J.err # Standard error

###Usefull for debuging 
set -e # Exit the script if any statement returns a non-true return value.
#set -x #Print each line of code being computed
#set -u #Exit the script if a variable was not initialized

usage() {
    echo "Usage: submit_Excalibur.sh [OPTIONS] projectdir"
    echo ""
    echo "-h : usage"
    echo "-A : path and name file of the cases data [if -I is set to 1, should be a file with cases ID]"
    echo "-B : path and name file of the controls data [if -I is set to 1, should be a file with controls ID]"
    echo "-E : email [mandatory]"
    echo "-G : path and name of GVCF file for the analysis [only used if importation method is set to 1]"
    echo "-0 : column name in GVCF that contain sample ID [by default is Indiv, only used if importation method is set to 1]"
    echo "-1 : column name in GVCF that contain chromosome info [by default is CHROM, only used if importation method is set to 1]"
    echo "-2 : column name in GVCF that contain position info [by default is POS, only used if importation method is set to 1]"
    echo "-3 : column name in GVCF that contain reference info [by default is REF, only used if importation method is set to 1]"
    echo "-4 : column name in GVCF that contain alternative info [by default is ALT, only used if importation method is set to 1]"
    echo "-5 : column name in GVCF that contain zyosity info [by default is will be computed based on gt_GT column from GCVF, only used if importation method is set to 1]"
    echo "-6 : column name in GVCF that contain gene_symbol info [by default is SYMBOL, only used if importation method is set to 1]"
    echo "-7 : column name in GVCF that contain allele frequency info [by default is gnomad40_AF, only used if importation method is set to 1]"
    echo "-8 : column name in GVCF that contain genotype info [by default is gt_GT, only used if importation method is set to 1]"
    echo "-9 : column name in GVCF that contain coverage info [by default is gt_DP, only used if importation method is set to 1]"
    echo "-I : Importation mehtod, if set to 0 the inputs are two tsv files [default], if set to 1 one VCF or GVCF"
    echo "-M : depending on your module on your cluster, you might need to do module load or module use, you can test that while installing the R package using the install_r_lib.sh function [default is module use PATHTOYOURRVERSION followed by module load r]"
    echo "-P : max number of jobs that can be launched [default 1 000]"
    echo "-Q partition : Slurm partition (queue) to use [default is parallel]"
    echo "-R : R module version to use [default is Extended-R-bundle-Bioconductor/3.15-foss-2022a-R-4.2.1]"
    echo "-L : which libraries of R to use [default is /storage/research/dduv/gehu-mv/rlib/4.2]"
    echo "-b : supply here full path of file and extension (.tsv) containing a gene list to be analyzed (must be HUGO annotation) [default is all genes present in patient and control data will be analyzed]"
    echo "-c : filtering to be applied for further analysis [TO BE DONE]"
    echo "-d : filtering level to be applied for further analysis [TO BE DONE]"
    echo "-e : filtering direction to be applied for further analysis [TO BE DONE]"
    echo "-f : Minimum number of unique variant to be present within a genetic region for the region to be analyzed [default = 2]"
    echo "-g : Maximum number of unique variant to be present within a genetic region for the region to be analyzed [default = 100000]"
    echo "-i : rare MAF threshold to be used to separate rare variants from common variant [default is set to 0.01]"
    echo "-j : submit full path of file and file extension of covariates to be included in the analysis. Be carefull that all individuals from the analysis are present in that file [TO BE DONE]"
    echo "-k : if set to 1, then argument -j will be ignored and the covariates will be the principal component computed based on the cohort data, otherwise set it to 0 [default is set to 1]"
    echo "-m : mapqv threshold to determine if a position is well covered or not [default mapqv_threshold = 20]"
    echo "-n : minimum coverage for a position to be well covered or not [default = 20]"
    echo "-o : remove individu that are badly covered [default is 1=yes]"
    echo "-p : remove position that are badly covered amongst the cohort [default is 1=yes]"
    echo "-q : perform data annotation and cleaning [default is yes]"
    echo "-r : path to table to annotate data with pathway and go terms [default is /storage/service/excalibur/data/Data/]"
    echo "-s : only use that argument if you only want to perform aggregation test on over represented pathways. By default, perform aggregation on all pathways that have at least two genes in your input data"
    echo "-u : version of Excalibur to use [default "optimal" verison, "all" is Excalibur_baseline, then you can put a list of method name to use specific methods]"
    echo "-v : path to Excalibur functions [default is ${EXCALIBUR}/data/functions/]"
    echo "-w : path to Excalibur Data [default is ${EXCALIBUR}/data/Data/]"
    echo "-x : remove private variant [default is no, private variant are taken into account]"
    echo ""
    echo "projectdir : directory of the project where to store the results. This directory must contains two tsv files named patient_data and control data (if all data for control and patient are already been retrieved) OR patient_list and control_list (if only the names of individuals should be use)((this option is mandatory)."
    echo ""
    echo "submit_Excalibur.sh submit the analysis of Excalibur on the given patients and controls data to the cluster scheduler."
    echo ""
}

## Formatting and colors (use echo -e)
RESET='\033[0m'
BOLD='\033[1m'
UL='\033[4m'
INV='\033[7m'
RED='\033[38;5;1m'
GREEN='\033[38;5;34m'
BLUE='\033[38;5;21m'
YELLOW='\033[38;5;226m'
MAGENTA='\033[38;5;165m'
PURPLE='\033[38;5;129m'
ORANGE='\033[38;5;208m'
PINK='\033[38;5;198m'
GREY='\033[38;5;246m'


if [ $# -eq 0 ]
then
    usage
    exit
fi

while getopts "hA:B:E:G:0:1:2:3:4:5:6:7:8:9:I:MP:Q:R:L:b:c:d:e:f:g:i:j:km:n:opqr:su:v:w:x" OPTION
do
    case $OPTION in
        h)
            usage
            exit 1
            ;;
        A)
            PATH_PATIENTS="$OPTARG"
            ;;
        B)
            PATH_CONTROLS="$OPTARG"
            ;;
        E)
            EMAIL="$OPTARG"
            ;;
        G)
            path_GVCF="$OPTARG"
            ;;
        0)
            sample_GVCF="$OPTARG"
            ;;
        1)
            chr_GVCF="$OPTARG"
            ;;
        2)
            pos_GVCF="$OPTARG"
            ;;
        3)
            ref_GVCF="$OPTARG"
            ;;
        4)
            alt_GVCF="$OPTARG"
            ;;
        5)
            zyosity_GVCF="$OPTARG"
            ;;
        6)
            gene_symbol_GVCF="$OPTARG"
            ;;
        7)
            frequency_GVCF="$OPTARG"
            ;;
        8)
            GT_GVCF="$OPTARG"
            ;;
        9)
            DP_GVCF="$OPTARG"
            ;;
        I)
            importation_method="$OPTARG"
            ;;
        M)
            load_moduleR=1
            ;;
        P)
            MAX_JOB="$OPTARG"
            ;;
        Q)
            QUEUE="$OPTARG"
            ;;
        R)
            RVERSION="$OPTARG"
            ;;
        L)
            RLIBRARIES="$OPTARG"
            ;;
        b)
            know_target_region="$OPTARG"
            ;;
        c)
            filtering_criteria="$OPTARG"
            ;;
        d)
            filtering_level="$OPTARG"
            ;;
        e)
            filtering_direction="$OPTARG"
            ;;
        f)
            min_var="$OPTARG"
            ;;
        g)
            max_var="$OPTARG"
            ;;
        i)
            rare_maf_threshold="$OPTARG"
            ;;
        j)
            var_covariate="$OPTARG"
            ;;
        k)
            var_covariate_pca=0
            ;;
        m)
            mapqv_threshold="$OPTARG"
            ;;
        n)
            min_coverage="$OPTARG"
            ;;
        o)
            remove_individual=0
            ;;
        p)
            remove_position=0
            ;;
        q)
            do_data_preparation=0
            ;;
        r)
            path_data_annotate="$OPTARG"
            ;;
        s)
            DO_over_representation=1
            ;;
        u)
            version="${OPTARG}"
            ;;
        v)
            path_Excalibur_function="${OPTARG}"
            ;;
        w)
            path_Excalibur_Data="${OPTARG}"
            ;;
        x)
            remove_private=1
            ;;
        ?)
            usage
            exit
            ;;
    esac
done

shift $(($OPTIND - 1))
if [ "$1" == "" ]
then
    usage
    exit
fi

if [ "$EMAIL" == "" ]
    then
        echo "email adress must be provided, please read usage -E"
        exit
fi
### TODO Change the name of the module to load or use for the right version of R on your cluster
if [ "$RVERSION" == "" ]
then
    RVERSION=/work/scitas-share/spack-r-gr-fe/share/spack/lmod/linux-rhel8-x86_64/Core
fi
### TODO Change here path to the right R libraries to be used
if [ "$RLIBRARIES" == "" ]
then
    RLIBRARIES="/work/gr-fe/R_4.3.1"
fi


### store arguments
PROJECTDIR=$(basename $1)
###Get date and time of run and add it to the name of the analysis
#DATE=`date +"%Y_%m_%d_%H_%M_%S"`
#PROJECTDIR="${DATE}_${PROJECTDIR}"

### paths
EXCALIBUR=/work/gr-fe/sboutry/excalibur
EXCALIBUR_SCRATCH=/scratch/sboutry/excalibur
LOG=${EXCALIBUR_SCRATCH}/logs/${PROJECTDIR}/
INPUT=${EXCALIBUR}/input/${PROJECTDIR}/
RESULTS=${EXCALIBUR_SCRATCH}/results/${PROJECTDIR}/
ZIP=${EXCALIBUR}/zip/

### We create a new folder within the logs folder for each analysis
#ulimit -c unlimited
#umask 002

mkdir -p ${LOG}
chmod 775 ${LOG} || true
mkdir -p ${RESULTS}
chmod 775 ${RESULTS} || true


### Retrieve error from this script
touch "${LOG}/submit_Excalibur.txt"
LOG_submit="${LOG}/submit_Excalibur.txt"
chmod 664 ${LOG_submit}
echo "assigning script variable" > ${LOG_submit}
if [ "$PATH_PATIENTS" == "" ]
then
    PATH_PATIENTS="${INPUT}/patient_data.tsv"
fi
if [ "$PATH_CONTROLS" == "" ]
then
    PATH_CONTROLS="${INPUT}/control_data.tsv"
fi
if [ "$path_GVCF" == "" ]
then
    path_GVCF=0
fi
if [ "$sample_GVCF" == "" ]
then
    sample_GVCF="Indiv"
fi
if [ "$chr_GVCF" == "" ]
then
    chr_GVCF="CHROM"
fi
if [ "$pos_GVCF" == "" ]
then
    pos_GVCF="POS"
fi
if [ "$ref_GVCF" == "" ]
then
    ref_GVCF="REF"
fi
if [ "$alt_GVCF" == "" ]
then
    alt_GVCF="ALT"
fi
if [ "$zyosity_GVCF" == "" ]
then
    zyosity_GVCF=0
fi
if [ "$gene_symbol_GVCF" == "" ]
then
    gene_symbol_GVCF="SYMBOL"
fi
if [ "$frequency_GVCF" == "" ]
then
    frequency_GVCF="gnomad40_AF"
fi
if [ "$GT_GVCF" == "" ]
then
    GT_GVCF="gt_GT"
fi
if [ "$DP_GVCF" == "" ]
then
    DP_GVCF="gt_DP"
fi
if [ "$load_moduleR" == "" ]
then
    load_moduleR=0
fi
if [ "$MAX_JOB" == "" ]
then
    MAX_JOB=1000
fi
if [ "$QUEUE" == "" ]
then
    #QUEUE=pipeline
    #QUEUE=parallel
    QUEUE=serial
fi
if [ "$know_target_region" == "" ]
then
    know_target_region=0
fi
if [ "$filtering_criteria" == "" ]
then
    filtering_criteria=0
fi
if [ "$filtering_level" == "" ]
then
    filtering_level=0
fi
if [ "$filtering_direction" == "" ]
then
    filtering_direction=0
fi
if [ "$min_var" == "" ]
then
    min_var=2
fi
if [ "$max_var" == "" ]
then
    max_var=100000
fi
if [ "$rare_maf_threshold" == "" ]
then
    rare_maf_threshold=0.01
fi
if [ "$var_covariate" == "" ]
then
    var_covariate=0
fi
if [ "$var_covariate_pca" == "" ]
then
    var_covariate_pca=1
fi
if [ "$importation_method" == "" ]
then
    importation_method=0
fi
if [ "$mapqv_threshold" == "" ]
then
    mapqv_threshold=20
fi
if [ "$min_coverage" == "" ]
then
    min_coverage=20
fi
if [ "$remove_individual" == "" ]
then
    remove_individual=1
fi
if [ "$remove_position" == "" ]
then
    remove_position=1
fi
if [ "$do_data_preparation" == "" ]
then
    do_data_preparation=1
fi
if [ "$path_data_annotate" == "" ]
then
    path_data_annotate="${EXCALIBUR}/data/Data/"
fi
if [ "$DO_over_representation" == "" ]
then
    DO_over_representation=0
fi
if [ "$version" == "" ]
then
    version="optimal"
fi
if [ "$path_Excalibur_function" == "" ]
then
    path_Excalibur_function="${EXCALIBUR}/data/functions/"
fi
if [ "$path_Excalibur_Data" == "" ]
then
    path_Excalibur_Data="${EXCALIBUR}/data/Data/"
fi
if [ "$remove_private" == "" ]
then
    remove_private=0
fi
echo -e "DONE assigning script variable\n" >> ${LOG_submit}


### Variables
PATH_SAMPLES=0


PROJECT=${PROJECTDIR}

#MAILING=/storage/service/administration/send_results.sh
ZIPPING=${path_Excalibur_function}do_zip.sh


###Gene list 
if [ "${know_target_region}" != "0" ]
then
    echo -e "Start importing gene_list\n" >> ${LOG_submit}
    LIST_GENE=""
    while IFS= read line || [ -n "${line}" ]
    do
        LIST_GENE+="'${line/$'\r'/}',"
    done < "${know_target_region}"
    LIST_GENE=${LIST_GENE::-1}
    echo -e "gene_list imported\n" >> ${LOG_submit}
#else
#    # Send message to standard error.
#    echo -e "gene_list file not found\n" >> ${LOG_submit}
#    # Exit with EX_DATAERR from sysexits.h.
#    exit 65
fi



### Submission
#put here variable mandatory for do_excalibur that are not in args of submit

#temporary --> when /storage/service/Excalibur exists on cluster, asign that value to PROGRAM
#mkdir -p ${PROGRAM}
#chmod 775 ${PROGRAM}


## Phase specific variables
#RUNNINGDIR=/tmp
RUNNINGDIR=/tmp
MEM_job=10g
THREADS_FINAL=12
MEM_FINAL=80g


#GROUP=$(/storage/service/administration/derive_group_path.sh "$PROJECT")
#echo -e "Group detected = '$GROUP'\n" >> ${LOG_submit}
#if [ "$GROUP" != "" ]
#then
#    ZIP=""
#    echo -e "## ${PINK}Creating symlinks of the results in service directory${RESET}"
#    sudo /storage/service/administration/link_data_to_service.sh $RESULTS
#fi


# Define proxy for all scripts (users on sss-intranet don't have those variables defined)
#export http_proxy=http://sss-proxy.icp.ucl.ac.be:3128
#export https_proxy=http://sss-proxy.icp.ucl.ac.be:3128
#export ftp_proxy=http://sss-proxy.icp.ucl.ac.be:3128

#Slurm
JOB_submit=`sbatch --parsable --job-name=${PROJECT} --qos=${QUEUE} --chdir=${RUNNINGDIR} --mail-user=${EMAIL} --mail-type=FAIL,REQUEUE --ntasks=1 --cpus-per-task=${THREADS_FINAL} --mem=${MEM_FINAL} --time=70:00:00 --output=${LOG}%j.%x --export=ALL,PROJECT=${PROJECT},MAX_JOB=${MAX_JOB},QUEUE=${QUEUE},RUNNINGDIR=${RUNNINGDIR},EMAIL=${EMAIL},THREADS_FINAL=${THREADS_FINAL},MEM_FINAL=${MEM_FINAL},MEM_job=${MEM_job},LOG=${LOG},PATH_SAMPLES=${PATH_SAMPLES},PATH_PATIENTS=${PATH_PATIENTS},PATH_CONTROLS=${PATH_CONTROLS},path_GVCF=${path_GVCF},sample_GVCF=${sample_GVCF},chr_GVCF=${chr_GVCF},pos_GVCF=${pos_GVCF},ref_GVCF=${ref_GVCF},alt_GVCF=${alt_GVCF},zyosity_GVCF=${zyosity_GVCF},gene_symbol_GVCF=${gene_symbol_GVCF},frequency_GVCF=${frequency_GVCF},GT_GVCF=${GT_GVCF},DP_GVCF=${DP_GVCF},know_target_region=${know_target_region},filtering_criteria=${filtering_criteria},filtering_level=${filtering_level},filtering_direction=${filtering_direction},min_var=${min_var},max_var=${max_var},rare_maf_threshold=${rare_maf_threshold},var_covariate=${var_covariate},var_covariate_pca=${var_covariate_pca},importation_method=${importation_method},mapqv_threshold=${mapqv_threshold},min_coverage=${min_coverage},remove_individual=${remove_individual},remove_position=${remove_position},INPUT=${INPUT},RESULTS=${RESULTS},do_data_preparation=${do_data_preparation},path_data_annotate=${path_data_annotate},DO_over_representation=${DO_over_representation},version=${version},path_Excalibur_function=${path_Excalibur_function},path_Excalibur_Data=${path_Excalibur_Data},RVERSION=${RVERSION},RLIBRARIES=${RLIBRARIES},load_moduleR=${load_moduleR},remove_private=${remove_private},ZIP=${ZIP},MAILING=${MAILING},ZIPPING=${MAILING} ${path_Excalibur_function}parallel_do_Excalibur.sh`
DEPENDNOTOK="afternotok:${JOB_submit}"


echo -e "Job submitted -- ALL DONE" >> ${LOG_submit}
