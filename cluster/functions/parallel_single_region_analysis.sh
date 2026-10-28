#!/bin/bash

# Exit the script if any statement returns a non-true return value.
set -e
#set -x

### The following mandatory variables are set in parallel_do_Excalibur.sh and are passed as environment variables by Slurm
#FOR R SCRIPT
##id
##PATH_DATA
##max_var
##rare_maf_threshold
##var_covariate
##var_covariate_pca
##version
##path_Excalibur_function
##RVERSION
##RLIBRARIES
##TMP_RESULT
##TMP_TIME
##load_moduleR


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


### Variables
if [ "${load_moduleR}" == 1 ]
then
    export ${RLIBRARIES}
fi
R_SCRIPT=${path_Excalibur_function}parallel_script_region_analysis.R

ulimit -c unlimited
umask 002

startStage "Parallel single region analysis"

startStep "Parallel single region analysis" $RVERSION

Rscript ${R_SCRIPT} ${id} ${PATH_DATA} ${max_var} ${rare_maf_threshold} ${var_covariate} ${var_covariate_pca} ${version} ${path_Excalibur_function} ${TMP_RESULT} ${TMP_TIME} ${RLIBRARIES}
finishStep "Parallel single region analysis"

finishStage "Parallel single region analysis"

