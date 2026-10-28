#!/bin/bash

# Exit the script if any statement returns a non-true return value.
set -e
#set -x

### The following mandatory variables are set in parallel_do_Excalibur.sh and are passed as environment variables by Slurm
#FOR R SCRIPT
##PATH_GLOBAL
##PATH_TO_RESULT
##PATH_OUTPUT
##path_Excalibur_function
##RVERSION
##RLIBRARIES
##load_moduleR
##SUB_LOG

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
R_SCRIPT=${path_Excalibur_function}parallel_script_finish.R

ulimit -c unlimited
umask 002

startStage "Parallel finish"

startStep "Parallel finish" $RVERSION
Rscript ${R_SCRIPT} ${PATH_GLOBAL} ${PATH_TO_RESULT} ${PATH_OUTPUT} ${RLIBRARIES}
finishStep "Parallel finish"


startStep "Cleaning temp folder"
#TODO rm region_log folder
rm -r ${PATH_TO_RESULT}/tmp
rm -r ${SUB_LOG}
finishStep "Cleaning temp folder"

finishStage "Parallel finish"


