## Parameters for Apocrita job submission
    # Will simply be commented out if not running on apocrita

#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -t 1-17
#$ -tc 1


####             Ch1_ReMap.3_trimmedQC              ####
# Written by: James B

# Original Date: 07.08.2024

# This code creates necessary folders for storing sample data and associated QC results of trimmed methylation reads
# before running QC analysis on each sample

# Tutorial for using fastQC: https://tinyurl.com/rut9jypw


####        Environment Preparation         ####

## Moving into top Ch1.ReMap directory from script directory (Assuming this script is being run from the script directory)
cd ..

## Loading sample list to create array
sample_list=$(pwd)/Metadata/Ch1.ReMap_SampleList.txt

## Pulling out sample IDs
    # cat =                 # sed =
ID=$( cat "$sample_list" | sed -n ${SGE_TASK_ID}p )

## Creating directory for each sample
mkdir $(pwd)/QC/${ID}

## Setting directory paths in and out of each of these folders
IN=$(pwd)/../../Ch1_dataStorage/Ch1.M_rawData/${ID}
OUT=$(pwd)/QC/${ID}

## Move into the directory of each sample in preparation for running QC reports
cd ${IN}

## If required, loading fastqc module
module load fastqc/0.11.9




####            Object Generation           ####

## Run fastQC on the current sample
fastqc -o ${OUT} -f fastq -t 1 ${ID}_1_trim.fq.gz ${ID}_2_trim.fq.gz
    #  -o   : Specifies output directory
    #  -f   : Specifies the format of the input files
    # fastq : The format of the input files
    #  -t   : Argument for the number of threads to be used
    #   1   : Specifies 1 thread to be used
    # (Two files provided as we have paired end read data)


## Unloading the fastQC module
module unload fastqc/0.11.9
