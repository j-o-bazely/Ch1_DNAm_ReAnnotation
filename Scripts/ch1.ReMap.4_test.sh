#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=10G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -m beas
#$ -t 1:2

####             Ch1_ReMap.4_MethCalling         ####

# Created by: James B
# Date: 12.08.2024

# Adapted from Charley  -  /data/SBCS-EizaguirreLab/Turtle_WGBS/00_Scripts/01_WGBS_Bismark_Pipeline
# Date: June 2023


# This script takes in the deduplicated bam files from the 17 samples being re-annotated
# These bam files undergo methylation calling chromosome by chromosome, are destranded and finally merged back into full genomes
# The final output is a file containing methylation calls across the whole genome of each sample
# These files can be loaded into R for analysis via the MethylKit package


####        Environment Preparation         ####

## Load required modules ##
module load samtools/1.9 # v.1.10 has incompatible gcc with bismark module
module load bismark # v.0.22.1

## Specify the number of cores to use
REPCORES=$((NSLOTS/2))

## Move into top Ch1.ReMap directory (Assuming this script is being run from the script directory)
cd ..

## Loading sample list to create array
sample_list=$(pwd)/Metadata/Ch1.ReMap_SampleList.txt

## Extract sample ID ##
sample_list=$(pwd)/Metadata/Ch1.ReMap_SampleList.txt
SAMPLE=$(sed -n "${SGE_TASK_ID}p" $sample_list)

## Set directories
BAM_DIR=/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_dataStorage/Ch1.ReMap_data/Ch1.ReMap.2_data/bams/$SAMPLE
echo "BAM_DIR: ${BAM_DIR}"
METH_CALL_DIR=/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_dataStorage/Ch1.ReMap_data/Ch1.ReMap.4_data/$SAMPLE
echo "METH_CALL_DIR: ${METH_CALL_DIR}"
GENOME_DIR=/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_dataStorage/Ch1_inputData/Ch1_ChangGenome
echo "GENOME_DIR: ${GENOME_DIR}"
OUT_DIR=/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_dataStorage/Ch1.ReMap_data/Ch1.ReMap.4_data/destrandedCalls/$SAMPLE
echo "OUT_DIR: ${OUT_DIR}"
