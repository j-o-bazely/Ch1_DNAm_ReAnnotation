## Parameters for Apocrita job submission

#!/bin/bash
#$ -pe smp 10
#$ -l h_vmem=13G
#$ -l h_rt=1:00:00
#$ -cwd
#$ -j y
#$ -m bea
#$ -l highmem


####             Ch1_ReMap.1_bismarkGenomePreparation              ####
# Written by: Charley 

# Adapted by: James B

# Original Date: 07.08.2024

# This code carries out preparation of the genome for alignment of methylation reads using bismark

# FAQ for GCA_ vs GCF_ files   : https://tinyurl.com/38n289x3

# NCBI glossary of genome terms: https://tinyurl.com/4jkevtxa


####        Environment Preparation         ####
## Load necessary modules
module load bismark # v.0.22.1

## Assign path to genome to variable GENOME_DIR
GENOME_DIR=/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_dataStorage/Ch1_inputData/Ch1_ChangGenome

## Assign path to output directory
OUTDIR=/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_dataStorage/Ch1.ReMap_data/Ch1.ReMap.1_data

## Specify the number of cores to use
REPCORES=$((NSLOTS/2))

## Run bismark_genome_preparation
bismark_genome_preparation \
--verbose --parallel $REPCORES \
$GENOME_DIR

## Move output files to correct folder
mv $GENOME_DIR/*.bt2 $OUTDIR