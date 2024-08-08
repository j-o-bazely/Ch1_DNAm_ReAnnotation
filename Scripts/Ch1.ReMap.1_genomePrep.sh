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


# Source for genome download   : https://tinyurl.com/56j75mcs
# Paper for genome assembly    : https://tinyurl.com/3ebf4d62

# FAQs for GCA_ vs GCF_ files   : https://tinyurl.com/38n289x3
# NCBI glossary of genome terms: https://tinyurl.com/4jkevtxa


####        Environment Preparation         ####
## Load necessary modules
module load bismark # v.0.22.1

## Assign path to genome to variable GENOME_DIR
GENOME_DIR=/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_dataStorage/Ch1_inputData/Ch1_ChangGenome

## Assign path to output directory
OUTDIR=/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_dataStorage/Ch1.ReMap_data/Ch1.ReMap.1_data

## N.B The genome is downloaded in .fna format, which Bismark doesn't recognise.
#      Luckily, .fna is synonymous with .fasta format, which bismark does recognise

## Changing the genome file extension to a suitable format for bismark
mv $GENOME_DIR/Chang23Genome.fna $GENOME_DIR/Chang23Genome.fasta

## Specify the number of cores to use
REPCORES=$((NSLOTS/2))

## Run bismark_genome_preparation
bismark_genome_preparation \
--verbose --parallel $REPCORES \
$GENOME_DIR

## Move output files to correct folder
mv $GENOME_DIR/Bisulfite_Genome $OUTDIR