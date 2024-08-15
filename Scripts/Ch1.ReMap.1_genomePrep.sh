## Parameters for Apocrita job submission

#!/bin/bash
#$ -pe smp 10
#$ -l h_vmem=13G
#$ -l h_rt=3:00:00
#$ -cwd
#$ -j y
#$ -m bea
#$ -l highmem


####             Ch1_ReMap.1_bismarkGenomePreparation              ####
# Written by: Charley 

# Adapted by: James B

# Original Date: 07.08.2024

# This code carries out preparation of the genome for alignment of methylation reads using bismark

# IMPORTANT: I tried this with the whole genome fasta file first, and it seems to have messed up somewhere in the annotation. 
#            I then tried running with each chromosome having it's own separate fasta file within the GENOME_DIR and it seems to have correctly indexed each chromosome
#            I found these files by diggining into teh directory on the NCBI genome website (Click on FTP -> ...
#            /genomes/all/GCF/023/653/815/GCF_023653815.1_GSC_CCare_1.0/GCF_023653815.1_GSC_CCare_1.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA


# Source for genome download   : https://tinyurl.com/56j75mcs
# Paper for genome assembly    : https://tinyurl.com/3ebf4d62

# FAQs for GCA_ vs GCF_ files   : https://tinyurl.com/38n289x3
# NCBI glossary of genome terms: https://tinyurl.com/4jkevtxa


####        Environment Preparation         ####
## Load necessary modules
module load bismark # v.0.22.1


## Assign path to genome to variable GENOME_DIR
GENOME_DIR=/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_dataStorage/Ch1.ReMap_data/Ch1.ReMap.1_data

## N.B The genome is downloaded in .fna format, which Bismark doesn't recognise.
#      Luckily, .fna is synonymous with .fasta format, which bismark does recognise.

## Changing the genome file extension to a suitable format for bismark
#mv $GENOME_DIR/Chang23Genome.fna $GENOME_DIR/Chang23Genome.fasta

## Specify the number of cores to use
REPCORES=$((NSLOTS/2))

## Run bismark_genome_preparation
bismark_genome_preparation \
--verbose --parallel $REPCORES \
$GENOME_DIR
