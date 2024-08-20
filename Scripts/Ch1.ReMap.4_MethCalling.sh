#!/bin/bash
#$ -j y
#$ -cwd
#$ -pe smp 10
#$ -l h_vmem=13G
#$ -l h_rt=36:00:00
#$ -t 1-3
#$ -l highmem
#$ -t 1:17

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

## Extract sample ID ##
sample_list=/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_DNAm_ReAnnotation/Metadata/Ch1.ReMap_SampleList.txt
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


####        Step 1: Methylation Calling         ####

## Load required modules ##
module load bismark # v.0.22.1

## Create directories for storing methylation calls if needed 
mkdir -p "$METH_CALL_DIR"

## Perform Bismark alignment ##
bismark_methylation_extractor \
--paired-end --comprehensive \
--split_by_chromosome --no_header --gzip \
--cytosine_report \
--output $METH_CALL_DIR \
--genome_folder $GENOME_DIR \
"$BAM_DIR/${SAMPLE}_deduplicated.sorted_by_name.bam"

echo -e "\n### Methylation calling finished for ${SAMPLE} ###\n"


### Clean up output ###
cd $METH_CALL_DIR
echo -e "### Entered directory: $PWD ###\n"
echo -e "\n### Cleaning up output files for ${SAMPLE} ###\n"


# Remove CHG and CHH context files to save space
rm CHG_context_*
rm CHH_context_*
echo -e "\n### Removed CHG and CHH context files ###\n"


# Rename cytosine report files to simple name with just SAMPLE_CHROM.CpG_report.txt.gz
for i in *deduplicated*.CpG_report.txt.gz; do
    mv $i ${i/"deduplicated.sorted_by_name.CpG_report.txt.chrSLK063_ragtag_"}
done

echo -e "\n### Renamed cytosine report files per chromosome ###\n"


# Create directory for cytosine reports
mkdir -p cytosineReports
mv *.CpG_report.txt.gz cytosineReports

echo -e "\n### Moved cytosine report files to new directory CytosineReports ###\n"


# Move reports to reports folder 
mkdir -p splitting_reports
mv *_splitting_report.txt splitting_reports

echo -e "\n### Moved splitting reports to splitting_reports directory ###\n"

# Unload modules to prevent clashes and save memory
module unload bismark
module unload samtools/1.9

####        Step 2: Destranding        ####
## Load necessary modules ##
module load python/3.10.7

### Set scripts ###
MERGE_CPG=/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_DNAm_ReAnnotation/Scripts/merge_CpG.py

### Navigate to directory containing meth calls
cd $METH_CALL_DIR

mkdir -p destranding_reports
mkdir -p destranding_reports/both_strands # Put both_strands files into diff directory, as not needed
mkdir -p destranding_reports/$SAMPLE

## Navigate to directory containing cytosine reports
cd cytosineReports

## Run loop for CytosineReports of all chromosomes ##
for i in *CpG_report.txt.gz; do
  # Run destranding script
  python $MERGE_CPG --input $i

  echo -e "\n### Finished destranding: $i ###"

  # Tidy output
  gzip *.cov # Zip output files
  chmod a+w *.cov.gz # Add permissions

  mv *_merged.cov.gz $OUT_DIR
  mv *_both_strands.cov.gz $METH_CALL_DIR/destranding_reports/both_strands
  mv *_merged.report.txt $METH_CALL_DIR/destranding_reports/$SAMPLE
done

echo -e "\n### DONE ###\n"

####        Step 3: Merging Destranded Chromosomes         ####
## Create directory if needed ##
mkdir -p $OUT_DIR/$SAMPLE

# Move into correct directory
cd $OUT_DIR/$SAMPLE

### Merge all chrom meth calls ###
# With chr0 at end -> doing it manually rather than with *

cat \
${SAMPLE}_chr1.CpG_merged.cov.gz \
${SAMPLE}_chr2.CpG_merged.cov.gz \
${SAMPLE}_chr3.CpG_merged.cov.gz \
${SAMPLE}_chr4.CpG_merged.cov.gz \
${SAMPLE}_chr5.CpG_merged.cov.gz \
${SAMPLE}_chr6.CpG_merged.cov.gz \
${SAMPLE}_chr7.CpG_merged.cov.gz \
${SAMPLE}_chr8.CpG_merged.cov.gz \
${SAMPLE}_chr9.CpG_merged.cov.gz \
${SAMPLE}_chr10.CpG_merged.cov.gz \
${SAMPLE}_chr11.CpG_merged.cov.gz \
${SAMPLE}_chr12.CpG_merged.cov.gz \
${SAMPLE}_chr13.CpG_merged.cov.gz \
${SAMPLE}_chr14.CpG_merged.cov.gz \
${SAMPLE}_chr15.CpG_merged.cov.gz \
${SAMPLE}_chr16.CpG_merged.cov.gz \
${SAMPLE}_chr17.CpG_merged.cov.gz \
${SAMPLE}_chr18.CpG_merged.cov.gz \
${SAMPLE}_chr19.CpG_merged.cov.gz \
${SAMPLE}_chr20.CpG_merged.cov.gz \
${SAMPLE}_chr21.CpG_merged.cov.gz \
${SAMPLE}_chr22.CpG_merged.cov.gz \
${SAMPLE}_chr23.CpG_merged.cov.gz \
${SAMPLE}_chr24.CpG_merged.cov.gz \
${SAMPLE}_chr25.CpG_merged.cov.gz \
${SAMPLE}_chr26.CpG_merged.cov.gz \
${SAMPLE}_chr27.CpG_merged.cov.gz \
${SAMPLE}_chr28.CpG_merged.cov.gz > $OUT_DIR/${SAMPLE}.CpG_merged.cov.gz

module unload python/3.10.7
