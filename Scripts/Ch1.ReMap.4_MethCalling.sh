#!/bin/bash
#$ -j y
#$ -cwd
#$ -pe smp 9
#$ -l h_vmem=11G
#$ -l h_rt=36:00:00
#$ -t 1-3
#$ -l highmem
#$ -t 10:17

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
METH_CALL_DIR=/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_dataStorage/Ch1.ReMap_data/Ch1.ReMap.4_data/$SAMPLE
GENOME_DIR=/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_dataStorage/Ch1_inputData/Ch1_ChangGenome
OUT_DIR=/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_dataStorage/Ch1.ReMap_data/Ch1.ReMap.4_data/destrandedCalls/$SAMPLE


####        Step 1: Methylation Calling         ####

## Load required modules ##
module load bismark # v.0.22.1

## Create directories for storing methylation calls if needed 
mkdir -p $METH_CALL_DIR

mkdir -p $OUT_DIR

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

                    # NAME CHANGE #
# Navigate to the directory containing the files
cd $METH_CALL_DIR/cytosineReports 

# Create directory for non-chromosome reads
mkdir NW_Files

# Loop through all reports
for file in *.CpG_report.txt.gz; do

    # If the file contains 'chrNW', move it to the NW_Files subdirectory 
    if [[ "$file" == *chrNW* ]]; then
        mv "$file" NW_Files
        echo "Moved $file to NW_Files/"
        continue  # Skip renaming since the file has been moved
    fi

    # Extract the chromosome part from the filename
    chr=$(echo "$file" | grep -oP 'chrNC_[0-9]+\.1')

    # Construct the new filename
    new_name="${SAMPLE}_${chr}.CpG_report.txt.gz"

    # Rename the file
    mv "$file" "$new_name"

    echo "Renamed $file to $new_name"
done

cd ../..

    #   NAME CHANGE DONE    #

mkdir -p $OUT_DIR

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
done

mv *_merged.cov.gz $OUT_DIR
mv *_both_strands.cov.gz $METH_CALL_DIR/destranding_reports/both_strands
mv *_merged.report.txt $METH_CALL_DIR/destranding_reports/$SAMPLE

echo -e "\n### DONE ###\n"


####        Step 3: Merging Destranded Chromosomes         ####
# Move into correct directory
echo 'Moving into $OUT_DIR'
cd $OUT_DIR

### Merge all chrom meth calls ###
# doing it manually rather than with *

echo 'Merging the chromosomes into one merged folder for the sample'
cat \
${SAMPLE}_chrNC_064473.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064474.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064475.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064476.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064477.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064478.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064479.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064480.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064481.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064482.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064483.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064484.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064485.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064486.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064487.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064488.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064489.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064490.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064491.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064492.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064493.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064494.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064495.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064496.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064497.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064498.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064499.CpG_merged.cov.gz \
${SAMPLE}_chrNC_064500.CpG_merged.cov.gz > $OUT_DIR/${SAMPLE}.CpG_merged.cov.gz

module unload python/3.10.7


# Generated by Job Script Builder on 2024-08-22
# For assistance, please email its-research-support@qmul.ac.uk
# Please cite use of Apocrita in your Research
# See https://docs.hpc.qmul.ac.uk/using/citing/ for details