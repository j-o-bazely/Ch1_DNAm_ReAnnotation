
## Parameters for Apocrita job submission
#!/bin/bash
#$ -j y
#$ -cwd
#$ -pe smp 8
#$ -l h_vmem=16G
#$ -l h_rt=240:00:00
#$ -t 1-17
#$ -l highmem

####             Ch1_ReMap.2_wholeGenomeAlignment              ####
# Written by: Charley 

# Adapted by: James B

# Original Date: 08.08.2024

# This is an array script that aligns the reads of each sample to the Chang reference genome 


# This script puts each sampe through one of two pathways:
    # 1: Runs from multiple 

####        Environment Preparation         ####
## Load necesasry modules
module load bismark # v.0.22.1
module load samtools/1.9 # v.1.10 has incompatible gcc with bismark module

## Specify the number of cores to use
REPCORES=$((NSLOTS/2))
echo "cores: ${REPCORES}"

## Move into top Ch1.ReMap directory (Assuming this script is being run from the script directory)
echo "Moving to top dir"
cd ..
echo "Currently in: $(pwd)"

## Loading sample list to create array
sample_list=$(pwd)/Metadata/Ch1.ReMap_SampleList.txt

## Extract sample ID ##
SAMPLE=$(sed -n "${SGE_TASK_ID}p" $sample_list)
echo "task ID: ${SGE_TASK_ID}"
echo "Running sample: ${SAMPLE}"

## Set directories
READS_DIR=/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_dataStorage/Ch1_inputData/Ch1.M_rawData/$SAMPLE
GENOME_DIR=/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_dataStorage/Ch1.ReMap_data/Ch1.ReMap.1_data
OUT_DIR=/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_dataStorage/Ch1.ReMap_data/Ch1.ReMap.2_data/bams/$SAMPLE

## Create outdir for bams
mkdir -p $OUT_DIR 


#################################
# Navigate to READS_DIR
echo "Moving to reads_dir"
cd $READS_DIR
echo "Currently in: $(pwd)"


####        Pipeline 1: Files are un-merged         ####
# Sequencing was performed across multiple flow cells and some samples have yet to be merged. In this case, we check for 3+ files in the sample data directory. 
#  If a sample was run across multiple lanes which haven't been merged, then this will be the case


# Count the number of files in the directory
FILE_COUNT=$(find . -maxdepth 1 -type f | wc -l)

echo "This sample has ${FILE_COUNT} files"

if [ "$FILE_COUNT" -gt 2 ]; then
    echo "There are more than 2 files for this sample, so merging is needed"
# Check if merged file exists yet -> only run alignment if it doesn't
    if [[ -f ${OUT_DIR}/${SAMPLE}_deduplicated.bam ]]; then
    echo -e "\n### Merged bam already exists for ${SAMPLE}. Not re-aligning ###\n"
    else
    echo -e "\n### Merged bam does not exist for ${SAMPLE}. Continue to alignment ###\n"
    
    # Save unique flowcell info from fastq filenames to variable
    FLOWCELLS=$(ls $READS_DIR | cut -d '_' -f 1,2,3 | uniq )
    echo "The flowcell IDs for this sample are: ${FLOWCELLS}"

    for FLOWCELL in $FLOWCELLS

    do
    # Read group ID in the format FlowCell.Lane (replace underscore in filename)
    rgid=$(echo $FLOWCELL | cut -d '_' -f 1,2 | sed 's/_/\./g')

    echo -e "\n### Starting alignment for flowcell: ${FLOWCELL}, sample: ${SAMPLE} ###\n"

    bismark \
    -p ${REPCORES} \
    --rg_sample $SAMPLE --rg_id $rgid \
    -o $OUT_DIR \
    --basename $FLOWCELL \
    $GENOME_DIR \
    -1 ${FLOWCELL}*_1_trim.fq.gz \
    -2 ${FLOWCELL}*_2_trim.fq.gz 

    echo -e "\n### Finished alignment for flowcell: ${FLOWCELL}, sample: ${SAMPLE} ###\n"
    
        ### Run Bismark deduplication ### 
    echo -e "\n### Starting deduplication for flowcell: ${FLOWCELL}, sample: ${SAMPLE} ###\n"

    deduplicate_bismark --paired --bam \
    --output_dir $OUT_DIR \
    --outfile ${FLOWCELL}_pe \
    ${OUT_DIR}/${FLOWCELL}_pe.bam

    echo -e "\n### Finished deduplication for flowcell: ${FLOWCELL}, sample: ${SAMPLE} ###\n"
    done

    ### Merge deduplicated bams from all flowcells ###
    
    echo -e "\n### Starting merging deduplicated flowcell bams for ${SAMPLE} ###\n"
    
    samtools merge -@ ${NSLOTS} \
    ${OUT_DIR}/${SAMPLE}_deduplicated.bam \
    ${OUT_DIR}/*_pe.deduplicated.bam
    
    echo -e "\n### Finished merging deduplicated flowcell bams for ${SAMPLE} ###\n"
    
    fi

    ### Delete all flowcell bams, only if merged bam exists ###
    
    if [[ -f ${OUT_DIR}/${SAMPLE}_deduplicated.bam ]]; then
    echo -e "\n### Merged bam for has been created for ${SAMPLE} ###"
    echo -e "\n### Deleting flowcell bams ###\n"
    
    # Remove flowcell bams whilst keeping all reports
    rm ${OUT_DIR}/V*.bam
    
    echo -e "### Deleted flowcell bams for ${SAMPLE} ###\n"
    
    else
    echo -e "\n### Merged bam does not exist for ${SAMPLE} ###"
    echo -e "### Will not delete flowcell bams ###\n"
    
    fi # End if statement for deleting flowcell bams


    ### Sort merged deduplicated bam by name ###

    # -n flag sorts by read name (required for meth calling)

    echo -e "### Sorting merged, deduplicated bam by name for ${SAMPLE} ###"

    samtools sort -n -@ ${NSLOTS} \
    -o ${OUT_DIR}/${SAMPLE}_deduplicated.sorted_by_name.bam \
    ${OUT_DIR}/${SAMPLE}_deduplicated.bam

    echo -e "### Sorted merged, deduplicated bam by name for ${SAMPLE} ###"


    if [[ -f ${OUT_DIR}/${SAMPLE}_deduplicated.sorted_by_name.bam ]]; then
    echo -e "\n### Sorted bam for has been created for ${SAMPLE} ###"
    echo -e "\n### Deleting unsorted bam ###\n"

    # Delete
    rm ${OUT_DIR}/${SAMPLE}_deduplicated.bam
    
    echo -e "### Deleted unsorted bam for ${SAMPLE} ###\n"
    
    else
    echo -e "\n### Sorted bam does not exist for ${SAMPLE} ###"
    echo -e "### Will not delete unsorted bams ###\n"
    fi


    echo -e "### All done ###"


else


####        Pipeline 2: Files are merged but not correctly named         ####

## Some samples have been pre-merged, but the flowcell IDs have been retained, rather than being replaced by sample ID
## This pipeline first renames the files to match the sample ID


## WORK IN PROGRESS


####        Pipeline 3: Files are merged AND correctly named            ####

## Now all merged files should be 

##

if [[ -f ${OUT_DIR}/${SAMPLE}_deduplicated.bam ]]; then
    echo -e "\n### Merged bam already exists for ${SAMPLE}. Not re-aligning ###\n"
    else
    echo -e "\n### Merged bam does not exist for ${SAMPLE}. Continue to alignment ###\n"

### Run Bismark alignment ###

echo -e "\n### Starting alignment for sample: ${SAMPLE} ###\n"

bismark \
-p ${REPCORES} \
--basename $SAMPLE \
-o $OUT_DIR \
$GENOME_DIR \
-1 ${SAMPLE}_1_trim.fq.gz \
-2 ${SAMPLE}_2_trim.fq.gz

echo -e "\n### Finished alignment for sample: ${SAMPLE} ###\n"

# NB. Automatically has _pe.bam added to end of basename by default

  
### Run Bismark deduplication ### 
echo -e "\n### Starting deduplication for sample: ${SAMPLE} ###\n"
  
deduplicate_bismark --paired --bam \
--output_dir $OUT_DIR \
--outfile ${SAMPLE} \
${OUT_DIR}/${SAMPLE}_pe.bam

# NB. basename has .deduplicated.bam added to end by default
  
echo -e "\n### Finished deduplication for sample: ${SAMPLE} ###\n"
  


### Sort deduplicated bam by name ###

# -n flag sorts by read name (required for meth calling)

echo -e "### Sorting deduplicated bam by name for ${SAMPLE} ###"

samtools sort -n -@ ${NSLOTS} \
-o ${OUT_DIR}/${SAMPLE}_deduplicated.sorted_by_name.bam \
${OUT_DIR}/${SAMPLE}.deduplicated.bam

echo -e "### Sorted deduplicated bam by name for ${SAMPLE} ###"


### Delete unsorted bam ###

if [[ -f ${OUT_DIR}/${SAMPLE}_deduplicated.sorted_by_name.bam ]]; 
  then
      echo -e "\n### Sorted bam for has been created for ${SAMPLE} ###"
      echo -e "\n### Deleting unsorted bam ###\n"

      # Delete
      rm ${OUT_DIR}/${SAMPLE}.deduplicated.bam
      rm ${OUT_DIR}/${SAMPLE}_pe.bam
     
      echo -e "### Deleted unsorted bam for ${SAMPLE} ###\n"
       
   else
       echo -e "\n### Sorted bam does not exist for ${SAMPLE} ###"
       echo -e "### Will not delete unsorted bams ###\n"
fi

echo -e "### All done ###"