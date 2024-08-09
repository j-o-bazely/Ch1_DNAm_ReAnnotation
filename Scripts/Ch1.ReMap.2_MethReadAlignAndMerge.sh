## Parameters for Apocrita job submission
#!/bin/bash
#$ -j y
#$ -cwd
#$ -pe smp 8
#$ -l h_vmem=16G
#$ -l h_rt=240:00:00
#$ -t 1-2
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

## Move into top Ch1.ReMap directory (Assuming this script is being run from the script directory)
cd ..

## Loading sample list to create array
sample_list=$(pwd)/Metadata/Ch1.ReMap_SampleList.txt

## Extract sample ID ##
sample_list=$(pwd)/Metadata/Ch1.ReMap_SampleList.txt
SAMPLE=$(sed -n "${SGE_TASK_ID}p" $Sample_File)



## Set directories
READS_DIR=/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_dataStorage/Ch1_inputData/Ch1.M_rawData/$SAMPLE
GENOME_DIR=/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_dataStorage/Ch1_inputData/Ch1_ChangGenome
OUT_DIR=/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_dataStorage/Ch1.ReMap_data/Ch1.ReMap.2_data/bams/$SAMPLE

## Create outdir for bams
mkdir -p $OUT_DIR 


##############################
# Navigate to READS_DIR
cd $READS_DIR
  
####        Checkpoint 1: Merged lanes         ####

# Some samples have yet to be merged. In this case, we check for 3+ files in the sample data directory. 
#  If a sample was run across multiple lanes which haven't been merged, then this will be the case


# Count the number of files in the directory
FILE_COUNT=$(find . -maxdepth 1 -type f | wc -l)

# Check if there are 3 or more files
if [ "$FILE_COUNT" -ge 3 ]; then

# Check if merged file exists yet -> only run alignment if it doesn't
    if [[ -f ${OUT_DIR}/${SAMPLE}_deduplicated.bam ]]; 
    then
    echo -e "\n### Merged bam already exists for ${SAMPLE}. Not re-aligning ###\n"
    else
    echo -e "\n### Merged bam does not exist for ${SAMPLE}. Continue to alignment ###\n"

    # Navigate to READS_DIR
    cd $READS_DIR
    
        # Save unique flowcell info from fastq filenames to variable
        FLOWCELLS=$(ls $READS_DIR | cut -d '_' -f 1,2,3 | uniq )

        # For each flowcell read pair, take read group ID (rgid) from file name
        # Sample name is specified (smid) and bismark is run
    
    for FLOWCELL in $FLOWCELLS

        do
        # Read group ID in the format FlowCell.Lane (replace underscore in filename)
        rgid=$(echo $FLOWCELL | cut -d '_' -f 1,2 | sed 's/_/\./g')
    
        ### Run Bismark alignment ###
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
    
    done # End flowcell loop
    

    ### Merge deduplicated bams from all flowcells ###
    
    echo -e "\n### Starting merging deduplicated flowcell bams for ${SAMPLE} ###\n"
    
    samtools merge -@ ${NSLOTS} \
    ${OUT_DIR}/${SAMPLE}_deduplicated.bam \
    ${OUT_DIR}/*_pe.deduplicated.bam
    
        echo -e "\n### Finished merging deduplicated flowcell bams for ${SAMPLE} ###\n"
    
    fi # End if statement for alignment and deduplication
    


    ### Delete all flowcell bams, only if merged bam exists ###
    
    if [[ -f ${OUT_DIR}/${SAMPLE}_deduplicated.bam ]]; 
        then
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


    ### Delete unsorted bam ###

    if [[ -f ${OUT_DIR}/${SAMPLE}_deduplicated.sorted_by_name.bam ]]; 
    then
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


fi


####        Checkpoint 2: File names         ####

# Some samples have been merged, but retain their IDs from the sequencing run
# We check this by seeing if the file names contain SAMPLE. 
# If they don't then we must change the name

if ! [[ -f "${READS_DIR}/${SAMPLE}_1_trim.fq.gz" && -f "${READS_DIR}/${SAMPLE}_2_trim.fq.gz" ]]; then
    for file in *; do
        # Extract the suffix (everything after the third underscore)
        suffix=$(echo "$file" | sed 's/^[^_]*_[^_]*_[^_]*_//')
        
        # Construct the new file name
        new_filename="${SAMPLE}_${suffix}"
        
        # Rename the file
        mv "$file" "$new_filename"
        
        # Output the new file name
        echo "Renamed: $file -> $new_filename"
    done
else
    #### Checkpoint 3: Alignment ####

    # If a sample's files satisfy both of these conditions, then the files are ready to undergo alignment

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
    echo -e "### Sorting deduplicated bam by name for ${SAMPLE} ###"

    samtools sort -n -@ ${NSLOTS} \
    -o ${OUT_DIR}/${SAMPLE}_deduplicated.sorted_by_name.bam \
    ${OUT_DIR}/${SAMPLE}.deduplicated.bam

    echo -e "### Sorted deduplicated bam by name for ${SAMPLE} ###"

    ### Delete unsorted bam ###
    if [[ -f ${OUT_DIR}/${SAMPLE}_deduplicated.sorted_by_name.bam ]]; then
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
fi # This closes the outer if-else block


