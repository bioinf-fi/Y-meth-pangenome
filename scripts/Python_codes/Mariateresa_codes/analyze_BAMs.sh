#!/bin/bash

# Define the directory containing your BAM files
BAM_DIR="./Work/Data/BAM/ont_nucflag_alignments_Pille_hallast/results_ceph"


# Set a minimum mapping quality threshold (e.g., 20 for basic filtering)
MIN_MAPQ=20

# Output file for summary results
SUMMARY_FILE="Work/Data/BAM/ont_nucflag_alignments_Pille_hallast/results_ceph/alignment_summary_ceph.txt"

# Print header to the summary file
echo -e "File_Name\tMean_Depth\tMean_MAPQ\tMean_Read_Length" > $SUMMARY_FILE

# Loop through all BAM files in the directory
for bam_file in "$BAM_DIR"/*.bam; do
    if [ -f "$bam_file" ]; then
        echo "Processing $bam_file..."

        # 1. Get coverage and MAPQ statistics
        # samtools coverage provides mean depth and mean MAPQ for the whole genome/each chromosome
        COV_STATS=$(samtools coverage -q $MIN_MAPQ "$bam_file" | awk -v OFS='\t' 'NR>1 {sum_depth+=$4*$7; sum_mapq+=$4*$9; sum_bases+=$4} END {print sum_depth/sum_bases, sum_mapq/sum_bases}')
        MEAN_DEPTH=$(echo $COV_STATS | cut -f 1 -d ' ')
        MEAN_MAPQ=$(echo $COV_STATS | cut -f 2 -d ' ')

        # 2. Get read length statistics
        # Extract read lengths and calculate the mean using awk
        MEAN_READ_LEN=$(samtools view -F 4 "$bam_file" | awk '{sum += length($10)} END {if(NR>0) print sum/NR; else print 0}')

        # Extract filename for the report
        FILENAME=$(basename "$bam_file")

        # Append results to the summary file
        echo -e "$FILENAME\t$MEAN_DEPTH\t$MEAN_MAPQ\t$MEAN_READ_LEN" >> $SUMMARY_FILE
    fi
done

echo "Done. Results saved to $SUMMARY_FILE"
