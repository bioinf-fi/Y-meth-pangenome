#!/bin/bash

# Define the directories for your input files
BAM_DIR="Work/Data/BAM/ont_nucflag_alignments_Pille_hallast/results_hgsvc" 
BED_DIR="Work/annotation/seq_classes_Pille_Hallast/2025-11_draft/t2tv2/BED"

echo "Starting first column modification for BED files in $BED_DIR..."
echo "---"

# Loop through all files ending with .bam in the BAM directory
for bam_file_path in "$BAM_DIR"/*.bam; do
    
    # Get the full filename with extension (e.g., "sample_chrY.bam")
    filename=$(basename -- "$bam_file_path")
    
    # 1. Remove the .bam extension: "sample_chrY"
    temp_name="${filename%.bam}"
    
    # 2. Remove the "_chrY" suffix: "sample"
    # This is your simple basename: "sample"
    base_name="${temp_name%_chrY}"

    # Define the *target* chromosome name for the BED file: "sample_chrY"
    TARGET_CHR_NAME="${base_name}_chrY"

    # Construct the full path for the corresponding BED file (e.g., "/path/to/your/bed_folder/sample.bed")
    bed_file_path="$BED_DIR/${base_name}.t2tv2.chrY-regions.err-strict.bed"

    # Check if the corresponding BED file exists before proceeding
    if [[ -f "$bed_file_path" ]]; then
        echo "Modifying $bed_file_path: changing first column from 'chrY' to '$TARGET_CHR_NAME'..."
        
        # Use awk to change the first column ($1) from "chrY" to the value of the 'TARGET_CHR_NAME' variable
        awk -i inplace -v OFS='\t' -v new_chr_name="$TARGET_CHR_NAME" '$1=="chrY" {$1=new_chr_name} 1' "$bed_file_path"

    else
        echo "Warning: Corresponding BED file not found for $filename at $bed_file_path"
    fi
done

echo "---"
echo "All applicable BED files processed."
