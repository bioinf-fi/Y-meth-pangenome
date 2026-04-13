#!/bin/bash

# Define the paths to your output directory
OUTPUT_DIR="./Work/Data/TSV"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all BAM files in the input directory
for bam_file in Work/Data/BAM/ont_nucflag_alignments_Pille_hallast/results_ceph/*_chrY.bam; do
    
    tag=$(basename "${bam_file}" _chrY.bam)
    fasta_ref="Work/Data/FASTA/pangenome_sequences_Jax/ceph/${tag}_chrY.fa"
    
    if [ -f "$bam_file" ]; then
        filename=$(basename -- "$bam_file") 
        output_full="$OUTPUT_DIR"/"${tag}_full_calls.tsv.gz"
        output_calls="$OUTPUT_DIR"/"${tag}_thresholded_calls.tsv.gz"

        echo "Processing $filename..."
        echo "Using reference: $fasta_ref"

        # --- FIX: Pass output path as argument. Remove | bgzip > redirect ---
        modkit extract full \
            --bgzf \
            --reference "$fasta_ref" \
            "$bam_file" \
            "$output_full" # This is the <OUT_PATH> argument

        # --- FIX: Pass output path as argument. Remove | bgzip > redirect ---
        modkit extract calls \
            --bgzf \
            --reference "$fasta_ref" \
            "$bam_file" \
            "$output_calls" # This is the <OUT_PATH> argument

        echo "Finished processing $filename."
    fi
done

echo "All files processed."
