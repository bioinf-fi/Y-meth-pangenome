#!/bin/bash

# Define the paths to your input BAM files folder, output directory, and reference FASTA location
OUTPUT_BED_DIR="./Work/Data/BED"
# Assumed location for BAM files based on your previous script
BAM_DIR="./Work/Data/BAM/ont_nucflag_alignments_Pille_hallast/results_ceph"
# Assumed location for per-individual FASTA references
FASTA_DIR="./Work/Data/FASTA/pangenome_sequences_Jax/ceph"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_BED_DIR"

echo "Starting modkit entropy processing..."

# Loop through all BAM files in the input directory
for bam_file in "$BAM_DIR"/*_chrY.bam; do
    
    # Check if the file actually exists and is a regular file
    if [ -f "$bam_file" ]; then
        # Extract the base tag needed for naming (e.g., NA12886)
        tag=$(basename "${bam_file}" _chrY.bam)
        
        # Define the specific reference FASTA path for this sample
        fasta_ref="$FASTA_DIR"/"${tag}_chrY.fa"
        
        # Define the output BED file path
        output_bed="$OUTPUT_BED_DIR"/"${tag}_entropy.bed"

        echo "Processing $bam_file..."
        echo "Using reference: $fasta_ref"
        echo "Outputting to: $output_bed"

        # Run modkit entropy
        # It requires the input BAM, output path, and reference FASTA as arguments.
        # Removed the unsupported --bgzf flag.
        modkit entropy --in-bam "$bam_file" -o "$output_bed" \
            --ref "$fasta_ref" \
            --cpg # Optional: use --cpg flag to focus on CpG sites, or remove for all modified bases

        echo "Finished processing $tag."
    fi
done

echo "All files processed."
