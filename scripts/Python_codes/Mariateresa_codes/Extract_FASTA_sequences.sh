#!/bin/bash

INPUT_DIR="Work/Data/BED/P_coordinates"
OUTPUT_DIR="Work/Data/FASTA/P_FASTA"
mkdir -p "$OUTPUT_DIR"

for file in "$INPUT_DIR"/*_P[1-8]_ranges.txt; do
    [ -e "$file" ] || continue
    
    filename=$(basename "$file")
    # Extract Sample (everything before the first _)
    sample=$(echo "$filename" | cut -d'_' -f1)
    p_class=$(echo "$filename" | grep -o "P[1-8]")
    
    input_fasta="Work/Data/FASTA/pangenome_sequences_Jax/ceph/${sample}_chrY.fa"
    output_fasta="${OUTPUT_DIR}/${sample}_${p_class}.fasta"
    
    if [[ -f "$input_fasta" ]]; then
        echo "Processing $sample ($p_class)..."
        
        # We construct the header format "Sample_chrY:Start-End" 
        # NR>1 skips the header row of your .txt file
        awk -v smp="$sample" 'NR>1 {print smp"_chrY:"$2"-"$3}' "$file" | \
        samtools faidx "$input_fasta" -r - > "$output_fasta"
    else
        echo "File missing: $input_fasta"
    fi
done
