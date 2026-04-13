#!/bin/bash

mkdir -p P2_fasta_outputs
echo "Individual,Identity_Percentage,Alignment_Length" > p2_ceph_identity_results.csv

# Loop through the CSV files created in R
for csv in Work/annotation/seq_classes_Pille_Hallast/P2_occurrences/*_P2.csv; do
    # Extract Individual Name from filename
    indiv=$(basename "$csv" _P2.csv)
    fasta="Work/Data/FASTA/pangenome_sequences_Jax/ceph/${indiv}_chrY.fa"
    
    # Check if Fasta exists
    if [[ ! -f "$fasta" ]]; then continue; fi

    # Get coordinates for the 1st occurrence (Row 2 because Row 1 is header)
    coords1=$(awk -F',' 'NR==2 {print $1":"$2"-"$3}' "$csv" | tr -d '"')
    # Get coordinates for the 2nd occurrence (Row 3)
    coords2=$(awk -F',' 'NR==3 {print $1":"$2"-"$3}' "$csv" | tr -d '"')

    # If coords2 is empty, skip (individual only has one P2)
    if [[ -z "$coords2" ]]; then continue; fi

    # 1. Extract Sequences
    samtools faidx "$fasta" "$coords1" > "Work/Data/FASTA/pangenome_sequences_Jax/P2_fasta_outputs/${indiv}_occ1.fa"
    samtools faidx "$fasta" "$coords2" > "Work/Data/FASTA/pangenome_sequences_Jax/P2_fasta_outputs/${indiv}_occ2.fa"

    # 2. Run BLAST (Local Alignment)
    # Using -outfmt 6 to get a simple tab-separated line
    identity_data=$(blastn -query "Work/Data/FASTA/pangenome_sequences_Jax/P2_fasta_outputs/${indiv}_occ2.fa" \
                           -subject "Work/Data/FASTA/pangenome_sequences_Jax/P2_fasta_outputs/${indiv}_occ1.fa" \
                           -outfmt "6 pident length" | head -n 1)

    # 3. Save result
    if [[ ! -z "$identity_data" ]]; then
        # Replace tabs with commas for CSV format
        formatted_data=$(echo "$identity_data" | tr '\t' ',')
        echo "${indiv},${formatted_data}" >> p2_ceph_identity_results.csv
    fi
done
