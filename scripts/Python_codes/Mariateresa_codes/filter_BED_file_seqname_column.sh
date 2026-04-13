#!/bin/bash

# Define the directory containing your BED files
BED_DIR="Work/annotation/seq_classes_Pille_Hallast/2025-11_draft/t2tv2/BED"

echo "Starting filtering process for BED files in $BED_DIR..."
echo "---"

# Loop through all files ending with .bed
for bed_file_path in "$BED_DIR"/*.bed; do
    
    if [[ -f "$bed_file_path" ]]; then
        echo "Filtering $bed_file_path..."
        
        # Use awk to filter lines: only print lines where the first column does NOT start with "chrY_random"
        # '^chrY_random' is a regular expression meaning "starts with chrY_random"
        awk -i inplace '$1 !~ /^chrY_random/' "$bed_file_path"

        # Alternative method using grep (might be simpler):
        # grep -v -e "^chrY_random" "$bed_file_path" > temp.bed && mv temp.bed "$bed_file_path"

    else
        echo "Warning: No BED files found in $BED_DIR"
    fi
done

echo "---"
echo "All applicable BED files filtered."
