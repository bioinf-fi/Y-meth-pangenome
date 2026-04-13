#### measure methylation heterogeneity
# Define the directories for your input files
BAM_DIR="Work/Data/BAM/ont_nucflag_alignments_Pille_hallast/results_hprc"
BED_DIR="Work/annotation/AmplGene_coords_Pille_Hallast/BED"
OUTPUT_DIR="Work/Pille_Hallast_methylation/modbamtools/AmplGenes" # Define an output directory

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

echo "Starting batch processing..."
echo "BAM source: $BAM_DIR"
echo "BED source: $BED_DIR"
echo "Output destination: $OUTPUT_DIR"
echo "---"
# Loop through all files ending with .bam in the BAM directory
for bam_file_path in "$BAM_DIR"/*.bam; do
    # Get the base filename with extension (e.g., "sample1.bam")
    filename=$(basename -- "$bam_file_path")
    
    # Get the base name without the extension (e.g., "sample1")
    # 1. Remove the .bam extension: "sample_chrY"
    temp_name="${filename%.bam}"
    # 2. Remove the "_chrY" suffix: "sample"
    # This specifically removes the pattern "_chrY" from the end of the string
    base_name="${temp_name%_chrY}"
    # Construct the full path for the corresponding BED file
    bed_file_path="$BED_DIR/${base_name}_ampl_genes_t2tv.bed"
    
    # Construct the full path for the output file
    output_file_path="$OUTPUT_DIR/${base_name}_ampl_genes.calcHet.txt"

    # Check if the corresponding BED file exists before proceeding
    if [[ -f "$bed_file_path" ]]; then
        echo "Processing $filename..."
        
        # Run the modbamtools calcHet command
        modbamtools calcHet -b "$bed_file_path" -o "$output_file_path" "$bam_file_path"

    else
        echo "Warning: Corresponding BED file not found for $filename at $bed_file_path"
    fi
done

echo "---"
echo "All files processed."