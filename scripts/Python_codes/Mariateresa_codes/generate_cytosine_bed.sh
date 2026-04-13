#!/bin/bash

# Configuration:
FASTA_FILE="Work/annotation/chm13v1.0_hg002Xv0.7_grch38Y.fasta"
START_COORD=1           # 1-based start coordinate of the range
END_COORD=46709983	# 1-based end coordinate of the range
OUTPUT_BED="Work/Simulations/fixed_range_methylation.bed"

# 1. Use samtools faidx to extract the specific sequence range
# The region format is "chr:start-end"
REGION="chr21:${START_COORD}-${END_COORD}"
echo "Extracting sequence for $REGION..." >&2

# samtools faidx outputs FASTA format, so we need to linearize it for processing
SEQUENCE=$(samtools faidx $FASTA_FILE $REGION | grep -v "^>" | tr -d "\n")

if [ -z "$SEQUENCE" ]; then
    echo "Error: Could not extract sequence for region $REGION."
    exit 1
fi

# 2. Use awk to find every cytosine in the extracted sequence
echo "$SEQUENCE" | awk -v chrom="chr21" -v start_offset="$START_COORD" '
BEGIN {
    OFS="\t";
    srand(); # Seed random number generator
}
{
    seq=$0;
    for (i=0; i<length(seq); i++) {
        if (substr(seq, i+1, 1) ~ /[Cc]/) {
            # Convert 0-based index (i) to 1-based genomic coordinates:
            genomic_start = start_offset + i - 1; 
            genomic_end = genomic_start + 1;
            
            # Generate random methylation percentage (0-100)
            methylation_percent = int(rand() * 101);
            
            print chrom, genomic_start, genomic_end, methylation_percent;
        }
    }
}' > $OUTPUT_BED

# 3. Sort the final BED file (optional, but good practice)
# Note: The output is already sorted by coordinate, but this ensures standard format
sort -k1,1 -k2,2n $OUTPUT_BED > ${OUTPUT_BED%.*}_sorted.bed

echo "Generated sorted BED file: ${OUTPUT_BED%.*}_sorted.bed"
