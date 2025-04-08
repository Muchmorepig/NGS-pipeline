#!/bin/bash

# Usage: ./bamtobw.sh <input_dir> <output_dir> [threads]
# Example: ./bamtobw.sh ~/wkdir/LiuK/data/merge_bam/chloroplast/ ./bigwig 20

# Check if required arguments are provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir> [threads]"
    echo "Example: $0 ~/wkdir/LiuK/data/merge_bam/chloroplast/ ./bigwig 20"
    exit 1
fi

input="$1"
output="$2"
threads=${3:-20}  # Default to 20 threads if not specified

# Remove trailing slash from input path if present
input=${input%/}

# Make the output directory if it doesn't exist
if [ ! -d "${output}" ]; then
    mkdir -p "${output}"
    echo "Created output directory: ${output}"
fi

echo ""
echo "Converting BAM files to BigWig format..."
echo "Input directory: ${input}"
echo "Output directory: ${output}"
echo "Using ${threads} threads"
echo "Normalization method: BPM"
echo ""
# Check if BAM files exist in the input directory
if ! ls "${input}"/*.bam 1> /dev/null 2>&1; then
    echo "No BAM files found in ${input}"
    exit 1
fi

# Process all BAM files in the input directory
for id in "${input}"/*.bam; do
    # Get base filename
    base=$(basename "$id")
    b=${base%%.bam}
    
    echo "Processing: ${base} -> ${b}.bw"
    
    # Convert BAM to BigWig (suppress logs)
    bamCoverage -b "$id" \
        --binSize 10 --numberOfProcessors ${threads} \
        -o "${output}/${b}.bw" --normalizeUsing BPM > /dev/null 2>&1
done

echo "Conversion complete! BigWig files are saved in ${output}"
