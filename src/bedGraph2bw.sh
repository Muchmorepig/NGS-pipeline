#!/bin/bash
set -euo pipefail

input=$1
output=$2
genome=$3

# Paths to utilities
bedSort="/data5/wmc_data/utils/ucsc-tools/bedSort"
bedGraphToBigWig="/data5/wmc_data/utils/ucsc-tools/bedGraphToBigWig"

# Check for necessary commands
for cmd in "$bedSort" "$bedGraphToBigWig" "parallel" "seqkit" "tail"; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "Error: $cmd is not installed or not in PATH." >&2
        exit 1
    fi
done

# Check if the input/output directories exist and are not the same
if [[ ! -d "$input" ]]; then
    echo "Error: Input directory does not exist." >&2
    exit 1
fi

if [[ "$input" == "$output" ]]; then
    echo "Error: Input and output directories must be different." >&2
    exit 1
fi

mkdir -p "$output"

# Generate the `chrom.size` file
seqkit fx2tab -nl "$genome" | awk '{print $1"\t"$2}' >"${output}/chrom.size"

convert_to_bigwig() {
    local bedGraphFile="$1"
    local bedGraphFileName="$(basename "$bedGraphFile")"
    local sortedBedGraphFile="${output}/tmp/${bedGraphFileName}.sorted"
    local bigWigFile="${output}/${bedGraphFileName%.bedGraph}.bw"

    echo "Converting $bedGraphFileName to BigWig..."

    # Create a sorted copy by skipping the first line of the `.bedGraph` file
    tail -n +2 "$bedGraphFile" | "$bedSort" stdin "${sortedBedGraphFile}"
    "$bedGraphToBigWig" "${sortedBedGraphFile}" "${output}/chrom.size" "$bigWigFile"
}

# Export the required variables and function
export bedSort bedGraphToBigWig output
export -f convert_to_bigwig

# Create temporary directory for intermediate sorted files
mkdir -p "${output}/tmp"

# Use GNU Parallel to convert `.bedGraph` files to `.bw` in parallel
find "$input" -name "*.bedGraph" -print0 | parallel --progress --will-cite -0 -j16 convert_to_bigwig "{}"

# Clean up the sorted temporary `.bedGraph` files
rm -rf "${output}/tmp"

echo "Conversion to .bw files is complete."
