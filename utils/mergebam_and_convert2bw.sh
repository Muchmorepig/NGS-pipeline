#!/bin/bash

# Function for displaying usage information
usage() {
    echo "Usage: $0 -i <input_dir> -o <output_dir> -s <suffix_pattern>"
    echo "  -i: Input directory containing BAM files"
    echo "  -o: Output directory for bigwig files"
    echo "  -s: Suffix pattern (e.g., '.rep1.bam' or '.bam')"
    echo "Example:"
    echo "  $0 -i /path/to/bams -o /path/to/output -s .bam"
    exit 1
}

# Function for logging messages
log_message() {
    echo -e "$1"
    echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >> "$log_file"
}

# Function to process each group
process_group() {
    local group=$1
    local input_dir=$2
    local output_dir=$3
    local suffix=$4
    local log_file=$5
    local temp_dir=$6

    cd "$input_dir"
    
    # Create group-specific log message function
    local_log_message() {
        echo -e "$1"
        echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] [$group] $1" >> "$log_file"
    }

    bam_files=$(find . -maxdepth 1 -name "${group}[-_\.][0-9]*${suffix}" | sed 's/\.\///')
    if [ -z "$bam_files" ]; then
        local_log_message "No files found for $group, skipping..."
        return 1
    fi
    local_log_message "Processing group: $group"
    local_log_message "Found files: \n$bam_files"

    # Create temporary merged BAM file
    temp_merged_bam="${temp_dir}/${group}_merged.bam"
    
    # Merge BAM files
    samtools merge -f "$temp_merged_bam" $bam_files -@ 10 
    if [ $? -ne 0 ]; then
        local_log_message "Error merging files for $group, skipping..."
        return 1
    fi

    # Index merged BAM file
    samtools index "$temp_merged_bam" -@ 10
    if [ $? -ne 0 ]; then
        local_log_message "Error indexing merged BAM for $group, skipping..."
        return 1
    fi

    # Convert to bigwig
    bamCoverage --bam "$temp_merged_bam" \
        --outFileName "${output_dir}/bigwig/${group}.bw" \
        --binSize 10 \
        --normalizeUsing BPM \
        --numberOfProcessors 10 > /dev/null 2>&1
    if [ $? -ne 0 ]; then
        local_log_message "Error converting to bigwig for $group"
    fi

    # Remove temporary merged BAM and its index
    rm -f "$temp_merged_bam" "${temp_merged_bam}.bai"

    local_log_message "Completed processing $group"
    local_log_message "------------------------"
}

# Parse command line arguments
while getopts "i:o:s:" opt; do
    case $opt in
        i) input_dir="$OPTARG";;
        o) output_dir="$OPTARG";;
        s) suffix="$OPTARG";;
        ?) usage;;
    esac
done

# Check if required arguments are provided
if [ -z "$input_dir" ] || [ -z "$output_dir" ] || [ -z "$suffix" ]; then
    usage
fi

# Check if input directory exists
if [ ! -d "$input_dir" ]; then
    echo "Error: Input directory $input_dir does not exist"
    exit 1
fi

# Create output directory and temp directory
mkdir -p "${output_dir}/bigwig"
temp_dir=$(mktemp -d)
trap 'rm -rf "$temp_dir"' EXIT

# Setup logging
timestamp=$(date +"%Y%m%d_%H%M%S")
log_file="${output_dir}/merge_bw_${timestamp}.log"

# Change to input directory
cd "$input_dir"

# Find unique groups based on provided suffix
escaped_suffix=$(echo "$suffix" | sed 's/\./\\\./g')
groups=$(find . -maxdepth 1 -name "*${suffix}" | sed -E "s/\.\/(.*)/\1/" | sed -E "s/([-_\.][0-9]+${escaped_suffix}$)//" | sort | uniq)

# Show grouping details
log_message "=== Starting BAM merge and bigwig conversion ==="
log_message "Using suffix pattern: $suffix"
echo -e "\nGroups and their files:"
echo "----------------------"
for group in $groups; do
    echo -e "\n[$group]:"
    find . -maxdepth 1 -name "${group}[-_\.][0-9]*${suffix}" | sed 's/\.\//  /'
done
echo -e "\nTotal groups: $(echo "$groups" | wc -l)"
echo -n "Continue with these groups? [y/n]: "
read answer
if [ "$answer" != "y" ]; then
    log_message "User aborted processing"
    echo "Aborting..."
    exit 1
fi

# Export functions for parallel
export -f process_group
export -f log_message

sleep 0.5
echo "Processing..."
# Process groups in parallel
echo "$groups" | parallel -j 4 process_group {} "$input_dir" "$output_dir" "$suffix" "$log_file" "$temp_dir"

log_message "=== Processing completed ==="
