#!/usr/bin/env bash
# C1 must succeed before C2 starts. Paths are resolved before changing directory.
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: bash $0 TRANSCRIPTOME RAWDATA_ROOT OUTPUT_DIR" >&2
    echo 'RAWDATA_ROOT must contain C1/ and C2/; FASTQ sample prefixes must be C1/C2.' >&2
    exit 2
fi
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
reference=$(cd -- "$1" && pwd -P)
rawdata=$(cd -- "$2" && pwd -P)
for sample in C1 C2; do
    [[ -d $rawdata/$sample ]] || { echo "Missing FASTQ directory: $rawdata/$sample" >&2; exit 2; }
done
mkdir -p -- "$3"
output=$(cd -- "$3" && pwd -P)
cd -- "$output"
# Unique log directory preserves previous attempts, including retries.
mkdir -p logs
log_dir=$(mktemp -d "$output/logs/run-$(date +%Y%m%dT%H%M%S)-XXXXXX")
printf 'Logs: %s\n' "$log_dir"
printf 'sample\tstate\texit_code\ttime\n' > "$log_dir/status.tsv"
for sample in C1 C2; do
    printf '%s\tRUNNING\t-\t%s\n' "$sample" "$(date -Is)" >> "$log_dir/status.tsv"
    if bash "$script_dir/count.sh" "$sample" "$sample" "$reference" "$rawdata/$sample" > "$log_dir/$sample.log" 2>&1; then
        printf '%s\tSUCCEEDED\t0\t%s\n' "$sample" "$(date -Is)" >> "$log_dir/status.tsv"
    else
        status=$?
        printf '%s\tFAILED\t%s\t%s\n' "$sample" "$status" "$(date -Is)" >> "$log_dir/status.tsv"
        echo "$sample failed (exit $status); see $log_dir/$sample.log" >&2
        exit "$status"
    fi
done
echo 'C1 and C2 commands both exited successfully; review output files and QC.'
