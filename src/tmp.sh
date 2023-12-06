#!/usr/bin/bash
set -o pipefail

extract_methylation() {
    # Usage and argument check
    if [ "$#" -lt 3 ]; then
        echo "Usage: extract_methylation <input_dir> <genome> <output_dir> [<threads>=10]"
        return 1
    fi

    # Variable assignment
    local input_dir="$1"
    local genome="$2"
    local output_dir="$3"
    local threads="${4:-10}"
    local log_dir="${output_dir}/logs"
    local methylKit_dir="${output_dir}/methylKit"
    local methylDackel_path="/data5/wmc_data/mamba/envs/wgbs/bin/MethylDackel"

    # Directory creation
    mkdir -p "${output_dir}" "${log_dir}" "${methylKit_dir}" || {
        echo "Failed to create directory $dir"
        return 1
    }

    # Confirm MethylDackel is accessible
    if ! command -v "${methylDackel_path}" &>/dev/null; then
        echo "MethylDackel cannot be found or is not executable."
        return 1
    fi

    # Processing
    echo ">>>>>> Begin to extract methylation using MethylDackel"
    date >&2
    echo

    # Retrieve an array of BAM files
    local bam_files=("${input_dir}"/*.bam)
    if [ "${#bam_files[@]}" -eq 0 ]; then
        echo "No BAM files found in ${input_dir}."
        return 1
    fi

    # Process each BAM file
    for bam_file in "${bam_files[@]}"; do
        local base
        base=$(basename "${bam_file}" .sorted.bam)
        local output_prefix="${output_dir}/${base}"
        local log_file="${log_dir}/${base}_methylDackel.log"
        local bias_log="${output_prefix}.bias"

        echo -e "Mbias generation for Sample: ${base}"
        "${methylDackel_path}" mbias --CHG --CHH -@ "${threads}" "${genome}" "${bam_file}" "${output_prefix}" >"${bias_log}" 2>&1

        echo -e "Extracting methylation for Sample: ${base}"
        if ! process_bam "${genome}" "${bam_file}" "${base}" "${output_prefix}" "${bias_log}" \
            "${log_file}" "${methylDackel_path}" "${threads}" "${methylKit_dir}"; then
            echo "Error occurred with ${base}, please check the log file."
            continue
        fi

        echo -e "${base} processing complete.\n"
    done

    echo "<<<<<< The results are in: ${output_dir}"
    date >&2
}

process_bam() {
    local genome="$1"
    local bam_file="$2"
    local base="$3"
    local output_prefix="$4"
    local bias_log="$5"
    local log_file="$6"
    local methylDackel_path="$7"
    local threads="$8"
    local methylKit_dir="$9"
    local ot ob
    {
        read -r ot
        read -r ob
    } < <(awk '/Suggested inclusion options:/ {print $4 " " $5; print $6 " " $7}' "$bias_log")

    local cmd_common=("$methylDackel_path" extract --CHG --CHH $ot $ob -@ "$threads" "$genome" "$bam_file" --opref "$output_prefix")
    # Extraction to MethylKit format
    local cmd_extra=("${cmd_common[@]}" --methylKit)
    echo -e "CMD:\n${cmd_extra[*]}" >>"$log_file"
    echo 'Generating files for methylKit...'
    if ! "${cmd_extra[@]}" >>"$log_file" 2>&1; then
        echo "Error during methylKit file generation for $base, check log file: $log_file."
        return 1
    fi

    # Move generated files to methylKit directory
    mv "${output_prefix}"*.methylKit "$methylKit_dir"

    # Extract cytosine report if the previous step succeeded
    local cmd_cyto=("${cmd_common[@]}" --cytosine_report)
    echo -e "CMD:\n${cmd_cyto[*]}" >>"$log_file"
    echo 'Generating cytosine report...'
    if ! "${cmd_cyto[@]}" >>"$log_file" 2>&1; then
        echo "Error during cytosine report generation for $base, check log file: $log_file."
        return 1
    fi
}

# This final line adds a check if the script is being sourced or executed directly
# and only calls the function if executed directly with the provided arguments.
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    extract_methylation "$@"
fi
