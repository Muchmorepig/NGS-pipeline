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
    local bedGraph_dir="${output_dir}/bedGraph"
    local methylDackel_path="/data5/wmc_data/mambaforge/envs/wgbs/bin/MethylDackel"

    # Directory creation
    mkdir -p "${output_dir}" "${log_dir}" "${methylKit_dir}" "${bedGraph_dir}" || {
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
        if ! process_bam "${genome}" "${bam_file}" "${output_prefix}" "${bias_log}" \
            "${log_file}" "${methylDackel_path}" "${threads}" "${methylKit_dir}" "${bedGraph_dir}"; then
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
    local output_prefix="$3"
    local bias_log="$4"
    local log_file="$5"
    local methylDackel_path="$6"
    local threads="$7"
    local methylKit_dir="$8"
    local bedGraph_dir="$9"
    local OT_OB=()
    local bias_line

    # Read the suggested options line differently depending on the shell
    if [ -n "$ZSH_VERSION" ]; then
        bias_line=$(awk '/Suggested inclusion options:/{print $4,$5,$6,$7}' "$bias_log")
        IFS=' ' read -r -A OT_OB <<<"$bias_line"
    else
        bias_line=$(awk '/Suggested inclusion options:/{print $4,$5,$6,$7}' "$bias_log")
        IFS=' ' read -r -a OT_OB <<<"$bias_line"
    fi

    echo "Suggested options: ${OT_OB[@]}"

    local cmd_common=("$methylDackel_path" extract --CHG --CHH ${OT_OB[@]} -@ $threads "$genome" "$bam_file" --opref "$output_prefix")

    local cmds_extras=(
        "--methylKit >> Generating files for methylKit..."
        "--cytosine_report >> Generating cytosine report..."
        "--fraction >> Generating fractional methylation bedGraph..."
        "--logit >> Generating logit(M/(M+U)) methylation bedGraph..."
        "--counts >> Generating methy-base-counts bedGraph..."
    )

    for cmd_extra in "${cmds_extras[@]}"; do
        local cmd_option="${cmd_extra%% *}"
        local description="${cmd_extra#*>> }"
        local cmd=("${cmd_common[@]}" "$cmd_option")
        execute_and_log "$description" "${cmd[@]}" "$log_file" || return 1

        if [[ "$cmd_option" == '--methylKit' ]]; then
            mv "${output_prefix}"*.methylKit "$methylKit_dir" 2>/dev/null || echo "Failed to move .methylKit files or none found."
        else
            mv "${output_prefix}"*.bedGraph "$bedGraph_dir" 2>/dev/null || echo " .bedGraph files Creating"
        fi
    done
}

# Function to perform a command and log it.
execute_and_log() {
    local description="$1"
    # Remove the first and the last parameter
    local -a cmd=("${@:2:$#-2}")
    local log_file="${@: -1:1}"

    echo -e "CMD:\n${cmd[*]}" >>"$log_file"
    echo "$description"
    if ! "${cmd[@]}" >>"$log_file" 2>&1; then
        echo "Error during $description, check log file: $log_file."
        return 1
    fi
}

# This final line adds a check if the script is being sourced or executed directly
# and only calls the function if executed directly with the provided arguments.
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    extract_methylation "$@"
fi
