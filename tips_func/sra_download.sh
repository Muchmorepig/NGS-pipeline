#!/bin/bash

max_jobs=2  # 最大并发任务数量
sra_list="sra.list"
output_dir="output"

while getopts 'j:l:o:' flag; do
    case "${flag}" in
        j) max_jobs="${OPTARG}" ;;
        l) sra_list="${OPTARG}" ;;
        o) output_dir="${OPTARG}" ;;
        *) 
            echo "Usage: $0 [-j max_jobs] [-l sra_list] [-o output_dir]"
            exit 1 
            ;;
    esac
done

temp_dir="${output_dir}/temp_$$"
mkdir -p "${output_dir}" "${temp_dir}" || { echo "Failed to create needed directories"; exit 1; }

job_count=0
declare -a job_pids

run_fastq_dump() {
    local sra="$1"
    local attempts=0
    local max_attempts=4  # 最多重试次数
    local log_file="${output_dir}/${sra}.log"

    while (( attempts < max_attempts )); do
        ((attempts++))

        if prefetch --progress --max-size 300G -q --output-directory "${output_dir}" "${sra}" && \
           fasterq-dump --mem 500MB --temp "${temp_dir}" --outdir "${output_dir}" --split-files "${output_dir}/${sra}/${sra}.sra" > "${log_file}" 2>&1; then
            rm -rf "${output_dir}/${sra}"
            echo "Successfully processed ${sra} on attempt ${attempts}." | tee -a "${log_file}"
            return 0
        else
            echo "Failed to process ${sra} on attempt ${attempts}. Retrying in 4 seconds..." | tee -a "${log_file}"
            sleep 4  # 等待4秒后重试
        fi
    done

    echo "Failed to process ${sra} after ${max_attempts} attempts." | tee -a "${log_file}"
    return 1
}

wait_for_jobs() {
    while (( job_count >= max_jobs )); do
        for pid in "${job_pids[@]}"; do
            if ! kill -0 "$pid" 2>/dev/null; then
                wait "$pid"
                job_count=$((job_count - 1))
                job_pids=("${job_pids[@]/$pid}")
                break
            fi
        done
        sleep 1
    done
}

while IFS= read -r sra; do
    [[ -z "${sra}" || "${sra}" =~ ^# ]] && continue

    wait_for_jobs
    run_fastq_dump "${sra}" &
    job_pids+=("$!")
    ((job_count++))
done < "${sra_list}"

for pid in "${job_pids[@]}"; do
    wait "$pid"
done

rm -rf "${temp_dir}"
echo "All tasks are complete."
