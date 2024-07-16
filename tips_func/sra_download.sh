#!/bin/bash

max_jobs=2  # 最大并发任务数量
sra_list="sra.list"
output_dir="output"
 
while getopts 'j:l:o:' flag; do
    case "${flag}" in
        j) max_jobs="${OPTARG}" ;;
        l) sra_list="${OPTARG}" ;;
        o) output_dir="${OPTARG}" ;;
        *) exit 1 ;;
    esac
done

job_count=0 

function run_fastq_dump {
    local sra=$1
    local attempt=0
    local max_attempts=4  # 最多重试次数

    while (( attempt < max_attempts )); do
        let attempt+=1
        fasterq-dump \
                --mem 500MB \
                --temp "${output_dir}/temp" \
                --outdir "${output_dir}" \
                --split-files "$sra" > "${output_dir}/${sra}.log" 2>&1

        if [ $? -eq 0 ]; then
            echo "Successfully processed $sra on attempt $attempt."
            return 0
        else
            echo "Failed to process $sra on attempt $attempt. Retrying..."
            sleep 4  # 等待5秒后重试
        fi
    done

    echo "Failed to process $sra after $max_attempts attempts."
    return 1
}

function wait_for_jobs {
    while (( job_count >= max_jobs )); do
        wait -n  # 等待任意一个后台任务完成
        ((job_count--))
    done
}

mkdir -p "${output_dir}/temp"

# 读取 SRA 列表文件并处理
while read sra; do
    wait_for_jobs
    run_fastq_dump "$sra" &
    ((job_count++))
done < "$sra_list"

wait

rm -rf "${output_dir}/temp"
