#!/usr/bin/env bash
run_bwa_mem2() {
    local ref_index=$1
    local fastq_1=$2
    local fastq_2=$3

    bwa-mem2 mem -t 30 -SP "$ref_index" "$fastq_1" "$fastq_2"
}

run_pairtools() {
    local sizes_file=$1
    local output_file=$2
    local threads=$3

    pairtools parse --nproc-in "$threads" --nproc-out "$threads" --drop-sam --drop-seq -c "$sizes_file" |
        pairtools sort --nproc "$threads" |
        pairtools dedup -p "$threads" --backend cython -o "$output_file"
}

# 一个函数，用于从输入目录获取所有唯一的样本名称
get_unique_sample_names() {
    local indir=$1
    ls "${indir}"/*_R[12].fq.gz | sed 's/_R[12]\..*//' | uniq
}

main_pipeline() {
    local indir=$1
    local ref_index=$2
    local sizes_file=$3
    local output_dir=$4
    local threads=$5
    local log_dir="${output_dir}/log"
    # 获取所有唯一的样本名称
    local sample_names=($(get_unique_sample_names "$indir"))

    mkdir -p "$output_dir"
    mkdir -p "$log_dir"

    # 遍历每个样本名并执行流程
    for base_path in "${sample_names[@]}"; do
        local fastq_1="${base_path}_R1.fq.gz"
        local fastq_2="${base_path}_R2.fq.gz"
        local base_name="$(basename ${base_path})"

        if [[ -f "$fastq_1" && -f "$fastq_2" ]]; then
            echo >&2 "Processing sample: $base_name"
            date >&2

            local output_file="${output_dir}/${base_name}.pairs"

            # 运行bwa-mem2并将输出传递给pairtools进行后续处理
            (
                run_bwa_mem2 "$ref_index" "$fastq_1" "$fastq_2" 2>>${log_dir}/${base_name}.mem.log |
                    run_pairtools "$sizes_file" "$output_file" "$threads"
            ) >${log_dir}/${base_name}.log 2>&1 &
            
            echo "Completed processing of sample: $base_name"

            while (($(jobs -r | wc -l) >= 2)); do wait -n; done

        else
            echo "Skipping sample ${base_name}: One or both FASTQ files are missing."
        fi
    done

    wait
    echo >&2 "All samples have been processed."
    date >&2
}
