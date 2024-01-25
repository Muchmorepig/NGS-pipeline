#!/bin/bash
run_bwa_mem2() {
    local ref_index=$1
    local fastq_1=$2
    local fastq_2=$3

    if [[ $# -ne 3 ]]; then
        echo "Usage: run_bwa_mem2 <ref_index> <fastq_r1> <fastq_r2>" >&2
        return 1
    fi

    bwa-mem2 mem -t 30 -SP "$ref_index" "$fastq_1" "$fastq_2"
}

parse_dedup_pairtools() {
    local sizes_file=$1
    local output_file=$2
    local threads=$3

    pairtools parse --nproc-in "$threads" --nproc-out "$threads" --drop-sam --drop-seq -c "$sizes_file" |
        pairtools sort --nproc "$threads" |
        pairtools dedup -p "$threads" --backend cython -o "$output_file"
}

bwa_dedup_pairs() {
    local indir=$1
    local ref_index=$2
    local sizes_file=$3
    local output_dir=$4
    local sample_names=$5
    local log_dir="${output_dir}/log"

    mkdir -p "$output_dir"
    mkdir -p "$log_dir"

    for base_name in "${sample_names[@]}"; do
        local fastq_1="${indir}/${base_name}_R1.fq.gz"
        local fastq_2="${indir}/${base_name}_R2.fq.gz"

        if [[ -f "$fastq_1" && -f "$fastq_2" ]]; then
            echo >&2 "Processing sample: $base_name"
            date >&2

            local output_file="${output_dir}/${base_name}.pairs"

            # 运行bwa-mem2并将输出传递给pairtools进行后续处理
            (
                run_bwa_mem2 "$ref_index" "$fastq_1" "$fastq_2" 2>>${log_dir}/${base_name}.mem.log |
                    parse_dedup_pairtools "$sizes_file" "$output_file" 16
            ) 2>&1 &

            while (($(jobs -r | wc -l) >= 2)); do wait -n; done
        else
            echo "Skipping sample ${base_name}: One or both FASTQ files are missing."
        fi
    done

    wait
    echo >&2 "All samples have been processed."
    date >&2
}

pairtools_select_type() {
    local indir=$(realpath "$1")
    local sample_names=$2
    echo >&2 "Selecting pairs from $indir"
    date >&2

    for sample in "${sample_names[@]}"; do
        local input_file="${indir}/${sample}.pairs"
        local output_file="${indir}/${sample}.selected.pairs"
        echo >&2 "Processing sample: $sample"
        pairtools select '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' \
            --nproc-in 12 --nproc-out 12 \
            -o "$output_file" "$input_file" && pigz -p 16 "$output_file" && pigz -p 16 "$input_file" &

        while (($(jobs -r | wc -l) >= 2)); do wait -n; done

    done

    wait
    date >&2
}

zoom() {
    local input_file=$1
    local output_file=$2

    cooler zoomify --nproc 12 \
        --balance \
        --resolutions 10000000,5000000,2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000 \
        --out "$output_file" \
        "$input_file"
}

pairs_cool() {
    local input_dir=$1
    local output_dir=$2
    local chromsizes_file=$3
    local assembly=$4
    local sample_names=$5
    local min_res=1000
    
    if [[ ! -d "$output_dir" ]]; then
        mkdir -p "$output_dir"
    fi
    date >&2

    for sample in "${sample_names[@]}"; do

        local input_file="${input_dir}/${sample}.selected.pairs.gz"
        local oo="${output_dir}/${sample}_1kb.cool"
        local ooo="${output_dir}/${sample}.mcool"
        echo >&2 "Processing file: $input_file"

        # 使用pigz解压缩输入文件，并通过管道传递给cooler
        cooler cload pairs \
            -c1 2 -p1 3 -c2 4 -p2 5 \
            --assembly "$assembly" \
            "${chromsizes_file}:${min_res}" <(pigz -dc -p 16 "$input_file") \
            "$oo" &>>"${output_dir}/${sample}.log" &&
            zoom "$oo" "$ooo" &>>"${output_dir}/${sample}.log" &

        while (($(jobs -r | wc -l) >= 2)); do wait -n; done
    done

    wait
    date >&2

}

#
# printf '%s\n' "${sample_names[@]}" | xargs -P 4 -I {} bash -c '
# mamba activate HiCpipe
# fastq_1=${indir}/{}_R1.fq.gz
# fastq_2=${indir}/{}_R2.fq.gz
# if [[ -f "$fastq_1" && -f "$fastq_2" ]]; then
#     run_bwa_mem2 "$ref_index" "$fastq_1" "$fastq_2" 2>>${log_dir}/{}.mem.log | \
#         parse_dedup_pairtools "$sizes_file" "${output_dir}/{}.pairs" "$threads"
# else
#     echo "Skipping sample {}: One or both FASTQ files are missing."
# fi
# '
