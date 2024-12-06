#!/bin/bash

input_dir=${1}
output_dir=${2:-output_bw}
threads=${3:-20}

if [ ! -d "${input_dir}" ]; then
    echo "Input directory ${input_dir} does not exist."
    exit 1
fi

if [ ! -d "${output_dir}" ]; then
    mkdir -p "${output_dir}"
fi

# 声明一个数组来存储所有独特前缀
declare -A prefix_map

# 遍历输入目录中的所有 BAM 文件
for file in "${input_dir}"/*.bam; do
    # 提取前缀，例如 F-0d-1.sorted.bam -> F-0d
    prefix=$(basename "$file" | sed -E 's/(-[0-9]+)?\.bam//')

    # 使用关联数组记录前缀
    prefix_map["$prefix"]=1
done

# 遍历每个唯一前缀进行合并操作和文件转换
for prefix in "${!prefix_map[@]}"; do
    # 创建一个包含该前缀所有文件名的列表
    files=("${input_dir}/${prefix}-"*.bam)

    # 检查是否存在匹配的文件
    if [ ${#files[@]} -gt 0 ]; then
        merged_bam="${input_dir}/${prefix}.merged.bam"
        
        samtools merge "$merged_bam" "${files[@]}" -@ 20
        
        # 列出哪些文件被合并
        echo "Merged files into $merged_bam:"
        for file in "${files[@]}"; do
            echo " - $file"
        done

        # 为合并后的 BAM 文件建立索引
        samtools index "$merged_bam"
        echo "Indexed merged BAM file: $merged_bam"

        # 转换合并后的 BAM 文件为 BW 文件
        bw_file="${output_dir}/${prefix}.bw"
        bamCoverage -b "$merged_bam" \
            --binSize 10 --numberOfProcessors ${threads} \
            -o "$bw_file" --normalizeUsing BPM > /dev/null 2>&1

        echo "Converted $merged_bam to $bw_file"

        rm "$merged_bam" "$merged_bam.bai"
        echo "Deleted merged BAM file and its index: $merged_bam"
    else
        echo "No files matched for prefix ${prefix}"
    fi
done
