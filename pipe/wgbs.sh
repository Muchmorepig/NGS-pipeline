#!/bin/bash

# 设置颜色
G="\033[0;32m"
R="\033[0;31m"
N="\033[0m"
B="\033[1m"

# 参数解析
index=$1
genome=$2

# 获取脚本的路径
script_path=$(dirname $(readlink -f $0))

# 引入所需的函数
source ${script_path}/../src/cmd_check.sh
#source ${script_path}/../src/fastp_cmd.sh
source ${script_path}/../src/trim_galore_cmd.sh
source ${script_path}/../src/bismark_pipe.sh
source ${script_path}/../src/sortbam.sh
source ${script_path}/../src/MethyDackel_pipe.sh

# 函数：确认工具是否存在
check_deps() {
    command_exists trim_galore
    command_exists bismark
    command_exists bowtie2
    command_exists samtools
}

# 函数：质量和接头的修剪
trim_adapters() {
    echo -e "${G}Trimming adapters...${N}"
    galore ./rawdata clean
    #trim_files ./rawdata clean 16
}

# 函数：Bismark比对
bismark_alignment() {
    echo -e "${G}Performing Bismark alignment...${N}"
    bismark_align clean ${index} align 10 true
}

# 函数：去重复
deduplicate_alignment() {
    echo -e "${G}Deduplicating alignments...${N}"
    deduplicate_bismark_paired align dedup
}

# 函数：排序bam文件
sort_bam_files() {
    echo -e "${G}Sorting BAM files...${N}"
    sort_bam dedup dedup/sorted 6
}

# 函数：提取甲基化信息
extract_methylation_info() {
    echo -e "${G}Extracting methylation information...${N}"
    extract_methylation dedup/sorted ${genome} MethyDackel 16
}

# 主流程函数
main() {
    check_deps               # 确认所有依赖存在
    trim_adapters            # 质量和接头修剪
    bismark_alignment        # Bismark比对
    deduplicate_alignment    # 去重复
    sort_bam_files           # 排序BAM文件
    extract_methylation_info # 提取甲基化信息
}

# 调用主流程函数
main "$@"

# optional
# methylation_extraction_aligned_bismark dedup ${index} methy_bismark 10
