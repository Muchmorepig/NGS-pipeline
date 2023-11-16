#!/usr/bin/bash

G="\033[0;32m"
R="\033[0;31m"
N="\033[0m"
B="\033[1m"

script_path=$(dirname $(readlink -f $0))
source ${script_path}/../src/cmd_check.sh

# command_exists trim_galore
# source ${script_path}/../src/trim_galore_cmd.sh
# galore ./rawdata clean

command_exists bismark
command_exists bowtie2
command_exists samtools

source ${script_path}/../src/bismark_pipe.sh

# bismark_align clean /data5/wmc_data/index/bismark/mp_tak2 align 10
# deduplicate_bismark_paired align dedup

methylation_extraction_aligned_bismark dedup /data5/wmc_data/index/bismark/mp_tak2 methy 10