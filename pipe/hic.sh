#!/bin/bash

source ~/NGS-pipeline/src/fastp_cmd2.sh
source ~/NGS-pipeline/src/bwa_mem2_pairtools_cmd.sh


trim_files rawdata clean


bwa_mem2_index="/data5/wmc_data/index/bwa_mem2/tak2_V6/tak2"
c_size="/data5/wmc_data/reference/Mpolymorpha/tak2_V6/Tak2_discard_short_scaffold.chrom.size"

sample_names=($(get_unique_sample_names "clean" ".R[12].fq.gz"))

bwa_dedup_pairs clean ${bwa_mem2_index} ${c_size} pair $sample_names

pairtools_select_type pair $sample_names

pairs_cool pair cool $c_size tak2_V6 $sample_names
