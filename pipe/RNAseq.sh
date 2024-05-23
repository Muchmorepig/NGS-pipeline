source ~/NGS-pipeline/src/fastp_cmd2.sh
trim_files rawdata clean

source ~/NGS-pipeline/src/hisat2_align.sh

hisat2idx="/data5/wmc_data/index/hisat2/Mpolymorpha/tak1_V6/tak1_V6"
hisat_align clean  alignment/map2Tak1V6 32 "normal"
