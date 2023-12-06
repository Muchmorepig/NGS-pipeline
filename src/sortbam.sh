#!/usr/bin/bash
# samtools sort
sort_bam(){
    local input=$1
    local output=$2
    local threads=$3
    
    if [ ! -d ${output} ]; then
        mkdir -p ${output}
    fi

    for i in ${input}/*.bam; do
        local tbam=${output}/$(basename $i .bam).tmp.bam
        local obam=${output}/$(basename $i .bam).sorted.bam
        local cmd1="samtools sort -@ ${threads} -o ${tbam} ${i}"
        local cmd2="samtools view -@ ${threads} -bF 2828 -q 20 ${tbam} -o ${obam}"
        local cmd3="samtools index ${obam}"
        eval ${cmd1}
        eval ${cmd2}
        eval ${cmd3}
        rm "$tbam"
    done

}