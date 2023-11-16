#!/bin/bash
bismark_align() {
    local indir=$1
    local ref=$2
    local odir=$3
    local th=$4
    local temp_dir='./temp'

    if [ ! -d "$indir" ]; then
        echo "The input directory $indir does not exist."
        exit 1
    fi

    mkdir -p $odir $odir/log $temp_dir

    echo "Starting alignment at $(date)"
    local file_arr=($(ls ${indir}/*gz | sed 's/_R[12]_val_[12]\..*//' | uniq))

    for i in "${file_arr[@]}"; do
        local fq_1="${i}_R1_val_1.fq.gz"
        local fq_2="${i}_R2_val_2.fq.gz"
        local base=$(basename $i)
        local output_file=${odir}/${base}

        if [ ! -f "$fq_1" ] || [ ! -f "$fq_2" ]; then
            echo "Input files do not exist." | tee -a ${odir}/log/bismark_align.log
            exit 2
        fi

        bismark --bowtie2 $ref \
            --parallel $th -p 4 --temp_dir $temp_dir \
            -o $odir \
            -1 $fq_1 \
            -2 $fq_2 \
            --bam 2>$odir/log/${i##*/}.log

        local tmp_bam=${output_file}*_bismark_bt2_pe.bam
        mv $tmp_bam ${output_file}.bam
    done

    echo "Finished alignment at $(date)"
}

deduplicate_bismark_paired() {
    local indir=$1
    local odir=$2

    if [ ! -d "$indir" ]; then
        echo -e "${R}The input directory $indir does not exist.${N}"
        return 1
    fi

    if [ ! -d "$odir" ]; then
        echo -e "${B}The output directory $odir does not exist. Creating it.${N}\n"
        mkdir -p $odir/reports
    fi

    echo "Starting deduplication at $(date)"
    local ff=$(ls ${indir}/*.bam)

    for i in $ff; do
        base=$(basename ${i})
        echo -e "  ${G}BAM: ${base%.*}${N}"

        deduplicate_bismark --bam --paired $i -o ${odir} >/dev/null 2>&1 &
    done

    wait
    mv ${odir}/*txt ${odir}/reports
    echo "Finished deduplication at $(date)"

    return 0
}

methylation_extraction_aligned_bismark() {
    local indir=$1
    local bismark_genome_folder=$2
    local odir=$3
    local th=$4

    if [ ! -d "$indir" ]; then
        echo "The input directory $indir does not exist."
        return 1
    fi

    if [ ! -d "$odir" ]; then
        echo "The output directory $odir does not exist. Creating it."
        mkdir -p $odir
    fi

    echo "Starting methylation extraction at $(date)"
    local ff=$(ls ${indir}/*.deduplicated.bam)

    for i in $ff; do
        base=$(basename ${i})
        echo -e "  ${G}BAM: ${base}${N}"
        bismark_methylation_extractor --paired-end --ignore_r2 2 --CX --comprehensive \
            --gzip --bedGraph --multicore $th --buffer_size 20G \
            --cytosine_report \
            --genome_folder $bismark_genome_folder \
            $i -o $odir >/dev/null 2>&1
    done

    echo "Finished methylation extraction at $(date)"

    return 0
}
