#!/bin/bash
bt2_options_normal="--phred33"
bt2_options_strict="--end-to-end --very-sensitive --no-mixed --no-discordan -I 10 -X 700"
bt2_options_frag_120="--dovetail"
bt2_options=${bt2_options_normal}

bowtie_align() {
    local indir=$1
    local odir=$2
    local th=$3
    local type=$4

    if [ ! -d "$indir" ]; then
        echo "The input directory $indir does not exist."
        return 1
    fi

    if [ ! -d "$odir" ]; then
        echo "The output directory $odir does not exist. Creating it."
        mkdir -p $odir
    fi

    if [ ! -d "${odir}/log" ]; then
        echo "The log directory ${odir}/log does not exist. Creating it."
        mkdir -p ${odir}/log
    fi

    echo "Starting alignment at $(date)" >>${odir}/log/bowtie_align.log

    case $type in
    "normal")
        bt2_options=${bt2_options_normal}
        ;;
    "strict")
        bt2_options=${bt2_options_strict}
        ;;
    "frag_120")
        bt2_options=${bt2_options_frag_120}
        ;;
    *)
        echo "Invalid type: $type" | tee -a ${odir}/log/bowtie_align.log
        return 1
        ;;
    esac

    local ff=$(ls ${indir}/*fq.gz | sed 's/_R[12]\..*//' | uniq)

    for i in $ff; do
        if [ ! -f "${i}_R1.fq.gz" -o ! -f "${i}_R2.fq.gz" ]; then
            echo "Input files do not exist." | tee -a ${odir}/log/bowtie_align.log
            return 1
        fi

        bowtie2 -p $th $bt2_options -x ${bt2idx} \
            -1 ${i}_R1.fq.gz -2 ${i}_R2.fq.gz \
            2>>${odir}/log/${i##*/}.bowtie2 |
            samtools sort -@ 20 -O bam -o ${odir}/${i##*/}.sorted.bam -
    done

    echo "Finished alignment at $(date)" >>${odir}/log/bowtie_align.log

    return 0
}
