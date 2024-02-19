#!/bin/bash
# Define default HISAT2 options
hisat2_options_normal=""
hisat2_options_strict="--no-softclip --no-mixed --no-discordant -I 10 -X 700"
# hisat2_options=${hisat2_options_normal}

hisat_align() {
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
        mkdir -p "$odir"
    fi

    if [ ! -d "${odir}/log" ]; then
        echo "The log directory ${odir}/log does not exist. Creating it."
        mkdir -p "${odir}/log"
    fi

    echo >&2 "Starting alignment at $(date)"

    # Select HISAT2 options based on type
    case $type in
        "normal")
            hisat2_options=${hisat2_options_normal}
            ;;
        "strict")
            hisat2_options=${hisat2_options_strict}
            ;;
        *)
            echo "Invalid type: $type"
            return 1
            ;;
    esac

    # Iterate over the fastq files and run HISAT2 alignment for each pair
    local ff=$(ls ${indir}/*fq.gz | sed 's/_R[12]\..*//' | uniq)

    for i in $ff; do
        if [ ! -f "${i}_R1.fq.gz" ] || [ ! -f "${i}_R2.fq.gz" ]; then
            echo "Input files do not exist."
            return 1
        fi
        echo >&2 "Processing $(basename $i)"
        hisat2 -p $th ${hisat2_options} -x ${hisat2idx} --new-summary \
            -1 "${i}_R1.fq.gz" -2 "${i}_R2.fq.gz" \
            2>>"${odir}/log/$(basename "$i").hisat2" | \
            samtools sort -@ $th -O bam -o "${odir}/$(basename "$i").sorted.bam" -

        samtools index -@ 8 "${odir}/$(basename "$i").sorted.bam"
    done

    echo >&2 "Finished alignment at $(date)"

    return 0
}