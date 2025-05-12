#!/bin/bash
trim_with_fastp() {
    local infile_r1=$1
    local infile_r2=$2
    local base=$3
    local odir=$4

    if [[ ! -f "${infile_r1}" || ! -f "${infile_r2}" ]]; then
        echo "Error: input files not found for sample ${base} ..."
        exit 1
    fi

    fastp --in1 "${infile_r1}" --in2 "${infile_r2}" \
        --out1 "${odir}/${base}_R1.fq.gz" --out2 "${odir}/${base}_R2.fq.gz" \
        --json "${odir}/log/${base}.fastp.json" --html "${odir}/log/${base}.fastp.html" \
        --thread 16 2>"${odir}/log/${base}.logs"
}

export -f trim_with_fastp

trim_files() {
    local indir=$1
    local odir=$2

    if [[ ! -d "${indir}" ]]; then
        echo "Error: input directory ${indir} does not exist."
        exit 1
    fi

    if [ ! -d "$odir" ]; then
        echo "The output directory $odir does not exist. Creating it."
        mkdir -p $odir
    fi

    if [ ! -d "${odir}/log" ]; then
        echo "The log directory ${odir}/log does not exist. Creating it."
        mkdir -p ${odir}/log
    fi

    date >&2
    echo >&2 "[info] Trimming files in ${indir} using 'fastp'"
    echo "       Output directory: $odir"
    
    # cite info
    parallel --citation
    find "${indir}" -name '*_R1.fq.gz' -print0 |
        sed -z 's/_R1\.fq\.gz$//' |
        parallel -0 -j 4 trim_with_fastp \
            "{}_R1.fq.gz" \
            "{}_R2.fq.gz" \
            "{= s:.*/:: =}" \
            "${odir}"

    date >&2
    echo >&2 "[info] Trimming has finished"
}

# 主函数调用
# trim_files "$@"
