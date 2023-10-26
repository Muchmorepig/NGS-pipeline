#!/usr/bin/bash
bowtie_align() {
    local indir=$1
    local odir=$2
    local th=$3

    # 判断是否存在工具
    if ! command -v bowtie2 &>/dev/null || ! command -v samtools &>/dev/null; then
        echo -e "${R}Need tools both 'bowtie2' and 'samtools' ${N}"
        echo -e "${R}Try to Run Comand: \n  module load bowtie/2.3.4.3 && module load samtools/1.10 ${N}"
        echo -e "${R}Then, Rerun Command${N} ${0} without ${B}'-t'${N} to skip trim"
        exit 1
    fi

    if [ -d "${indir}" ] && ls ${indir}/*fq.gz &>/dev/null; then
        date >&2
        ff=$(ls ${indir}/*fq.gz | sed 's/_R[12]\..*//' | uniq)

        echo >&2 -e "[info] Index: ${B}${bt2idx}${N}\n"

        # align
        if [ ! -d "${odir}" ]; then
            echo -e "  ${G}Created directory: ${odir}${N}"
            mkdir -p ${odir}
        fi
        sleep 0.2

        bowtie2_log=${odir}/log/bowtie2_log
        if [ ! -d "${bowtie2_log}" ]; then
            echo -e "  ${G}Created directory: ${bowtie2_log} ${N} \n"
            mkdir -p ${bowtie2_log}
        fi
        sleep 0.2

        if [ "$frag_120" == "TRUE" ]; then
            echo "[info] use Bowtie2 command: --dovetail --phred33"
            echo "[info] The dovetail mode is enabled [as parameter '-f|--frag_120' is on]"

            for i in $ff; do
                base=$(basename ${i})
                echo -e "  ${G}Sample: ${base}${N}"
                # echo ${i}_R1.fq.gz
                # echo ${odir}/${base}.bam

                bowtie2 -p $th --dovetail --end-to-end --very-sensitive --no-mixed --no-discordan -x ${bt2idx} \
                    -1 ${i}_R1.fq.gz -2 ${i}_R2.fq.gz \
                    2>$logdir/"$base".bowtie2 | samtools sort -@ 20 -O bam -o ${odir}/${base}.sorted.bam -

                samtools index ${odir}/${base}.sorted.bam

            done
        else
            echo "[info] use Bowtie2 command: --end-to-end --very-sensitive --no-mixed --no-discordan -I 10 -X 700"
            echo "[info] The dovetail mode is off [as parameter frag_120 is off]"
            # --very-sensitive-local
            for i in $ff; do
                base=$(basename ${i})
                echo -e "  ${G}Sample: ${base}${N}"
                #
                (bowtie2 -p $th --end-to-end --very-sensitive --no-mixed --no-discordan -I 10 -X 700 -x ${bt2idx} \
                    -1 ${i}_R1.fq.gz -2 ${i}_R2.fq.gz) \
                    2>${bowtie2_log}/${base}.bowtie2 |
                    samtools sort -@ 20 -O bam -o ${odir}/${base}.sorted.bam -
                # | samtools view -bS - > ${odir}/${base}.bam
                samtools index ${odir}/${base}.sorted.bam

            done
        fi
    fi
}
