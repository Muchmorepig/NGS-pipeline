#!/usr/bin/bash
odir=cuttag
indir=$wkdir/rawdata
genome=tair10
th=8
trim=false
frag_120=false
align=false
shift=false
filter=false
call=false
help=false

G="\033[0;32m"
R="\033[0;31m"
N="\033[0m"
B="\033[1m"

# 使用方法
usage() {
  echo "需要工具 fastp bowtie2 samtools sambamba deeptools MACS2;
    请先运行 module load ..." >&2
  echo "Usage: bash $0 [-i fastq/path] [-o output/dir] -taslc -g tair10" >&2
  echo -e "
    -i | --input    双端测序文件位置 (eg: fq.gz 格式为 sample_R1.fq.gz sample_R2.fq.gz)
    -o | --output   项目路径, 后续结果文件都会保存于此 (Defult: cuttag)
    -g | --genome   基因组, 默认是tair10, 其他可选: mpTak1_V5、mpTak1_V6、mpTak2、irgsp
    -t | --trim     是否去接头 (Defult: false)
    -a | --align    是否比对到基因组 (Defult: false)
    -l | --filter   是否去除PCR重复等 (Defult: false)
    -s | --shift    是否执行'deeptools'的ATAC-shift (Defult: false) 
    -c | --callpeak 是否执行 Call Peaks (Defult: false)
    -p | --parallel 比对等步骤用的线程数 (Defult: 8)" >&2
  exit 1
}

if [ $# -eq 0 ]; then
  usage
  exit 0
fi

# Parse command line arguments
# 使用 getopt 整理参数
ARGS=$(getopt \
  -o 'i:o:g:p:tafslchL::' \
  -l 'input:,output:,genome:,parallel:,trim,align,frag_120,shift,callpeak,help,log-level::' -- "$@")

if [ $? != 0 ]; then
  echo "Parse error! Terminating..." >&2
  exit 1
fi
# 将参数设置为 getopt 整理后的参数
# $ARGS 需要用引号包围
eval set -- "$ARGS"

while true; do
  case "$1" in
  -i | --input)
    indir="$2"
    shift 2
    ;;
  -o | --output)
    odir="$2"
    shift 2
    ;;
  -g | --genome)
    genome="$2"
    shift 2
    ;;
  -p | --parallel)
    th="$2"
    shift 2
    ;;
  -t | --trim)
    trim=true
    shift
    ;;
  -a | --align)
    align=true
    shift
    ;;
  -f | --frag_120)
    frag_120=true
    shift
    ;;
  -s | --shift)
    shift=true
    shift
    ;;
  -l | --filter)
    filter=true
    shift
    ;;
  -c | --callpeak)
    call=true
    shift
    ;;
  -h | --help)
    help=true
    shift
    ;;
  -L | --log-level)
    case "$2" in
    "")
      CONN_LOG_LEVEL=1
      shift 2
      ;;
    *)
      CONN_LOG_LEVEL="$2"
      shift 2
      ;;
    esac
    ;;
  --)
    shift
    break
    ;;
  *)
    echo "Internal error!"
    exit 1
    ;;
  esac
done

if [ "$help" == "true" ]; then
  usage
  exit 0
fi

indir=$(realpath $indir)
odir=$(realpath $odir)

echo >&2 -e "[info] Input Fastq Folder: ${B}$indir${N}"
echo >&2 -e "[info] Result will be saved in: ${B}$odir${N} \n"
sleep 0.2

file_array=()
for file in "${indir}"/*fq.gz; do
  if [[ -f "$file" ]]; then
    filename=$(basename "$file")
    file_array+=("$filename")
  fi
done

# 输出样本信息
sn=$(printf "%s\n" "${file_array[@]}" | sed 's/_R[12]\..*//' | uniq)

echo -e "${B}Sample:${N}"
for i in ${sn}; do printf "  %s" "${i}"; done

# 设置对应索引和基因组大小
script_path=$(dirname $(readlink -f $0))
source ${script_path}/../src/cmd_check.sh
source ${script_path}/../src/get_info.sh

get_info $genome

echo -e "\n${B}Genome:${N}\n  $genome\n${B}Size:${N}\n  $gs\n${B}Bowtie2-index:${N}\n  $bt2idx"
echo " "
sleep 0.5

echo >&2 -e "[info] Trim: ${B}$trim${N}"

if [ "$trim" == "true" ]; then
  command_exists "fastp"
  source ${script_path}/../src/fastp_cmd2.sh
  trim_files $indir $odir/clean
else
  date >&2
  echo >&2 -e "[info] Skipping Trim ... \n"
  sleep 0.2
fi

echo >&2 -e "[info] Align: ${B}$align${N}"

if [ "$align" == "true" ]; then
  command_exists "bowtie2"
  command_exists "samtools"
  source ${script_path}/../src/bowtie2_align.sh
  bowtie_align $odir/clean $odir/align $th "normal"

  date >&2
  echo >&2 -e "[info] Aligning Finished...\n"
  sleep 0.2
else
  date >&2
  echo -e "[info] Skip Align...\n"
fi

echo >&2 -e "[info] Filter: ${B}$filter${N}"

if [ "$filter" == "true" ]; then
  # mark PCR-dup
  echo >&2 "[info] Mark PCR-Duplicate"
  source ${script_path}/../src/filter.sh
  filter_bam ${odir}/align 16

  # rm ${odir}/align/*sorted.bam*
else
  date >&2
  echo -e "[info] Skip Filter...\n"
fi

# use --ATACshift
echo >&2 -e "[info] Shift: ${B}$shift${N}"

if [ "$shift" == "true" ]; then

  date >&2
  echo >&2 "[info] Reads shifting "
  echo >&2 -e "${G}[info]${N} ${B}All reads aligning to the + strand were offset by +4 bp,
       and all reads aligning to the - strand were offset -5 bp, 
       since the CUT&Tag approach uses Tn5 transposase which has 
       been shown to bind as a dimer and insert two adaptors separated by 9 bp...${N}\n"

  echo >&2 "[info] The shifted BAM files will be used in the following analysis "
  echo >&2 "[info] Shifting and indexing "

  for dup in ${odir}/align/*flt.bam; do
    # echo $dup
    base=$(basename $dup .flt.bam)
    echo -e "\t${G}${base}${N}"
    alignmentSieve --numberOfProcessors $th --ATACshift \
      --bam $dup -o ${odir}/align/${base}_shift.bam
    samtools sort -@ 20 -O bam -o ${odir}/align/${base}_shift.sorted.bam ${odir}/align/${base}_shift.bam 2>${odir}/log/tmp.log
    samtools index ${odir}/align/${base}_shift.sorted.bam
    rm ${odir}/align/${base}_shift.bam
  done

  bam_file=${odir}/align/*shift.sorted.bam
else
  date >&2
  echo >&2 -e "[info] Skip ATAC-Shift...\n"
  bam_file=${odir}/align/*flt.sorted.bam
fi

# bam_file=${odir}/align/*shift.sorted.bam
echo >&2 -e "[info] Call Peak: ${B}$call${N}"

if [ "$call" == "true" ]; then
  date >&2
  echo >&2 "[info] Peaks Calling using MACS2... "

  if [ ! -d "${odir}/peak" ]; then
    echo -e "  ${G}Created directory: ${odir}/peak/narrow & peak/broad ${N}"
    mkdir -p ${odir}/peak/narrow
    mkdir -p ${odir}/peak/broad
  fi
  sleep 0.2

  macs2_log=${odir}/log/macs2
  if [ ! -d "${macs2_log}" ]; then
    echo -e "  ${G}Created directory: ${macs2_log} to store 'MACS2' logs${N}\n"
    mkdir -p ${macs2_log}
  fi
  sleep 0.2

  bedopsbin="/data3/wanglab/wmc/tools/bin"

  for shift in ${bam_file}; do
    base=$(basename ${shift} _shift.sorted.bam)
    echo -e "  ${G}Sample: ${base}${N}"
    macs2 callpeak -t ${shift} -f BAMPE \
      -n ${base} \
      -g $gs \
      --outdir ${odir}/peak/narrow 2>${macs2_log}/${base}.narrow

    macs2 callpeak -t ${shift} -f BAMPE \
      -g $gs \
      -n ${base} \
      --outdir ${odir}/peak/broad \
      --broad --broad-cutoff 0.1 2>${macs2_log}/${base}.broad

    python /data5/wmc_data/cuttag/scr/get_summits_broadPeak.py peak/broad/${base}_peaks.broadPeak |
      $bedopsbin/sort-bed - >peak/broad/${base}_summits.bed
  done
else
  date >&2
  echo >&2 -e "[info] Skip Peaks Calling...\n"
fi

echo >&2 "#########################"
echo >&2 "#     All Finished      #"
echo >&2 "#         ^ _ ^         #"
echo >&2 "#########################"

# 将已经处理过的选项(参数)移除掉
shift $((OPTIND - 1))
