#!/usr/bin/bash
trim_files() {
  local indir=$1
  local odir=$2
  local th=$3

  if ! command -v fastp &>/dev/null; then
    echo -e "${R}Need tool: fastp ...${N}"
    exit
  fi

  if [ ! -d "${odir}/clean" ]; then
    echo -e "  ${G}Created directory: ${odir}/clean to store trimmed data${N}"
    mkdir -p ${odir}/clean
  fi
  sleep 0.1

  fplog=${odir}/log/fastp
  if [ ! -d "$fplog" ]; then
    echo -e "  ${G}Created directory: ${fplog} to store 'fastp' logs${N} \n"
    mkdir -p ${fplog}
  fi
  sleep 0.3

  date >&2
  echo >&2 -e "[info] Trimming files in $indir using 'fastp'\n"

  ff=$(ls ${indir}/*.gz | sed 's/_R[12]\..*//' | uniq)

  for i in $ff; do
    base=$(basename "${i}")
    echo -e "  ${G}Trimming Sample: ${base}${N}"
    fastp --thread $th \
      -i ${i}_R1.fq.gz -I ${i}_R2.fq.gz \
      -o ${odir}/clean/${base}_R1.fq.gz -O ${odir}/clean/${base}_R2.fq.gz \
      -j ${fplog}/${base}.fastp.json -h /dev/null \
      2> ${fplog}/${base}.logs
  done

  date >&2
  echo >&2 -e "[info] Trimming has finished \n"
}