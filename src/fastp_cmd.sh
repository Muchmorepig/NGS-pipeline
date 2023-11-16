#!/bin/bash
trim_files() {
  local indir=$1
  local odir=$2
  local th=$3
  
  date >&2
  echo >&2 -e "[info] Trimming files in $indir using 'fastp'\n"

  ff=$(find ${indir} -name '*.gz' | sed -E 's/_R[12]\..*//' | sort | uniq)

  for i in $ff; do
    base=$(basename "${i}")
    if [[ ! -f "${i}_R1.fq.gz" || ! -f "${i}_R2.fq.gz" ]]; then
      echo -e "${R}Error: input files not found for sample ${base} ...${N}"
      exit 1
    fi
    echo -e "  ${G}Trimming Sample: ${base}${N}"

    mkdir -p "${odir}"
    mkdir -p "${odir}/log"

    fastp --thread $th \
      -i ${i}_R1.fq.gz -I ${i}_R2.fq.gz \
      -o ${odir}/${base}_R1.fq.gz -O ${odir}/${base}_R2.fq.gz \
      -j ${odir}/log/${base}.fastp.json -h /dev/null \
      2>${odir}/log/${base}.logs
  done

  date >&2
  echo >&2 -e "[info] Trimming has finished \n"
}
