#!/bin/bash
galore() {
  local indir=$1
  local odir=$2

  date >&2
  echo >&2 -e "[info] Trimming files in $indir using 'trim_galore'\n"

  ff=$(find ${indir} -name '*.gz' | sed -E 's/_R[12]\..*//' | sort | uniq)

  for i in $ff; do
    base=$(basename "${i}")
    if [[ ! -f "${i}_R1.fq.gz" || ! -f "${i}_R2.fq.gz" ]]; then
      echo -e "${R}Error: input files not found for sample ${base} ...${N}"
      exit 1
    fi
    echo -e "  ${G}Trimming Sample: ${base}${N}"

    mkdir -p "${odir}/log"

    trim_galore --paired \
      -o ${odir} \
      ${i}_R1.fq.gz ${i}_R2.fq.gz \
      2>${odir}/log/${base}.log &

  done

  wait

  date >&2
  echo >&2 -e "[info] Trimming has finished \n"
}