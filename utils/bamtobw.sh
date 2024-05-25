#!/bin/bash

input=${1}
output=${2}
threads=${3-20}

# make the output file
if [ ! -d ${output} ]; then
	mkdir -p ${output}
fi

# change bam to bw
for id in ./alignment/map2Tak1V6/*.bam; do
	base=$(basename "${id}")
	b=${base%%.*}
	bamCoverage -b ${id} \
		--binSize 10 --numberOfProcessors ${threads} \
		-o ${output}/${b}.bw --normalizeUsing BPM
done
