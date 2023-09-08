#!/usr/bin/bash
input=$1

mkdir bw


# change bam to bw
ls ${input}/*flt.bam | while read id;
do
	base=$(basename ${id} .flt.bam)

	bamCoverage -b ${id} \
	--binSize 10 --numberOfProcessors 20 \
	-o bw/${base}.bw --normalizeUsing BPM
done
