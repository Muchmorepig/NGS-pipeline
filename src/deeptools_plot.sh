#!/bin/bash
input=${1}
output=${2}
prefix=${3}
bed_file=${4}

threads=20

if [ ! -f "${bed_file}" ]; then
	echo "BED file does not exist."
	exit 1
fi

if [ ! -d "${output}" ]; then
	mkdir -p ${output}
fi

computeMatrix scale-regions -S ${input} -R ${bed_file} \
	--beforeRegionStartLength 2000 --regionBodyLength 3000 --afterRegionStartLength 2000 \
	--numberOfProcessors ${threads} \
	--skipZeros -o ${output}/matrix.mat.gz


computeMatrix reference-point \
       --referencePoint center -b 2000 -a 2000 \
       -R ${bed_file} \
       -S ${input} \
       --skipZeros --numberOfProcessors ${threads} \
       -o ${output}/matrix_center.gz

# plotHeatmap -m ${output}/matrix_center.gz -out ${output}/Heatmap1.pdf
plotProfile -m ${output}/matrix_center.gz -out ${output}/${prefix}_profile_center.png
plotProfile -m ${output}/matrix.mat.gz -out ${output}/${prefix}_profile_tsstes.png

rm ${output}/matrix_center.gz
rm ${output}/matrix.mat.gz
