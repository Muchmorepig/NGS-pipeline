#!bin/bash
# set up the software environment
module load deeptools/2.0

# make the output file
mkdir -p result/03_bamtobw

# set up file names
input=result/02_alignment
bw_out=result/03_bamtobw

# change bam to bw
ls ${input}/*.bam | while read id;
do
	base=$(basename ${id} .bam)

	bamCoverage -b ${id} \
	--binSize 10 --numberOfProcessors 16 \
	-o ${bw_out}/${base}.bw --normalizeUsing RPKM
done

