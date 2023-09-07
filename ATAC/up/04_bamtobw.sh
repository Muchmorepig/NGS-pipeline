# set up the software environment
module load deeptools/2.0

# set threads
threads=${1-20}

# make the output file
mkdir -p result/05_bamtobw

# set up file names
input=result/03_filter_alignment
bw_out=result/05_bamtobw

# change bam to bw
ls ${input}/*.bam | while read id;
do
	base=$(basename ${id} .bam)

	bamCoverage -b ${id} \
	--binSize 10 --numberOfProcessors ${threads} \
	-o ${bw_out}/${base}.bw --normalizeUsing BPM
done

