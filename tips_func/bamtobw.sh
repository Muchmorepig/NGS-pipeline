# set up the software environment
# module load deeptools/2.0

input=${1}
# set threads
threads=${2-20}

# make the output file
mkdir bw

# set up file names
bw_out=bw

# change bam to bw
ls ${input}/*.bam | while read id;
do
	base=$(basename ${id} .bam)

	bamCoverage -b ${id} \
	--binSize 10 --numberOfProcessors ${threads} \
	-o ${bw_out}/${base}.bw --normalizeUsing BPM
done

