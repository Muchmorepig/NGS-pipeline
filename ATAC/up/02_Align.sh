# set up the software environment
module load multiqc/1.8
module load bowtie/2.3.4.3
module load bedtools/2.25.0
module load samtools/1.10

index=${1}

# set threads
threads=${2}

echo "your index: ${index}"

# make the sample_for_align.txt
realpath result/01_cleandata/*.gz | sed 's/_[12]\..*//' | uniq > tmp/sample_for_align.txt

# make all of the output directories
mkdir -p result/02_alignment
mkdir -p result/03_filter_alignment
mkdir -p logs/bowtie2_alignment
mkdir -p logs/sambamba_markdup

# set up file names
alignment_out=result/02_alignment
filter_alignment_out=result/03_filter_alignment
bowtie2_log=logs/bowtie2_alignment
sambamba_log=logs/sambamba_markdup

# alignment using bowtie2
echo ========== Bowtie2 Begins O_O ==========

cat tmp/sample_for_align.txt | while read id;
do
	base=$(basename ${id})
	echo The sample name is ${base}

    bowtie2 -p ${threads} -x ${index} \
    -1 ${id}_1.clean.fq.gz \
    -2 ${id}_2.clean.fq.gz 2> ${bowtie2_log}/${base}.log \
    | samtools sort -@ 20 -O bam -o ${alignment_out}/${base}.sorted.bam - 
    
    samtools index ${alignment_out}/${base}.sorted.bam

    echo The sample ${base} is done
done

echo ========== Bowtie2 Ends ^_^ ==========

# markdup with sambamba
echo ---------- Sambamba Begins ----------

for i in ${alignment_out}/*.sorted.bam
do
	base=$(basename ${i} .sorted.bam)
	echo The sample name is ${base}

	sambamba markdup -t 5 ${i} ${alignment_out}/${base}.sorted.markdup.bam \
	2> ${sambamba_log}/${base}.log

	echo The sample ${base} is done
done

echo ---------- Sambamba Ends ^_^ ----------


# filter duplcation、multi-mappers、low_quality reads with samtools
echo ========== Filter begins ==========

for j in ${alignment_out}/*.sorted.markdup.bam
do
	base=$(basename ${j} .sorted.markdup.bam)
	echo The sample name is ${base}

	samtools view -@ 20 -bF 1804 -q 20 ${j} -o ${filter_alignment_out}/${base}.flt.bam 

	samtools index -@ 20 ${filter_alignment_out}/${base}.flt.bam
	echo The sample ${base} is done
done

echo ========== Filter Ends ^_^ ==========

# filter ChrM ChrCh with bedtools
#echo bedtools filter organelle begins

#for i in ${filter_alignment_out}/*.flt.bam
#do
#	base=$(basename ${i} .flt.bam)
#	echo The sample name is ${base}
#
#    bedtools intersect -abam ${i} -b filter.bed -v > ${filter_alignment_out}/${base}.rm_organelle.bam
#    samtools index ${filter_alignment_out}/${base}.rm_organelle.bam
#    
#    echo The sample ${base} is done
#done

#echo bedtools filter is done

echo Begin remove tmp files
rm -f ${alignment_out}/*.sorted.markdup.bam
rm -f ${alignment_out}/*.sorted.markdup.bam.bai

#rm -f ${filter_alignment_out}/*.flt.bam
#rm -f ${filter_alignment_out}/*.flt.bam.bai
echo -------- Removment Complete --------

# multiqc 
module load multiqc/1.8
multiqc -o ${bowtie2_log} ${bowtie2_log}
