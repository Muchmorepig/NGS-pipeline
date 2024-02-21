#!bin/bash
# set up the software environment
module load multiqc/1.8
module load HISAT/2.2.0
module load samtools/1.9

# set directory with hisat2 genome index
index=${1}

# make the sample_for_align.txt
realpath result/01_cleandata/*.gz | sed 's/_[12]\..*//' | uniq > tmp/sample_for_align.txt

# make the output dir
mkdir -p result/02_alignment
mkdir -p logs/hisat2_alignment

# set up file names
alignment_out=result/02_alignment
hisat2_log=logs/hisat2_alignment
echo ${alignment_out}
echo ${hisat2_log}
# alignment using hisat2
echo hisat2 begins

cat tmp/sample_for_align.txt | while read id;
do
    base=$(basename ${id})
    echo The sample name is ${base}

    hisat2 -t -p 16 --new-summary -x ${index} \
    -1 ${id}_1.clean.fq.gz \
    -2 ${id}_2.clean.fq.gz \
    2> ${hisat2_log}/${base}.log \
    | samtools sort -O bam -@ 8 -o ${alignment_out}/${base}.sort.bam -

    samtools index -@ 8 ${alignment_out}/${base}.sort.bam
done

echo hisat2 ends

# multiqc
multiqc -o ${hisat2_log} ${hisat2_log}


