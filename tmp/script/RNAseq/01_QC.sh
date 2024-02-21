#!bin/bash
threads=$1

module load multiqc/1.8
module load FastQC/0.11.7

# make the sample_for_QC.txt
mkdir tmp
realpath rawdata/*.gz | sed 's/_R[12]\..*//' | uniq > tmp/sample_for_QC.txt

# make the output file
mkdir -p result/01_cleandata
mkdir -p logs/fastqc_v1
mkdir -p logs/fastp

# set up file names
cleandata_out=result/01_cleandata
fastqc_v1_log=logs/fastqc_v1
fastp_log=logs/fastp

# fastp
echo fastp begins

cat tmp/sample_for_QC.txt | while read id; 
do
    base=$(basename ${id})
    echo Sample name is ${base}

    fastp --thread ${threads} \
	--length_required 21 \
    -i ${id}_R1.fastq.gz \
    -I ${id}_R2.fastq.gz \
    -o ${cleandata_out}/${base}_1.clean.fq.gz \
    -O ${cleandata_out}/${base}_2.clean.fq.gz \
    -j ${fastp_log}/${base}.fastp.json \
    -h ${fastp_log}/${base}.fastp.html \
    2> ${fastp_log}/${base}.logs

    echo Sample ${base} is done
done

echo fastp ends

# fastqc_v1
fastqc -q -t ${threads+4} ${cleandata_out}/*.fq.gz -o ${fastqc_v1_log}

# multiqc
multiqc -o ${fastqc_v1_log} ${fastqc_v1_log}
multiqc -o ${fastp_log} ${fastp_log}
