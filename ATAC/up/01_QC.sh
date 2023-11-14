module load FastQC/0.11.7

if [ ! -d "rawdata" ]; then
 	echo please make sure : there is a dir named rawdata
  	exit 1
else
	echo *********Let Us Begin ^_^**********
fi

mkdir -p result logs tmp

# make the sample_for_QC.txt
realpath rawdata/*.gz | sed 's/_[12]\..*//' | uniq > tmp/sample_for_QC.txt

# make the output file
mkdir -p result/01_cleandata
mkdir -p logs/fastqc_v1
mkdir -p logs/fastp

# set up file names
cleandata_out=result/01_cleandata
fastqc_v1_log=logs/fastqc_v1
fastp_log=logs/fastp

# fastp
echo *********Fastp Begins ^_^*********

cat tmp/sample_for_QC.txt | while read id; 
do
    base=$(basename ${id})
    echo Sample name is ${base}

    #-a CTGTCTCTTATACACATCT \
    fastp --thread 16 \
    -i ${id}_1.fq.gz \
    -I ${id}_2.fq.gz \
    -o ${cleandata_out}/${base}_1.clean.fq.gz \
    -O ${cleandata_out}/${base}_2.clean.fq.gz \
    -j ${fastp_log}/${base}.fastp.json \
    -h ${fastp_log}/${base}.fastp.html \
    2> ${fastp_log}/${base}.logs

    echo Sample ${base} is done
done

echo *********Fastp Ends ^_^*********

# fastqc
fastqc -q -t 20 ${cleandata_out}/*.fq.gz -o ${fastqc_v1_log}


# multiqc
module load multiqc/1.8
multiqc -o ${fastqc_v1_log} ${fastqc_v1_log}
multiqc -o ${fastp_log} ${fastp_log}

