index=$1
gtf=$2
th=$3

echo "********************QC**************************"
module load multiqc/1.8
module load FastQC/0.11.7

bin/01_QC.sh ${th}

echo "*******************ALIGN************************"
module load multiqc/1.8
module load HISAT/2.2.0
module load samtools/1.9

bin/02_align.sh ${index}

echo "******************BAMTOBW***********************"
module load deeptools/2.0

bin/03_bamtobw.sh

echo "***************FEATURECOUNTS********************"
module load subread/1.6.2

bin/04_FeatureCount.sh ${gtf}

echo "********************END*************************"
