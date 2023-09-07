module load sratoolkit/2.9.2 
sra=$1
dir=$2
cat $sra|while read i
do
	fastq-dump --gzip --split-3 \
	--defline-qual '+' --defline-seq '@$ac-$si' sra/${i}.sra -O $dir
done
