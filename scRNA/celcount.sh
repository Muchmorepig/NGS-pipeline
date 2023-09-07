module load cellranger/3.1.0
id=$1
fq=$2
cel=$3
cellranger count --id=$id \
	--sample=$id \
	--transcriptome=Araport11V3 \
	--fastqs=$fq \
	--nosecondary \
	--localcores=70 \
	--expect-cells=$cel
