# module load cellranger/7.0.0
id=$1
sample=$2
tr=$3
fq=$4

cellranger count \
	--id=$id \
	--sample=$sample \
	--transcriptome=$tr \
	--fastqs=$fq \
	--nosecondary --localcores=45
