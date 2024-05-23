# module load cellranger/7.0.0
id=$1
sample=$2
tr=$3
fq=$4

cellranger count --id=$id \
	--sample=$sample \
	--transcriptome=$tr \
	--fastqs=$fq \
	--nosecondary \
	--localcores=30 \

# nohup bash celcount.sh cim22d SC20210915S1 ~/reference/Oryza_sativa/forSCRNA/IRGSP_V1 ./rawdata/cim22d 2>&1 > cim22d.log &
# nohup time bash celcount.sh sim12d SC20200629S1 ~/reference/Oryza_sativa/forSCRNA/IRGSP_V1 ./rawdata/sim12d 2>&1 > sim12d.log &
# nohup time bash celcount.sh sim4d SC20211014S1 ~/reference/Oryza_sativa/forSCRNA/IRGSP_V1 ./rawdata/sim4d 2>&1 > sim4d.log &
