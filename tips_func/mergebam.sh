txt=$1

cat ${txt}|while read id;
do 
	arr=($id)
        out=${arr[0]}
        bam1=${arr[1]} 
	bam2=${arr[2]}
	echo ${out}
	samtools merge ${out}.bam ${bam1}.flt.bam ${bam2}.flt.bam -@ 10 
done
