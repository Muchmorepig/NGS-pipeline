txt=$1
outdir=$2
module load EMBOSS/6.60

#mkdir temp_fa

if [ ! -d "$outdir" ]; then
  mkdir $outdir
fi


cat $txt|while read id;
do
	arr=($id)
	v6=${arr[0]%%.*}
	v5=${arr[1]%%.*}
	echo $v6 $v5
	
	cat Mpolymorpha/tak2v6_gene.bed|grep $v6 > temp_v6.bed
	cat Mpolymorpha/tak1v5_gene.bed|grep $v5 > temp_v5.bed
	seqkit subseq --bed temp_v6.bed -u 5000 -d 5000 Mpolymorpha/tak2v6_male.fasta > ${outdir}/${v6}.fa
	seqkit subseq --bed temp_v5.bed -u 5000 -d 5000 Mpolymorpha/tak1v5.1_with_female_chr.fasta > ${outdir}/${v5}.fa
	needle -asequence ${outdir}/${v6}.fa -bsequence ${outdir}/${v5}.fa -auto -outfile ${outdir}/${v6}_${v5}.needle
	water -asequence ${outdir}/${v6}.fa -bsequence ${outdir}/${v5}.fa -auto -outfile ${outdir}/${v6}_${v5}.water
done

rm temp*.bed

cd $outdir
ls *.water|while read id
do
	i=${id%%.*}
	echo $i
	cat ${id}|grep "# "|sed '1,3d'|sed '3,12d'|sed 's/tity:/tity: /g'|sed 's/ \+/ /g'|sed 's/-//g'|sed 's/temp_fa\///g'|sed 's/.fa//g'|sed 's/# //g'|sed '$d'|sed 's/ ( / (/g'|sed 's/ (/(/g'\
 	> ${i}.log
done

paste *.log | sed -r 's/\t[a-zA-Z:]+ / /g' | sed -r 's/[0-9]{0,5}\/[0-9]{5}//g'|sed 's/(//g'|sed 's/)//g'|sed 's/://g'|sed 's/ /\t/g'|sed 's/%//g'|sed -r 's/ity/ity(%)/g' > allwater.txt

awk 'BEGIN{OFS="\t"}{for(i=1;i<=NF;i++)a[NR,i]=$i}END{for(j=1;j<=NF;j++)for(k=1;k<=NR;k++)printf k==NR?a[k,j]RS:a[k,j] OFS}' allwater.txt > all_water.txt

mkdir {log,fa,needle,water}
mv *.fa fa/
rm fa -r
mv *.needle needle/
mv *.water water/
mv *.log log/

