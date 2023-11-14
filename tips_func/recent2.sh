
ls *.water|while read id
do
	i=${id%%.*}
	echo $i
	cat ${id}|grep "# "|sed '1,3d'|sed '3,12d'|sed 's/tity:/tity: /g'|sed 's/ \+/ /g'|sed 's/-//g'|sed 's/temp_fa\///g'|sed 's/.fa//g'|sed 's/# //g'|sed '$d'|sed 's/ ( / (/g'|sed 's/ (/(/g'\
 	> ${i}.log
done

paste *.log | sed -r 's/\t[a-zA-Z:]+ / /g' | sed -r 's/[0-9]{0,5}\/[0-9]{5}//g'|sed 's/(//g'|sed 's/)//g'|sed 's/://g'|sed 's/ /\t/g' > allwater.txt
