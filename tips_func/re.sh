cat tmp.txt|while read id
do
    arr=($id)
    or=${arr[0]}
    re=${arr[1]}

    mv ${or}_R1.fq.gz ${re}_R1.fq.gz
    mv ${or}_R2.fq.gz ${re}_R2.fq.gz
done
