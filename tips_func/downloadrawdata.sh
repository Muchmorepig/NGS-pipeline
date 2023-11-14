cat wget.txt|while read id;
do
	arr=($id)
	name=${arr[0]}
	url=${arr[1]}
	echo $name
	echo ${url}
	wget -c ${url} -O ./${name}
done		
