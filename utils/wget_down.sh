cat md5_link.txt|while read id; do
	base=${id##*Data/}
	name=${base%%\?*}
	n=$(echo ${name} | sed 's/\//_/g')
	echo ${n}
	wget -O ${n} ${id} &
done
