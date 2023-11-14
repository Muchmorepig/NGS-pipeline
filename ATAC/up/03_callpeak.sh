module load MACS/2.1.2

genome_size=${1}

mkdir -p result/04_callpeak/merge_peak
mkdir -p result/04_callpeak/div_peak
#mkdir -p logs/macs2/merge_peak
mkdir -p logs/macs2/div_peak

input_file=result/03_filter_alignment
div_callpeak_out=result/04_callpeak/div_peak
div_callpeak_log=logs/macs2/div_peak

# MACS2
echo Macs2 div_callpeak will begin 

ls ${input_file}/*.bam | while read id;
do
    base=$(basename ${id} .flt.bam)
	echo ${base}

    macs2 callpeak -t ${id} -f BAMPE \
    -n ${base} \
    -g ${genome_size} \
    --outdir ${div_callpeak_out} 2> ${div_callpeak_log}/${base}.log 
done
