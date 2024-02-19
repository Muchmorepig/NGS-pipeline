#!/bin/bash

indir=$1
info=$2

cat $info | while read line; do
    echo $line
    arr=($line)
    control=${arr[0]}
    treatment=${arr[1]}
    # mamba run -n so bamCompare \
    bamCompare \
        -b1 "$indir/$treatment.flt.bam" -b2 "$indir/$control.flt.bam" \
        -o "$indir/${treatment}.bw" \
        --binSize 10 \
        --operation subtract \
        --numberOfProcessors 32
    # --normalizeUsing BPM \
done
