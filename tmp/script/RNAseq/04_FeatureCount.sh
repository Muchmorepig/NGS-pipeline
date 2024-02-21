#!bin/bash
module load subread/1.6.2

gtf=$1

cd result/02_alignment
ls -lh *.bam

featureCounts -T 16 -p -a ${gtf} -o ../AllCounts.txt *.bam 1> ../../logs/featureCounts.log
