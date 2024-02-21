index=$1
gtf=$2
th=$3

echo "********************QC**************************"
bash bin/01_QC.sh ${th}
echo "*******************ALIGN************************"
bash bin/02_align.sh ${index}
echo "******************BAMTOBW***********************"
bash bin/03_bamtobw.sh
echo "***************FEATURECOUNTS********************"
bash bin/04_FeatureCount.sh ${gtf}
echo "********************END*************************"
