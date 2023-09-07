#bash bin/01_QC.sh
bash bin/02_Align.sh ~/index/bowtie2/Oryza_sativa/Oryza_sativa 16
bash bin/03_callpeak.sh 3.7e+8
bash bin/04_bamtobw.sh
