#!/bin/bash

gene_list=$1

gtf_file="/data5/wmc_data/reference/Mpolymorpha/tak1_V6/MpTak1_v6.1r2.agat.gtf"
genome="/data5/wmc_data/reference/Mpolymorpha/tak1_V6/MpTak1_v6.1r2.genome.fasta"

while IFS= read -r gene; do
    awk -v gene="$gene" '$0 ~ gene && ($3 == "gene" || $3 == "CDS")' "$gtf_file" > temp.gtf
    
    if [ ! -s temp.gtf ]; then
        echo "Warning: Gene or CDS for $gene not found in $gtf_file"
        continue
    fi

    gffread temp.gtf -g "$genome" -w "${gene}.fasta" && echo "Extracted whole gene sequence for $gene"

    gffread temp.gtf -g "$genome" -x "${gene}_cds.fasta" && echo "Extracted CDS for $gene"
done < "$gene_list"

rm temp.gtf

echo "Finish"

