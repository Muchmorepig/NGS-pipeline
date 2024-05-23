#!/bin/bash

# 输入的基因列表文件
gene_list=$1

# 参考文件路径
gtf_file="/data5/wmc_data/reference/Mpolymorpha/tak1_V6/MpTak1_v6.1r2.agat.gtf"
genome="/data5/wmc_data/reference/Mpolymorpha/tak1_V6/MpTak1_v6.1r2.genome.fasta"

# 检查输入文件是否存在
if [ ! -f "$gene_list" ]; then
    echo "Gene list file not found!"
    exit 1
fi

# 创建主目录用于保存所有生成的文件
main_dir="gene_sequences"
mkdir -p "$main_dir"

# 循环读取每个基因
while IFS= read -r gene; do
    # 创建基因名称对应的目录
    gene_dir="$main_dir/$gene"
    mkdir -p "$gene_dir"
    
    # 使用临时文件更安全
    temp_gtf=$(mktemp)
    
    # 提取基因和CDS的信息到临时文件
    awk -v gene="$gene" '$0 ~ gene && ($3 == "gene" || $3 == "CDS")' "$gtf_file" > "$temp_gtf"
    
    if [ ! -s "$temp_gtf" ]; then
        echo "Warning: Gene or CDS for $gene not found in $gtf_file"
        rm "$temp_gtf"
        continue
    fi

    # 提取基因序列
    seqkit subseq --quiet \
        --gtf "$temp_gtf" \
        --feature gene \
        -u 10000 -d 10000 "$genome" > "${gene_dir}/${gene}.fasta"
        
    if [ $? -eq 0 ] && [ -s "${gene_dir}/${gene}.fasta" ]; then
        echo "Extracted Gene for $gene"
    else
        echo "Failed to extract Gene for $gene"
    fi

    # 提取CDS序列
    seqkit subseq --quiet \
        --gtf "$temp_gtf" \
        --feature CDS "$genome" > "${gene_dir}/${gene}_cds.fasta"

    if [ $? -eq 0 ] && [ -s "${gene_dir}/${gene}_cds.fasta" ]; then
        echo "Extracted CDS for $gene"
        
        # 分割每个CDS条目到单独文件
        awk -v gene_dir="$gene_dir" -v gene="$gene" '
        /^>/ {
            if (seq) {
                filename = gene_dir "/" gene "_" substr(header, 2) ".fasta"
                print ">" header "\n" seq > filename
                seq = ""
            }
            header = $0
            next
        }
        {
            seq = seq $0
        }
        END {
            if (seq) {
                filename = gene_dir "/" gene "_" substr(header, 2) ".fasta"
                print ">" header "\n" seq > filename
            }
        }' "${gene_dir}/${gene}_cds.fasta"

    else
        echo "Failed to extract CDS for $gene"
    fi

    # 删除临时文件
    rm "$temp_gtf"
done < "$gene_list"

echo "Finish"
