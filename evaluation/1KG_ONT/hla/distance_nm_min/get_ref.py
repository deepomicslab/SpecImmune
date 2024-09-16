#!/usr/bin/env python3

import os
from Bio import SeqIO

def extract_first_fasta(input_file, output_dir):
    gene_seen = set()  # 用于记录已经处理的基因

    # 确保输出目录存在，如果不存在则创建
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 使用SeqIO解析输入的FASTA文件
    for record in SeqIO.parse(input_file, "fasta"):
        # 提取基因名，假设格式是例如 '>HLA-A*01:01:01:01'
        gene_name = record.id.split('*')[0]
        
        # 如果该基因未被处理过
        if gene_name not in gene_seen:
            gene_seen.add(gene_name)
            
            # 创建输出文件路径，以基因名命名
            output_file = os.path.join(output_dir, f"{gene_name}.fasta")
            
            # 将序列写入单独的输出文件
            with open(output_file, 'w') as outfile:
                SeqIO.write(record, outfile, "fasta")

if __name__ == "__main__":
    # 输入和输出文件路径
    input_fasta = "input.fasta"   # 替换为你的输入文件路径
    output_dir = "output_dir"     # 替换为你的输出目录路径
    
    # 调用函数进行提取
    extract_first_fasta(input_fasta, output_dir)