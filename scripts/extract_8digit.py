import pysam
import argparse

# 定义基因列表
genes = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]

def extract_gene_and_allele(comment):
    # 从描述字段中提取基因信息并检查是否有四个冒号
    parts = comment.split()
    if len(parts) > 0:
        #print(parts)
        gene_allele = parts[0]
        if gene_allele.count(":") == 3:
            gene = gene_allele.split("*")[0]
            allele = gene_allele.split("*")[1]
            return gene, allele
    return None, None

def main(input_fasta):
    # 创建一个字典来存储每个基因的序列
    gene_sequences = {gene: [] for gene in genes}

    # 使用pysam读取FASTA文件
    with pysam.FastxFile(input_fasta) as fasta_file:
        for entry in fasta_file:
            # 提取contig name和序列
            contig_name = entry.name
            sequence = entry.sequence
            description = entry.comment  # 获取描述信息

            # 提取基因和等位基因信息
            gene, allele = extract_gene_and_allele(description) if description else (None, None)
            print(gene, allele)
            if gene and allele and gene in genes:
                print(gene, allele)
                new_contig_name = f">{gene}*{allele}"
                gene_sequences[gene].append((new_contig_name, sequence))
            else:
                pass
                print(f"Skipping contig with invalid gene or allele format: {description}")

    # 写入新的FASTA文件
    for gene, sequences in gene_sequences.items():
        if sequences:
            output_fasta = f"/mnt/d/HLAPro_backup/Nanopore_optimize/SpecHLA/db2/HLA_{gene}.fasta"
            with open(output_fasta, "w") as outfile:
                for new_contig_name, sequence in sequences:
                    outfile.write(f"{new_contig_name}\n")
                    outfile.write(f"{sequence}\n")

    print("FASTA files have been successfully generated.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process and split a FASTA file by gene.")
    parser.add_argument("input_fasta", help="Input FASTA file")
    args = parser.parse_args()
    
    main(args.input_fasta)
