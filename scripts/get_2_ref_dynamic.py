import csv
import sys
import os
from Bio import SeqIO

def read_hla_file(filename):
    # 使用字典来存储每个样本的HLA类型
    hla_data = {}

    with open(filename, 'r') as file:
        # 创建一个csv阅读器，假设字段是通过制表符分隔的
        
        # 遍历文件中的每一行
        for idx, row in enumerate(file):
            if idx!=2:
                continue
            items=row.strip().split("\t")
            for iti, it in enumerate(items):
                if iti==0:
                    continue
                allele=it.split(",")[0]
                gene_tag=allele.split("*")[0]
                gene_dict[gene_tag].append(allele)

def replace_single_contig_name(input_file, output_file, new_contig_name):
    # 读取FASTA文件
    record = next(SeqIO.parse(input_file, "fasta"))
    
    # 更改contig名字
    record.id = new_contig_name
    record.description = new_contig_name
    
    # 写入新的FASTA文件
    SeqIO.write([record], output_file, "fasta")
    print("Contig name replaced successfully.")

def build_ref():
    # Write the reference file
    # shell code
    for gene, alleles in gene_dict.items():
        gene_dir=f"{HLA_dir}/{gene}"
        if not os.path.exists(gene_dir):
            os.makedirs(gene_dir)
        db_ref=db_ref_dir+f"/HLA_{gene}.fasta"
        for allele_idx, allele in enumerate(alleles):
            allele_idx+=1
            print(f"processing {gene} {allele} {allele_idx}")
            cmd=f"""
                samtools faidx {db_ref} {allele}>{gene_dir}/HLA_{gene}.raw.{allele_idx}.fa 
            """
            os.system(cmd)
            replace_single_contig_name(f"{gene_dir}/HLA_{gene}.raw.{allele_idx}.fa", f"{gene_dir}/HLA_{gene}.{allele_idx}.fa", f"HLA_{gene}_ref{allele_idx}")
            
            os.system(cmd)
        merged_ref=f"{gene_dir}/HLA_{gene}.2allele.fa"


        cmd=f"""
            cat {gene_dir}/HLA_{gene}.1.fa  {gene_dir}/HLA_{gene}.2.fa > {merged_ref}
        """
        if len(set(alleles))==1:
            cmd=f"""
            cp {gene_dir}/HLA_{gene}.1.fa {merged_ref}
            """
        os.system(cmd)
        cmd=f"""
        samtools faidx {merged_ref}
        bwa index {merged_ref}
        faToTwoBit {merged_ref} {gene_dir}/HLA_{gene}.2allele.2bit
        makeblastdb -in {merged_ref} -dbtype nucl -parse_seqids -out {gene_dir}/HLA_{gene}.2allele
        """
        os.system(cmd)
        ref_count_file_f.write(f"{gene}\t{len(set(alleles))}\n")

def main():
    filename=f"{filedir}/hla.new.result.txt"
    read_hla_file(filename)
    print(gene_dict)
    build_ref()

if __name__ == "__main__":
    if len(sys.argv) <2:
        print("Usage: python script.py <filename>")
        sys.exit(1)

    filedir = sys.argv[1]
    sample=sys.argv[2]
    db_ref_dir=sys.argv[3]
    db_build_dir=sys.argv[4]
    ref_count_file=sys.argv[5]
    HLA_dir = f"{db_build_dir}/HLA"
    if not os.path.exists(HLA_dir):
        os.makedirs(HLA_dir)
    genes=["A","B","C","DPA1","DPB1","DQA1","DQB1", "DRB1"]
    gene_dict={}
    for gene in genes:
        gene_dict[gene]=[]
    ref_count_file_f=open(ref_count_file,"w")
    main()
    ref_count_file_f.close()