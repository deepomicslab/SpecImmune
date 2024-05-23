import sys
import os
from Bio import SeqIO


def find_highest_depth_alleles(input_file):
    gene_allele_depth = {}
    # Read the file and parse the data
    with open(input_file, 'r') as file:
        next(file)  # Skip the header row
        for line in file:
            parts = line.strip().split()
            if len(parts) != 4:
                continue

            gene, allele, depth, _ = parts

            # Attempt to convert depth to float, skip line if conversion fails
            try:
                depth = float(depth)
            except ValueError:
                continue  # Skip entries where depth is not a float

            # Check for the highest depth allele for each gene
            if gene not in gene_allele_depth or gene_allele_depth[gene][1] < depth:
                gene_allele_depth[gene] = (allele, depth)

    # Output the highest depth allele for each gene
    for gene, (allele, depth) in gene_allele_depth.items():

        print(f'Gene: {gene}, Highest Depth Allele: {allele}, Depth: {depth}')
    return gene_allele_depth


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
    for gene, (allele, depth) in gene_allele_depth.items():
        gene_dir=f"{outdir}/HLA/{gene}"
        if not os.path.exists(gene_dir):
            os.makedirs(gene_dir)

        print(f"processing {gene} {allele} {depth}")
        cmd=f"""
            samtools faidx {ref_db} {allele}>{gene_dir}/HLA_{gene}.raw.fa 
        """
        os.system(cmd)
        replace_single_contig_name(f"{gene_dir}/HLA_{gene}.raw.fa", f"{gene_dir}/HLA_{gene}.fa", f"HLA_{gene}")
        cmd=f"""
        samtools faidx {gene_dir}/HLA_{gene}.fa
        bwa index {gene_dir}/HLA_{gene}.fa
        faToTwoBit {gene_dir}/HLA_{gene}.fa {gene_dir}/{gene}.2bit
        makeblastdb -in {gene_dir}/HLA_{gene}.fa -dbtype nucl -parse_seqids -out {gene_dir}/HLA_{gene}
        """
        os.system(cmd)


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python script.py input_file")
        sys.exit(1)

    input_file = sys.argv[1]
    ref_db= sys.argv[2] 
    outdir= sys.argv[3]
    HLA_dir = f"{outdir}/HLA"
    if not os.path.exists(HLA_dir):
        os.makedirs(HLA_dir)
    
    gene_allele_depth = find_highest_depth_alleles(input_file)
    build_ref()
