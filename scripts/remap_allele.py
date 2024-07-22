
import sys
import os
import pysam
import gzip
import argparse
from db_objects import My_db
import csv
import sys
import os
from Bio import SeqIO
from determine_gene import get_folder_list
from alignment_modules import Read_Type
from folder_objects import My_folder

def remove_characters(s):
    """
    Remove '*' and ':' characters from the input string.

    Parameters:
    s (str): The input string.

    Returns:
    str: The string with '*' and ':' characters removed.
    """
    s = s.replace('*', '')
    s = s.replace(':', '')
    return s


def read_hla_file(filename, some_dict):
    with open(filename, 'r') as file:
        next(file)
        next(file)
        for idx, row in enumerate(file):
            # print(row.strip())
            items=row.strip().split("\t")
            if len(items)<3:
                # remove gene_tag from some_dict
                gene_tag = items[0]
                some_dict[gene_tag].append("-")
                continue
            elif items[2] == "":
                continue
            alleles_str=items[2].split(";")
            print("alleles_str", alleles_str)
            for iti, it in enumerate(alleles_str):
                if iti>0:
                    continue
                allele=it.split(",")[0]
                # gene_tag=allele.split("*")[0]
                gene_tag = items[0]
                # print("result :", gene_tag, allele) 
                some_dict[gene_tag].append(allele)

def parse_full_allele(filename, full_allele_dict):
    with open(filename, 'r') as file:
        next(file)
        next(file)
        for idx, row in enumerate(file):
            print(row.strip())
            items=row.strip().split("\t")
            if len(items)<3:
                # remove gene_tag from some_dict
                gene_tag = items[0]
                continue
            elif items[2] == "":
                continue
        
            alleles_str=items[2].split(";")
            # print("alleles_str", alleles_str)
            for iti, it in enumerate(alleles_str):
                allele=it.split(",")[0]
                # gene_tag=allele.split("*")[0]
                gene_tag = items[0]
                gene_idx=int(items[1])-1
                print("result :", gene_tag, allele) 
                full_allele_dict[gene_tag][gene_idx].append(allele)

def extract_allele_fasta(alleles, gene_dir, allele_idx):
    # extract seq in db_ref to allele fasta
    for allele in alleles:
        fmt_allele=remove_characters(allele)
        allele_seq = db_ref_dict[fmt_allele]
        fa=f"{gene_dir}/{fmt_allele}.{allele_idx}.fasta"
        with open(fa, "w") as f:
            f.write(f">{fmt_allele}\n{allele_seq}\n")
    
def bwa_map_all_allele(alleles, gene_dir, ref, allele_idx):
    for allele in alleles:
        fmt_allele=remove_characters(allele)
        fq=f"{gene_dir}/{fmt_allele}.{allele_idx}.fasta"

        bam=f"{gene_dir}/{fmt_allele}.{allele_idx}.bam"
        # index ref
        cmd=f"bwa index {ref}"
        os.system(cmd)
        # perform bwa long reads mapping
        # cmd=f"bwa mem -x pacbio -t {threads} {ref} {fq} | samtools view -bS -F 0x804 -| samtools sort - >{bam}"
        # os.system(cmd)
        # perform minimap2 long reads mapping
        cmd=f"minimap2 -a -t {threads} {ref} {fq} | samtools view -bS -F 0x804 -| samtools sort - >{bam}"
        os.system(cmd)
        # index bam
        cmd=f"samtools index {bam}"
        os.system(cmd)

    

def allele_remap():
    for gene, alleles in step1_res_dict.items():
        gene_dir=f"{remap_allele_dir}/{gene}"
        if not os.path.exists(gene_dir):
            os.mkdir(gene_dir)
        if len(alleles)==0:
            continue
        else:
            if '-' == alleles[0]:
                continue
    # for step1 hom
        if alleles[0] == alleles[1]:
            # for step2 hom
            if len(step2_res_dict[gene]) == 0:
                continue
            else:
                if "-" == step2_res_dict[gene][0]:
                    continue
            if step2_res_dict[gene][0] == step2_res_dict[gene][1]:
                ref=f"{my_folder.sequence_dir}/{gene_class}.allele.1.{gene}.fasta"
                # ref=f"{outdir}/{sample}/{gene_class}.allele.1.{gene}.fasta"
                alleles=full_allele_dict[gene][0]
                # extract allele fasta to single file
                if os.path.exists(ref):
                    extract_allele_fasta(alleles, gene_dir, 0)
                    bwa_map_all_allele(alleles, gene_dir, ref, 0)               
            else:
                # for step2 het
                ref1=f"{my_folder.sequence_dir}/{gene_class}.allele.1.{gene}.fasta"
                ref2=f"{my_folder.sequence_dir}/{gene_class}.allele.2.{gene}.fasta"
                # ref1=f"{outdir}/{sample}/{gene_class}.allele.1.{gene}.fasta"
                # ref2=f"{outdir}/{sample}/{gene_class}.allele.2.{gene}.fasta"
                if os.path.exists(ref1):
                    alleles=full_allele_dict[gene][0]
                    extract_allele_fasta(alleles, gene_dir, 0)
                    bwa_map_all_allele(alleles, gene_dir, ref1, 0)
                if os.path.exists(ref2):
                    alleles=full_allele_dict[gene][1]
                    extract_allele_fasta(alleles, gene_dir, 1)
                    bwa_map_all_allele(alleles, gene_dir, ref2, 1)

        else:
            # for step1 het
            for allele_idx, allele in enumerate(alleles):
                # for step2 hom
                if len(step2_res_dict[gene]) == 0:
                    continue
                else:
                    if "-" == step2_res_dict[gene][0]:
                        continue
                if step2_res_dict[gene][0] == step2_res_dict[gene][1]:
                    ref=f"{my_folder.sequence_dir}/{gene_class}.allele.1.{gene}.fasta"
                    # ref=f"{outdir}/{sample}/{gene_class}.allele.1.{gene}.fasta"
                    if os.path.exists(ref):
                        alleles=full_allele_dict[gene][allele_idx]
                        extract_allele_fasta(alleles, gene_dir, allele_idx)
                        bwa_map_all_allele(alleles, gene_dir, ref, allele_idx)
                else:
                    # for step2 het
                    ref1=f"{my_folder.sequence_dir}/{gene_class}.allele.1.{gene}.fasta"
                    ref2=f"{my_folder.sequence_dir}/{gene_class}.allele.2.{gene}.fasta"
                    # ref1=f"{outdir}/{sample}/{gene_class}.allele.1.{gene}.fasta"
                    # ref2=f"{outdir}/{sample}/{gene_class}.allele.2.{gene}.fasta"
                    if os.path.exists(ref1):
                        alleles=full_allele_dict[gene][0]
                        extract_allele_fasta(alleles, gene_dir, 0)
                        bwa_map_all_allele(alleles, gene_dir, ref1, 0)
                    if os.path.exists(ref2):
                        alleles=full_allele_dict[gene][1]
                        extract_allele_fasta(alleles, gene_dir, 1)
                        bwa_map_all_allele(alleles, gene_dir, ref2, 1)
                    break


def main():
    # parse step 2 res
    read_hla_file(step1_result, step1_res_dict)
    read_hla_file(step2_result, step2_res_dict)
    parse_full_allele(step2_result, full_allele_dict)
    allele_remap()


if __name__ == "__main__":   
    if len(sys.argv) <2:
        print("Usage: python script.py <filename>")
        sys.exit(1)
    
    sample = sys.argv[1]
    gene_class = sys.argv[2]
    outdir = sys.argv[3]
    data_type = sys.argv[4]
    seq_tech = sys.argv[5]
    RNA_type = sys.argv[6]
    threads = sys.argv[7]
    db_ref = sys.argv[8]
    my_folder = My_folder({"o": outdir, "n":sample})
    read_type = Read_Type(seq_tech, data_type, RNA_type)

    step1_result = f"{my_folder.sample_prefix}.{gene_class}.type.result.txt"
    step2_result = f"{my_folder.sample_prefix}.{gene_class}.final.type.result.txt"
    
    minimap_para = read_type.get_minimap2_param()
    step1_res_dict = {}
    step2_res_dict = {}
    full_allele_dict = {}
    db_folder=os.path.dirname(db_ref)
    gene_list = get_folder_list(db_folder)
    for gene in gene_list:
        step1_res_dict[gene] = []
        step2_res_dict[gene] = []
        full_allele_dict[gene] = [[],[]]
    remap_allele_dir=my_folder.for_viz_dir + "/remap_allele"
    if not os.path.exists(remap_allele_dir):
        os.mkdir(remap_allele_dir)
    # open db ref here
    db_ref_fasta = SeqIO.parse(db_ref, "fasta")
    # to dict
    db_ref_dict = {}
    for record in db_ref_fasta:
        db_ref_dict[remove_characters(record.id)] = record.seq
    main()
        



