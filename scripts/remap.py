
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
from determine_gene import get_focus_gene, get_folder_list
from alignment_modules import Read_Type


def call_variants(bam, ref, ovcf):
            # longshot -F -c 2 -C 100000 -P {args["strand_bias_pvalue_cutoff"]} -r {interval_dict[gene]} --bam {bam} --ref {hla_ref} --out {parameter.outdir}/{parameter.sample}.{gene}.longshot.vcf 
    cmd=f"""
    longshot -F -c 2 -C 100000 -P 0.01 --bam {bam} --ref {ref} --out {ovcf} 
    """
    os.system(cmd)

def read_hla_file(filename, some_dict):
    with open(filename, 'r') as file:
        next(file)
        next(file)
        for idx, row in enumerate(file):
            print(row.strip())
            items=row.strip().split("\t")
            if len(items)<3:
                # remove gene_tag from some_dict
                gene_tag = items[0]
                some_dict[gene_tag].append("-")
                continue
            alleles_str=items[2].split(";")
            print("alleles_str", alleles_str)
            for iti, it in enumerate(alleles_str):
                if iti>0:
                    continue
                allele=it.split(",")[0]
                # gene_tag=allele.split("*")[0]
                gene_tag = items[0]
                print("result :", gene_tag, allele) 
                some_dict[gene_tag].append(allele)
def remap():
    for gene, alleles in step1_res_dict.items():
        if len(alleles)==0:
            continue
        else:
            if '-' == alleles[0]:
                continue
    # for step1 hom
        if alleles[0] == alleles[1]:
            fq=f"{outdir}/{sample}/{gene}.long_read.fq.gz"
            # for step2 hom
            if len(step2_res_dict[gene]) == 0:
                continue
            else:
                if "-" == step2_res_dict[gene][0]:
                    continue
            if step2_res_dict[gene][0] == step2_res_dict[gene][1]:
                ref=f"{outdir}/{sample}/{gene_class}.allele.1.{gene}.fasta"
                bam=f"{outdir}/{sample}/{gene}.remap2res.bam"
                cmd=f"""
                    minimap2 -t {threads} -a {ref} {minimap_para} {fq} | samtools view -bS -F 0x804 -| samtools sort - >{bam}
                    samtools index {bam}
                """
                # check ref exists
                if os.path.exists(ref):
                    os.system(cmd)
                    call_variants(bam, ref, f"{outdir}/{sample}/{gene}.remap2res.vcf")
                
            else:
                # for step2 het
                ref1=f"{outdir}/{sample}/{gene_class}.allele.1.{gene}.fasta"
                ref2=f"{outdir}/{sample}/{gene_class}.allele.2.{gene}.fasta"
                bam1=f"{outdir}/{sample}/{gene}.remap2res1.bam"
                bam2=f"{outdir}/{sample}/{gene}.remap2res2.bam"
                
                cmd1=f"""
                    minimap2 -t {threads} -a {ref1} {minimap_para} {fq} | samtools view -bS -F 0x804 -| samtools sort - >{bam1}
                    samtools index {bam1}
                """
                cmd2=f"""
                    minimap2 -t {threads} -a {ref2} {minimap_para} {fq} | samtools view -bS -F 0x804 -| samtools sort - >{bam2}
                    samtools index {bam2}
                """
                # check ref1 exists
                if os.path.exists(ref1):
                    os.system(cmd1)
                    call_variants(bam1, ref1, f"{outdir}/{sample}/{gene}.remap2res1.vcf")
                # check ref2 exists
                if os.path.exists(ref2):
                    os.system(cmd2)
                    call_variants(bam2, ref2, f"{outdir}/{sample}/{gene}.remap2res2.vcf")
        else:
            # for step1 het
            for allele_idx, allele in enumerate(alleles):
                fq=f"{outdir}/{sample}/{allele}.fq.gz"
                # for step2 hom
                if len(step2_res_dict[gene]) == 0:
                    continue
                else:
                    if "-" == step2_res_dict[gene][0]:
                        continue
                if step2_res_dict[gene][0] == step2_res_dict[gene][1]:
                    ref=f"{outdir}/{sample}/{gene_class}.allele.1.{gene}.fasta"
                    bam=f"{outdir}/{sample}/{gene}.remap2res{allele_idx+1}.bam"
                    cmd=f"""
                        minimap2 -t {threads} -a {ref} {minimap_para} {fq} | samtools view -bS -F 0x804 -| samtools sort - >{bam}
                        samtools index {bam}
                    """
                    # check ref exists
                    if os.path.exists(ref):
                        os.system(cmd)
                        call_variants(bam, ref, f"{outdir}/{sample}/{gene}.remap2res{allele_idx+1}.vcf")
                else:
                    # for step2 het
                    ref1=f"{outdir}/{sample}/{gene_class}.allele.1.{gene}.fasta"
                    ref2=f"{outdir}/{sample}/{gene_class}.allele.2.{gene}.fasta"
                    bam1=f"{outdir}/{sample}/{gene}.remap2res1.bam"
                    bam2=f"{outdir}/{sample}/{gene}.remap2res2.bam"
                    cmd1=f"""
                        minimap2 -t {threads} -a {ref1} {minimap_para} {fq} | samtools view -bS -F 0x804 -| samtools sort - >{bam1}
                        samtools index {bam1}
                    """
                    cmd2=f"""
                        minimap2 -t {threads} -a {ref2} {minimap_para} {fq} | samtools view -bS -F 0x804 -| samtools sort - >{bam2}
                        samtools index {bam2}
                    """
                    # check ref1 exists
                    if os.path.exists(ref1):
                        os.system(cmd1)
                        call_variants(bam1, ref1, f"{outdir}/{sample}/{gene}.remap2res1.vcf")
                    # check ref2 exists
                        os.system(cmd2)
                        call_variants(bam2, ref2, f"{outdir}/{sample}/{gene}.remap2res2.vcf")
                    break




def main():
    # parse step one res
    read_hla_file(step1_result, step1_res_dict)
    read_hla_file(step2_result, step2_res_dict)
    remap()



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
    step1_result = f"{outdir}/{sample}/{sample}.{gene_class}.type.result.txt"
    step2_result = f"{outdir}/{sample}/hlala.like.results.txt"

    read_type = Read_Type(seq_tech, data_type, RNA_type)
    minimap_para = read_type.get_minimap2_param()
    step1_res_dict = {}
    step2_res_dict = {}
    db_folder=os.path.dirname(db_ref)
    gene_list = get_folder_list(db_folder)
    for gene in gene_list:
        step1_res_dict[gene] = []
        step2_res_dict[gene] = []
    main()
        



