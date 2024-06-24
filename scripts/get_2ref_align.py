import csv
import sys
import os
from Bio import SeqIO
from determine_gene import get_focus_gene
from alignment_modules import Read_Type

def read_hla_file(filename):
    with open(filename, 'r') as file:
        next(file)
        next(file)
        for idx, row in enumerate(file):
            print(row.strip())
            items=row.strip().split("\t")
            alleles_str=items[2].split(";")
            print("alleles_str", alleles_str)
            for iti, it in enumerate(alleles_str):
                if iti>0:
                    continue
                allele=it.split(",")[0]
                # gene_tag=allele.split("*")[0]
                gene_tag = items[0]
                print("result :", gene_tag, allele) 
                gene_ref_dict[gene_tag].append(allele)

            
def replace_single_contig_name(input_file, output_file, new_contig_name):
    record = next(SeqIO.parse(input_file, "fasta"))
    
    record.id = new_contig_name
    record.description = new_contig_name
    
    SeqIO.write([record], output_file, "fasta")
    print("Contig name replaced successfully.")

def build_HLA_ref():
    # Write the reference file
    # shell code
    for gene, alleles in gene_ref_dict.items():
        if len(alleles)==0:
            continue
        else:
            if '-' == alleles[0]:
                continue
        allele_dir=f"{db_build_dir}/{gene}"
        if not os.path.exists(allele_dir):
            os.makedirs(allele_dir)

        # for hom
        if alleles[0] == alleles[1]:
            allele=alleles[0]
            print(f"processing {gene} {allele} homo")
            cmd=f"""
                samtools faidx {db_ref} {allele}>{allele_dir}/{gene}.raw.fasta
            """
            os.system(cmd)
            replace_single_contig_name(f"{allele_dir}/{gene}.raw.fasta", f"{allele_dir}/{gene}.fasta", f"{gene}")
            index_cmd=f"""
                samtools faidx {allele_dir}/{gene}.fasta
                bwa index {allele_dir}/{gene}.fasta
                makeblastdb -in {allele_dir}/{gene}.fasta -dbtype nucl -parse_seqids -out {allele_dir}/{gene}
            """
            os.system(index_cmd)
        else:
            # for het
            for allele_idx, allele in enumerate(alleles):
                allele_idx+=1
                print(f"processing {gene} {allele} {allele_idx}")
                cmd=f"""
                    samtools faidx {db_ref} {allele}>{allele_dir}/{gene}.raw.{allele_idx}.fasta
                """
                os.system(cmd)
                replace_single_contig_name(f"{allele_dir}/{gene}.raw.{allele_idx}.fasta", f"{allele_dir}/{gene}.{allele_idx}.fasta", f"{gene}_ref{allele_idx}")
                index_cmd=f"""
                    samtools faidx {allele_dir}/{gene}.{allele_idx}.fasta
                    bwa index {allele_dir}/{gene}.{allele_idx}.fasta
                    makeblastdb -in {allele_dir}/{gene}.{allele_idx}.fasta -dbtype nucl -parse_seqids -out {allele_dir}/{gene}.{allele_idx}
                """
                os.system(index_cmd)

def map_phased_reads_2_ref():
    for gene, alleles in gene_ref_dict.items():
        if len(alleles)==0:
            continue
        else:
            if '-' == alleles[0]:
                continue
        print (f"processing alignment for {gene}")
        # for hom
        if alleles[0] == alleles[1]:
            fq=f"{outdir}/{sample}/{gene}.long_read.fq.gz"
            ref=f"{db_build_dir}/{gene}/{gene}.fasta"
            bam=f"{outdir}/{sample}/{gene}.bam"
            depth_file=f"{outdir}/{sample}/{gene}.depth"
            # minimap
            # minimap2 -t %s %s -a $hla_ref $outdir/$hla.fq.gz | samtools view -bS -F 0x800 -| samtools sort - >$outdir/$hla.bam
            cmd=f"""
                minimap2 -t {threads} -a {ref} {minimap_para} {fq} | samtools view -bS -F 0x800 -| samtools sort - >{bam}
                samtools index {bam}
                samtools depth -d 1000000 -aa {bam} > {depth_file}
            """
            os.system(cmd)
        else:
            # for het
            for allele_idx, allele in enumerate(alleles):
                fq=f"{outdir}/{sample}/{allele}.fq.gz"
                ref=f"{db_build_dir}/{gene}/{gene}.{allele_idx+1}.fasta"
                bam=f"{outdir}/{sample}/{gene}.{allele_idx}.bam"
                depth_file=f"{outdir}/{sample}/{gene}.{allele_idx}.depth"
                # minimap
                # minimap2 -t %s %s -a $hla_ref $outdir/$hla.%s.fq.gz | samtools view -bS -F 0x800 -| samtools sort - >$outdir/$hla.bam
                cmd=f"""
                    minimap2 -t {threads} -a {ref} {minimap_para} {fq} | samtools view -bS -F 0x800 -| samtools sort - >{bam}
                    samtools index {bam} 
                    samtools depth -d 1000000 -aa {bam} > {depth_file}
                """
                os.system(cmd)
            

def main():
    read_hla_file(ref_file)
    print(gene_ref_dict)
    build_HLA_ref()
    map_phased_reads_2_ref()


if __name__ == "__main__":
    if len(sys.argv) <2:
        print("Usage: python script.py <filename>")
        sys.exit(1)

    sample=sys.argv[1]
    db_ref=sys.argv[2]
    db_build_dir=sys.argv[3]
    outdir=sys.argv[4]
    data_type=sys.argv[5]
    threads=sys.argv[6]
    gene_class = sys.argv[7]
    seq_tech = sys.argv[8]
    RNA_type = sys.argv[9]
    ref_file = f"{outdir}/{sample}/{sample}.{gene_class}.type.result.txt"

    read_type = Read_Type(seq_tech, data_type, RNA_type)
    minimap_para = read_type.get_minimap2_param()

    if not os.path.exists(db_build_dir):
        os.makedirs(db_build_dir)
    # for HLA allele
    gene_list, interval_dict =  get_focus_gene(gene_class)
    gene_ref_dict={}
    for gene in gene_list:
        gene_ref_dict[gene]=[]
    main()