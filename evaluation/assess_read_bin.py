import gzip
from Bio import SeqIO
import sys


def extract_read_names(file_path):
    read_names = {}
    with open(file_path, 'r') as file:
        for line in file:
            read_name, gene_name = line.strip().split('\t')
            if gene_name not in read_names:
                read_names[gene_name] = []
            read_names[gene_name].append(read_name)
    return read_names


def extract_read_allele(hap_fasta):
    # print (hap_fasta)
    read_allele_dict = {}
    with open(hap_fasta, 'r') as file:
        index = 1
        for line in file:
            
            if line[0] != ">":
                continue
            
            gene_name = line[1:].strip().split("*")[0]
            index_name = f"S{index}"

            index += 1

            read_allele_dict[index_name] = gene_name
    # print (read_allele_dict)
    return read_allele_dict


def extract_read_names_from_fastq(gzipped_fastq_file, read_allele_dict):
    read_names = {}

    with gzip.open(gzipped_fastq_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            # read_names.append(record.id)

            field = record.id.split("_")
            index_name = field[0]
            gene_name = read_allele_dict[index_name]
            read_name = record.id
            # print (field)

            if gene_name not in read_names:
                read_names[gene_name] = []
            read_names[gene_name].append(read_name)
    # print (read_names)
    return read_names

def compare(true_read_names, infer_read_names):
    for gene in true_read_names:
        if gene not in infer_read_names:
            print ("error")
            sys.exit(0)
        true_reads = true_read_names[gene]
        infer_reads = infer_read_names[gene]

        recall, accuracy = calculate_recall_and_accuracy(true_reads, infer_reads)
        print (gene, recall, accuracy)


def calculate_recall_and_accuracy(true_list, inferred_list):
    # Calculate the number of true positives
    true_positives = len(set(true_list) & set(inferred_list))

    # Calculate the total number of positive samples
    total_positives = len(true_list)

    # Calculate the total number of samples
    total_samples = len(inferred_list)

    # Calculate recall
    recall = true_positives / total_positives

    # Calculate accuracy
    accuracy = true_positives / total_samples

    return recall, accuracy


def cal_bin_accuracy(hap_fasta, fastq, read_assign):
    read_allele_dict = extract_read_allele(hap_fasta)
    true_read_names = extract_read_names_from_fastq(fastq, read_allele_dict)

    infer_read_names = extract_read_names(read_assign)

    return compare(true_read_names, infer_read_names)

if __name__ == "__main__":  

    # read_assign = "/home/wangshuai/softwares/SpecLong/test/test/test_HLA/test_HLA.assign.txt"
    # hap_fasta = "/mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/test.HLA.sep.fa"
    # fastq = "/mnt/d/HLAPro_backup/Nanopore_optimize/data/sim_hap/test.fastq.gz"

    # print (read_names)

    # python3 ../evaluation/assess_read_bin.py $outdir/$sample/${sample}.assign.txt $outdir/$sample/$sample.HLA.sep.fa $outdir/$sample/${sample}.fastq.gz 

    read_assign = sys.argv[1]
    hap_fasta = sys.argv[2]
    fastq = sys.argv[3]

    cal_bin_accuracy(hap_fasta, fastq, read_assign)
