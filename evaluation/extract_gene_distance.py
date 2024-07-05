import re
import sys
import pickle
from collections import defaultdict

sys.path.insert(0, sys.path[0]+'/../scripts/')
from determine_gene import get_focus_gene
from read_objects import get_interval_distance


## given a gene GTF file, extract the chr and interval of each gene, then given a gene list, return the distance matrix of these genes
def extract_gene_distance(gtf_file):
    gene_dict = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip().split("\t")
            # if line[2] == "transcript":
            if line[2] == "gene":
                chr = line[0]
                start = int(line[3])
                end = int(line[4])
                # print (line[-1].split(";"))
                gene_name = re.search('gene_name "(.*?)"', line[-1]).group(1)
                # print (gene_name, chr, start, end)
                gene_dict[gene_name] = (chr, start, end)
    return gene_dict

## given a gene list, return the distance matrix of these genes
def cal_distance(gene_list, gene_dict):
    distance_matrix = {}
    for i in range(len(gene_list)):
        for j in range(i+1, len(gene_list)):
            gene1 = gene_list[i]
            gene2 = gene_list[j]

            if gene1 not in gene_dict:
                
                print (gene1, "not in the gene dict")
                distance = float("inf")
            elif gene2 not in gene_dict:
                print (gene2, "not in the gene dict")
                distance = float("inf")
            else:
                chr1, start1, end1 = gene_dict[gene1]
                chr2, start2, end2 = gene_dict[gene2]
                if chr1 != chr2:
                    distance = float("inf")
                else:
                    distance = get_interval_distance([start1, end1], [start2, end2])

            distance_matrix[(gene1, gene2)] = distance
            distance_matrix[(gene2, gene1)] = distance
    return distance_matrix

def get_VDJ_names():
    fai = "../db/IG_TR/IG_TR.fasta.fai"
    VDJ_gene_list = []
    with open(fai, 'r') as f:
        for line in f:
            line = line.strip().split("\t")
            gene = line[0].split("*")[0]
            if gene not in VDJ_gene_list:
                VDJ_gene_list.append(gene)
    return VDJ_gene_list

def get_VDJ_bed():
    gene_dict = extract_gene_distance(gtf_file)
    VDJ_gene_list = get_VDJ_names()

    gene_list = []
    f = open("../gene_dist/" + 'IG_TR.gene.bed', 'w')
    interval_dict = defaultdict(dict)
    for gene in VDJ_gene_list:
        if gene in gene_dict:
            chr, start, end = gene_dict[gene]
            interval_dict[chr][gene] = (start, end)
        # else:
            # print(gene, "-", "-", "-") 
    order_gene_list = sort_gene(interval_dict) 
    for gene in order_gene_list:
        chr, start, end = gene_dict[gene]
        print(gene, chr, start, end, file = f)  
        gene_list.append(gene)
    print (len(gene_list))
    print (gene_list)

def sort_gene(interval_dict):
    order_gene_list = []
    for genome in interval_dict:
        # sort the gene on it by their start position
        interval_dict[genome] = dict(sorted(interval_dict[genome].items(), key=lambda item: item[1][0]))
        for gene in interval_dict[genome]:
            order_gene_list.append(gene)
    return order_gene_list

if __name__ == "__main__":
    #### https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz
    # gtf_file = "/mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38.ncbiRefSeq.gtf"

    ### https://www.gencodegenes.org/human/
    gtf_file = "/mnt/d/HLAPro_backup/Nanopore_optimize/data/gencode.v46.annotation.gtf"

    # args = {}
    # # args["i"] = "HLA"
    # # args["i"] = "CYP"
    # args["i"] = "KIR"
    # gene_list, interval_dict =  get_focus_gene(args)
    # gene_dict = extract_gene_distance(gtf_file)
    # distance_matrix = cal_distance(gene_list, gene_dict)
    # # store the distance matrix with pickle

    # with open("../gene_dist/" + args["i"] + '.distance_matrix.pkl', 'wb') as f:
    #     pickle.dump(distance_matrix, f)
    # # print (distance_matrix)
    # with open("../gene_dist/" + args["i"] + '.gene.bed', 'w') as f:
    #     for gene in gene_list:
    #         if gene in gene_dict:
    #             chr, start, end = gene_dict[gene]
    #             print(gene, chr, start, end, file = f)
    #         else:
    #             print(gene, "-", "-", "-", file=f)

    get_VDJ_bed()
