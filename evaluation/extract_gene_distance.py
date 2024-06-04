import re
import sys
import pickle

sys.path.insert(0, sys.path[0]+'/../scripts/')
from determine_gene import get_focus_gene


## given a gene GTF file, extract the chr and interval of each gene, then given a gene list, return the distance matrix of these genes
def extract_gene_distance(gtf_file):
    gene_dict = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip().split("\t")
            if line[2] == "transcript":
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

            if gene1 not in gene_dict or gene2 not in gene_dict:
                print (gene1, gene2, "not in the gene dict")
                distance = float("inf")
            else:
                chr1, start1, end1 = gene_dict[gene1]
                chr2, start2, end2 = gene_dict[gene2]
                if chr1 != chr2:
                    distance = float("inf")
                else:
                    distance = abs(start1 - start2)

            distance_matrix[(gene1, gene2)] = distance
            distance_matrix[(gene2, gene1)] = distance
    return distance_matrix

gtf_file = "/mnt/d/HLAPro_backup/Nanopore_optimize/data/hg38.ncbiRefSeq.gtf"

args = {}
args["i"] = "HLA"
# args["i"] = "CYP"
gene_list, interval_dict =  get_focus_gene(args)
gene_dict = extract_gene_distance(gtf_file)
distance_matrix = cal_distance(gene_list, gene_dict)
# store the distance matrix with pickle

with open(args["i"] + '.distance_matrix.pkl', 'wb') as f:
    pickle.dump(distance_matrix, f)
print (distance_matrix)
