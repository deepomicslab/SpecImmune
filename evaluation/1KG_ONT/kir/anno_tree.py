# use pandas to load a csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
from ete3 import Tree

sys.path.insert(0, sys.path[0]+'/../../../scripts/')

from get_lite_db import convert_field_for_allele
from determine_gene import get_focus_gene
from db_objects import My_db

HLA_gene_list, interval_dict =  get_focus_gene('HLA')
KIR_gene_list, interval_dict =  get_focus_gene('KIR')
CYP_gene_list, interval_dict =  get_focus_gene('CYP')
VDJ_gene_list, interval_dict =  get_focus_gene('IG_TR')

class1_genes = ['A', 'B', 'C', 'E', 'F', 'G']
class1_pseudogenes = ['H', 'J', 'K', 'L', 'N', 'P', 'S', 'T', 'U', 'V', 'W', 'Y']
class2_genes = ['DRA', 'DQA1', 'DQA2', 'DQB1', 'DQB2', 'DPA1', 'DPA2', 'DPB1', 'DPB2', 'DMA', 'DMB', 'DOA', 'DOB', 
                'DRB1', 'DRB3', 'DRB4', 'DRB5']
non_hla_genes = ['MICA', 'MICB', 'TAP1', 'TAP2']

def group_genes(locus):
    if len(locus.split("-")) == 2:
        gene = locus.split("-")[1]
    else:
        gene = locus
    if gene in class1_genes:
        return "Class I"
    if gene in class1_pseudogenes:
        return "Pseudogenes"
    if gene in class2_genes:
        return "Class II"
    if gene in non_hla_genes:
        return "Non-HLA"
    return "other"

## classify genes
def classify_genes(gene):
    if gene in HLA_gene_list:
        return 'HLA'
    if gene in KIR_gene_list:
        return 'KIR'
    if gene in CYP_gene_list:
        return 'CYP'
    if gene in VDJ_gene_list:
        if gene[:2] == 'TR':
            return "TCR"
        elif gene[:2] == 'IG':
            return "IG"
    return 'other'

def classify_genes2(gene):
    if gene in HLA_gene_list:
        return 'HLA'
    if gene in KIR_gene_list:
        return 'KIR'
    if gene in CYP_gene_list:
        return 'CYP'
    if gene in VDJ_gene_list:
        return "IG/TCR"
    return 'other'

def define_tree_color():
    # load the csv
    df = pd.read_csv("hla_kir_cyp_vdj_LD_values.csv")
    # enumerate the df
    data = []
    gene_dict = {}
    for index, row in df.iterrows():
        gene_dict[row['gene1']] = 1
        gene_dict[row['gene2']] = 1

    colors = {"TCR":"#d9e6eb", "IG":"#9fc3d5", "KIR":"#8f96bd", "CYP":"#2a347a", "HLA":"#d6d69b"}


    f = open("tree/gene_node_color.txt", "w")
    print ("TREE_COLORS\nSEPARATOR SPACE\nDATA\n", file=f)
    for gene in gene_dict:
        type = classify_genes(gene)
        color = colors[type]
        print (f"{gene} range {color} {type}", file=f)
            
## load a tree and compute the Page Rank with networkx
def analyze_tree():
    import networkx as nx
    import matplotlib.pyplot as plt
    from networkx.drawing.nx_agraph import graphviz_layout

    df = pd.read_csv("hla_kir_cyp_vdj_LD_values.csv")
    # enumerate the df
    data = []
    depth_dict = {}
    G = nx.DiGraph()
    for index, row in df.iterrows():
        G.add_edge(row['gene1'], row['gene2'], weight=row['Wab'])

    # compute the Page Rank
    # pr = nx.pagerank(G, alpha=0.9)
    ## compute Centrality 
    centrality = nx.eigenvector_centrality_numpy(G)
    # sort the node by centrality
    sorted_pr = sorted(centrality.items(), key=lambda x: x[1], reverse=True)
    print (sorted_pr)

    ## use Louvain Method to find the community
    import community as community_louvain
    partition = community_louvain.best_partition(G)
    print (partition)
    # save the Page Rank
    ## sort the node by PR
    # sorted_pr = sorted(pr.items(), key=lambda x: x[1], reverse=True)
    f = open("tree/node_PR.txt", "w")
    print ("DATASET_SIMPLEBAR\nSEPARATOR COMMA\nDATASET_LABEL,my_data\nCOLOR,#45818e\nDATA", file=f)
    for node, pr in sorted_pr:
        print (f"{node},{pr},label1", file=f)
    f.close()

def compare():
    df = pd.read_csv("hla_kir_cyp_vdj_LD_values.csv")
    # enumerate the df
    data = []
    depth_dict = {}
    for index, row in df.iterrows():
        gene1_t  = classify_genes(row['gene1'])
        gene2_t  = classify_genes(row['gene2'])
        if gene1_t == "KIR" and gene2_t == "HLA":
            HLA_class = group_genes(row['gene2'])

            data.append([HLA_class, row['Wab'], "Wab"])
        elif gene2_t == "KIR" and gene1_t == "HLA":
            HLA_class = group_genes(row['gene1'])

            data.append([HLA_class, row['Wab'], "Wba"])
        elif gene1_t == "KIR" and gene2_t not in  ["HLA", "KIR"]:
            data.append([gene2_t, row['Wab'], "Wab"])
        elif gene2_t == "KIR" and gene1_t not in ["HLA", "KIR"]:
            data.append([gene1_t, row['Wab'], "Wba"])

    df = pd.DataFrame(data, columns = ['Class', 'ALD', 'type'])
    df.to_csv("tree/hla_class_kir.csv", index=False)

def compare2():
    df = pd.read_csv("hla_kir_cyp_LD_values.csv")
    # enumerate the df
    data = []
    count_dict = {}
    for index, row in df.iterrows():
        

        gene1_t  = classify_genes(row['gene1'])
        gene2_t  = classify_genes(row['gene2'])

        if gene1_t not in count_dict:
            count_dict[gene1_t] = set()
        count_dict[gene1_t].add(row['gene1'])
        if gene2_t not in count_dict:
            count_dict[gene2_t] = set()
        count_dict[gene2_t].add(row['gene2'])

        if row['gene1'] in ['HLA-DRB1','HLA-DPB1', 'HLA-DQB1'] and gene2_t == "KIR":
            data.append([row['gene2'], row['gene1'], row['min_w'], "Wab", 'class II'])
        if row['gene1'] in ['HLA-A','HLA-B', 'HLA-C'] and gene2_t == "KIR":
            data.append([row['gene2'], row['gene1'], row['min_w'], "Wab", 'class I'])

        # elif gene1_t == "KIR" and gene2_t not in  ["HLA", "KIR"]:
        #     data.append([gene2_t, row['min_w'], "Wab"])
    for gene_class in count_dict:
        print (gene_class, len(count_dict[gene_class]))
            
    df = pd.DataFrame(data, columns = ['KIR', 'gene', 'ALD', 'type', 'class'])
    df.to_csv("tree/hla_class_kir.csv", index=False)

def compare3():
    df = pd.read_csv("hla_kir_cyp_LD_values.csv")
    # enumerate the df
    data = []
    depth_dict = {}
    for index, row in df.iterrows():
        

        gene1_t  = classify_genes(row['gene1'])
        gene2_t  = classify_genes(row['gene2'])

        if gene1_t == "CYP" and gene2_t == "KIR":
            data.append([row['gene2'], row['gene1'], row['min_w'], "CYP2D6 & KIR"])

        if gene1_t == "CYP" and gene2_t == "HLA":
            data.append([row['gene2'], row['gene1'], row['min_w'], "CYP2D6 & HLA"])

        # elif gene1_t == "KIR" and gene2_t not in  ["HLA", "KIR"]:
        #     data.append([gene2_t, row['min_w'], "Wab"])

            
    df = pd.DataFrame(data, columns = ['CYP', 'gene', 'ALD', 'type'])
    df.to_csv("tree/cyp_hla_kir.csv", index=False)

def compare4():
    df = pd.read_csv("hla_kir_cyp_vdj_LD_values.csv")
    # enumerate the df
    data = []
    count_dict = {}
    for index, row in df.iterrows():
        
        gene1_t  = classify_genes2(row['gene1'])
        gene2_t  = classify_genes2(row['gene2'])

        if row['gene1'] == row['gene2']:
            continue

        if gene1_t not in count_dict:
            count_dict[gene1_t] = set()
        count_dict[gene1_t].add(row['gene1'])
        if gene2_t not in count_dict:
            count_dict[gene2_t] = set()
        count_dict[gene2_t].add(row['gene2'])

        # if gene1_t == "IG/TCR" and gene2_t != "IG/TCR":
        #     data.append([row['gene1'], row['gene2'], row['min_w'], gene2_t])
        if gene1_t == "IG/TCR" :
            data.append([row['gene1'], row['gene2'], row['min_w'], gene2_t])

            
    df = pd.DataFrame(data, columns = ['vdj', 'gene', 'ALD', 'type'])
    df.to_csv("tree/vdj_cyp_hla_kir.csv", index=False)
    for gene_class in count_dict:
        print (gene_class, len(count_dict[gene_class]))


def compare5():
    df = pd.read_csv("hla_kir_cyp_vdj_LD_values.csv")
    # enumerate the df
    data = []
    count_dict = {}
    save_past = {}
    for index, row in df.iterrows():
        
        gene1_t  = classify_genes(row['gene1'])
        gene2_t  = classify_genes(row['gene2'])

        if row['gene1'] == row['gene2']:
            continue
        gene_pair = "xx".join(sorted([row['gene1'], row['gene2']]))

        if gene_pair in save_past:
            continue
        
        save_past[gene_pair] = 1
        # if gene1_t == "IG/TCR" and gene2_t != "IG/TCR":
        #     data.append([row['gene1'], row['gene2'], row['min_w'], gene2_t])
        # if gene1_t == "IG/TCR" :
        #     data.append([row['gene1'], row['gene2'], row['min_w'], gene2_t])
        
        if gene1_t == gene2_t:
            compare_type = gene1_t
        else:
            sort_type = sorted([gene1_t, gene2_t])
            compare_type = "(" + sort_type[0] + "," + sort_type[1] + ")"

        data.append([row['gene1'], row['gene2'], row['min_w'], compare_type])
        

            
    df = pd.DataFrame(data, columns = ['vdj', 'gene', 'ALD', 'type'])
    df.to_csv("tree/vdj_cyp_hla_kir.csv", index=False)
    for gene_class in count_dict:
        print (gene_class, len(count_dict[gene_class]))


def just_count():
    df = pd.read_csv("hla_kir_cyp_vdj_LD_values.csv")
    # enumerate the df
    data = []
    count_dict = {}
    for index, row in df.iterrows():
        
        gene1_t  = classify_genes(row['gene1'])
        gene2_t  = classify_genes(row['gene2'])

        if gene1_t not in count_dict:
            count_dict[gene1_t] = set()
        count_dict[gene1_t].add(row['gene1'])
        if gene2_t not in count_dict:
            count_dict[gene2_t] = set()
        count_dict[gene2_t].add(row['gene2'])
     

    for gene_class in count_dict:
        print (gene_class, len(count_dict[gene_class]))

# analyze_tree()
compare2()
compare3()
compare4()
compare5()
just_count()