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

analyze_tree()