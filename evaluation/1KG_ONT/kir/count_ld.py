# use pandas to load a csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys

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
        return 'VDJ'
    return 'other'

# load the csv
df = pd.read_csv("hla_kir_cyp_vdj_LD_values.csv")
# enumerate the df
data = []
depth_dict = {}
for index, row in df.iterrows():
    gene1_t  = classify_genes(row['gene1'])
    gene2_t  = classify_genes(row['gene2'])
    if gene1_t == 'other' or gene2_t == 'other':
        continue
    if gene1_t == 'CYP' or gene2_t == 'CYP':
        continue
    if gene1_t == gene2_t:
        data.append(["intra", row['Wn'], gene1_t])
    else:
        data.append(["inter", row['Wn'], gene1_t ])
        data.append(["inter", row['Wn'], gene2_t ])
## data to df
df = pd.DataFrame(data, columns = ['type', 'Wn', 'Class'])
## save the df
df.to_csv("intra_inter_LD.csv", index=False)
        
