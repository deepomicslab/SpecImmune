import pandas as pd
from collections import defaultdict
from scipy.stats import fisher_exact
from sklearn.metrics import mutual_info_score
import numpy as np
import re
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
## use seaborn for better aesthetics
import seaborn as sns
from scipy.stats import pearsonr


def read_tab(csv, sample_pop_dict, family="HLA"):
    sample_set = set()
    allele_sample_dict = defaultdict(dict)
    
    df = pd.read_csv(csv)
    for index, row in df.iterrows():
        sample = row['Sample']
        ## if row["One_guess"] is empty or not a string, skip
        if not isinstance(row["One_guess"], str) or row["One_guess"].strip()== "":
            continue
        if sample not in sample_pop_dict:
            continue
        if row["Dataset"] != "1KGP":
            continue
        pop = sample_pop_dict[sample]
        
        gene = row["One_guess"].split("*")[0]
        ## get 1-field HLA
        if family == "HLA":
            field = row["One_guess"].split(":")
            # if field[0] != gene + "*" + "01":
            #     field[0] = gene + "*" + "other"

            # if len(field) < 2:
            #     allele = row["One_guess"]
            # else:
            #     allele = field[0] + ":" + field[1]
            allele = field[0]
            if gene not in ["HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1"]:
                continue
        if family == "KIR":
            field = row["One_guess"].split("*")
            allele = field[0] + "*" + field[1][:3]
        if allele.split("*")[0] in ["MICA", "MICB", "TAP2", "TAP1"]:
            allele = "HLA-" + allele
        if pop not in allele_sample_dict[allele]:
            allele_sample_dict[allele][pop] = 0
        # print (pop)
        allele_sample_dict[allele][pop] += 1
        sample_set.add(sample)
        # if re.search("HLA-A\*01", allele):
        #     print (allele_sample_dict["HLA-A*01"])
    
    pop_num_dict = defaultdict(int)
    for sample in sample_set:
        pop = sample_pop_dict[sample]
        pop_num_dict[pop] += 1
    pop_list = sorted(list(pop_num_dict.keys()))
    print (pop_list)
    ## get frequency of each allele in each population
    allele_pop_freq_dict = defaultdict(dict)
    for allele in allele_sample_dict:
        allele_pop_freq_dict[allele] = {}
        for pop in pop_list:
            if pop not in allele_sample_dict[allele]:
                allele_pop_freq_dict[allele][pop] = 0
            else:
                allele_pop_freq_dict[allele][pop] = allele_sample_dict[allele][pop] / (pop_num_dict[pop]*2)

    return allele_pop_freq_dict, pop_list

def read_cyp_csv(csv, sample_pop_dict):
    sample_set = set()
    allele_sample_dict = defaultdict(dict)
    
    df = pd.read_csv(csv)
    for index, row in df.iterrows():
        sample = row['Sample']
        ## if row["One_guess"] is empty or not a string, skip
        if not isinstance(row["One_guess"], str) or row["One_guess"].strip()== "":
            continue
        if sample not in sample_pop_dict:
            continue
        if row["Dataset"] != "1KGP":
            continue
        pop = sample_pop_dict[sample]
        
        gene = row["Locus"]
        allele = row["One_guess"].split(".")[0]

        if pop not in allele_sample_dict[allele]:
            allele_sample_dict[allele][pop] = 0
        allele_sample_dict[allele][pop] += 1
        sample_set.add(sample)
    
    pop_num_dict = defaultdict(int)
    for sample in sample_set:
        pop = sample_pop_dict[sample]
        pop_num_dict[pop] += 1
    pop_list = sorted(list(pop_num_dict.keys()))
    print (pop_list)
    ## get frequency of each allele in each population
    allele_pop_freq_dict = defaultdict(dict)
    for allele in allele_sample_dict:
        allele_pop_freq_dict[allele] = {}
        for pop in pop_list:
            if pop not in allele_sample_dict[allele]:
                allele_pop_freq_dict[allele][pop] = 0
            else:
                allele_pop_freq_dict[allele][pop] = allele_sample_dict[allele][pop] / (pop_num_dict[pop]*2)

    return allele_pop_freq_dict, pop_list

def read_ig_tcr_csv(csv, sample_pop_dict, family="TR"):
    sample_set = set()
    allele_sample_dict = defaultdict(dict)
    
    df = pd.read_csv(csv)
    for index, row in df.iterrows():
        sample = row['Sample']
        if sample not in sample_pop_dict:
            continue
        if row["Dataset"] != "1KG":
            continue
        pop = sample_pop_dict[sample]
        
        gene = row["gene"]
        if gene[:2] != family:
            continue
        allele_list = [row["allele_1"], row["allele_2"]]
        for allele in allele_list:
            if pop not in allele_sample_dict[allele]:
                allele_sample_dict[allele][pop] = 0
            allele_sample_dict[allele][pop] += 1
        sample_set.add(sample)
    
    pop_num_dict = defaultdict(int)
    for sample in sample_set:
        pop = sample_pop_dict[sample]
        pop_num_dict[pop] += 1
    pop_list = sorted(list(pop_num_dict.keys()))
    print (pop_list)
    ## get frequency of each allele in each population
    allele_pop_freq_dict = defaultdict(dict)
    for allele in allele_sample_dict:
        allele_pop_freq_dict[allele] = {}
        for pop in pop_list:
            if pop not in allele_sample_dict[allele]:
                allele_pop_freq_dict[allele][pop] = 0
            else:
                allele_pop_freq_dict[allele][pop] = allele_sample_dict[allele][pop] / (pop_num_dict[pop]*2)

    return allele_pop_freq_dict, pop_list



## read the file using pandas, sep is \t
def get_super_pop(super_pop_file):
    df = pd.read_csv(super_pop_file, sep='\t')
    super_pop_dict = {}
    for index, row in df.iterrows():
        if  pd.isna(row['Population Code']) or pd.isna(row['Super Population']):
            continue
        super_pop_dict[row['Population Code']] = row['Super Population']
    # print (super_pop_dict)
    return super_pop_dict

def get_sample_pop(sample_pop_file, super_pop_dict):
    pop_set = set()
    super_pop_set = set()
    # read excel file into df
    df = pd.read_excel(sample_pop_file)
    sample_pop_dict = {}
    sample_super_pop_dict = {}
    for index, row in df.iterrows():
        sample_pop_dict[row['Sample']] = row['Population']
        super_pop = super_pop_dict.get(row['Population'], None)
        pop_set.add(row['Population'])
        super_pop_set.add(super_pop)
        sample_super_pop_dict[row['Sample']] = super_pop
        # print (f"{row['Sample']} -> {row['Population']}")
    # print (sample_pop_dict)
    return sample_pop_dict, pop_set, sample_super_pop_dict, super_pop_set

def plot(final_freq_dict, allele_1, allele_2, pop_list, super_pop_dict):
    data2 = []
    for i in range(len(final_freq_dict[allele_1])):
        data2.append((final_freq_dict[allele_1][i], final_freq_dict[allele_2][i], pop_list[i], super_pop_dict[pop_list[i]]))
    df2 = pd.DataFrame(data2, columns=[f"{allele_1}_freq", f"{allele_2}_freq", "Population", "Super_Population"])
    ## plot the data

    sns.scatterplot(data=df2, x=f"{allele_1}_freq", y=f"{allele_2}_freq", hue="Super_Population")
    sns.regplot(data=df2, x=f"{allele_1}_freq", y=f"{allele_2}_freq", scatter=False, color='black')
    plt.xlabel(f"{allele_1} Frequency")
    plt.ylabel(f"{allele_2} Frequency")
    ## save the plot
    plt.savefig(f"./corr_{allele_1}_{allele_2}.png")

def correlation_analysis(HLA_allele_pop_freq_dict, KIR_allele_pop_freq_dict, pop_list):
    corr_num = 0
    data = []
    final_freq_dict = {}
    for HLA_allele in HLA_allele_pop_freq_dict:
        for KIR_allele in KIR_allele_pop_freq_dict:
            hla_freqs = [HLA_allele_pop_freq_dict[HLA_allele][pop] for pop in pop_list]
            kir_freqs = [KIR_allele_pop_freq_dict[KIR_allele][pop] for pop in pop_list]
            ## the freq in each pop should be >0
            if any(freq ==0 for freq in hla_freqs) or any(freq == 0 for freq in kir_freqs):
                continue
            final_freq_dict[HLA_allele] = hla_freqs
            final_freq_dict[KIR_allele] = kir_freqs
            ## get pearson correlation
            corr, pval = pearsonr(hla_freqs, kir_freqs)
            data.append((HLA_allele, KIR_allele, corr, pval))
            if pval < 0.05 :
                print(f"{HLA_allele} vs {KIR_allele}: Pearson r={corr:.3f}, p={pval:.3g}")
                # print (f"HLA allele frequencies: {hla_freqs}")
                # print (f"KIR allele frequencies: {kir_freqs}")
                # print ("##############################################")
                corr_num += 1

    print (f"Total number of significant correlations: {corr_num}")
    ## to df
    df = pd.DataFrame(data, columns=["HLA_allele", "KIR_allele", "Pearson_r", "p_value"])
    ## correct p-values

    p_values = df["p_value"].values
    rejected, pvals_corrected, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')  ##multipletests(p_values, method='bonferroni')
    df["p_value_corrected"] = pvals_corrected
    ## keep only significant associations
    df = df[df["p_value_corrected"] < 0.05]
    print(f"Number of significant associations after correction: {len(df)}")
    ## sort by p_value_corrected
    df = df.sort_values(by="p_value_corrected", ascending=True)

    print (df)
    plot(final_freq_dict, "IGKV2-29*01", "HLA-DPA1*01", pop_list, super_pop_dict)



hla_csv = "/home/shuaiw/methylation/data/hla/genotypes/data/HLA_Table_S6.csv"
kir_csv = "/home/shuaiw/methylation/data/hla/genotypes/data/KIR_Table_S7.csv"
cyp_csv = "/home/shuaiw/methylation/data/hla/genotypes/data/CYP_Table_S10.csv"
ig_tcr_csv = "/home/shuaiw/methylation/data/hla/genotypes/data/IG_TR_Table_S9.csv"
super_pop_file = "../1KG_ONT/hla/20131219.populations.tsv"
sample_pop_file = "../1KG_ONT/hla/20130606_sample_info.xlsx"

super_pop_dict = get_super_pop(super_pop_file)
sample_pop_dict, pop_set, sample_super_pop_dict, super_pop_set = get_sample_pop(sample_pop_file, super_pop_dict)
allele_count_dict = defaultdict(int)
print (super_pop_set)


all_allele_sample_dict = {}
all_sample_set = set()

TCR_pop_freq_dict, pop_list = read_ig_tcr_csv(ig_tcr_csv, sample_pop_dict, family="TR")
IG_pop_freq_dict, pop_list = read_ig_tcr_csv(ig_tcr_csv, sample_pop_dict, family="IG")
# CYP_allele_pop_freq_dict, pop_list = read_cyp_csv(cyp_csv, sample_pop_dict)
HLA_allele_pop_freq_dict, pop_list = read_tab(hla_csv, sample_pop_dict, family="HLA")
# KIR_allele_pop_freq_dict, pop_list = read_tab(kir_csv, sample_pop_dict, family="KIR")


correlation_analysis(IG_pop_freq_dict, HLA_allele_pop_freq_dict, pop_list)
# correlation_analysis(KIR_allele_pop_freq_dict, CYP_allele_pop_freq_dict, pop_list)

