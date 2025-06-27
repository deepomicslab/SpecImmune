import pandas as pd
from collections import defaultdict
from scipy.stats import fisher_exact
from sklearn.metrics import mutual_info_score
import numpy as np


def read_tab(csv, sample_pop_dict, pop, family="HLA"):
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
        if sample_pop_dict[sample] != pop:
            continue
        # if sample_pop_dict[sample] not in  ['SAS', 'AFR', 'EUR', 'AMR']: 
        #     continue
        # print (sample, row["One_guess"])
        ## get 1-field HLA
        if family == "HLA":
            allele = row["One_guess"].split(":")[0]
            if allele.split("*")[0] not in ["HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1"]:
                continue
        if family == "KIR":
            allele = row["One_guess"][:11]
        if allele.split("*")[0] in ["MICA", "MICB", "TAP2", "TAP1"]:
            allele = "HLA-" + allele
        allele_sample_dict[allele][sample] = 1
        allele_count_dict[allele] += 1
        sample_set.add(sample)
    # print (len(allele_sample_dict), "alleles in", family)
    return allele_sample_dict, sample_set


def get_mutual_info(allele_sample_dict, allele1, allele2, sample_set):
    allele1_binary = [1 if sample in allele_sample_dict[allele1] else 0 for sample in sample_set]
    allele2_binary = [1 if sample in allele_sample_dict[allele2] else 0 for sample in sample_set]
    

    mi, p = permutation_test_mi(allele1_binary, allele2_binary)
    return mi, p


def permutation_test_mi(x, y, n_permutations=1000):
    observed_mi = mutual_info_score(x, y)
    permuted_mis = [mutual_info_score(np.random.permutation(x), y) for _ in range(n_permutations)]
    pval = np.mean([mi >= observed_mi for mi in permuted_mis])
    return observed_mi, pval




def association_test(all_allele_sample_dict, all_sample_set):
    allele_list = list(all_allele_sample_dict.keys())
    # allele_list = allele_list[:5]
    data = []
    z = 0
    for i in range(len(allele_list)):
        for j in range(i + 1, len(allele_list)):
            allele1 = allele_list[i]
            allele2 = allele_list[j]
            if allele1.split("*")[0] == allele2.split("*")[0]:  ## remove same gene
                continue
            
            if allele1[:3] == allele2[:3]:  ## remove same family
                continue
            if allele1 not in all_allele_sample_dict or allele2 not in all_allele_sample_dict:
                continue
            # cutoff = 50
            # if allele_count_dict[allele1] < cutoff or allele_count_dict[allele2] < cutoff:
            #     continue
            # print (allele_count_dict[allele1], allele_count_dict[allele2], allele1, allele2)
            mi, p = get_mutual_info(all_allele_sample_dict, allele1, allele2, all_sample_set)
            data.append((allele1, allele2, mi, p,allele_count_dict[allele1], allele_count_dict[allele2]))
            print (f"{z}, {allele1} and {allele2}: mi={mi}, p={p}, Allele1_count={allele_count_dict[allele1]}, Allele2_count={allele_count_dict[allele2]}")
            z += 1
            # if p_value < 0.01 and oddsratio > 1:
            #     print(f"{allele1} and {allele2}: a={a}, b={b}, c={c}, d={d}, oddsratio={oddsratio}, p_value={p_value}")
    print (f"start correction..")
    ## perform multiple testing correction
    from statsmodels.stats.multitest import multipletests
    df = pd.DataFrame(data, columns=["Allele1", "Allele2", "mi", "p_value", "Allele1_count", "Allele2_count"])
    ## use rejected, pvals_corrected, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
    p_values = df["p_value"].values
    rejected, pvals_corrected, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
    df["p_value_corrected"] = pvals_corrected
    ## sort by p_value_corrected
    df = df.sort_values(by="p_value_corrected", ascending=True)
    ## print top 10 associations
    print (len(df), "associations tested")
    print(df.head(100))
    ## filter by p_value_corrected < 0.05
    # df_filtered = df[df["p_value_corrected"] < 0.05]
    # print(f"Number of associations with p_value_corrected < 0.05: {len(df_filtered)}")
    # ## save to file


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

hla_csv = "/home/shuaiw/methylation/data/hla/genotypes/data/HLA_Table_S6.csv"
kir_csv = "/home/shuaiw/methylation/data/hla/genotypes/data/KIR_Table_S7.csv"
super_pop_file = "../1KG_ONT/hla/20131219.populations.tsv"
sample_pop_file = "../1KG_ONT/hla/20130606_sample_info.xlsx"

super_pop_dict = get_super_pop(super_pop_file)
sample_pop_dict, pop_set, sample_super_pop_dict, super_pop_set = get_sample_pop(sample_pop_file, super_pop_dict)
allele_count_dict = defaultdict(int)
print (super_pop_set)
for pop in super_pop_set:
    all_allele_sample_dict = {}
    all_sample_set = set()
    HLA_allele_sample_dict, HLA_sample_set = read_tab(hla_csv, sample_super_pop_dict, pop, family="HLA")
    KIR_allele_sample_dict, KIR_sample_set = read_tab(kir_csv, sample_super_pop_dict, pop, family="KIR")
    all_allele_sample_dict.update(HLA_allele_sample_dict)
    all_allele_sample_dict.update(KIR_allele_sample_dict)
    all_sample_set.update(HLA_sample_set)
    all_sample_set.update(KIR_sample_set)
    print (pop, "pop...", len(all_allele_sample_dict), "alleles,",
           len(all_sample_set), "samples")
    association_test(all_allele_sample_dict, all_sample_set)
    break