import pandas as pd
from collections import defaultdict
from scipy.stats import fisher_exact
from sklearn.metrics import mutual_info_score
import numpy as np
import re


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
            cutoff = 50
            if allele_count_dict[allele1] < cutoff or allele_count_dict[allele2] < cutoff:
                continue
            # print (allele_count_dict[allele1], allele_count_dict[allele2], allele1, allele2)
            mi, p = get_mutual_info(all_allele_sample_dict, allele1, allele2, all_sample_set)
            data.append((allele1, allele2, mi, p,allele_count_dict[allele1], allele_count_dict[allele2]))
            print (f"{z}, {allele1} and {allele2}: mi={mi}, p={p}, Allele1_count={allele_count_dict[allele1]}, Allele2_count={allele_count_dict[allele2]}")
            z += 1
        #     break
        # break
    print (f"start correction..")
    ## perform multiple testing correction
    from statsmodels.stats.multitest import multipletests
    df = pd.DataFrame(data, columns=["Allele1", "Allele2", "mi", "p_value", "Allele1_count", "Allele2_count"])
    ## use rejected, pvals_corrected, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
    p_values = df["p_value"].values
    # rejected, pvals_corrected, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
    rejected, pvals_corrected, _, _ = multipletests(p_values, method='bonferroni')
    df["p_value_corrected"] = pvals_corrected
    ## sort by p_value_corrected
    df = df.sort_values(by="p_value_corrected", ascending=True)
    ## print top 10 associations
    print (len(df), "associations tested")
    print(df.head(100))
    ## filter by p_value_corrected < 0.05
    df_filtered = df[df["p_value_corrected"] < 100]
    # print(f"Number of associations with p_value_corrected < 0.05: {len(df_filtered)}")
    # ## save to file
    df_filtered.to_csv(f"/home/shuaiw/methylation/data/hla/genotypes/data/co_evo_mutual_2_field_{pop}.csv", index=False)


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


all_allele_sample_dict = {}
all_sample_set = set()
HLA_allele_pop_freq_dict, pop_list = read_tab(hla_csv, sample_pop_dict, family="HLA")
KIR_allele_pop_freq_dict, pop_list = read_tab(kir_csv, sample_pop_dict, family="KIR")
# print (HLA_allele_pop_freq_dict["HLA-A*01"])
from scipy.stats import pearsonr
corr_num = 0
data = []
for HLA_allele in HLA_allele_pop_freq_dict:
    for KIR_allele in KIR_allele_pop_freq_dict:
        hla_freqs = [HLA_allele_pop_freq_dict[HLA_allele][pop] for pop in pop_list]
        kir_freqs = [KIR_allele_pop_freq_dict[KIR_allele][pop] for pop in pop_list]
        ## the freq in each pop should be >0
        if any(freq ==0 for freq in hla_freqs) or any(freq == 0 for freq in kir_freqs):
            continue
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
from statsmodels.stats.multitest import multipletests
p_values = df["p_value"].values
rejected, pvals_corrected, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')  ##multipletests(p_values, method='bonferroni')
df["p_value_corrected"] = pvals_corrected
## keep only significant associations
df = df[df["p_value_corrected"] < 0.05]
print(f"Number of significant associations after correction: {len(df)}")
## sort by p_value_corrected
df = df.sort_values(by="p_value_corrected", ascending=True)

print (df)

    

# all_allele_sample_dict.update(HLA_allele_sample_dict)
# all_allele_sample_dict.update(KIR_allele_sample_dict)
# all_sample_set.update(HLA_sample_set)
# all_sample_set.update(KIR_sample_set)
# print (pop, "pop...", len(all_allele_sample_dict), "alleles,",
#         len(all_sample_set), "samples")
# association_test(all_allele_sample_dict, all_sample_set)
# # break